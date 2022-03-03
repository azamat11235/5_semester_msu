#include <math.h>
#include <cblas.h>
#include "qr.h"
#include "parameters.h"


void compute_params(double aii, double aji, double* c, double* s) {
    *c = aii / sqrt(aii*aii + aji*aji);
    *s = -aji / sqrt(aii*aii + aji*aji);
}

void rotate(double* xi, double* xj, double c, double s) {
    double xi_ = (*xi)*c - (*xj)*s;
    double xj_ = (*xi)*s + (*xj)*c;
    *xi = xi_;
    *xj = xj_;
}

void bcache(double* a, int na, double* cache, int i, int j, int k) {
    for (int ii = 0; ii < _b; ++ii) {
            int row_abs = i + ii;
            for (int jj = 0; jj < _b; ++jj) {
                int col_abs = j + jj;
                cache[k*_b*_b + ii*_b + jj] = a[row_abs*na + col_abs];
            }
        }
}

void bflush(double* a, int na, double* cache, int i, int j, int k) {
    for (int ii = 0; ii < _b; ++ii) {
            int row_abs = i + ii;
            for (int jj = 0; jj < _b; ++jj) {
                int col_abs = j + jj;
                a[row_abs*na + col_abs] = cache[k*_b*_b + ii*_b + jj];
            }
        }
}

void qr(double* a, double* q, int n) {
    double cache[4*_b*_b] = {0};
    for (int jb = 0; jb < n; jb += _b) {
        bcache(a, n, cache, jb, jb, 2); // кешируем диаг. блок
        // вращаем диаг. блок
        for (int j = 0; j < _b-1; ++j) {
            for (int i = j+1; i < _b; ++i) {
                int jj = 2*_b*_b + j*_b + j;
                int ij = 2*_b*_b + i*_b + j;
                double ajj = cache[jj];
                double aij = cache[ij];
                double c;
                double s;
                cblas_drotg(&ajj, &aij, &c, &s);
                cache[i*_b + j] = c;
                cache[j*_b + i] = s;
                cblas_drot(_b-j, &cache[jj], 1, &cache[ij], 1, c, s);
            }
        }
        bflush(q, n, cache, jb, jb, 0); // sin, cos
        // вращаем столбец (поддиаг. блоки)
        for (int ib = jb+_b; ib < n; ib += _b) {
            bcache(a, n, cache, ib, jb, 3); // поддиаг. блок
            for (int j = 0; j < _b; ++j) {
                for (int i = 0; i < _b; ++i) {
                    int jj = 2*_b*_b + j*_b + j;
                    int ij = 3*_b*_b + i*_b + j;
                    double ajj = cache[jj];
                    double aij = cache[ij];
                    double c;
                    double s;
                    cblas_drotg(&ajj, &aij, &c, &s);
                    cache[i*_b + j] = c;
                    cache[j*_b + i] = s;
                    cblas_drot(_b-j, &cache[jj], 1, &cache[ij], 1, c, s);
                }
            }
            bflush(q, n, cache, ib, jb, 0); // cos
            bflush(q, n, cache, jb, ib, 1); // sin
            bflush(a, n, cache, ib, jb, 3); // поддиаг. блок
        }
        bflush(a, n, cache, jb, jb, 2); // диаг. блок
        // обновляем строку (блоки справа от диаг.)
        for (int jb2 = jb+_b; jb2 < n; jb2 += _b) {
            bcache(a, n, cache, jb, jb2, 2); // внедиаг. блок
            bcache(q, n, cache, jb, jb, 0);  // cos, sin диаг. блока
            // вращения диаг. блока
            for (int j = 0; j < _b-1; ++j) {
                for (int i = j+1; i < _b; ++i) {
                    double c;
                    double s;
                    c = cache[i*_b + j];
                    s = cache[j*_b + i];
                    int jrow = 2*_b*_b + j*_b;
                    int irow = 2*_b*_b + i*_b;
                    cblas_drot(_b, &cache[jrow], 1, &cache[irow], 1, c, s);
                }
            }
            // вращения поддиаг. блоков
            for (int ib = jb+_b; ib < n; ib += _b) {
                bcache(a, n, cache, ib, jb2, 3); // внедиаг. блок (нижний)
                bcache(q, n, cache, ib, jb, 0);  // cos поддиаг. блока
                bcache(q, n, cache, jb, ib, 1);  // sin поддиаг. блока
                for (int j = 0; j < _b; ++j) {
                    for (int i = 0; i < _b; ++i) {
                        double c;
                        double s;
                        c = cache[i*_b + j];
                        s = cache[_b*_b + j*_b + i];
                        int jrow = 2*_b*_b + j*_b;
                        int irow = 3*_b*_b + i*_b;
                        cblas_drot(_b, &cache[jrow], 1, &cache[irow], 1, c, s);
                    }
                }
                bflush(a, n, cache, ib, jb2, 3); // внедиаг. блок (нижний)
            }
            bflush(a, n, cache, jb, jb2, 2); // недиаг. блок
        }
    }
 }

void restore_q(double* q, double* restored_q, int n) {
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            restored_q[n * i + j] = (double)(i == j);
     for (int jb = n - 1; jb >= 0; jb -= _b) {
        for (int ib = n - 1; ib >= 0; ib -= _b) {
            for (int j = 0; j < _b; ++j) {
                for (int i = 0; i < _b; ++i) {
                    int row_abs = ib - i;
                    int col_abs = jb - j;
                    if (row_abs > col_abs) {
                        double c = q[row_abs*n + col_abs];
                        double s = -q[col_abs*n + row_abs];
                        for (int k = 0; k < n; ++k)
                            rotate(&restored_q[col_abs*n + k], &restored_q[row_abs*n + k], c, s);
                    }
                }
            }
         }
     }
}

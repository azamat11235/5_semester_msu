#include <stdio.h>
/*
#include "qr.h"
#include "qr_batch.h"
#include "qr_cblas.h"
#include "qr_lapack.h"
#include "debug.h"
*/

#include <time.h>
#include <stdlib.h>
#include <string.h>
#include "routine.h"
#include "parameters.h"

typedef void qr_func_type(double*, double*, int);
typedef void restore_q_func_type(double*, double*, int);

void compute_params(double aii, double aji, double* c, double* s) {
    *c = aii / sqrt(aii * aii + aji * aji);
    *s = -aji / sqrt(aii * aii + aji * aji);
}

void rotate(double* xi, double* xj, double c, double s) {
    double xi_ = (*xi) * c - (*xj) * s;
    double xj_ = (*xi) * s + (*xj) * c;
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
                double ajj = cache[2*_b*_b + j*_b + j];
                double aij = cache[2*_b*_b + i*_b + j];
                double c;
                double s;
                compute_params(ajj, aij, &c, &s);
                cache[i*_b + j] = c;
                cache[j*_b + i] = s;
                for (int k = j; k < _b; ++k) {
                    int jk = 2*_b*_b + j*_b + k;
                    int ik = 2*_b*_b + i*_b + k;
                    rotate(&cache[jk], &cache[ik], c, s);
                }

            }
        }
        bflush(q, n, cache, jb, jb, 0); // sin, cos
        // вращаем столбец (поддиаг. блоки)
        for (int ib = jb+_b; ib < n; ib += _b) {
            bcache(a, n, cache, ib, jb, 3); // поддиаг. блок
            for (int j = 0; j < _b; ++j) {
                for (int i = 0; i < _b; ++i) {
                    double ajj = cache[2*_b*_b + j*_b + j];
                    double aij = cache[3*_b*_b + i*_b + j];
                    double c;
                    double s;
                    compute_params(ajj, aij, &c, &s);
                    cache[i*_b + j] = c;
                    cache[_b*_b + j*_b + i] = s;
                    for (int k = j; k < _b; ++k) {
                       int jk = 2*_b*_b + j*_b + k;
                       int ik = 3*_b*_b + i*_b + k;
                       rotate(&cache[jk], &cache[ik], c, s);
                    }
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
                    for (int k = 0; k < _b; ++k) {
                        int jk = 2*_b*_b + j*_b + k;
                        int ik = 2*_b*_b + i*_b + k;
                        rotate(&cache[jk], &cache[ik], c, s);
                    }
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
                        for (int k = 0; k < _b; ++k) {
                            int jk = 2*_b*_b + j*_b + k;
                            int ik = 3*_b*_b + i*_b + k;
                            rotate(&cache[jk], &cache[ik], c, s);
                        }
                    }
                }
                bflush(a, n, cache, ib, jb2, 3); // внедиаг. блок (нижний)
            }
            bflush(a, n, cache, jb, jb2, 2); // недиаг. блок
        }
    }
 }

double compute_time(qr_func_type* qr_func, int size) {
    const int num_of_iters = 5;
    double t = 0;
    double* matrix;
    double* q;
    allocMatrix(&q, size);
    allocMatrix(&matrix, size);
    for (int i = 0; i < num_of_iters; ++i) {
        clock_t start;
        clock_t end;
        fillMatrix(matrix, size);
        start = clock();
        (*qr_func)(matrix, q, size);
        end = clock();
        t += end - start;
    }
    free(matrix);
    free(q);
    return t / num_of_iters / CLOCKS_PER_SEC;
}

int main(int argv, char** argc) {
    for (int n = 256; n <= 2048; n *= 2)
        printf("size = %4d, time (qr_batch):   %f s.\n", n,  compute_time(&qr, n));
    /*
    printf("-------------------------------------------\n");
    printf("tests:\n");
    if (test(&qr, &restore_q) == 0)
        printf("qr: OK!\n");
    if (test(&qr_cblas, &restore_q) == 0)
        printf("qr_cblas: OK!\n");
    if (test(&qr_batch, &restore_q_batch) == 0)
        printf("qr_batch: OK!\n");
    printf("-------------------------------------------\n");

    for (int n = 256; n <= 2048; n *= 2)
        printf("size = %4d, time (qr_batch):   %f s.\n", n,  compute_time(&qr_batch, n));
    printf("-------------------------------\n");  
    for (int n = 256; n <= 2048; n *= 2)
        printf("size = %4d, time (qr):         %f s.\n", n,  compute_time(&qr, n));
    printf("-------------------------------------------\n");
    for (int n = 256; n <= 2048; n *= 2)
        printf("size = %4d, time (qr_lapack):  %f s.\n", n,  compute_time(&qr_lapack, n));
    printf("-------------------------------\n");
    for (int n = 256; n <= 2048; n *= 2)
        printf("size = %4d, time (qr_cblas):   %f s.\n", n, compute_time(&qr_cblas, n));
    printf("-------------------------------\n");
       */
    return 0;
}

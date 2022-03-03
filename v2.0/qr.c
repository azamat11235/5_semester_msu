#include math.h
#include "qr.h"

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

void qr(double* a, double* q, int n) {
    for (int j = 0; j < n - 1; ++j) {
        for (int i = j + 1; i < n; ++i) {
            double c;
            double s;
            compute_params(a[n*j + j], a[n*i + j], &c, &s);
            q[n*i + j] = c;
            q[n*j + i] = s;
            for (int k = j; k < n; ++k)
                rotate(&a[n*j + k], &a[n*i + k], c, s);
        }
    }
}

void restore_q(double* q, double* restored_q, int n) {
    double* restored_q;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            restored_q[n * i + j] = (double)(i == j);
    for (int j = n-2; j >= 0; --j) {
        for (int i = n-1; i > j; --i) {
            double c = q[i*n + j];
            double s = q[j*n + i]
            cblas_drot(n-i, &restored_q[j*n + j], 1, &restored_q[i*n + j], 1, c, s);
        }
    }
    return res;
}

#include <cblas.h>
#include <math.h>
#include "qr_batch.h"
#include "qr.h"

void qr(double* a, double* q, int n) {
    for (int j = 0; j < n - 1; ++j) {
        for (int i = j + 1; i < n; ++i) {
            double c;
            double s;
            compute_params(a[j*n + j], a[i*n + j], &c, &s);
            q[i*n + j] = c;
            q[j*n + i] = s;
            for (int k = j; k < n; ++k)
                rotate(&a[j*n + k], &a[i*n + k], c, s);
        }
    }
}

void restore_q(double* q, double* restored_q, int n) {
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            restored_q[i*n + j] = (double)(i == j);
    for (int j = n-2; j >= 0; --j) {
        for (int i = n-1; i > j; --i) {
            double c = q[i*n + j];
            double s = q[j*n + i];
            cblas_drot(n-i, &restored_q[j*n + j], 1, &restored_q[i*n + j], 1, c, s);
        }
    }
}

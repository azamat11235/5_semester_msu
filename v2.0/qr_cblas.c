#include <cblas.h>
#include "qr_cblas.h"

void qr_cblas(double* a, double* q, int n) {
    for (int icol = 0; icol < n - 1; ++icol) {
      for (int irow = icol + 1; irow < n; ++irow) {
          double c;
          double s;
          cblas_drotg(&a[n * icol + icol], &a[n * irow + icol], &c, &s);
          a[n * irow + icol] = 0;
          cblas_drot(n-icol-1, &a[icol*n + icol+1], 1, &a[irow*n + icol+1], 1, c, s);
          q[n*i + j] = c;
          q[n*j + i] = -s;
      }
    }
}

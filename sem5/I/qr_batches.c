#include <limits.h>
#include <time.h>
#include "qr_batches.h"

#define MAX(x, y) ((x > y) ? x : y)

void qr_batches(ELEMENT_TYPE* a, ELEMENT_TYPE* q, int n) {
    int q_size = n*(n-1);
    for (int irow_abs = 0; irow_abs < n; irow_abs += batchHeight) {
        for (int icol_abs = 0; icol_abs < n; icol_abs += batchWidth) {
            
            for (int icol_abs2 = 0; icol_abs2 < icol_abs; icol_abs2 += batchWidth) {
                for (int icol = 0; icol < batchWidth; ++icol) {
                    int j = icol_abs2 + icol;
                    for (int irow = MAX(j-irow_abs+1, 0); irow < batchHeight; ++irow) {
                        int i = irow_abs + irow;
                        ELEMENT_TYPE c = q[j*(2*n-1-j) + 2*(i-j-1)];
                        ELEMENT_TYPE s = -q[j*(2*n-1-j) + 2*(i-j-1)+1];
                        for (int k = 0; k < batchWidth && k+icol_abs < n; ++k)
                            rotate(&a[j*n + (icol_abs+k)], &a[i*n + (icol_abs+k)], c, s);
                    }
                }
            }
            
            for (int icol = 0; icol < batchWidth && icol+icol_abs < n; ++icol) {
                int j = icol_abs + icol;
                for (int irow =  MAX(j-irow_abs+1, 0); irow < batchHeight; ++irow) {
                    int i = irow_abs + irow;
                    ELEMENT_TYPE c;
                    ELEMENT_TYPE s;
                    compute_params(a[j*n + j], a[i*n + j], &c, &s);
                    q[j*(2*n-1-j) + 2*(i-j-1)] = c;
                    q[j*(2*n-1-j) + 2*(i-j-1)+1] = -s;
                    for (int k = 0; k < batchWidth && k+icol_abs < n; ++k)
                        rotate(&a[j*n + (icol_abs+k)], &a[i*n + (icol_abs+k)], c, s);
                }
            }
        }
    }
}

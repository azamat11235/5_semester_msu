#include <stdio.h>
/*
#include "qr.h"
#include "qr_batch.h"
#include "qr_cblas.h"
#include "qr_lapack.h"
#include "debug.h"
*/

#include "qr_batch.h"

#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include "routine.h"
#include "parameters.h"

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
        printf("size = %4d, time (qr_batch):   %f s.\n", n,  compute_time(&qr_batch, n));
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

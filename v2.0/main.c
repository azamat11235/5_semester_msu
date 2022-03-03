#include <stdio.h>
#include "qr.h"
#include "qr_batch.h"
#include "qr_cblas.h"
#include "qr_lapack.h"
#include "debug.h"

int main(int argv, char** argc) {
    printf("-------------------------------------------\n");
    if (test(&qr, &restore_q) == 0)
        printf("qr: all tests passed!\n");
    if (test(&qr_cblas, &restore_q) == 0)
        printf("qr_cblas: all tests passed!\n");
    if (test(&qr_batch, &restore_q_batch) == 0)
        printf("qr_batches: all tests passed!\n");
    printf("-------------------------------------------\n");

    for (int n = 256; n <= 2048; n *= 2)
        printf("size = %4d, time (qr_lapack):  %f s.\n", n,  compute_time(&qr_lapack, n));
    printf("-------------------------------\n");
    for (int n = 256; n <= 2048; n *= 2)
        printf("size = %4d, time (qr_cblas):   %f s.\n", n, compute_time(&qr_cblas, n));
    printf("-------------------------------\n");
    for (int n = 256; n <= 2048; n *= 2)
        printf("size = %4d, time (qr_batches): %f s.\n", n,  compute_time(&qr_batch, n));
    printf("-------------------------------\n");  
    for (int n = 256; n <= 2048; n *= 2)
        printf("size = %4d, time (qr):         %f s.\n", n,  compute_time(&qr, n));
    printf("-------------------------------------------\n");
        
    return 0;
}

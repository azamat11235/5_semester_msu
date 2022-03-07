#include <stdio.h>
#include "debug.h"
#include "qr.h"
#include "qr_omp.h"


#define LINE "------------------------------"

int main() {
    if (test(qr) == 0)
        printf("test (qr):     OK\n");
    if (test(qr_omp) == 0)
        printf("test (qr_omp): OK\n");
    printf(LINE);
    
    for (int n = 256; n <= 2048; n *= 2)
        printf("size = %4d, time = %f\n", n, compute_time(qr_omp, n));
    prinf(LINE);
    for (int n = 256; n <= 2048; n *= 2)
        printf("size = %4d, time = %f\n", n, compute_time(qr, n));
    return 0;
}

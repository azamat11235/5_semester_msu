#include <stdio.h>
#include "debug.h"
#include "qr.h"
#include "qr_omp.h"

int main() {
    if (test(qr) == 0)
        printf("test: OK\n");
    for (int n = 256; n <= 2048; n *= 2)
        printf("size = %4d, time = %f\n", n, compute_time(qr_omp, n));
    return 0;
}

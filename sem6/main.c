#include <stdio.h>
#include "debug.h"
#include "qr.h"

int main() {
    if (test(qr) == 0)
        printf("test: OK\n");
    for (int n = 256; n <= 2048; n *= 2)
        printf("size = %4d, time = %f\n", n, compute_time(qr, n));
    return 0;
}

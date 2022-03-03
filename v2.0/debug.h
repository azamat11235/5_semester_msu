#ifndef _DEBUG_
#define _DEBUG_

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include "routine.h"

typedef void qr_func_type(double*, double*, int);
typedef void restore_q_func_type(double*, double*, int);
double compute_time(qr_func_type* qr_func, int size);
int test(qr_func_type* qr_func, restore_q_func_type* restore_q_func);

#endif

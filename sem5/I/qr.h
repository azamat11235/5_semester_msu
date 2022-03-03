#ifndef _QR_
#define _QR_

void compute_params(double aii, double aji, double* c, double* s);
void rotate(double* xi, double* xj, double c, double s);
void qr(double* a, double* q, int n);
void restore_q(double* q, double* restored_q, int n);

#endif

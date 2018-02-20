#ifndef RHS_H
#define RHS_H

double rhs_f0 (double t, double x, double y);
double rhs_f1 (double t, double x, double y);

double rhs_test_f0 (double t, double x, double y);
double rhs_test_f1 (double t, double x, double y);

double p_drv_linear (double z);
double p_drv_polynomial (double z);

#endif /* RHS_H */

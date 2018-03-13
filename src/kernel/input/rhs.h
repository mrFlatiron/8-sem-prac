#ifndef RHS_H
#define RHS_H

#include "kernel/kernel_typedefs.h"

double rhs_f0 (double t, double x, double y, double mu, pressure_func_t pr);
double rhs_f1 (double t, double x, double y, double mu, pressure_func_t pr);
double rhs_f2 (double t, double x, double y, double mu, pressure_func_t pr);

double rhs_test_f0 (double t, double x, double y, double mu, pressure_func_t pr);
double rhs_test_f1 (double t, double x, double y, double mu, pressure_func_t pr);
double rhs_test_f2 (double t, double x, double y, double mu, pressure_func_t pr);

double rhs_test_f0_zero (double t, double x, double y, double mu, pressure_func_t pr);
double rhs_test_f1_zero (double t, double x, double y, double mu, pressure_func_t pr);
double rhs_test_f2_zero (double t, double x, double y, double mu, pressure_func_t pr);

double p_drv_linear (double z);
double p_drv_polynomial (double z);

#endif /* RHS_H */

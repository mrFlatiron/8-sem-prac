#ifndef T0_FUNCTIONS_H
#define T0_FUNCTIONS_H

double t0_vx  (double x, double y, double border_omega);
double t0_vy  (double x, double y, double border_omega);
double t0_g (double x, double y, double border_omega);
double t0_h (double x, double y, double border_omega);

double t0_vx_test  (double x, double y, double border_omega);
double t0_vy_test  (double x, double y, double border_omega);
double t0_g_test (double x, double y, double border_omega);
double t0_h_test (double x, double y, double border_omega);

double t0_vx_test_zero  (double x, double y, double border_omega);
double t0_vy_test_zero  (double x, double y, double border_omega);
double t0_g_test_zero (double x, double y, double border_omega);

#endif /* T0_FUNCTIONS_H */

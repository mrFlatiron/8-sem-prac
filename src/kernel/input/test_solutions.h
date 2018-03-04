#ifndef TEST_SOLUTIONS_H
#define TEST_SOLUTIONS_H

double test_g (double t, double x, double y);

double test_dg_dt (double t, double x, double y);
double test_dg_dx (double t, double x, double y);
double test_dg_dy (double t, double x, double y);

double test_vx (double t, double x, double y);

double test_dvx_dt (double t, double x, double y);
double test_dvx_dx (double t, double x, double y);
double test_dvx_dy (double t, double x, double y);

double test_dvx_dxdx (double t, double x, double y);
double test_dvx_dxdy (double t, double x, double y);
double test_dvx_dydy (double t, double x, double y);

double test_vy (double t, double x, double y);

double test_dvy_dt (double t, double x, double y);
double test_dvy_dx (double t, double x, double y);
double test_dvy_dy (double t, double x, double y);

double test_dvy_dxdx (double t, double x, double y);
double test_dvy_dxdy (double t, double x, double y);
double test_dvy_dydy (double t, double x, double y);
#endif /* TEST_SOLUTIONS_H */

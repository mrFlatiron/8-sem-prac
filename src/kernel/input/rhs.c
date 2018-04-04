#include "rhs.h"
#include "test_solutions.h"
#include "common/debug_utils.h"

#include <math.h>
#include "common/math_utils.h"

double rhs_f0 (double t, double x, double y, double mu, pressure_func_t pr)
{
  FIX_UNUSED (t);
  FIX_UNUSED (x);
  FIX_UNUSED (y);
  FIX_UNUSED (mu);
  FIX_UNUSED (pr);
  return 0;
}

double rhs_f1 (double t, double x, double y, double mu, pressure_func_t pr)
{
  FIX_UNUSED (t);
  FIX_UNUSED (x);
  FIX_UNUSED (y);
  FIX_UNUSED (mu);
  FIX_UNUSED (pr);
  return 0;
}

double rhs_f2 (double t, double x, double y, double mu, pressure_func_t pr)
{
  FIX_UNUSED (t);
  FIX_UNUSED (x);
  FIX_UNUSED (y);
  FIX_UNUSED (mu);
  FIX_UNUSED (pr);
  return 0;
}

double p_drv_linear (double z)
{
  FIX_UNUSED (z);
  return 10;
}

double p_drv_polynomial (double z)
{
  return 1.4 * pow (z, 0.4);
}

double p_drv (double z, pressure_func_t pr)
{
  switch (pr)
    {
    case pressure_linear:
      return p_drv_linear (z);
    case pressure_polynomial:
      return p_drv_polynomial (z);
    }

  ASSERT_RETURN (0 , 0);
}

double rhs_test_f0 (double t, double x, double y, double mu, pressure_func_t pr)
{
  FIX_UNUSED (pr);
  FIX_UNUSED (mu);

  return
      + test_dg_dt (t, x, y)
      + 0.5 * (
        + test_vx (t, x, y) * test_dg_dx (t, x, y)
        + test_dvx_dx (t, x, y) * test_g (t, x, y)
        + test_dg_dx (t, x, y) * test_vx (t, x, y)
        + (2 - test_g (t, x, y)) * test_dvx_dx (t, x, y))
      + 0.5 * (
        + test_vy (t, x, y) * test_dg_dy (t, x, y)
        + test_dvy_dy (t, x, y) * test_g (t, x, y)
        + test_dg_dy (t, x, y) * test_vy (t, x, y)
        + (2 - test_g (t, x, y)) * test_dvy_dy (t, x, y));
}

double rhs_test_f1 (double t, double x, double y, double mu, pressure_func_t pr)
{

  return
      + test_dvx_dt (t, x, y)
      + (1. / 3.) * (test_vx (t, x, y) * test_dvx_dx (t, x, y)
                     + 2 * test_vx (t, x,y) * test_dvx_dx (t, x ,y))
      + (1. / 2.) * (test_vy (t, x, y) * test_dvx_dy (t, x, y)
                     + test_dvy_dy (t, x, y) * test_vx (t, x, y) + test_dvx_dy (t, x, y) * test_vy (t, x, y)
                     - test_vx (t, x, y) * test_dvy_dy (t, x, y))
      + p_drv (exp (test_g (t, x, y)), pr) * test_dg_dx (t, x, y)
      - (mu / exp (test_g (t, x, y))) * ((4. / 3.) * test_dvx_dxdx (t, x, y)
                                         + test_dvx_dydy (t, x, y)
                                         + (1. / 3.) * (test_dvy_dxdy (t, x, y)));

  return
      + test_dvx_dt (t, x, y)
      + test_vx (t, x, y) * test_dg_dt (t, x, y)
      + 2 * test_vx (t, x, y) * test_dvx_dx (t, x, y)
      + test_vx (t, x, y) * test_vx (t, x, y) * test_dg_dx (t, x, y)
      + test_vy (t, x, y) * test_dvx_dy (t, x, y)
      + test_vx (t, x, y) * test_dvy_dy (t, x, y)
      + test_vx (t, x, y) * test_vy (t, x, y) * test_dg_dy (t, x ,y)
      + p_drv (exp (test_g (t, x, y)), pr) * test_dg_dx (t, x, y)
      - (4. / 3.) * mu / exp (test_g (t, x ,y)) * test_dvx_dxdx (t,  x, y)
      - mu / exp (test_g (t, x, y)) * test_dvx_dydy (t, x, y)
      - (1. / 3.) * mu / exp (test_g (t, x, y)) * test_dvy_dxdy (t, x, y);

}



double rhs_test_f2 (double t, double x, double y, double mu, pressure_func_t pr)
{
  return
      + test_dvy_dt (t, x, y)
      + (1. / 3.) * (test_vy (t, x, y) * test_dvy_dy (t, x, y)
                     + 2 * test_vy (t, x,y) * test_dvy_dy (t, x ,y))
      + (1. / 2.) * (test_vx (t, x, y) * test_dvy_dx (t, x, y)
                     + test_dvx_dx (t, x, y) * test_vy (t, x, y) + test_dvy_dx (t, x, y) * test_vx (t, x, y)
                     - test_vy (t, x, y) * test_dvx_dx (t, x, y))
      + p_drv (exp (test_g (t, x, y)), pr) * test_dg_dy (t, x, y)
      - (mu / exp (test_g (t, x, y))) * ((4. / 3.) * test_dvy_dydy (t, x, y)
                                         + test_dvy_dxdx (t, x, y)
                                         + (1. / 3.) * (test_dvx_dxdy (t, x, y)));
}

double rhs_test_f0_zero (double t, double x, double y, double mu, pressure_func_t pr)
{
  FIX_UNUSED (t);
  FIX_UNUSED (x);
  FIX_UNUSED (y);
  FIX_UNUSED (mu);
  FIX_UNUSED (pr);
  return 0;
}

double rhs_test_f1_zero (double t, double x, double y, double mu, pressure_func_t pr)
{
  FIX_UNUSED (t);
  FIX_UNUSED (x);
  FIX_UNUSED (y);
  FIX_UNUSED (mu);
  FIX_UNUSED (pr);
  /*return 0;*/
  /*return cos (x) * p_drv (exp (test_g_zero (t, x, y)), pr);*/
  return sin (y) * cos (x) * p_drv (exp (test_g_zero (t, x ,y)), pr);
}

double rhs_test_f2_zero (double t, double x, double y, double mu, pressure_func_t pr)
{
  FIX_UNUSED (t);
  FIX_UNUSED (x);
  FIX_UNUSED (y);
  FIX_UNUSED (mu);
  FIX_UNUSED (pr);
  /*return 0;*/
  /*return cos (y) * p_drv (exp (test_g_zero (t, x, y)), pr);*/
  return sin (x) * cos (y) * p_drv (exp (test_g_zero (t, x ,y)), pr);
}

double rhs_test_sok_f0 (double t, double x, double y, double mu, pressure_func_t pr)
{
  FIX_UNUSED (mu);
  FIX_UNUSED (pr);
  return
      + test_dg_dt (t, x, y)
      + (
        + test_dg_dx (t, x, y) * test_vx (t, x, y)
        + test_dvx_dx (t, x, y) * test_g (t, x, y)
        )
      + (
        + test_dg_dy (t, x, y) * test_vy (t, x, y)
        + test_dvy_dy (t, x, y) * test_g (t, x, y)
        );
}

double rhs_test_sok_f1 (double t, double x, double y, double mu, pressure_func_t pr)
{
  FIX_UNUSED (pr);
  ASSERT_RETURN (!math_is_null (test_g (t, x, y)), 0);
  return
      + 1 / test_g (t, x, y) * (
        + (
          + test_dg_dt (t, x, y) * test_vx (t, x, y)
          + test_dvx_dt (t, x, y) * test_g (t, x, y)
          )
        + (
          + test_dg_dx (t, x, y) * test_vx (t, x, y) * test_vx (t, x, y)
          + 2 * test_g (t, x, y) * test_dvx_dx (t, x, y) * test_vx (t, x, y)
          )
        + (
          + test_dg_dy (t, x, y) * test_vx (t, x, y) * test_vy (t, x, y)
          + test_g (t, x, y) * test_dvx_dy (t, x, y) * test_vy (t, x, y)
          + test_g (t, x, y) * test_vx (t, x, y) * test_dvy_dy (t, x, y)
          )
        + (
          + test_dg_dx (t, x, y)
          )
        - mu * (
          + (4. / 3.) * test_dvx_dxdx (t, x, y)
          + test_dvx_dydy (t, x, y)
          + (1. / 3.) * test_dvy_dxdy (t, x, y)
          )
        );
}

double rhs_test_sok_f2 (double t, double x, double y, double mu, pressure_func_t pr)
{
  FIX_UNUSED (pr);
  ASSERT_RETURN (!math_is_null (test_g (t, x, y)), 0);
  return
      + 1 / test_g (t, x, y) * (
        + (
          + test_dg_dt (t, x, y) * test_vy (t, x, y)
          + test_dvy_dt (t, x, y) * test_g (t, x, y)
          )
        + (
          + test_dg_dy (t, x, y) * test_vy (t, x, y) * test_vy (t, x, y)
          + 2 * test_g (t, x, y) * test_dvy_dy (t, x, y) * test_vy (t, x, y)
          )
        + (
          + test_dg_dx (t, x, y) * test_vx (t, x, y) * test_vy (t, x, y)
          + test_g (t, x, y) * test_dvx_dx (t, x, y) * test_vy (t, x, y)
          + test_g (t, x, y) * test_vx (t, x, y) * test_dvy_dx (t, x, y)
          )
        + (
          + test_dg_dy (t, x, y)
          )
        - mu * (
          + (4. / 3.) * test_dvy_dydy (t, x, y)
          + test_dvy_dxdx (t, x, y)
          + (1. / 3.) * test_dvx_dxdy (t, x, y)
          )
        );
}

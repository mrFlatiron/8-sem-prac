#include "test_solutions.h"

#include "common/debug_utils.h"
#include <math.h>

double test_g (double t, double x, double y)
{
  return cos (2 * x) * cos (2 * y) * exp (t);
}

double test_dg_dt (double t, double x, double y)
{
  return test_g (t, x, y);
}

double test_dg_dx (double t, double x, double y)
{
  return -2 * sin (2 * x) * cos (2 * y) * exp (t);
}

double test_dg_dy (double t, double x, double y)
{
  return -2 * sin (2 * y) * cos (2 * x) * exp (t);
}

double test_vx (double t, double x, double y)
{
  return sin (x) * sin (y) * exp (t);
}

double test_dvx_dt (double t, double x, double y)
{
  return test_vx (t, x, y);
}

double test_dvx_dx (double t, double x, double y)
{
  return cos (x) * sin (y) * exp (t);
}

double test_dvx_dxdx (double t, double x, double y)
{
  return -sin (x) * sin (y) * exp (t);
}

double test_dvx_dxdy(double t, double x, double y)
{
  return cos (x) * cos (y) * exp (t);
}

double test_dvx_dy (double t, double x, double y)
{
  return sin (x) * cos (y) * exp (t);
}

double test_dvx_dydy (double t, double x, double y)
{
  return - sin (x) * sin (y) * exp (t);
}

double test_vy (double t, double x, double y)
{
  return sin (x) * sin (y) * exp (-t);
}

double test_dvy_dt (double t, double x, double y)
{
  return -test_vy (t, x, y);
}

double test_dvy_dx (double t, double x, double y)
{
  return cos (x) * sin (y) * exp (-t);
}

double test_dvy_dxdx (double t, double x, double y)
{
  return -sin (x) * sin (y) * exp (-t);
}

double test_dvy_dxdy (double t, double x, double y)
{
  return cos (x) * cos (y) * exp (-t);
}

double test_dvy_dy (double t, double x, double y)
{
  return sin (x) * cos (y) * exp (-t);
}

double test_dvy_dydy (double t, double x, double y)
{
  return - sin (x) * sin (y) * exp (-t);
}










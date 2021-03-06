#include "t0_functions.h"
#include "common/debug_utils.h"
#include "kernel/kernel_typedefs.h"
#include "test_solutions.h"
#include "common/math_utils.h"


double t0_vx (double x, double y, double border_omega)
{
  FIX_UNUSED (x);
  FIX_UNUSED (y);

  if (math_is_null (x))
    return border_omega;

  return 0;

}

double t0_vy (double x, double y, double border_omega)
{
  FIX_UNUSED (x);
  FIX_UNUSED (y);
  FIX_UNUSED (border_omega);
  return 0;
}

double t0_g (double x, double y, double border_omega)
{
  FIX_UNUSED (x);
  FIX_UNUSED (y);
  FIX_UNUSED (border_omega);
  return log (t0_h (x, y, border_omega));
}

double t0_vx_test (double x, double y, double border_omega)
{
  FIX_UNUSED (border_omega);
  return test_vx (0, x, y);
}

double t0_vy_test (double x, double y, double border_omega)
{
  FIX_UNUSED (border_omega);
  return test_vy (0, x, y);
}

double t0_g_test (double x, double y, double border_omega)
{
  FIX_UNUSED (border_omega);
  return test_g (0, x, y);
}

double t0_h_test (double x, double y, double border_omega)
{
  FIX_UNUSED (border_omega);
  return test_h (0, x, y);
}

double t0_vx_test_zero (double x, double y, double border_omega)
{
  FIX_UNUSED (border_omega);
  return test_vx_zero (0, x, y);
}

double t0_vy_test_zero (double x, double y, double border_omega)
{
  FIX_UNUSED (border_omega);
  return test_vy_zero (0, x, y);
}

double t0_g_test_zero (double x, double y, double border_omega)
{
  FIX_UNUSED (border_omega);
  return test_g_zero (0, x, y);
}

double t0_h (double x, double y, double border_omega)
{
  FIX_UNUSED (x);
  FIX_UNUSED (y);
  FIX_UNUSED (border_omega);
  if (x < MY_PI / 5)
    return RHO_LEFTMOST;
  else
    return 0.1;
}

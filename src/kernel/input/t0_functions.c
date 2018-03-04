#include "t0_functions.h"
#include "common/debug_utils.h"
#include "kernel/kernel_typedefs.h"
#include "test_solutions.h"
#include "common/math_utils.h"


double t0_vx (double x, double y)
{
  if (x <= M_PI / 4)
    return BORDER_OMEGA * cos (x);

  return BORDER_OMEGA * sin (x) * sin (y) / sqrt (2) * 2.;

  return 0;
}

double t0_vy (double x, double y)
{
  return 0;
}

double t0_g (double x, double y)
{
  return log (RHO_LEFTMOST);
}

double t0_vx_test (double x, double y)
{
  return test_vx (0, x, y);
}

double t0_vy_test (double x, double y)
{
  return test_vy (0, x, y);
}

double t0_g_test (double x, double y)
{
  return test_g (0, x, y);
}

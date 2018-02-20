#include "rhs.h"
#include "common/debug_utils.h"
#include <math.h>

double rhs_f0 (double t, double x, double y)
{
  FIX_UNUSED (t);
  FIX_UNUSED (x);
  FIX_UNUSED (y);
  return 0;
}

double rhs_f1 (double t, double x, double y)
{
  FIX_UNUSED (t);
  FIX_UNUSED (x);
  FIX_UNUSED (y);
  return 0;
}

double p_drv_linear (double z)
{
  FIX_UNUSED (z);
  return 1.4;
}

double p_drv_polynomial (double z)
{
  return 1.4 * pow (z, 0.4);
}

#include "tests.h"
#include "common/math_utils.h"
#include "common/debug_utils.h"

int test_eq_double (double actual, double must, double prec)
{
  if (!math_eq_w_prec (actual, must, prec))
    {
      double dif = fabs (actual - must);
      FIX_UNUSED (dif);
      DEBUG_PAUSE ("Test Failed");
      return 0;
    }
  return 1;
}

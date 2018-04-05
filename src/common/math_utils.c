#include "math_utils.h"

#include <math.h>

int math_fuzzy_eq (double lhs, double rhs)
{
  return (fabs (lhs - rhs) <= EQ_BOUNDARY);
}

int math_is_null (double a)
{
  return math_fuzzy_eq (a, 0);
}

int math_sign (double x)
{
  if (math_is_null (x))
    return 0;

  if (x < 0)
    return -1;

  return 1;
}

int math_is_pos (double x)
{
  if (math_is_null (x))
    return 0;

  if (x > 0)
    return 1;

  return 0;
}

int math_is_neg (double x)
{
  if (math_is_null (x))
    return 0;

  if (x < 0)
    return 1;

  return 0;
}

int math_eq_w_prec (double lhs, double rhs, double prec)
{
  if (fabs (lhs - rhs) <= prec)
    return 1;

  return 0;
}

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

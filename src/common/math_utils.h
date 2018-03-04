#ifndef MATH_UTILS_H
#define MATH_UTILS_H
#include <math.h>


#define EQ_BOUNDARY 1e-6

int math_fuzzy_eq (double lhs, double rhs);

int math_is_null (double a);

#endif /* MATH_UTILS_H */

#ifndef MATH_UTILS_H
#define MATH_UTILS_H
#include <math.h>


#define EQ_BOUNDARY 1e-6

int math_fuzzy_eq (double lhs, double rhs);

int math_is_null (double a);

int math_sign (double x);

int math_is_pos (double x);

int math_is_neg (double x);

#endif /* MATH_UTILS_H */

#ifndef VECTOR_OPS_H
#define VECTOR_OPS_H

#include "common/vectors_fwd.h"

double dot_product (vector_double_t a, vector_double_t b, int size);

void linear_combination_1 (const vector_double_t a,
                           double coef, const vector_double_t b,
                           vector_double_t out,
                           int size);

void linear_combination_w_override_1 (vector_double_t a,
                           double coef, const vector_double_t b,
                           int size);

double l2_norm (const vector_double_t a, int size);

double c_norm (const vector_double_t a, int size);

double c_norm_w_index (const vector_double_t a, int size, int *index);

#endif /* VECTOR_OPS_H */

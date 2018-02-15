#ifndef MSR_MATRIX_H
#define MSR_MATRIX_H

#include "common/vectors_fwd.h"

typedef struct
{
  vector_double_t AA;
  vector_int_t JA;

  int matrix_size;
  int array_size;

} msr_matrix;

int msr_init (msr_matrix *matrix,
              int N,
              int non_zero_el /*>= N*/);

void msr_destroy (msr_matrix *matrix);

double msr_ij (const msr_matrix *matrix, int i, int j);

/* Below are private methods */

void msr_get_ja_row_bounds (const msr_matrix *matrix, int i, int *begin, int *end);



#endif /* MSR_MATRIX_H */

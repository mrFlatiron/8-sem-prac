#include "msr_matrix.h"

#include "common/vectors.h"

int msr_init (msr_matrix *matrix, int N, int non_zero_el)
{
  matrix->matrix_size = N;
  matrix->array_size = non_zero_el + 1;
  matrix->AA = VECTOR_CREATE (double, matrix->array_size);
  matrix->JA = VECTOR_CREATE (int, matrix->array_size);

  if (!matrix->AA || !matrix->JA)
    return 1;

  return 0;
}

void msr_destroy (msr_matrix *matrix)
{
  VECTOR_DESTROY (matrix->AA);
  VECTOR_DESTROY (matrix->JA);
}

double msr_ij (const msr_matrix *matrix, int i, int j)
{
  int begin;
  int end;
  int ja_iter;

  msr_get_ja_row_bounds (matrix, i, &begin, &end);

  for (ja_iter = begin; ja_iter < end; ja_iter++)
    if (matrix->JA[ja_iter] == j)
      return matrix->AA[ja_iter];

  return 0;
}

void msr_get_ja_row_bounds (const msr_matrix *matrix, int i, int *begin, int *end)
{
  *begin = matrix->JA[i];
  *end = matrix->JA[i + 1];
}

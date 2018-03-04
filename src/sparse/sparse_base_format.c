#include "sparse_base_format.h"

#include "common/vectors.h"

int sparse_base_init (sparse_base_format *matr, int N, int max_row_nz)
{
  matr->N = N;
  matr->nnz_in_rows     = VECTOR_CREATE (int, N);
  matr->column_indecies = VECTOR_CREATE (int, max_row_nz * N);
  matr->values          = VECTOR_CREATE (double, max_row_nz * N);

  return 0;
}

void sparse_base_destroy (sparse_base_format *matr)
{
  VECTOR_DESTROY (matr->nnz_in_rows);
  VECTOR_DESTROY (matr->column_indecies);
  VECTOR_DESTROY (matr->values);
}

void sparse_base_add_row (sparse_base_format *matr, int row, int cols[], double values[], int nnz)
{
  int start = 0;
  int i;

  for (i = 0; i < row; i++)
    start += matr->nnz_in_rows[i];

  matr->nnz_in_rows[row] = nnz;

  for (i = 0; i < nnz; i++)
    {
      matr->column_indecies[i] = cols[i];
      matr->values[i] = values[i];
    }
}

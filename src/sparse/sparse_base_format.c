#include "sparse_base_format.h"

#include "common/vectors.h"

int sparse_base_init (sparse_base_format *matr, int N, int max_row_nz)
{
  matr->N = N;
  matr->nnz_total = 0;
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
  int start = matr->nnz_total;
  int i;

  matr->nnz_in_rows[row] = nnz;

  for (i = 0; i < nnz; i++)
    {
      matr->column_indecies[start + i] = cols[i];
      matr->values[start + i] = values[i];
    }

  matr->nnz_total += nnz;
}

int sparse_base_nnz_total (const sparse_base_format *matr)
{
  return matr->nnz_total;
}

void sparse_base_to_init_state (sparse_base_format *matr)
{
  matr->nnz_total = 0;
}

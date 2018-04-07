#include "sparse_base_format.h"
#include "common/math_utils.h"
#include "common/vectors.h"
#include "common/debug_utils.h"

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

void sparse_base_fill_nz (double *vals, int *cols, int *nnz, double val, int col)
{
  if (math_is_null (val))
    return;

  vals[*nnz] = val;
  cols[*nnz] = col;
  (*nnz)++;
}

void sparse_base_fill_nz_s (nz_row_t *nz_row, double val, double col)
{
  sparse_base_fill_nz (nz_row->vals, nz_row->cols, &nz_row->nnz, val, col);
}

void sparse_base_add_row_s (sparse_base_format *mat, const nz_row_t *nz_row)
{
  sparse_base_add_row (mat, nz_row->row, nz_row->cols, nz_row->vals, nz_row->nnz);
}

void nz_row_init (nz_row_t *nz_row, int max_row_nz)
{
  nz_row->vals = VECTOR_CREATE (double, max_row_nz);
  nz_row->cols = VECTOR_CREATE (int, max_row_nz);
  nz_row->nnz = 0;
  nz_row->row = 0;
}

void nz_row_destroy (nz_row_t *nz_row)
{
  VECTOR_DESTROY (nz_row->vals);
  VECTOR_DESTROY (nz_row->cols);
}

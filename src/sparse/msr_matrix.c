#include "msr_matrix.h"

#include "common/vectors.h"
#include "common/math_utils.h"

#include "common/debug_utils.h"

int msr_init (msr_matrix *matrix, int N, int non_zero_el)
{
  matrix->matrix_size = N;
  matrix->array_size = non_zero_el + 1;
  matrix->AA = VECTOR_CREATE (double, matrix->array_size);
  matrix->JA = VECTOR_CREATE (int, matrix->array_size);

  if (!matrix->AA || !matrix->JA)
    return 1;

  matrix->AA[N] = 0;

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

  if (i == j)
    return matrix->AA[i];

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

void msr_mult_vector (const msr_matrix *matrix, const vector_double_t src, vector_double_t dest)
{
  int i;
  int ja_iter;
  int begin_row;
  int end_row;
  int N;

  N = matrix->matrix_size;
  VECTOR_SET (double, dest, 0, N);

  for (i = 0; i < N; i++)
    dest[i] = matrix->AA[i] * src[i];

  for (i = 0; i < N; i++)
    {
      msr_get_ja_row_bounds (matrix, i, &begin_row, &end_row);

      for (ja_iter = begin_row; ja_iter < end_row; ja_iter++)
        dest[i] += matrix->AA[ja_iter] * src[matrix->JA[ja_iter]];
    }
}

int msr_init_from_vector (msr_matrix *matrix, const vector_double_t dense_matrix, int N)
{
  int i;
  int j;
  int k;
  int nnz = 0;
  int ja_iter = 0;
  int zero_row_begin;
  int zero_row_end;

  for (i = 0; i < N; i++)
    {
      for (j = 0; j < N; j++)
        {
          if (!math_is_null (dense_matrix[j + i * N]))
              nnz++;
        }
    }

  if (msr_init (matrix, N, nnz))
    return 1;

  for (i = 0; i < N; i++)
    {
      for (j = 0; j < N; j++)
        {
          double value = dense_matrix[j + i * N];

          if (math_is_null (value) && (i == j))
            return 2;

          if (math_is_null (value))
            continue;

          if (i == j)
            matrix->AA[i] = value;
          else
            {
              matrix->AA[N + 1 + ja_iter] = value;
              matrix->JA[N + 1 + ja_iter] = j;
              ja_iter++;
            }
        }
    }

  ja_iter = N + 1;

  zero_row_begin = 0;
  zero_row_end = 0;
  for (i = 0; i < N; i++)
    {
      int nz_found = 0;

      for (j = 0; j < N; j++)
        {
          if (i == j)
            continue;

          if (!math_is_null (dense_matrix[j + i * N]))
            nz_found++;

        }

      if (nz_found)
        {
          matrix->JA[i] = ja_iter;
          matrix->JA[i + 1] = ja_iter + nz_found;
          ja_iter += nz_found;

          for (k = zero_row_begin; k < zero_row_end; k++)
            matrix->JA[k] = matrix->JA[i];

          zero_row_begin = i + 1;
          zero_row_end = i + 1;
        }
      else
        zero_row_end++;
    }

  for (k = zero_row_begin; k < zero_row_end; k++)
    matrix->JA[k] = matrix->array_size;

  matrix->JA[N] = matrix->array_size;

  return 0;
}

void msr_dump (const msr_matrix *matrix, FILE *fout)
{
  int i;
  int j;
  fprintf (fout, "MSR matrix dump:\n"
                 "N          = %d\n"
                 "array_size = %d\n",
           matrix->matrix_size, matrix->array_size);

  fprintf (fout, "AA = ");

  for (i = 0; i < matrix->array_size; i++)
    fprintf (fout, "%.3f\t", matrix->AA[i]);

  fprintf (fout, "\nJA = ");

  for (i = 0; i < matrix->array_size; i++)
    fprintf (fout, "%d\t", matrix->JA[i]);

  fprintf (fout, "\ndata:\n");

  for (i = 0; i < matrix->matrix_size; i++)
    {
      for (j = 0; j < matrix->matrix_size; j++)
        fprintf (fout, "%.3f\t", msr_ij (matrix, i, j));

      fprintf (fout, "\n");
    }

  fflush (fout);
}

int msr_reinit (msr_matrix *matrix, int N, int non_zero_el)
{
  if (non_zero_el + 1 > matrix->array_size)
    {
      msr_destroy (matrix);
      msr_init (matrix, N, non_zero_el);
    }
  else
    {
      matrix->matrix_size = N;
      matrix->array_size = non_zero_el + 1;
    }

  return 0;
}

int msr_fill_from_sparse_base (msr_matrix *matrix, const sparse_base_format *matrix_base)
{
  int nnz = 0;
  int i;
  int ja_iter = 0;
  int k;
  int N = matrix_base->N;
  int zero_row_begin;
  int zero_row_end;

  msr_reinit (matrix, N, sparse_base_nnz_total (matrix_base));

  nnz = 0;
  for (i = 0; i < N; i++)
    {
      int j;

      for (j = 0; j < matrix_base->nnz_in_rows[i]; j++)
        {
          int col;
          double value = matrix_base->values[j + nnz];

          if (math_is_null (value))
            return 2;

          col = matrix_base->column_indecies[nnz + j];
          if (i == col)
            matrix->AA[i] = value;
          else
            {
              matrix->AA[N + 1 + ja_iter] = value;
              matrix->JA[N + 1 + ja_iter] = col;
              ja_iter++;
            }
        }
      nnz += matrix_base->nnz_in_rows[i];
    }

  ja_iter = N + 1;

  zero_row_begin = 0;
  zero_row_end = 0;
  nnz = 0;
  for (i = 0; i < N; i++)
    {
      int nz_found = matrix_base->nnz_in_rows[i] - 1;

      if (nz_found)
        {
          matrix->JA[i] = ja_iter;
          matrix->JA[i + 1] = ja_iter + nz_found;
          ja_iter += nz_found;

          for (k = zero_row_begin; k < zero_row_end; k++)
            matrix->JA[k] = matrix->JA[i];

          zero_row_begin = i + 1;
          zero_row_end = i + 1;
        }
      else
        zero_row_end++;

      nnz += matrix_base->nnz_in_rows[i];
    }

  for (k = zero_row_begin; k < zero_row_end; k++)
    matrix->JA[k] = matrix->array_size;

  matrix->JA[N] = matrix->array_size;

  return 0;
}

int msr_init_empty (msr_matrix *matrix)
{
  matrix->AA = NULL;
  matrix->JA = NULL;

  matrix->array_size = 0;
  matrix->matrix_size = 0;

  return 0;
}

#ifndef SPARSE_BASE_FORMAT_H
#define SPARSE_BASE_FORMAT_H

#include "common/vectors_fwd.h"

typedef struct
{
  int N;
  int nnz_total;
  vector_int_t nnz_in_rows; /*N size*/
  vector_double_t values;   /*nnz size*/
  vector_int_t column_indecies; /*nnz size*/
} sparse_base_format;

typedef struct
{
  double *vals;
  int *cols;
  int nnz;
  int row;
} nz_row_t;

int sparse_base_init (sparse_base_format *matr,
                      int N,
                      int max_row_nz);

void sparse_base_to_init_state (sparse_base_format *matr);

void sparse_base_add_row (sparse_base_format *matr,
                          int row,
                          int cols[],
                          double values[],
                          int nnz);

int sparse_base_nnz_total (const sparse_base_format *matr);

void sparse_base_destroy (sparse_base_format *matr);

void sparse_base_fill_nz (double *vals, int *cols, int *nnz, double val, int col);

void sparse_base_fill_nz_s (nz_row_t *nz_row, double val, double col);
void sparse_base_add_row_s (sparse_base_format *mat,
                            const nz_row_t *nz_row);

void nz_row_init (nz_row_t *nz_row, int max_row_nz);
void nz_row_destroy (nz_row_t *nz_row);


#endif /* SPARSE_BASE_FORMAT_H */

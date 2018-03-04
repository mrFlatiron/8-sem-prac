#ifndef SPARSE_BASE_FORMAT_H
#define SPARSE_BASE_FORMAT_H

#include "common/vectors_fwd.h"

typedef struct
{
  int N;
  vector_int_t nnz_in_rows; /*N size*/
  vector_double_t values;   /*nnz size*/
  vector_int_t column_indecies; /*nnz size*/
} sparse_base_format;

int sparse_base_init (sparse_base_format *matr,
                      int N,
                      int max_row_nz);

void sparse_base_add_row (sparse_base_format *matr,
                          int row,
                          int cols[],
                          double values[],
                          int nnz);

void sparse_base_destroy (sparse_base_format *matr);

#endif /* SPARSE_BASE_FORMAT_H */

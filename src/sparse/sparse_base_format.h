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

#endif /* SPARSE_BASE_FORMAT_H */

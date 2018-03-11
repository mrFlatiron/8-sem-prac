#ifndef LASPACK_MATRIX_H
#define LASPACK_MATRIX_H

#include "3rd_party/laspack/qmatrix.h"
#include "sparse_base_format.h"
#include <stdio.h>

typedef struct
{
  QMatrix_L raw;
} laspack_matrix;

int laspack_matrix_init (laspack_matrix *mat, const sparse_base_format *base);
void laspack_matrix_destroy (laspack_matrix *mat);

void laspack_matrix_dump (laspack_matrix *mat, FILE *fout);

#endif /* LASPACK_MATRIX_H */

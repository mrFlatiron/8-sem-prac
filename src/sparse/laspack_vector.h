#ifndef LASPACK_VECTOR_H
#define LASPACK_VECTOR_H

#include "3rd_party/laspack/vector.h"
#include <stdio.h>
#include "common/vectors_fwd.h"

typedef struct
{
  Vector raw;
} laspack_vector;

int laspack_vector_init (laspack_vector *vec, int size);
void laspack_vector_destroy (laspack_vector *vec);

void laspack_vector_fill (laspack_vector *vec, const vector_double_t from);

void laspack_vector_dump (laspack_vector *vec, FILE *fout);

#endif /* LASPACK_VECTOR_H */

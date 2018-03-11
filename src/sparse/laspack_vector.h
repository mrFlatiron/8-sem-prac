#ifndef LASPACK_VECTOR_H
#define LASPACK_VECTOR_H

#include "3rd_party/laspack/vector.h"
#include "common/vectors_fwd.h"

typedef struct
{
  Vector raw;
} laspack_vector;

int laspack_vector_init (laspack_vector *vec, int size);
void laspack_vector_destroy (laspack_vector *vec);

void laspack_vector_fill (laspack_vector *vec, const vector_double_t from);

#endif /* LASPACK_VECTOR_H */

#include "laspack_vector.h"


int laspack_vector_init (laspack_vector *vec, int size)
{
  V_Constr (&vec->raw, "Vector", size, Normal, True);
  V_SetAllCmp (&vec->raw, 0.0);
  return 1;
}

void laspack_vector_fill (laspack_vector *vec, const vector_double_t from)
{
  size_t i;
  for (i = 0; i < vec->raw.Dim; i++)
    V_SetCmp (&vec->raw, i + 1, from[i]);
}

void laspack_vector_destroy (laspack_vector *vec)
{
  V_Destr (&vec->raw);
}

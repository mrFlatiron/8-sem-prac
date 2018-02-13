#include "vectors.h"

#include <stdlib.h>
#include <string.h>

vector_int_t vector_int_create (int size)
{
  vector_int_t retval;
  retval = (vector_int_t)malloc (size * sizeof (int));
  memset (retval, 0, size * sizeof (int));
  return retval;
}

void vector_int_copy (const vector_int_t from, vector_int_t to, int size)
{
  memcpy (to, from, size * sizeof (int));
}

vector_double_t vector_double_create (int size)
{
  vector_double_t retval;
  retval = (vector_double_t)malloc (size * sizeof (double));
  memset (retval, 0, size * sizeof (double));
  return retval;
}

void vector_double_copy (const vector_double_t from, vector_double_t to, int size)
{
  memcpy (to, from, size * sizeof (double));
}

void vector_move (void **from, void **to)
{
  *to = *from;
  *from = NULL;
}

void vector_destroy (void **vector)
{
  if (!*vector)
    return;

  free (*vector);

  *vector = NULL;
}

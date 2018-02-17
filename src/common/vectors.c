#include "vectors.h"

#include <stdlib.h>
#include <string.h>

vector_int_t vector_int_create (int size)
{
  vector_int_t retval;
  retval = (vector_int_t)malloc (size * sizeof (int));
  vector_int_set (retval, 0, size);
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
  vector_double_set (retval, 0, size);
  return retval;
}

void vector_double_copy (const vector_double_t from, vector_double_t to, int size)
{
  memcpy (to, from, size * sizeof (double));
}

void vector_destroy (void **vector)
{
  if (!*vector)
    return;

  free (*vector);

  *vector = NULL;
}

void vector_int_set (vector_int_t vector, int value, int size)
{
  memset (vector, value, size * sizeof (int));
}

void vector_double_set (vector_double_t vector, double value, int size)
{
  int i;

  for (i = 0; i < size; i++)
    vector[i] = value;
}

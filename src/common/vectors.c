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

string_t vector_char_create (int size)
{
  string_t retval;
  retval = (string_t)malloc ((size + 1) * sizeof (char));
  vector_char_set (retval, 0, size);
  return retval;
}

void vector_char_copy (const string_t from, string_t to, int size)
{
  memcpy (to, from, (size + 1) * sizeof (char));
}

void vector_char_set (string_t vector, char value, int size)
{
  memset (vector, value, size * sizeof (char));
  vector[size] = 0;
}

string_t *vector_string_t_create (int size)
{
  string_t *retval;
  retval = (string_t *)malloc (size * sizeof (string_t));
  vector_string_t_set (retval, NULL, size);
  return retval;
}

void vector_string_t_copy (const string_t *from, string_t *to, int size)
{
  int i;
  for (i = 0; i < size; i++)
    strcpy (to[i], from[i]);

}

void vector_string_t_set (string_t *vector, string_t value, int size)
{
  int i;
  for (i = 0; i < size; i++)
    {
      if (value == NULL)
        {
          vector[i] = NULL;
          continue;
        }

      strcpy (vector[i], value);
    }
}

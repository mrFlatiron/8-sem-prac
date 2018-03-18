#ifndef VECTORS_H
#define VECTORS_H

/*
 * Supported types are
 *   int
 *   double
 */


#include "vectors_fwd.h"

vector_int_t vector_int_create (int size);
void vector_int_copy   (const vector_int_t from, vector_int_t to, int size);
void vector_int_set (vector_int_t vector, int value, int size);



vector_double_t vector_double_create (int size);
void vector_double_copy   (const vector_double_t from, vector_double_t to, int size);
void vector_double_set (vector_double_t vector, double value, int size);

string_t vector_char_create (int size);
void vector_char_copy   (const string_t from, string_t to, int size);
void vector_char_set (string_t vector, char value, int size);

string_t *vector_string_t_create (int size);
void vector_string_t_copy (const string_t *from, string_t *to, int size);
void vector_string_t_set (string_t *vector, string_t value, int size);

void vector_destroy (void **vector);

#define PP_CAST (void**)

#define VECTOR_CREATE(one_word_type, size) vector_ ## one_word_type ## _create (size)
#define VECTOR_SET(one_word_type, vector, value, size) vector_ ## one_word_type ## _set (vector, value, size)
#define VECTOR_COPY(one_word_type, from, to, size) vector_ ## one_word_type ## _copy (from, to, size)
#define VECTOR_DESTROY(vector) vector_destroy (PP_CAST &vector)

#endif /* VECTORS_H */

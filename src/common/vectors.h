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



vector_double_t vector_double_create (int size);
void vector_double_copy   (const vector_double_t from, vector_double_t to, int size);

void vector_move (void **from, void **to);
void vector_destroy (void **vector);

#define PP_CAST (void**)

#define VECTOR_CREATE(one_word_type, size) vector_ ## one_word_type ## _create (size)
#define VECTOR_COPY(one_word_type, from, to, size) vector_ ## one_word_type ## _copy (from, to, size)
#define VECTOR_MOVE(from, to) vector_move (PP_CAST &from, PP_CAST &to)
#define VECTOR_DESTROY(vector) vector_destroy (PP_CAST &vector)

#endif /* VECTORS_H */

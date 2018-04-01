#ifndef HALF_NODES_VALUES_H
#define HALF_NODES_VALUES_H

#include "common/vectors_fwd.h"
#include "kernel/kernel_typedefs.h"

typedef struct
{
  int MX;
  int MY;
  int N;

  int vector_size;
  int layer_size;

  /*vector_size*/
  vector_double_t vals;
} half_nodes_values;

typedef struct
{
  int mx;
  int my;
  const half_nodes_values *vs;
  grid_area_t area;
  int is_end;
} hn_border_iter;

int hn_values_init (half_nodes_values *vs,
                int MX,
                int MY,
                int N);

void hn_values_destroy (half_nodes_values *vs);

int hn_values_index (const half_nodes_values *vs, int n, int mx, int my);
void hn_values_get_mx_my (const half_nodes_values *vs, int loc_layer_index, int *mx_ptr, int *my_ptr);

double hn_values_val_by_index (const half_nodes_values *vs, int n, int lli);
double hn_values_mx_my_val (const half_nodes_values *vs, int n, int mx, int my);

int hn_values_is_border (const half_nodes_values *vs, int lli);

void hn_border_iter_init (hn_border_iter *iter, const half_nodes_values *vs);
void hn_border_iter_next (hn_border_iter *iter);

#endif /* HALF_NODES_VALUES_H */

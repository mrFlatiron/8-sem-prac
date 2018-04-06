#ifndef NODES_VALUES_H
#define NODES_VALUES_H

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
} nodes_values;

typedef struct
{
  int mx;
  int my;
  const nodes_values *vs;
  int is_end;
  grid_area_t area;

} nodes_border_iter;

int nodes_values_init (nodes_values *vs,
                int MX,
                int MY,
                int N);

void nodes_values_destroy (nodes_values *vs);

int nodes_values_index (const nodes_values *vs, int n, int mx, int my);
void nodes_values_get_mx_my(const nodes_values *vs, int loc_layer_index, int *mx_ptr, int *my_ptr);
double nodes_values_mx_my_val (const nodes_values *vs, int n, int mx, int my);
double nodes_values_val (const nodes_values *vs, int n, int lli);

grid_area_t nodes_values_get_area (const nodes_values *vs, int mx, int my);

void nodes_border_iter_init (nodes_border_iter *iter, const nodes_values *vs);
void nodes_border_iter_next (nodes_border_iter *iter);

double nodes_avg_fwd_x (const nodes_values *vs, int n, int mx, int my);
double nodes_avg_fwd_y (const nodes_values *vs, int n, int mx, int my);


#endif /* NODES_VALUES_H */

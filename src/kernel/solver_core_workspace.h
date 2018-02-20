#ifndef SOLVER_CORE_WORKSPACE_H
#define SOLVER_CORE_WORKSPACE_H

#include "common/vectors_fwd.h"
#include "kernel_typedefs.h"
#include "sparse/msr_matrix.h"

#define SQUARES_COUNT           4
#define DIMENSIONS              2
#define DIMENSIONS_W_TIME       DIMENSIONS + 1
#define UNKNOWN_FUNCTIONS_COUNT 3

typedef struct
{
  int M1;
  int M2;
  int N;
  double X;
  double Y;
  double T;
  double h1;
  double h2;
  double tau;

  int vectors_size;
  int layer_size;

  vector_double_t vx;
  vector_double_t vy;
  vector_double_t g;

  int matrix_size;


  vector_double_t vector_to_compute;
  vector_double_t rhs_vector;
  msr_matrix matrix;

} solver_core_workspace;


int solver_workspace_data_init (solver_core_workspace *solver,
                          int M1,
                          int M2,
                          int N,
                          double X,
                          double Y,
                          double T);

void solver_workspace_data_destroy (solver_core_workspace *solver);


double solver_workspace_grid_vx (const solver_core_workspace *data, int n, int mx, int my);
double solver_workspace_grid_vy (const solver_core_workspace *data, int n, int mx, int my);
double solver_workspace_grid_g (const solver_core_workspace *data, int n, int mx, int my);
void   solver_workspace_fill_layer (solver_core_workspace *data, int n);

#endif /* SOLVER_CORE_WORKSPACE_H */

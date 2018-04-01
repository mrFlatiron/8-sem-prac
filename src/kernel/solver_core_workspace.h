#ifndef SOLVER_CORE_WORKSPACE_H
#define SOLVER_CORE_WORKSPACE_H

#include "common/vectors_fwd.h"
#include "kernel_typedefs.h"
#include "sparse/msr_matrix.h"
#include "sparse/sparse_base_format.h"
#include "cgs_solver.h"
#include "sparse/laspack_matrix.h"
#include "sparse/laspack_vector.h"






typedef struct
{
  solver_mode_t mode;
  linear_solver_t linear_solver;

  int MX;
  int MY;
  int N;
  double X;
  double Y;
  double T;
  double hx;
  double hy;
  double tau;

  double border_omega;

  int vectors_size;
  int layer_size;

  vector_double_t vx;
  vector_double_t vy;
  vector_double_t g;

  int matrix_size;


  /*matrix_size*/
  vector_double_t vector_to_compute;
  vector_double_t rhs_vector;
  laspack_vector  rhs_vector_l;
  laspack_vector  vector_to_compute_l;

  sparse_base_format matrix_base;
  msr_matrix matrix;
  laspack_matrix matrix_l;

  cgs_solver cgs_linear_solver;

  FILE *log_file;

} solver_core_workspace;



typedef enum
{
  grid_g,
  grid_vx,
  grid_vy
} grid_func_t;


int solver_workspace_data_init (solver_core_workspace *solver,
                                solver_mode_t mode,
                                int M1,
                                int M2,
                                int N,
                                double X,
                                double Y,
                                double T,
                                double border_omega,
                                double solver_prec,
                                int solver_max_iter,
                                preconditioner_t precond,
                                linear_solver_t linear_solver);

void solver_workspace_data_destroy (solver_core_workspace *solver);


double        solver_workspace_grid_vx             (const solver_core_workspace *data, int n, int mx, int my);
double        solver_workspace_grid_vy             (const solver_core_workspace *data, int n, int mx, int my);
double        solver_workspace_grid_g              (const solver_core_workspace *data, int n, int mx, int my);
double        solver_workspace_grid_val            (const solver_core_workspace *data, int n, int mx, int my, grid_func_t func);
int           solver_workspace_func_col             (const solver_core_workspace *data, int loc_layer_index, grid_func_t func);

int           solver_workspace_final_index         (const solver_core_workspace *data, int n, int mx, int my);

int           solver_workspace_layer_begin_index            (const solver_core_workspace *data, int n);

void          solver_workspace_fill_layer          (solver_core_workspace *data, int n);
void          solver_workspace_get_mx_my           (const solver_core_workspace *data, int loc_layer_index, int *mx_ptr, int *my_ptr);
grid_area_t   solver_workspace_get_area            (const solver_core_workspace *data, int mx, int my);


#endif /* SOLVER_CORE_WORKSPACE_H */

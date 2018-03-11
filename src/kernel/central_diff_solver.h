#ifndef CENTRAL_DIFF_SOLVER_H
#define CENTRAL_DIFF_SOLVER_H

#include "solver_core_workspace.h"

typedef struct
{
  solver_core_workspace ws;

  double mu;

  int layer;

  double (*p_drv) (double);
  pressure_func_t p_drv_type;

  rhs_func_t f0;
  rhs_func_t f1;
  rhs_func_t f2;

  layer_func_t start_vx;
  layer_func_t start_vy;
  layer_func_t start_g;

  time_layer_func_t test_solution_g;
  time_layer_func_t test_solution_vx;
  time_layer_func_t test_solution_vy;

} central_diff_solver;

int cdiff_solver_init (central_diff_solver *solver,
                       solver_mode_t mode,
                       int M1,
                       int M2,
                       int N,
                       double X,
                       double Y,
                       double T,
                       double border_omega,
                       linear_solver_t linear_solver,
                       time_layer_func_t test_solution_g,
                       time_layer_func_t test_solution_vx,
                       time_layer_func_t test_solution_vy);

void cdiff_solver_destroy (central_diff_solver *solver);


int cdiff_solver_compute (central_diff_solver *solver, pressure_func_t p_func,
                          double mu,
                          rhs_func_t f0,
                          rhs_func_t f1,
                          rhs_func_t f2,
                          layer_func_t start_vx,
                          layer_func_t start_vy,
                          layer_func_t start_g);

#endif /* CENTRAL_DIFF_SOLVER_H */

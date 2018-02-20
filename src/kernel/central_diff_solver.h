#ifndef CENTRAL_DIFF_SOLVER_H
#define CENTRAL_DIFF_SOLVER_H

#include "solver_core_workspace.h"

typedef struct
{
  solver_core_workspace ws;

  double mu;

  double (*p_drv) (double);

  area_func_t f0;
  area_func_t f1;

  layer_func_t start_vx;
  layer_func_t start_vy;
  layer_func_t start_g;

} central_diff_solver;

int cdiff_solver_init (central_diff_solver *solver,
                       int M1,
                       int M2,
                       int N,
                       double X,
                       double Y,
                       double T);

void cdiff_solver_destroy (central_diff_solver *solver);


int cdiff_solver_compute (central_diff_solver *solver, pressure_func_t p_func,
                          double mu,
                          area_func_t f0,
                          area_func_t f1,
                          layer_func_t start_vx,
                          layer_func_t start_vy,
                          layer_func_t start_g);

#endif /* CENTRAL_DIFF_SOLVER_H */

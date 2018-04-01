#ifndef SOKOLOV_SOLVER_H
#define SOKOLOV_SOLVER_H

#include "kernel_typedefs.h"
#include "kernel/mesh_info.h"
#include "half_nodes_values.h"
#include "nodes_values.h"
#include "linear_system_composer.h"

typedef struct
{
  solver_mode_t mode;

  mesh_info_t mesh_info;

  half_nodes_values h_obj;
  half_nodes_values *h;

  nodes_values vx_obj;
  nodes_values *vx;

  nodes_values vy_obj;
  nodes_values *vy;

  linear_system_composer composer_h_obj;
  linear_system_composer *composer_h;

  linear_system_composer composer_v_obj;
  linear_system_composer *composer_v;

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

} sokolov_solver;

int sokolov_solver_init (sokolov_solver *solver,
                         mesh_info_t mesh_info,
                         solver_mode_t mode,
                         pressure_func_t p_func,
                         double solver_prec,
                         int solver_max_iter,
                         preconditioner_t precond,
                         linear_solver_t linear_solver);

void sokolov_solver_destroy (sokolov_solver *solver);

int sokolov_solver_compute (sokolov_solver *solver,
                           time_layer_func_t test_solution_g,
                           time_layer_func_t test_solution_vx,
                           time_layer_func_t test_solution_vy,
                           rhs_func_t f0,
                           rhs_func_t f1,
                           rhs_func_t f2,
                           layer_func_t start_vx,
                           layer_func_t start_vy,
                           layer_func_t start_g);

#endif /* SOKOLOV_SOLVER_H */

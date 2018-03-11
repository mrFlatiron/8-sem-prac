#ifndef SOLVER_TESTER_H
#define SOLVER_TESTER_H

#include "solver_core_workspace.h"

typedef struct
{
  const solver_core_workspace *ws;

  double g_l2_norm;
  double vx_l2_norm;
  double vy_l2_norm;
  double g_c_norm;
  double vx_c_norm;
  double vy_c_norm;
  double g_w21_norm;
  double vx_w21_norm;
  double vy_w21_norm;

  time_layer_func_t g_func;
  time_layer_func_t vx_func;
  time_layer_func_t vy_func;
} solver_tester;

void solver_tester_init (solver_tester *tester, const solver_core_workspace *ws,
                         time_layer_func_t g_func,
                         time_layer_func_t vx_func,
                         time_layer_func_t vy_func);

double grid_l2_norm (const solver_tester *tester, vector_double_t vec);
double grid_c_norm (const solver_tester *tester, vector_double_t vec);
double grid_w21_norm (const solver_tester *tester, vector_double_t vec);


#endif /* SOLVER_TESTER_H */

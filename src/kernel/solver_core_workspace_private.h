#ifndef SOLVER_CORE_WORKSPACE_PRIVATE_H
#define SOLVER_CORE_WORKSPACE_PRIVATE_H

#include "solver_core_workspace.h"

int solver_workspace_data_check_input (int M1,
                           int M2,
                           int N,
                           double X1,
                           double X2,
                           double T);


int solver_workspace_top_square_begin_index (const solver_core_workspace *data, int n);

void solver_workspace_fill_zero_layer (solver_core_workspace *data);

void solver_workspace_check_layer (const solver_core_workspace *data, int n);

int  solver_workspace_check_border_value  (const solver_core_workspace *data, grid_area_t area, grid_func_t func, double actual, double *must_be);


#endif /* SOLVER_CORE_WORKSPACE_PRIVATE_H */

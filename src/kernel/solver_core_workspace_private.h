#ifndef SOLVER_CORE_WORKSPACE_PRIVATE_H
#define SOLVER_CORE_WORKSPACE_PRIVATE_H

#include "solver_core_workspace.h"

int solver_workspace_data_check_input (int M1,
                           int M2,
                           int N,
                           double X1,
                           double X2,
                           double T);

int solver_workspace_layer_begin_index (const solver_core_workspace *data, int n);

int solver_workspace_final_index (const solver_core_workspace *data, int n, int mx, int my);

#endif /* SOLVER_CORE_WORKSPACE_PRIVATE_H */

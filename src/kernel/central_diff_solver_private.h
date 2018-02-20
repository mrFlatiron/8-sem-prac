#ifndef CENTRAL_DIFF_SOLVER_PRIVATE_H
#define CENTRAL_DIFF_SOLVER_PRIVATE_H

#include "central_diff_solver.h"

void cdiff_solver_init_first_layer (central_diff_solver *solver);

void cdiff_solver_fill_matrix (central_diff_solver *solver);

void cdiff_solver_fill_rhs (central_diff_solver *solver);

void cdiff_solver_solve_system (central_diff_solver *solver);

#endif /* CENTRAL_DIFF_SOLVER_PRIVATE_H */

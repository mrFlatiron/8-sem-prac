#ifndef SOKOLOV_SOLVER_PRIVATE_H
#define SOKOLOV_SOLVER_PRIVATE_H

#include "sokolov_solver.h"

void sokolov_solver_fill_zero_layer (sokolov_solver *solver);
void sokolov_solver_fill_borders (sokolov_solver *solver);

void sokolov_solver_fill_h_matrix_w_rhs (sokolov_solver *solver);
void sokolov_solver_solve_h_system (sokolov_solver *solver);
void sokolov_solver_fill_h_layer (sokolov_solver *solver);

void sokolov_solver_fill_v_matrix_w_rhs (sokolov_solver *solver);
void sokolov_solver_solve_v_system (sokolov_solver *solver);
void sokolov_solver_fill_v_layer (sokolov_solver *solver);

void sokolov_solver_fill_v_values_from_functions (sokolov_solver *solver, int layer, int mx, int my);
void sokolov_solver_fill_h_values_from_functions (sokolov_solver *solver, int layer, int mx, int my);
void sokolov_solver_fill_borders_on_layer (sokolov_solver *solver, int layer);

void sokolov_solver_fill_x_init_w_real_values (sokolov_solver *solver);

void sokolov_solver_set_hn_x_y (const sokolov_solver *solver, int mx, int my, double *x_ptr, double *y_ptr);

#endif /* SOKOLOV_SOLVER_PRIVATE_H */

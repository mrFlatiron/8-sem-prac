#ifndef CENTRAL_DIFF_SOLVER_PRIVATE_H
#define CENTRAL_DIFF_SOLVER_PRIVATE_H

#include "central_diff_solver.h"

typedef struct
{
  solver_core_workspace *ws;
  central_diff_solver *svr;
  double (*g_val)  (const solver_core_workspace *, int, int, int);
  double (*vx_val) (const solver_core_workspace *, int, int, int);
  double (*vy_val) (const solver_core_workspace *, int, int, int);
  double tau;
  double hx;
  double hy;

  double mu_wave;

  int loc_layer_index;
  int mx;
  int my;
  int n;
  int row;

  double nz_values[MAX_ROW_NZ];
  int    nz_cols[MAX_ROW_NZ];
  int    nnz;

  int g_curr;
  int g_left;
  int g_right;
  int g_top;
  int g_bot;

  int vx_curr;
  int vx_left;
  int vx_right;
  int vx_top;
  int vx_bot;

  int vy_curr;
  int vy_left;
  int vy_right;
  int vy_top;
  int vy_bot;

} eq_filler_t;

void eq_filler_init (eq_filler_t *ef, central_diff_solver *solver);

void eq_filler_inc_layer_index (eq_filler_t *ef);

void cdiff_solver_init_first_layer (central_diff_solver *solver);

void cdiff_solver_init_borders (central_diff_solver *solver);

void cdiff_solver_init_layer_borders (central_diff_solver *solver, int layer);

void cdiff_solver_fill_values_from_functions (central_diff_solver *solver, int layer, int mx, int my);

void cdiff_solver_fill_matrix_w_rhs (central_diff_solver *solver);

void cdiff_solver_solve_system (central_diff_solver *solver);

double cdiff_solver_mu_wave (const central_diff_solver *solver);

void cdiff_solver_eq_5_1 (eq_filler_t *ef);
void cdiff_solver_eq_5_2 (eq_filler_t *ef);
void cdiff_solver_eq_5_3 (eq_filler_t *ef);
void cdiff_solver_eq_5_4 (eq_filler_t *ef);
void cdiff_solver_eq_5_5 (eq_filler_t *ef);
void cdiff_solver_eq_5_6 (eq_filler_t *ef);
void cdiff_solver_eq_5_7 (eq_filler_t *ef);
void cdiff_solver_eq_5_8 (eq_filler_t *ef);
void cdiff_solver_eq_trivial (eq_filler_t *ef, grid_func_t func);
void cdiff_solver_eq_border (eq_filler_t *ef, grid_area_t border);

#endif /* CENTRAL_DIFF_SOLVER_PRIVATE_H */

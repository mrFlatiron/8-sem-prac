#ifndef SOLVER_TESTER_H
#define SOLVER_TESTER_H

#include "central_diff_solver.h"

#define MAX_DOUBLE_LENGTH 50

typedef double (*grid_ws_func) (const solver_core_workspace *, int/*n*/, int/*mx*/, int/*my*/);

typedef struct
{
  int N_start;
  int MX_start;
  int MY_start;
  int N_mult;
  int MXY_mult;

  int N_mult_count;
  int MXY_mult_count;

  int cases_count;
  int vector_size;

  vector_double_t c_norms;
  string_t *cg_norms_text;
  string_t *cvx_norms_text;
  string_t *cvy_norms_text;
  vector_double_t l2_norms;
  string_t *l2g_norms_text;
  string_t *l2vx_norms_text;
  string_t *l2vy_norms_text;
  vector_double_t w21_norms;
  string_t *w21g_norms_text;
  string_t *w21vx_norms_text;
  string_t *w21vy_norms_text;

  central_diff_solver cds;

  double border_omega;
  double mu;
  double X;
  double Y;
  double T;

  time_layer_func_t g_func;
  time_layer_func_t vx_func;
  time_layer_func_t vy_func;

  rhs_func_t f0;
  rhs_func_t f1;
  rhs_func_t f2;

  layer_func_t start_vx;
  layer_func_t start_vy;
  layer_func_t start_g;

} solver_tester;

void solver_tester_init (solver_tester *tester,
                         int N_start,
                         int MX_start,
                         int MY_start,
                         int N_mult,
                         int MXY_mult,
                         int N_mult_count,
                         int MXY_mult_count,
                         double border_omega,
                         double mu,
                         double X,
                         double Y,
                         double T,
                         time_layer_func_t g_func,
                         time_layer_func_t vx_func,
                         time_layer_func_t vy_func,
                         rhs_func_t f0,
                         rhs_func_t f1,
                         rhs_func_t f2,
                         layer_func_t start_vx,
                         layer_func_t start_vy,
                         layer_func_t start_g);

void solver_tester_test (solver_tester *tester,
                         double solver_prec,
                         int solver_max_iter,
                         preconditioner_t precond);

void solver_tester_destroy (solver_tester *tester);

double tester_grid_dif_c_norm (const solver_tester *tester, grid_func_t f);
double tester_grid_dif_l2_norm (const solver_tester *tester, grid_func_t f);
double tester_grid_dif_w21_norm (const solver_tester *tester, grid_func_t f);

double tester_grid_true_val (const solver_tester *tester, grid_func_t f, int loc_layer_index);


#endif /* SOLVER_TESTER_H */

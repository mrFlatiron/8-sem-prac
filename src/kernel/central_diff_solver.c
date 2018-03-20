#include "central_diff_solver_private.h"
#include "input/rhs.h"
#include "cgs_solver.h"
#include "common/debug_utils.h"
#include "sparse/sparse_base_format.h"
#include "common/math_utils.h"
#include "3rd_party/laspack/itersolv.h"

/*TEMPCODE_BEGIN*/
#include "common/vectors.h"
#include "linear_ops/vector_ops.h"
/*TEMPCODE_END*/
/*
 * TODO
 *
 *  (DONE) fix max_iter_exceeded :(
 *  (DONE) use f0
 *  (DONE) fill_matrix
 *  (DONE) solve_matrix
 *  (DONE) compute difnorms
 *
 *  (DONE) solve system without new allocations
 *  (DONE) use the previous solution as x_init
 */

void fill_nz (double *vals, int *cols, int *nnz, double val, int col)
{
/*  DEBUG_ASSERT (!math_is_null (val)); */
  if (math_is_null (val))
    return;

  /*DEBUG_ASSERT (val < 1e6);*/

  vals[*nnz] = val;
  cols[*nnz] = col;
  (*nnz)++;
}

int cdiff_solver_init (central_diff_solver *solver,
                       solver_mode_t mode,
                       int MX,
                       int MY,
                       int N,
                       double X,
                       double Y,
                       double T,
                       double border_omega,
                       double mu,
                       double solver_prec,
                       int solver_max_iter,
                       preconditioner_t precond,
                       linear_solver_t linear_solver)
{
  solver_workspace_data_init (&solver->ws, mode, MX, MY, N, X, Y, T, border_omega,
                              solver_prec, solver_max_iter,
                              precond, linear_solver);
  solver->mu = mu;

  return 0;
}



void cdiff_solver_destroy (central_diff_solver *solver)
{
  solver_workspace_data_destroy (&solver->ws);
}



int cdiff_solver_compute (central_diff_solver *solver,
                          pressure_func_t p_func,
                          time_layer_func_t test_solution_g,
                          time_layer_func_t test_solution_vx,
                          time_layer_func_t test_solution_vy,
                          rhs_func_t f0,
                          rhs_func_t f1,
                          rhs_func_t f2,
                          layer_func_t start_vx,
                          layer_func_t start_vy,
                          layer_func_t start_g)
{
  if (solver->ws.mode == test_mode)
    {
      DEBUG_ASSERT (test_solution_g && test_solution_vx && test_solution_vy);
      solver->test_solution_g = test_solution_g;
      solver->test_solution_vx = test_solution_vx;
      solver->test_solution_vy = test_solution_vy;
    }
  else
    {
      solver->test_solution_g = NULL;
      solver->test_solution_vx = NULL;
      solver->test_solution_vy = NULL;
    }

  solver->f0 = f0;
  solver->f1 = f1;
  solver->f2 = f2;
  solver->start_vx = start_vx;
  solver->start_vy = start_vy;
  solver->start_g  = start_g;

  solver->layer = 0;

  solver->p_drv_type = p_func;

  switch (p_func)
    {
    case pressure_linear:
      solver->p_drv = &p_drv_linear;
      break;
    case pressure_polynomial:
      solver->p_drv = &p_drv_polynomial;
      break;
    }



  cdiff_solver_init_first_layer (solver);
  cdiff_solver_init_borders (solver);


  for (solver->layer = 1; solver->layer <= solver->ws.N; solver->layer++)
    {
      cdiff_solver_fill_matrix_w_rhs (solver);
      cdiff_solver_solve_system (solver);
      solver_workspace_fill_layer (&solver->ws, solver->layer);
    }

  return 0;
}

void cdiff_solver_init_first_layer (central_diff_solver *solver)
{
  int mx;
  int my;
  double x;
  double y;
  double border_omega = solver->ws.border_omega;
  int i = 0;

  for (my = 0; my <= solver->ws.MY; my++)
    {
      for (mx = 0; mx <= solver->ws.MX * BOT_ROW_SQUARES_COUNT; mx++)
        {
          x = mx * solver->ws.hx;
          y = my * solver->ws.hy;
          solver->ws.vector_to_compute[UNKNOWN_FUNCTIONS_COUNT * i]     = solver->start_g  (x, y, border_omega);
          solver->ws.vector_to_compute[UNKNOWN_FUNCTIONS_COUNT * i + 1] = solver->start_vx (x, y, border_omega);
          solver->ws.vector_to_compute[UNKNOWN_FUNCTIONS_COUNT * i + 2] = solver->start_vy (x, y, border_omega);
          i++;
        }
    }

  for (my = solver->ws.MY + 1; my <= 2 * solver->ws.MY; my++)
    {
      for (mx = solver->ws.MX; mx <= 2 * solver->ws.MX; mx++)
        {
          x = mx * solver->ws.hx;
          y = my * solver->ws.hy;
          solver->ws.vector_to_compute[UNKNOWN_FUNCTIONS_COUNT * i]     = solver->start_g  (x, y, border_omega);
          solver->ws.vector_to_compute[UNKNOWN_FUNCTIONS_COUNT * i + 1] = solver->start_vx (x, y, border_omega);
          solver->ws.vector_to_compute[UNKNOWN_FUNCTIONS_COUNT * i + 2] = solver->start_vy (x, y, border_omega);
          i++;
        }
    }

  solver_workspace_fill_layer (&solver->ws, 0);
}

void cdiff_solver_fill_values_from_functions (central_diff_solver *solver, int layer, int mx, int my)
{
  double x = solver->ws.hx * mx;
  double y = solver->ws.hy * my;
  double t = solver->ws.tau * layer;
  int index = solver_workspace_final_index (&solver->ws, layer, mx, my);
  solver->ws.vx[index] = solver->test_solution_vx (t, x, y);
  solver->ws.vy[index] = solver->test_solution_vy (t, x, y);
  solver->ws.g[index]  = solver->test_solution_g (t, x, y);
}
void cdiff_solver_init_layer_borders (central_diff_solver *solver, int layer)
{
  if (solver->ws.mode == solve_mode)
    {
      int mx;
      int my;
      int index;

      /*border_leftmost*/
      mx = 0;
      for (my = 0; my <= solver->ws.MY; my++)
        {
          index = solver_workspace_final_index (&solver->ws, layer, mx, my);
          solver->ws.vx[index] = solver->ws.border_omega;
          solver->ws.g[index] = log (RHO_LEFTMOST);
        }

      /*
       * Other borders are
       * already zero
       */
    }
  else
    {
      int mx;
      int my;

      /*leftmost*/
      mx = 0;
      for (my = 0; my <= solver->ws.MY; my++)
        cdiff_solver_fill_values_from_functions (solver, layer, mx, my);

      /*botmost*/
      my = 0;
      for (mx = 0; mx <= BOT_ROW_SQUARES_COUNT *  solver->ws.MX; mx++)
        cdiff_solver_fill_values_from_functions (solver, layer, mx, my);

      /*rightmost*/
      mx = BOT_ROW_SQUARES_COUNT * solver->ws.MX;
      for (my = 0; my <= solver->ws.MY; my++)
        cdiff_solver_fill_values_from_functions (solver, layer, mx, my);

      /*topmost*/
      my = 2 * solver->ws.MY;
      for (mx = solver->ws.MX; mx <= 2 * solver->ws.MX; mx++)
        cdiff_solver_fill_values_from_functions (solver, layer, mx, my);

      if (layer == solver->ws.N)
        DEBUG_ASSERT (solver_workspace_final_index (&solver->ws, layer, 2 * solver->ws.MX, my) ==
                      solver->ws.vectors_size - 1);

      /*left_hor*/
      my = solver->ws.MY;
      for (mx = 1; mx <= solver->ws.MX; mx++)
        cdiff_solver_fill_values_from_functions (solver, layer, mx, my);

      /*right_hor*/
      for (mx = 2 * solver->ws.MX; mx < BOT_ROW_SQUARES_COUNT * solver->ws.MX; mx++)
        cdiff_solver_fill_values_from_functions (solver, layer, mx, my);

      /*left_top*/
      mx = solver->ws.MX;
      for (my = solver->ws.MY + 1; my < 2 * solver->ws.MY; my++)
        cdiff_solver_fill_values_from_functions (solver, layer, mx, my);

      /*right_top*/
      mx = 2 * solver->ws.MX;
      for (my = solver->ws.MY + 1; my < 2 * solver->ws.MY; my++)
        cdiff_solver_fill_values_from_functions (solver, layer, mx, my);
    }
}

void cdiff_solver_init_borders (central_diff_solver *solver)
{
  int layer;
  for (layer = 1; layer <= solver->ws.N; layer++)
    cdiff_solver_init_layer_borders (solver, layer);

}

void eq_filler_init (eq_filler_t *ef, central_diff_solver *solver)
{
  ef->ws = &solver->ws;
  ef->svr = solver;

  ef->tau = solver->ws.tau;
  ef->hx = solver->ws.hx;
  ef->hy = solver->ws.hy;

  ef->mu_wave = cdiff_solver_mu_wave (solver);

  ef->g_val = solver_workspace_grid_g;
  ef->vx_val = solver_workspace_grid_vx;
  ef->vy_val = solver_workspace_grid_vy;

  ef->loc_layer_index = 0;
  ef->mx = 0;
  ef->my= 0;
  ef->n = solver->layer - 1;
  ef->row = 0;

  ef->nnz = 0;

  ef->g_curr = 0;
  ef->g_left = -UNKNOWN_FUNCTIONS_COUNT;
  ef->g_right = UNKNOWN_FUNCTIONS_COUNT;
  ef->g_top = UNKNOWN_FUNCTIONS_COUNT * (BOT_ROW_SQUARES_COUNT * ef->ws->MX + 1);
  ef->g_bot = -ef->g_top;

  ef->vx_left = -UNKNOWN_FUNCTIONS_COUNT + 1;
  ef->vx_curr = 1;
  ef->vx_right = 1 + UNKNOWN_FUNCTIONS_COUNT;
  ef->vx_top = ef->vx_curr + UNKNOWN_FUNCTIONS_COUNT * (BOT_ROW_SQUARES_COUNT * ef->ws->MX + 1);
  ef->vx_bot = ef->vx_curr - UNKNOWN_FUNCTIONS_COUNT * (BOT_ROW_SQUARES_COUNT * ef->ws->MX + 1);

  ef->vy_left = -UNKNOWN_FUNCTIONS_COUNT + 2;
  ef->vy_curr = 2;
  ef->vy_right = 2 + UNKNOWN_FUNCTIONS_COUNT;
  ef->vy_top = ef->vy_curr + UNKNOWN_FUNCTIONS_COUNT * (BOT_ROW_SQUARES_COUNT * ef->ws->MX + 1);
  ef->vy_bot = ef->vy_curr - UNKNOWN_FUNCTIONS_COUNT * (BOT_ROW_SQUARES_COUNT * ef->ws->MX + 1);
}

void eq_filler_inc_layer_index (eq_filler_t *ef)
{
  int loc_top = 0;
  int loc_bot = 0;
  int mx, my;

  ef->loc_layer_index++;

  solver_workspace_get_mx_my (ef->ws, ef->loc_layer_index, &mx, &my);

  if (my < ef->ws->MY
      || (
        mx >= ef->ws->MX
        && mx <= 2 * ef->ws->MX
        && my < 2 * ef->ws->MY))
    loc_top = solver_workspace_final_index (ef->ws, 0, mx, my + 1);

  if (my > 0)
    loc_bot = solver_workspace_final_index (ef->ws, 0, mx, my - 1);



  ef->g_curr += UNKNOWN_FUNCTIONS_COUNT;
  ef->g_left = ef->g_curr - UNKNOWN_FUNCTIONS_COUNT;
  ef->g_right = ef->g_curr + UNKNOWN_FUNCTIONS_COUNT;
  ef->g_bot = UNKNOWN_FUNCTIONS_COUNT * loc_bot;
  ef->g_top = UNKNOWN_FUNCTIONS_COUNT * loc_top;


  ef->vx_curr += UNKNOWN_FUNCTIONS_COUNT;
  ef->vx_left = ef->vx_curr - UNKNOWN_FUNCTIONS_COUNT;
  ef->vx_right = ef->vx_curr + UNKNOWN_FUNCTIONS_COUNT;
  ef->vx_bot = UNKNOWN_FUNCTIONS_COUNT * loc_bot + 1;
  ef->vx_top = UNKNOWN_FUNCTIONS_COUNT * loc_top + 1;

  ef->vy_curr += UNKNOWN_FUNCTIONS_COUNT;
  ef->vy_left = ef->vy_curr - UNKNOWN_FUNCTIONS_COUNT;
  ef->vy_right = ef->vy_curr + UNKNOWN_FUNCTIONS_COUNT;
  ef->vy_bot = UNKNOWN_FUNCTIONS_COUNT * loc_bot + 2;
  ef->vy_top = UNKNOWN_FUNCTIONS_COUNT * loc_top + 2;
}

void cdiff_solver_eq_5_2 (eq_filler_t *ef)
{
  int my = ef->my;
  int mx = ef->mx;

  double x = mx * ef->ws->hx;
  double y = my * ef->ws->hy;
  double t = ef->n * ef->ws->tau;
  int n = ef->n;

  ef->nnz = 0;

  fill_nz (ef->nz_values, ef->nz_cols, &ef->nnz,
           4,
           ef->g_curr);
  fill_nz (ef->nz_values, ef->nz_cols, &ef->nnz,
           - ef->tau / ef->hx * (
             + ef->vx_val (ef->ws, n, mx - 1, my)
             + ef->vx_val (ef->ws, n, mx, my)),
           ef->g_left);
  fill_nz (ef->nz_values, ef->nz_cols, &ef->nnz,
           - ef->tau / ef->hy * (
             + ef->vy_val (ef->ws, n, mx, my - 1)
             + ef->vy_val (ef->ws, n, mx, my)),
           ef->g_bot);
  fill_nz (ef->nz_values, ef->nz_cols, &ef->nnz,
           ef->tau / ef->hx * (
             + ef->vx_val (ef->ws, n, mx + 1, my)
             + ef->vx_val (ef->ws, n, mx, my)),
           ef->g_right);
  fill_nz (ef->nz_values, ef->nz_cols, &ef->nnz,
           ef->tau / ef->hy * (
             + ef->vy_val (ef->ws, n, mx, my + 1)
             + ef->vy_val (ef->ws, n, mx, my)),
           ef->g_top);
  fill_nz (ef->nz_values, ef->nz_cols, &ef->nnz,
           - 2 * ef->tau / ef->hx,
           ef->vx_left);
  fill_nz (ef->nz_values, ef->nz_cols, &ef->nnz,
           - 2 * ef->tau / ef->hy,
           ef->vy_bot);
  fill_nz (ef->nz_values, ef->nz_cols, &ef->nnz,
           2 * ef->tau / ef->hx,
           ef->vx_right);
  fill_nz (ef->nz_values, ef->nz_cols, &ef->nnz,
           2 * ef->tau / ef->hy,
           ef->vy_top);

  sparse_base_add_row (&ef->ws->matrix_base, ef->row, ef->nz_cols, ef->nz_values, ef->nnz);

  ef->ws->rhs_vector[ef->row] =
      + 4 * ef->g_val (ef->ws, n, mx, my)
      + ef->tau * ef->g_val (ef->ws, n, mx, my) * (
        + (ef->vx_val (ef->ws, n, mx + 1, my) - ef->vx_val (ef->ws, n, mx - 1, my)) / ef->hx
        + (ef->vy_val (ef->ws, n, mx, my + 1) - ef->vy_val (ef->ws, n, mx, my - 1))/ ef->hy)
      + 4 * ef->tau * ef->svr->f0 (t, x, y, ef->svr->mu, ef->svr->p_drv_type);

  ef->row++;
}

void cdiff_solver_eq_5_3 (eq_filler_t *ef)
{
  int my = ef->my;
  int mx = ef->mx;
  double x = mx * ef->ws->hx;
  double y = my * ef->ws->hy;
  double t = ef->n * ef->ws->tau;
  int n = ef->n;
  double full_rhs_val;

  DEBUG_ASSERT (mx == 0 || mx == ef->ws->MX);

  ef->nnz = 0;

  fill_nz (ef->nz_values, ef->nz_cols, &ef->nnz,
            2 - ef->tau / ef->hx * ef->vx_val (ef->ws, n, mx, my),
            ef->g_curr);

  fill_nz (ef->nz_values, ef->nz_cols, &ef->nnz,
            ef->tau / ef->hx * ef->vx_val (ef->ws, n, mx + 1, my),
            ef->g_right);

  fill_nz (ef->nz_values, ef->nz_cols, &ef->nnz,
           - 2 * ef->tau / ef->hx,
           ef->vx_curr);

  fill_nz (ef->nz_values, ef->nz_cols, &ef->nnz,
           2 * ef->tau / ef->hx,
           ef->vx_right);

  sparse_base_add_row (&ef->ws->matrix_base, ef->row, ef->nz_cols, ef->nz_values, ef->nnz);

  full_rhs_val =
      + 2 * ef->g_val (ef->ws, n, mx, my)
      + ef->tau / ef->hx * ef->g_val (ef->ws, n, mx, my) * (+ ef->vx_val (ef->ws, n, mx + 1, my)
                                                            - ef->vx_val (ef->ws, n, mx, my))
      + ef->tau / ef->hx * (
        + ef->g_val (ef->ws, n, mx, my) * ef->vx_val (ef->ws, n, mx, my)
        - 2 * ef->g_val (ef->ws, n, mx + 1, my) * ef->vx_val (ef->ws, n, mx + 1, my)
        + ef->g_val (ef->ws, n, mx + 2, my) * ef->vx_val (ef->ws, n, mx + 2, my)
        - 0.5 * (+ ef->g_val (ef->ws, n, mx + 1, my) * ef->vx_val (ef->ws, n, mx + 1, my)
                 - 2 * ef->g_val (ef->ws, n, mx + 2, my) * ef->vx_val (ef->ws, n, mx + 2, my)
                 + ef->g_val (ef->ws, n, mx + 3, my) * ef->vx_val (ef->ws, n, mx + 3, my))
        + (2 - ef->g_val (ef->ws, n, mx, my)) * (
          + ef->vx_val (ef->ws, n, mx, my)
          - 2 * ef->vx_val (ef->ws, n, mx + 1, my)
          + ef->vx_val (ef->ws, n, mx + 2, my)
          - 0.5 * (+ ef->vx_val (ef->ws, n, mx + 1, my)
                   - 2 * ef->vx_val (ef->ws, n, mx + 2, my)
                   + ef->vx_val (ef->ws, n, mx + 3, my))))
      + 2 * ef->tau * ef->svr->f0 (t, x, y, ef->svr->mu, ef->svr->p_drv_type);

  ef->ws->rhs_vector[ef->row] = full_rhs_val;

  ef->row++;
}

void cdiff_solver_eq_5_4 (eq_filler_t *ef)
{
  int mx = ef->mx;
  int my = ef->my;
  double x = mx * ef->ws->hx;
  double y = my * ef->ws->hy;
  double t = ef->n * ef->ws->tau;
  int n = ef->n;
  double full_rhs_val;

  DEBUG_ASSERT (mx == 2 * ef->ws->MX || mx == BOT_ROW_SQUARES_COUNT * ef->ws->MX);

  ef->nnz = 0;

  fill_nz (ef->nz_values, ef->nz_cols, &ef->nnz,
            2 + ef->tau / ef->hx * ef->vx_val (ef->ws, n, mx, my),
            ef->g_curr);

  fill_nz (ef->nz_values, ef->nz_cols, &ef->nnz,
            -ef->tau / ef->hx * ef->vx_val (ef->ws, n, mx - 1, my),
            ef->g_left);

  fill_nz (ef->nz_values, ef->nz_cols, &ef->nnz,
           - 2 * ef->tau / ef->hx,
           ef->vx_left);

  fill_nz (ef->nz_values, ef->nz_cols, &ef->nnz,
           + 2 * ef->tau / ef->hx,
           ef->vx_curr);

  sparse_base_add_row (&ef->ws->matrix_base, ef->row, ef->nz_cols, ef->nz_values, ef->nnz);

  full_rhs_val =
      + 2 * ef->g_val (ef->ws, n, mx, my)
      + ef->tau / ef->hx * ef->g_val (ef->ws, n, mx, my) * (+ ef->vx_val (ef->ws, n, mx, my)
                                                            - ef->vx_val (ef->ws, n, mx - 1, my))
      - ef->tau / ef->hx * (
        + ef->g_val (ef->ws, n, mx - 2, my) * ef->vx_val (ef->ws, n, mx - 2, my)
        - 2 * ef->g_val (ef->ws, n, mx - 1, my) * ef->vx_val (ef->ws, n, mx - 1, my)
        + ef->g_val (ef->ws, n, mx, my) * ef->vx_val (ef->ws, n, mx, my)
        - 0.5 * (+ ef->g_val (ef->ws, n, mx - 3, my) * ef->vx_val (ef->ws, n, mx - 3, my)
                 - 2 * ef->g_val (ef->ws, n, mx - 2, my) * ef->vx_val (ef->ws, n, mx - 2, my)
                 + ef->g_val (ef->ws, n, mx - 1, my) * ef->vx_val (ef->ws, n, mx - 1, my))
        + (2 - ef->g_val (ef->ws, n, mx, my)) * (
          + ef->vx_val (ef->ws, n, mx - 2, my)
          - 2 * ef->vx_val (ef->ws, n, mx - 1, my)
          + ef->vx_val (ef->ws, n, mx, my)
          - 0.5 * (+ ef->vx_val (ef->ws, n, mx - 3, my)
                   - 2 * ef->vx_val (ef->ws, n, mx - 2, my)
                   + ef->vx_val (ef->ws, n, mx - 1, my))))
      + 2 * ef->tau * ef->svr->f0 (t, x, y, ef->svr->mu, ef->svr->p_drv_type);

  FIX_UNUSED (full_rhs_val);
  ef->ws->rhs_vector[ef->row] = full_rhs_val;

  ef->row++;
}

void cdiff_solver_eq_5_5 (eq_filler_t *ef)
{
  int mx = ef->mx;
  int my = ef->my;
  double x = mx * ef->ws->hx;
  double y = my * ef->ws->hy;
  double t = ef->n * ef->ws->tau;
  int n = ef->n;
  double full_rhs_val;

  DEBUG_ASSERT (my == 0);

  ef->nnz = 0;

  fill_nz (ef->nz_values, ef->nz_cols, &ef->nnz,
            2 - ef->tau / ef->hy * ef->vy_val (ef->ws, n, mx, my),
            ef->g_curr);

  fill_nz (ef->nz_values, ef->nz_cols, &ef->nnz,
            + ef->tau / ef->hy * ef->vy_val (ef->ws, n, mx, my + 1),
            ef->g_top);

  fill_nz (ef->nz_values, ef->nz_cols, &ef->nnz,
           + 2 * ef->tau / ef->hy,
           ef->vy_top);

  fill_nz (ef->nz_values, ef->nz_cols, &ef->nnz,
           - 2 * ef->tau / ef->hy,
           ef->vy_curr);

  sparse_base_add_row (&ef->ws->matrix_base, ef->row, ef->nz_cols, ef->nz_values, ef->nnz);

  full_rhs_val =
      + 2 * ef->g_val (ef->ws, n, mx, my)
      + ef->tau / ef->hy * ef->g_val (ef->ws, n, mx, my) * (+ ef->vy_val (ef->ws, n, mx, my + 1)
                                                            - ef->vy_val (ef->ws, n, mx, my))
      + ef->tau / ef->hy * (
        + ef->g_val (ef->ws, n, mx, my) * ef->vy_val (ef->ws, n, mx, my)
        - 2 * ef->g_val (ef->ws, n, mx, my + 1) * ef->vy_val (ef->ws, n, mx, my + 1)
        + ef->g_val (ef->ws, n, mx, my + 2) * ef->vy_val (ef->ws, n, mx, my + 2)
        - 0.5 * (+ ef->g_val (ef->ws, n, mx + 1, my) * ef->vy_val (ef->ws, n, mx, my + 1)
                 - 2 * ef->g_val (ef->ws, n, mx, my + 2) * ef->vy_val (ef->ws, n, mx, my + 2)
                 + ef->g_val (ef->ws, n, mx, my + 3) * ef->vy_val (ef->ws, n, mx, my + 3))
        + (2 - ef->g_val (ef->ws, n, mx, my)) * (
          + ef->vy_val (ef->ws, n, mx, my)
          - 2 * ef->vy_val (ef->ws, n, mx, my + 1)
          + ef->vy_val (ef->ws, n, mx, my + 2)
          - 0.5 * (+ ef->vy_val (ef->ws, n, mx, my + 1)
                   - 2 * ef->vy_val (ef->ws, n, mx, my + 2)
                   + ef->vy_val (ef->ws, n, mx, my + 3))))
      + 2 * ef->tau * ef->svr->f0 (t, x, y, ef->svr->mu, ef->svr->p_drv_type);

  ef->ws->rhs_vector[ef->row] = full_rhs_val;

  ef->row++;
}

void cdiff_solver_eq_5_6 (eq_filler_t *ef)
{
  int mx = ef->mx;
  int my = ef->my;
  double x = mx * ef->ws->hx;
  double y = my * ef->ws->hy;
  double t = ef->n * ef->ws->tau;
  int n = ef->n;
  double full_rhs_val;

  ef->nnz = 0;

  DEBUG_ASSERT (my == ef->ws->MY || my == 2 * ef->ws->MY);

  fill_nz (ef->nz_values, ef->nz_cols, &ef->nnz,
            2 + ef->tau / ef->hy * ef->vy_val (ef->ws, n, mx, my),
            ef->g_curr);

  fill_nz (ef->nz_values, ef->nz_cols, &ef->nnz,
            - ef->tau / ef->hy * ef->vy_val (ef->ws, n, mx, my - 1),
            ef->g_bot);

  fill_nz (ef->nz_values, ef->nz_cols, &ef->nnz,
           - 2 * ef->tau / ef->hy,
           ef->vy_bot);

  fill_nz (ef->nz_values, ef->nz_cols, &ef->nnz,
           + 2 * ef->tau / ef->hy,
           ef->vy_curr);

  sparse_base_add_row (&ef->ws->matrix_base, ef->row, ef->nz_cols, ef->nz_values, ef->nnz);

  full_rhs_val =
      + 2 * ef->g_val (ef->ws, n, mx, my)
      + ef->tau / ef->hy * ef->g_val (ef->ws, n, mx, my) * (+ ef->vy_val (ef->ws, n, mx, my)
                                                            - ef->vy_val (ef->ws, n, mx, my - 1))
      - ef->tau / ef->hy * (
        + ef->g_val (ef->ws, n, mx, my - 2) * ef->vy_val (ef->ws, n, mx, my - 2)
        - 2 * ef->g_val (ef->ws, n, mx, my - 1) * ef->vy_val (ef->ws, n, mx, my - 1)
        + ef->g_val (ef->ws, n, mx, my) * ef->vy_val (ef->ws, n, mx, my)
        - 0.5 * (+ ef->g_val (ef->ws, n, mx, my - 3) * ef->vy_val (ef->ws, n, mx, my - 3)
                 - 2 * ef->g_val (ef->ws, n, mx, my - 2) * ef->vy_val (ef->ws, n, mx, my - 2)
                 + ef->g_val (ef->ws, n, mx, my - 1) * ef->vy_val (ef->ws, n, mx, my - 1))
        + (2 - ef->g_val (ef->ws, n, mx, my)) * (
          + ef->vy_val (ef->ws, n, mx, my - 2)
          - 2 * ef->vy_val (ef->ws, n, mx, my - 1)
          + ef->vy_val (ef->ws, n, mx, my)
          - 0.5 * (+ ef->vy_val (ef->ws, n, mx, my - 3)
                   - 2 * ef->vy_val (ef->ws, n, mx, my - 2)
                   + ef->vy_val (ef->ws, n, mx, my - 1))))
      + 2 * ef->tau * ef->svr->f0 (t, x, y, ef->svr->mu, ef->svr->p_drv_type);

  ef->ws->rhs_vector[ef->row] = full_rhs_val;

  ef->row++;
}

void cdiff_solver_eq_5_7 (eq_filler_t *ef)
{
  int mx = ef->mx;
  int my = ef->my;
  double x = mx * ef->ws->hx;
  double y = my * ef->ws->hy;
  double t = ef->n * ef->ws->tau;
  int n = ef->n;
  double mu_wave = ef->mu_wave;

  ef->nnz = 0;


  fill_nz (ef->nz_values, ef->nz_cols, &ef->nnz,
           + 6
           + 4 * ef->tau * mu_wave * (
             + 4 / (ef->hx * ef->hx)
             + 3 / (ef->hy * ef->hy)),
           ef->vx_curr);

  fill_nz (ef->nz_values, ef->nz_cols, &ef->nnz,
           - (+ ef->tau / ef->hx * (
                + ef->vx_val (ef->ws, n, mx - 1, my)
                + ef->vx_val (ef->ws, n, mx, my))
              + 8 * ef->tau * mu_wave / (ef->hx * ef->hx)),
           ef->vx_left);

  fill_nz (ef->nz_values, ef->nz_cols, &ef->nnz,
           - (+ 3 * ef->tau / (2 * ef->hy) * (
                + ef->vy_val (ef->ws, n, mx, my - 1)
                + ef->vy_val (ef->ws, n, mx, my))
              + 6 * ef->tau * mu_wave / (ef->hy * ef->hy)),
           ef->vx_bot);

  fill_nz (ef->nz_values, ef->nz_cols, &ef->nnz,
           + ef->tau / ef->hx * (
             + ef->vx_val (ef->ws, n, mx + 1, my)
             + ef->vx_val (ef->ws, n, mx, my))
           - 8 * ef->tau * mu_wave / (ef->hx * ef->hx),
           ef->vx_right);

  fill_nz (ef->nz_values, ef->nz_cols, &ef->nnz,
           + 3 * ef->tau / (2 * ef->hy) * (
             + ef->vy_val (ef->ws, n, mx, my + 1)
             + ef->vy_val (ef->ws, n, mx, my))
           - 6 * ef->tau * mu_wave / (ef->hy * ef->hy),
           ef->vx_top);

  fill_nz (ef->nz_values, ef->nz_cols, &ef->nnz,
           - 3 * ef->tau * ef->svr->p_drv (exp (ef->g_val (ef->ws, n, mx, my))) / ef->hx,
           ef->g_left);

  fill_nz (ef->nz_values, ef->nz_cols, &ef->nnz,
           3 * ef->tau * ef->svr->p_drv (exp (ef->g_val (ef->ws, n, mx, my))) / ef->hx,
           ef->g_right);

  sparse_base_add_row (&ef->ws->matrix_base, ef->row, ef->nz_cols, ef->nz_values, ef->nnz);

  ef->ws->rhs_vector[ef->row] =
      + 6 * ef->vx_val (ef->ws, n, mx, my)
      + 3 * ef->tau / (2 * ef->hy) * ef->vx_val (ef->ws, n, mx, my) * (
        + ef->vy_val (ef->ws, n, mx, my + 1)
        - ef->vy_val (ef->ws, n, mx, my - 1))
      + 6 * ef->tau * (ef->svr->mu / exp (ef->g_val (ef->ws, n, mx, my)) - mu_wave) * (+ 4. / (3. * ef->hx * ef->hx ) * (
                                                                                         + ef->vx_val (ef->ws, n, mx + 1, my)
                                                                                         - 2 * ef->vx_val (ef->ws, n, mx, my)
                                                                                         + ef->vx_val (ef->ws, n, mx - 1, my))
                                                                                       + 1. / (ef->hy * ef->hy) * (
                                                                                         + ef->vx_val (ef->ws, n, mx, my + 1)
                                                                                         - 2. * ef->vx_val (ef->ws, n, mx, my)
                                                                                         + ef->vx_val (ef->ws, n, mx, my - 1)))
      + ef->tau * ef->svr->mu / (2 * exp (ef->g_val (ef->ws, n, mx, my)) * ef->hx * ef->hy) * (
        + ef->vy_val (ef->ws, n, mx + 1, my + 1)
        - ef->vy_val (ef->ws, n, mx - 1, my + 1)
        - ef->vy_val (ef->ws, n, mx + 1, my - 1)
        + ef->vy_val (ef->ws, n, mx - 1, my - 1))
      + 6 * ef->tau * ef->svr->f1 (t, x, y, ef->svr->mu, ef->svr->p_drv_type);

  ef->row++;
}

void cdiff_solver_eq_5_8 (eq_filler_t *ef)
{
  int mx = ef->mx;
  int my = ef->my;
  double x = mx * ef->ws->hx;
  double y = my * ef->ws->hy;
  double t = ef->n * ef->ws->tau;
  int n = ef->n;
  double mu_wave = ef->mu_wave;

  ef->nnz = 0;

  fill_nz (ef->nz_values, ef->nz_cols, &ef->nnz,
           + 6
           + 4 * ef->tau * mu_wave * (
             + 3 / (ef->hx * ef->hx)
             + 4 / (ef->hy * ef->hy)),
           ef->vy_curr);

  fill_nz (ef->nz_values, ef->nz_cols, &ef->nnz,
           - (+ 3 * ef->tau / (2 * ef->hx) * (
                + ef->vx_val (ef->ws, n, mx - 1, my)
                + ef->vx_val (ef->ws, n, mx, my))
              + 6 * ef->tau * mu_wave / (ef->hx * ef->hx)),
           ef->vy_left);

  fill_nz (ef->nz_values, ef->nz_cols, &ef->nnz,
           - (ef->tau / (ef->hy) * (
                + ef->vy_val (ef->ws, n, mx, my - 1)
                + ef->vy_val (ef->ws, n, mx, my))
              + 8 * ef->tau * mu_wave / (ef->hy * ef->hy)),
           ef->vy_bot);

  fill_nz (ef->nz_values, ef->nz_cols, &ef->nnz,
           + 3 * ef->tau / (2 * ef->hx) * (
             + ef->vx_val (ef->ws, n, mx + 1, my)
             + ef->vx_val (ef->ws, n, mx, my))
           - 6 * ef->tau * mu_wave / (ef->hx * ef->hx),
           ef->vy_right);

  fill_nz (ef->nz_values, ef->nz_cols, &ef->nnz,
           + ef->tau / (ef->hy) * (
             + ef->vy_val (ef->ws, n, mx, my + 1)
             + ef->vy_val (ef->ws, n, mx, my))
           - 8 * ef->tau * mu_wave / (ef->hy * ef->hy),
           ef->vy_top);

  fill_nz (ef->nz_values, ef->nz_cols, &ef->nnz,
           - 3 * ef->tau * ef->svr->p_drv (exp (ef->g_val (ef->ws, n, mx, my))) / ef->hy,
           ef->g_bot);

  fill_nz (ef->nz_values, ef->nz_cols, &ef->nnz,
           3 * ef->tau * ef->svr->p_drv (exp (ef->g_val (ef->ws, n, mx, my))) / ef->hy,
           ef->g_top);

  sparse_base_add_row (&ef->ws->matrix_base, ef->row, ef->nz_cols, ef->nz_values, ef->nnz);

  ef->ws->rhs_vector[ef->row] =
      + 6 * ef->vy_val (ef->ws, n, mx, my)
      + 3 * ef->tau / (2 * ef->hx) * ef->vy_val (ef->ws, n, mx, my) * (
        + ef->vx_val (ef->ws, n, mx + 1, my)
        - ef->vx_val (ef->ws, n, mx - 1, my))
      + 6 * ef->tau * (ef->svr->mu / exp (ef->g_val (ef->ws, n, mx, my)) - mu_wave) * (1 / (ef->hx * ef->hx ) * (
                                                                                         + ef->vy_val (ef->ws, n, mx + 1, my)
                                                                                         - 2 * ef->vy_val (ef->ws, n, mx, my)
                                                                                         + ef->vy_val (ef->ws, n, mx - 1, my))
                                                                                       + 4. / (3 * ef->hy * ef->hy) * (
                                                                                         + ef->vy_val (ef->ws, n, mx, my + 1)
                                                                                         - 2 * ef->vy_val (ef->ws, n, mx, my)
                                                                                         + ef->vy_val (ef->ws, n, mx, my - 1)))
      + ef->tau * ef->svr->mu / (2 * exp (ef->g_val (ef->ws, n, mx, my)) * ef->hx * ef->hy) * (
        + ef->vx_val (ef->ws, n, mx + 1, my + 1)
        - ef->vx_val (ef->ws, n, mx - 1, my + 1)
        - ef->vx_val (ef->ws, n, mx + 1, my - 1)
        + ef->vx_val (ef->ws, n, mx - 1, my - 1))
      + 6 * ef->tau * ef->svr->f2 (t, x, y, ef->svr->mu, ef->svr->p_drv_type);

  ef->row++;
}

void cdiff_solver_eq_trivial (eq_filler_t *ef, grid_func_t func)
{
  int mx = ef->mx;
  int my = ef->my;
  int n = ef->n;

  ef->nnz = 0;
  fill_nz (ef->nz_values, ef->nz_cols, &ef->nnz,
           1,
           solver_workspace_func_col (ef->ws, ef->loc_layer_index, func));

  sparse_base_add_row (&ef->ws->matrix_base, ef->row, ef->nz_cols, ef->nz_values, ef->nnz);

  ef->ws->rhs_vector[ef->row] = solver_workspace_grid_val (ef->ws, n + 1, mx, my, func);

  /*if (func == grid_vx || func == grid_vy)
    DEBUG_ASSERT (math_is_null (ef->ws->rhs_vector[ef->row]));*/

  ef->row++;
}

void cdiff_solver_eq_border (eq_filler_t *ef, grid_area_t border)
{
  grid_func_t func;
  int prev_index = 0;
  int mx = ef->mx;
  int my = ef->my;

  DEBUG_ASSERT (border == border_rightmost || border == border_topmost);

  if (border == border_rightmost)
    {
      func = grid_vx;
      prev_index = ef->loc_layer_index - 1;
    }

  if (border == border_topmost)
    {
      func = grid_vy;
      prev_index = solver_workspace_final_index (ef->ws, 0, mx, my - 1);
    }


  ef->nnz = 0;

  fill_nz (ef->nz_values, ef->nz_cols, &ef->nnz,
           1,
           solver_workspace_func_col (ef->ws, ef->loc_layer_index, func));

  fill_nz (ef->nz_values, ef->nz_cols, &ef->nnz,
           -1,
           solver_workspace_func_col (ef->ws, prev_index, func));

  sparse_base_add_row (&ef->ws->matrix_base, ef->row, ef->nz_cols, ef->nz_values, ef->nnz);

  ef->ws->rhs_vector[ef->row] = 0;

  ef->row++;
}

double cdiff_solver_mu_wave (const central_diff_solver *solver)
{
  double max = 0;

  int layer_index = solver_workspace_layer_begin_index (&solver->ws, solver->layer - 1);
  int i;
  int size = solver->ws.layer_size;

  for (i = 0; i < size; i++)
    {
      if (i == 0)
        max = -solver->ws.g[layer_index];
      else
        {
          double val = -solver->ws.g[layer_index];

          max = (max > val) ? max : val;
        }

      layer_index++;
    }

  return exp (max)  * solver->mu;
}

void cdiff_solver_fill_matrix_w_rhs (central_diff_solver *solver)
{

  eq_filler_t eqfill_obj;
  eq_filler_t *ef = &eqfill_obj;
  int mx_max;

  eq_filler_init (ef, solver);
  sparse_base_to_init_state (&ef->ws->matrix_base);

  mx_max = BOT_ROW_SQUARES_COUNT * ef->ws->MX;

  ef->my = 0;

  ef->mx = 0;
  /*left_bot*/
  if (solver->ws.mode == solve_mode)
    cdiff_solver_eq_trivial (ef, grid_g);
  else
    cdiff_solver_eq_5_3 (ef);

  cdiff_solver_eq_trivial (ef, grid_vx);
  cdiff_solver_eq_trivial (ef, grid_vy);
  eq_filler_inc_layer_index (ef);

  for (ef->mx = 1; ef->mx < mx_max; ef->mx++)
    {
      /*botmost*/
      cdiff_solver_eq_5_5 (ef);
      cdiff_solver_eq_trivial (ef, grid_vx);
      cdiff_solver_eq_trivial (ef, grid_vy);
      eq_filler_inc_layer_index (ef);
    }

  ef->mx = mx_max;
  /*right_bot*/
  cdiff_solver_eq_5_4 (ef);
  if (solver->ws.mode == solve_mode)
    cdiff_solver_eq_border (ef, border_rightmost);
  else
    cdiff_solver_eq_trivial (ef, grid_vx);

  cdiff_solver_eq_trivial (ef, grid_vy);
  eq_filler_inc_layer_index (ef);

  for (ef->my = 1; ef->my < ef->ws->MY; ef->my++)
    {
      ef->mx = 0;
      /*left_most*/
      if (solver->ws.mode == solve_mode)
        cdiff_solver_eq_trivial (ef, grid_g);
      else
        cdiff_solver_eq_5_3 (ef);

      cdiff_solver_eq_trivial (ef, grid_vx);
      cdiff_solver_eq_trivial (ef, grid_vy);
      eq_filler_inc_layer_index (ef);

      for (ef->mx = 1; ef->mx < mx_max; ef->mx++)
        {
          /*internal*/
          cdiff_solver_eq_5_2 (ef);
          cdiff_solver_eq_5_7 (ef);
          cdiff_solver_eq_5_8 (ef);
          eq_filler_inc_layer_index (ef);
        }

      ef->mx = mx_max;
      /*rightmost*/
      cdiff_solver_eq_5_4 (ef);
      if (solver->ws.mode == solve_mode)
        cdiff_solver_eq_border (ef, border_rightmost);
      else
        cdiff_solver_eq_trivial (ef, grid_vx);

      cdiff_solver_eq_trivial (ef, grid_vy);
      eq_filler_inc_layer_index (ef);
    }

  ef->my = ef->ws->MY;
  ef->mx = 0;
  /*left_top*/
  if (solver->ws.mode == solve_mode)
    cdiff_solver_eq_trivial (ef, grid_g);
  else
    cdiff_solver_eq_5_3 (ef);

  cdiff_solver_eq_trivial (ef, grid_vx);
  cdiff_solver_eq_trivial (ef, grid_vy);
  eq_filler_inc_layer_index (ef);

  for (ef->mx = 1; ef->mx < ef->ws->MX; ef->mx++)
    {
      /*left_hor*/
      cdiff_solver_eq_5_6 (ef);
      cdiff_solver_eq_trivial (ef, grid_vx);
      cdiff_solver_eq_trivial (ef, grid_vy);
      eq_filler_inc_layer_index (ef);
    }

  ef->mx = ef->ws->MX;
  /*left_hor_angle*/
  cdiff_solver_eq_5_6 (ef);
  cdiff_solver_eq_trivial (ef, grid_vx);
  cdiff_solver_eq_trivial (ef, grid_vy);
  eq_filler_inc_layer_index (ef);

  for (ef->mx = ef->ws->MX + 1; ef->mx < 2 * ef->ws->MX; ef->mx++)
    {
      /*internal*/
      cdiff_solver_eq_5_2 (ef);
      cdiff_solver_eq_5_7 (ef);
      cdiff_solver_eq_5_8 (ef);
      eq_filler_inc_layer_index (ef);
    }

  ef->mx = 2 * ef->ws->MX;
  /*right_hor_angle*/
  cdiff_solver_eq_5_6 (ef);
  cdiff_solver_eq_trivial (ef, grid_vx);
  cdiff_solver_eq_trivial (ef, grid_vy);
  eq_filler_inc_layer_index (ef);

  for (ef->mx = 2 * ef->ws->MX + 1; ef->mx < mx_max; ef->mx++)
    {
      /*right_hor*/
      cdiff_solver_eq_5_6 (ef);
      cdiff_solver_eq_trivial (ef, grid_vx);
      cdiff_solver_eq_trivial (ef, grid_vy);
      eq_filler_inc_layer_index (ef);
    }

  ef->mx = mx_max;
  /*right_top*/
  cdiff_solver_eq_5_4 (ef);
  if (solver->ws.mode == solve_mode)
    cdiff_solver_eq_border (ef, border_rightmost);
  else
    cdiff_solver_eq_trivial (ef, grid_vx);
  cdiff_solver_eq_trivial (ef, grid_vy);
  eq_filler_inc_layer_index (ef);

  for (ef->my = ef->ws->MY + 1; ef->my < 2 * ef->ws->MY; ef->my++)
    {
      ef->mx = ef->ws->MX;
      /*left_vert*/
      cdiff_solver_eq_5_3 (ef);
      cdiff_solver_eq_trivial (ef, grid_vx);
      cdiff_solver_eq_trivial (ef, grid_vy);
      eq_filler_inc_layer_index (ef);

      for (ef->mx = ef->ws->MX + 1; ef->mx < 2 * ef->ws->MX; ef->mx++)
        {
          /*internal*/
          cdiff_solver_eq_5_2 (ef);
          cdiff_solver_eq_5_7 (ef);
          cdiff_solver_eq_5_8 (ef);
          eq_filler_inc_layer_index (ef);
        }

      ef->mx = 2 * ef->ws->MX;
      /*right_vert*/
      cdiff_solver_eq_5_4 (ef);
      cdiff_solver_eq_trivial (ef, grid_vx);
      cdiff_solver_eq_trivial (ef, grid_vy);
      eq_filler_inc_layer_index (ef);
    }

  ef->my = 2 * ef->ws->MY;

  ef->mx = ef->ws->MX;
  /*top_left*/
  cdiff_solver_eq_5_6 (ef);
  cdiff_solver_eq_trivial (ef, grid_vx);
  if (solver->ws.mode == solve_mode)
    cdiff_solver_eq_border (ef, border_topmost);
  else
    cdiff_solver_eq_trivial (ef, grid_vy);
  eq_filler_inc_layer_index (ef);

  for (ef->mx = ef->ws->MX + 1; ef->mx < 2 * ef->ws->MX; ef->mx++)
    {
      /*topmost*/
      cdiff_solver_eq_5_6 (ef);
      cdiff_solver_eq_trivial (ef, grid_vx);
      if (solver->ws.mode == solve_mode)
        cdiff_solver_eq_border (ef, border_topmost);
      else
        cdiff_solver_eq_trivial (ef, grid_vy);
      eq_filler_inc_layer_index (ef);
    }

  ef->mx = 2 * ef->ws->MX;
  /*top_right*/
  cdiff_solver_eq_5_6 (ef);
  cdiff_solver_eq_trivial (ef, grid_vx);
  if (solver->ws.mode == solve_mode)
    cdiff_solver_eq_border (ef, border_topmost);
  else
    cdiff_solver_eq_trivial (ef, grid_vy);
}

void cdiff_solver_solve_system (central_diff_solver *solver)
{
  if (solver->ws.linear_solver == custom_cgs)
    {
      cgs_solver *linear_solver = &solver->ws.cgs_linear_solver;
      cgs_solver_error_t error;
      vector_double_t x_init = NULL;
      vector_double_t DELETE_LATER = VECTOR_CREATE (double, solver->ws.matrix_size);
      double c_res;
#ifdef DEBUG
      int i;
#endif

      if (solver->layer > 1)
        x_init = solver->ws.vector_to_compute;

      msr_fill_from_sparse_base (&solver->ws.matrix, &solver->ws.matrix_base);

      /*TEMCODE_BEGIN*/
#ifdef DEBUG
      if (solver->ws.mode == test_mode)
        {
          x_init = VECTOR_CREATE (double, solver->ws.matrix_size);
          for (i = 0; i < solver->ws.matrix_size; i += 3)
            {int mx; int my;
              double x;
              double y;
              double t;

              solver_workspace_get_mx_my (&solver->ws, i / 3, &mx, &my);

              x = solver->ws.hx * mx;
              y = solver->ws.hy * my;
              t = solver->layer * solver->ws.tau;

              x_init[i] = solver->test_solution_g (t, x, y);
              x_init[i + 1] = solver->test_solution_vx (t, x, y);
              x_init[i + 2] = solver->test_solution_vy (t, x, y);
            }
        }
#endif

      /*TEMCODE_END*/

      if (solver->ws.MX == 3 && solver->ws.MY == 3)
        msr_dump (&solver->ws.matrix, solver->ws.log_file);

      error = cgs_solver_solve (linear_solver, &solver->ws.matrix, solver->ws.rhs_vector, x_init, solver->ws.vector_to_compute);

      fprintf (stdout, "Layer: %d, MX = %d, MY = %d, N = %d\n", solver->layer, solver->ws.MX, solver->ws.MY, solver->ws.N);
      if (error)
        {
          fprintf (stdout, "Failed to converge\n");
          DEBUG_PAUSE ("failed");
        }

      /*msr_mult_vector (&solver->ws.matrix, solver->ws.vector_to_compute, DELETE_LATER);*/

      /*linear_combination_w_override_1 (x_init, -1, solver->ws.vector_to_compute, solver->ws.matrix_size);*/
      VECTOR_DESTROY (DELETE_LATER);
      /*c_res = c_norm (x_init, solver->ws.matrix_size);*/
      /*fprintf (stdout, "Layer: %d, C discrepancy: %f\n", solver->layer, c_res);*/
      FIX_UNUSED (c_res);
      if (error)
        return;
    }
  else
    {
      laspack_matrix_init (&solver->ws.matrix_l, &solver->ws.matrix_base);
      laspack_vector_fill (&solver->ws.rhs_vector_l, solver->ws.rhs_vector);

      BiCGSTABIter (&solver->ws.matrix_l.raw, &solver->ws.vector_to_compute_l.raw, &solver->ws.rhs_vector_l.raw, 2000, JacobiPrecond, 1.2);

      laspack_matrix_destroy (&solver->ws.matrix_l);
    }

  solver_workspace_fill_layer (&solver->ws, solver->layer);
}

#include "central_diff_solver_private.h"

#include "central_diff_solver_private.h"
#include "input/rhs.h"
#include "cgs_solver.h"
#include "common/debug_utils.h"
#include "sparse/sparse_base_format.h"
#include "common/math_utils.h"

/*
 * TODO
 *  memory alloc arguments in init, others in compute
 *  use f0
 *  fill_matrix
 *  solve_matrix
 *  compute difnorms
 */

void fill_nz (double *vals, int *cols, int *nnz, double val, int col)
{
  if (math_is_null (val))
    return;

  vals[*nnz] = val;
  cols[*nnz] = col;
  (*nnz)++;
}

int cdiff_solver_init (central_diff_solver *solver,
                       solver_mode_t mode,
                       int M1,
                       int M2,
                       int N,
                       double X,
                       double Y,
                       double T,
                       double border_omega,
                       time_layer_func_t test_solution_g,
                       time_layer_func_t test_solution_vx,
                       time_layer_func_t test_solution_vy)
{
  solver_workspace_data_init (&solver->ws, mode, M1, M2, N, X, Y, T, border_omega);

  if (mode == test_mode)
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

  return 0;
}



void cdiff_solver_destroy (central_diff_solver *solver)
{
  solver_workspace_data_destroy (&solver->ws);
}



int cdiff_solver_compute (central_diff_solver *solver,
                          pressure_func_t p_func,
                          double mu,
                          rhs_func_t f0,
                          rhs_func_t f1,
                          rhs_func_t f2,
                          layer_func_t start_vx,
                          layer_func_t start_vy,
                          layer_func_t start_g)
{
  int i;

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

  solver->mu = mu;

  cdiff_solver_init_first_layer (solver);
  cdiff_solver_init_borders (solver);


  for (solver->layer = 1; solver->layer <= solver->ws.N; solver->layer++)
    {
      cdiff_solver_fill_matrix_w_rhs (solver);
      cdiff_solver_fill_rhs (solver);
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
  int i = 0;

  for (my = 0; my <= solver->ws.MY; my++)
    {
      for (mx = 0; mx <= solver->ws.MX * BOT_ROW_SQUARES_COUNT; mx++)
        {
          x = mx * solver->ws.hx;
          y = my * solver->ws.hy;
          solver->ws.vector_to_compute[UNKNOWN_FUNCTIONS_COUNT * i]     = solver->start_g  (x, y);
          solver->ws.vector_to_compute[UNKNOWN_FUNCTIONS_COUNT * i + 1] = solver->start_vx (x, y);
          solver->ws.vector_to_compute[UNKNOWN_FUNCTIONS_COUNT * i + 2] = solver->start_vy (x, y);
          i++;
        }
    }

  for (my = solver->ws.MY + 1; my <= 2 * solver->ws.MY; my++)
    {
      for (mx = solver->ws.MX; mx <= 2 * solver->ws.MX; mx++)
        {
          x = mx * solver->ws.hx;
          y = my * solver->ws.hy;
          solver->ws.vector_to_compute[UNKNOWN_FUNCTIONS_COUNT * i]     = solver->start_g  (x, y);
          solver->ws.vector_to_compute[UNKNOWN_FUNCTIONS_COUNT * i + 1] = solver->start_vx (x, y);
          solver->ws.vector_to_compute[UNKNOWN_FUNCTIONS_COUNT * i + 2] = solver->start_vy (x, y);
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

void eq_filler_init (eq_filler_t *ef, const central_diff_solver *solver)
{
  ef->ws = &solver->ws;
  ef->svr = solver;
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
  int mx, my;
  int loc_top = 0;
  int loc_bot = 0;
  int layer_begin_index = solver_workspace_layer_begin_index (ef->ws, ef->n);

  ef->loc_layer_index++;

  solver_workspace_get_mx_my (ef->ws, ef->loc_layer_index, &mx, &my);

  if (my < ef->ws->MY
      || (
        mx >= ef->ws->MX
        && mx <= 2 * ef->ws->MX
        && my < 2 * ef->ws->MY))
    loc_top = solver_workspace_final_index (ef->ws, ef->n, mx, my + 1) - layer_begin_index;

  if (my > 0 && my <= ef->ws->MY
      || (
        mx >= ef->ws->MX
        && mx <= 2 * ef->ws->MX
        && my <= 2 * ef->ws->MY))
    loc_bot = solver_workspace_final_index (ef->ws, ef->n, mx, my - 1) - layer_begin_index;



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
  ef->nnz = 0;

  fill_nz (ef->nz_values, ef->nz_cols, &ef->nnz,
           4,
           ef->g_curr);
  fill_nz (ef->nz_values, ef->nz_cols, &ef->nnz,
           - ef->tau / ef->hx * (
             + ef->vx_val (ef->n, ef->mx - 1, ef->my)
             + ef->vx_val (ef->n, ef->mx, ef->my)),
           ef->g_left);
  fill_nz (ef->nz_values, ef->nz_cols, &ef->nnz,
           - ef->tau / ef->hy * (
             + ef->vy_val (ef->n, ef->mx, ef->my - 1)
             + ef->vy_val (ef->n, ef->mx, ef->my)),
           ef->g_bot);
  fill_nz (ef->nz_values, ef->nz_cols, &ef->nnz,
           ef->tau / ef->hx * (
             + ef->vx_val (ef->n, ef->mx + 1, ef->my)
             + ef->vx_val (ef->n, ef->mx, ef->my)),
           ef->g_right);
  fill_nz (ef->nz_values, ef->nz_cols, &ef->nnz,
           ef->tau / ef->hy * (
             + ef->vy_val (ef->n, ef->mx, ef->my + 1)
             + ef->vy_val (ef->n, ef->mx, ef->my)),
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
      + 4 * ef->g_val (ef->n, ef->mx, ef->my)
      + ef->tau * ef->g_val (ef->n, ef->mx, ef->my) * (
        + (ef->vx_val (ef->n, ef->mx + 1, ef->my) - ef->vx_val (ef->n, ef->mx - 1, ef->my)) / ef->hx
        + (ef->vy_val (ef->n, ef->mx, ef->my + 1) - ef->vy_val (ef->n, ef->mx, ef->my - 1))/ ef->hy);

  ef->row++;
}

void cdiff_solver_eq_5_3 (eq_filler_t *ef)
{
  ef->nnz = 0;

  fill_nz (ef->nz_values, ef->nz_cols, &ef->nnz,
            2,
            ef->g_curr);

  fill_nz (ef->nz_values, ef->nz_cols, &ef->nnz,
            ef->tau / ef->hx * ef->vx_val (ef->ws, ef->n, ef->mx + 1, ef->my),
            ef->g_right);

  fill_nz (ef->nz_values, ef->nz_cols, &ef->nnz,
           2 * ef->tau / ef->hx,
           ef->vx_right);

  sparse_base_add_row (&ef->ws->matrix_base, ef->row, ef->nz_cols, ef->nz_values, ef->nnz);

  ef->ws->rhs_vector[ef->row] =
      + 2 * ef->g_val (n, ef->mx, ef->my)
      + ef->tau / ef->hx * ef->g_val (ef->n, ef->mx, ef->my) * ef->vx_val (ef->n, ef->mx + 1, ef->my)
      + 2 * ef->tau / ef->hx * (
        - 2.5 * ef->g_val (ef->n, ef->mx + 1, ef->my) * ef->vx_val (ef->n, ef->mx + 1, ef->my)
        + 2 * ef->g_val (ef->n, ef->mx + 2, ef->my) * ef->vx_val (ef->n, ef->mx + 2, ef->my)
        - 0.5 * ef->g_val (ef->n, ef->mx + 3, ef->my) * ef->vx_val (ef->n, ef->mx + 3, ef->my)
        + (2 - ef->g_val (ef->n, ef->mx, ef->my)) * (
          - 2.5 * ef->vx_val (ef->n, ef->mx + 1, ef->my)
          + 2 * ef->vx_val (ef->n, ef->mx + 2, ef->my)
          - 0.5 * ef->vx_val (ef->n, ef->mx + 3, ef->my)));

  ef->row++;
}

void cdiff_solver_eq_5_4 (eq_filler_t *ef)
{
  ef->nnz = 0;

  fill_nz (ef->nz_values, ef->nz_cols, &ef->nnz,
            2,
            ef->g_curr);

  fill_nz (ef->nz_values, ef->nz_cols, &ef->nnz,
            -ef->tau / ef->hx * ef->vx_val (ef->ws, ef->n, ef->mx - 1, ef->my),
            ef->g_left);

  fill_nz (ef->nz_values, ef->nz_cols, &ef->nnz,
           - 2 * ef->tau / ef->hx,
           ef->vx_left);

  sparse_base_add_row (&ef->ws->matrix_base, ef->row, ef->nz_cols, ef->nz_values, ef->nnz);

  ef->ws->rhs_vector[ef->row] =
      + 2 * ef->g_val (n, ef->mx, ef->my)
      - ef->tau / ef->hx * ef->g_val (ef->n, ef->mx, ef->my) * ef->vx_val (ef->n, ef->mx - 1, ef->my)
      - 2 * ef->tau / ef->hx * (
        - 2.5 * ef->g_val (ef->n, ef->mx - 1, ef->my) * ef->vx_val (ef->n, ef->mx - 1, ef->my)
        + 2 * ef->g_val (ef->n, ef->mx - 2, ef->my) * ef->vx_val (ef->n, ef->mx - 2, ef->my)
        - /*0.5*/ ef->g_val (ef->n, ef->mx - 3, ef->my) * ef->vx_val (ef->n, ef->mx - 3, ef->my)
        + (2 - ef->g_val (ef->n, ef->mx, ef->my)) * (
          - 2.5 * ef->vx_val (ef->n, ef->mx - 1, ef->my)
          + 2 * ef->vx_val (ef->n, ef->mx - 2, ef->my)
          - 0.5 * ef->vx_val (ef->n, ef->mx - 3, ef->my)));

  ef->row++;
}

void cdiff_solver_eq_5_5 (eq_filler_t *ef)
{
  ef->nnz = 0;

  fill_nz (ef->nz_values, ef->nz_cols, &ef->nnz,
            2,
            ef->g_curr);

  fill_nz (ef->nz_values, ef->nz_cols, &ef->nnz,
            + /*2*/ ef->tau / ef->hy * ef->vy_val (ef->ws, ef->n, ef->mx, ef->my + 1),
            ef->g_top);

  fill_nz (ef->nz_values, ef->nz_cols, &ef->nnz,
           + 2 * ef->tau / ef->hy,
           ef->vy_top);

  sparse_base_add_row (&ef->ws->matrix_base, ef->row, ef->nz_cols, ef->nz_values, ef->nnz);

  ef->ws->rhs_vector[ef->row] =
      + 2 * ef->g_val (n, ef->mx, ef->my)
      + ef->tau / ef->hy * ef->g_val (ef->n, ef->mx, ef->my) * ef->vy_val (ef->n, ef->mx, ef->my + 1)
      + 2 * ef->tau / ef->hy * (
        - 2.5 * ef->g_val (ef->n, ef->mx, ef->my + 1) * ef->vy_val (ef->n, ef->mx, ef->my + 1)
        + 2 * ef->g_val (ef->n, ef->mx, ef->my + 2) * ef->vy_val (ef->n, ef->mx, ef->my + 2)
        - 0.5  * ef->g_val (ef->n, ef->mx, ef->my + 3) * ef->vy_val (ef->n, ef->mx, ef->my + 3)
        + (2 - ef->g_val (ef->n, ef->mx, ef->my)) * (
          - 2.5 * ef->vy_val (ef->n, ef->mx, ef->my + 1)
          + 2 * ef->vy_val (ef->n, ef->mx, ef->my + 2)
          - 0.5 * ef->vy_val (ef->n, ef->mx, ef->my + 3)));

  ef->row++;
}

void cdiff_solver_eq_5_6 (eq_filler_t *ef)
{
  ef->nnz = 0;

  fill_nz (ef->nz_values, ef->nz_cols, &ef->nnz,
            2,
            ef->g_curr);

  fill_nz (ef->nz_values, ef->nz_cols, &ef->nnz,
            - 2 * ef->tau / ef->hy * ef->vy_val (ef->ws, ef->n, ef->mx, ef->my - 1),
            ef->g_bot);

  fill_nz (ef->nz_values, ef->nz_cols, &ef->nnz,
           - 2 * ef->tau / ef->hy,
           ef->vy_bot);

  sparse_base_add_row (&ef->ws->matrix_base, ef->row, ef->nz_cols, ef->nz_values, ef->nnz);

  ef->ws->rhs_vector[ef->row] =
      + 2 * ef->g_val (n, ef->mx, ef->my)
      - ef->tau / ef->hy * ef->g_val (ef->n, ef->mx, ef->my) * ef->vy_val (ef->n, ef->mx, ef->my - 1)
      - /*2*/ ef->tau / ef->hy * (
        - 2.5 * ef->g_val (ef->n, ef->mx, ef->my - 1) * ef->vy_val (ef->n, ef->mx, ef->my - 1)
        + 2 * ef->g_val (ef->n, ef->mx, ef->my - 2) * ef->vy_val (ef->n, ef->mx, ef->my - 2)
        - 0.5  * ef->g_val (ef->n, ef->mx, ef->my - 3) * ef->vy_val (ef->n, ef->mx, ef->my - 3)
        + (2 - ef->g_val (ef->n, ef->mx, ef->my)) * (
          - 2.5 * ef->vy_val (ef->n, ef->mx, ef->my - 1)
          + 2 * ef->vy_val (ef->n, ef->mx, ef->my - 2)
          - 0.5 * ef->vy_val (ef->n, ef->mx, ef->my - 3)));

  ef->row++;
}

void cdiff_solver_eq_5_7 (eq_filler_t *ef)
{
  ef->nnz = 0;

  fill_nz (ef->nz_values, ef->nz_cols, &ef->nnz,
           + 6
           + 4 * ef->tau * cdiff_solver_mu_wave (ef->svr) * (
             + 4 / (ef->hx * ef->hx)
             + 3 / (ef->hy * ef->hy)),
           ef->vx_curr);

  fill_nz (ef->nz_values, ef->nz_cols, &ef->nnz,
           - (+ ef->tau / ef->hx * (
                + ef->vx_val (ef->n, ef->mx - 1, ef->my)
                + ef->vx_val (ef->n, ef->mx, ef->my))
           + 8 * ef->tau * cdiff_solver_mu_wave (ef->svr) / (ef->hx * ef->hx)),
           ef->vx_left);

  fill_nz (ef->nz_values, ef->nz_cols, &ef->nnz,
           - (3 * ef->tau / (2 * ef->hy) * (
                + ef->vy_val (ef->n, ef->mx, ef->my - 1)
                + ef->vy_val (ef->n, ef->mx, ef->my))
              + 6 * ef->tau * cdiff_solver_mu_wave (ef->svr) / (ef->hy * ef->hy)),
           ef->vx_bot);

  fill_nz (ef->nz_values, ef->nz_cols, &ef->nnz,
           + ef->tau / ef->hx * (
             + ef->vx_val (ef->n, ef->mx + 1, ef->my)
             + ef->vx_val (ef->n, ef->mx, ef->my))
           - 8 * ef->tau * cdiff_solver_mu_wave (ef->svr) / (ef->hx * ef->hx),
           ef->vx_right);

  fill_nz (ef->nz_values, ef->nz_cols, &ef->nnz,
           + 3 * ef->tau / (2 * ef->hy) * (
             + ef->vy_val (ef->n, ef->mx, ef->my + 1)
             + ef->vy_val (ef->n, ef->mx, ef->my))
           - 6 * ef->tau * cdiff_solver_mu_wave (ef->svr) / (ef->hy * ef->hy),
           ef->vx_top);

  fill_nz (ef->nz_values, ef->nz_cols, &ef->nnz,
           - 3 * ef->tau * ef->svr->p_drv (exp (ef->g_val (ef->n, ef->mx, ef->my))) / ef->hx,
           ef->g_left);

  fill_nz (ef->nz_values, ef->nz_cols, &ef->nnz,
           3 * ef->tau * ef->svr->p_drv (exp (ef->g_val (ef->n, ef->mx, ef->my))) / ef->hx,
           ef->g_right);

  sparse_base_add_row (&ef->ws->matrix_base, ef->row, ef->nz_cols, ef->nz_values, ef->nnz);

  ef->ws->rhs_vector[ef->row] =
      + 6 * ef->vx_val (ef->n, ef->mx, ef->my)
      + 3 * ef->tau / (2 * ef->hy) * ef->vx_val (ef->n, ef->mx, ef->my) * (
        + ef->vy_val (ef->n, ef->mx, ef->my + 1)
        - ef->vy_val (ef->n, ef->mx, ef->my - 1))
      + 6 * ef->tau * (ef->svr->mu / exp (ef->g_val (ef->n, ef->mx, ef->my)) - cdiff_solver_mu_wave (ev->svr)) * 4 / (3 * ef->hx * ef->hx ) * (
        + ef->vx_val (ef->n, ef->mx + 1, ef->my)
        - 2 * ef->vx_val (ef->n, ef->mx, ef->my)
        + ef->vx_val (ef->n, ef->mx - 1, ef->my))
      + 1. / (ef->hy * ef->hy) * (
        + ef->vx_val (ef->n, ef->mx, ef->my + 1)
        - 2 * ef->vx_val (ef->n, ef->mx, ef->my)
        + ef->vx_val (ef->n, ef->mx, ef->my - 1))
      + ef->tau * ef->svr->mu / (2 * exp (ef->g_val (ef->n, ef->mx, ef->my)) * ef->hx * ef->hy) * (
        + ef->vy_val (ef->n, ef->mx + 1, ef->my + 1)
        - ef->vy_val (ef->n, ef->mx - 1, ef->my + 1)
        - ef->vy_val (ef->n, ef->mx + 1, ef->my - 1)
        + ef->vy_val (ef->n, ef->mx - 1, ef->my - 1))
      + 6 * ef->tau * ef->svr->f1 (ef->n, ef->mx, ef->my, ef->svr->mu, ef->svr->p_drv_type);

  ef->row++;
}

void cdiff_solver_eq_5_8 (eq_filler_t *ef)
{
  ef->nnz = 0;

  fill_nz (ef->nz_values, ef->nz_cols, &ef->nnz,
           + 6
           + 4 * ef->tau * cdiff_solver_mu_wave (ef->svr) * (
             + 3 / (ef->hx * ef->hx)
             + 4 / (ef->hy * ef->hy)),
           ef->vy_curr);

  fill_nz (ef->nz_values, ef->nz_cols, &ef->nnz,
           - (+ 3 * ef->tau / (2 * ef->hx) * (
                + ef->vx_val (ef->n, ef->mx - 1, ef->my)
                + ef->vx_val (ef->n, ef->mx, ef->my))
           + 6 * ef->tau * cdiff_solver_mu_wave (ef->svr) / (ef->hx * ef->hx)),
           ef->vy_left);

  fill_nz (ef->nz_values, ef->nz_cols, &ef->nnz,
           - (ef->tau / (ef->hy) * (
                + ef->vy_val (ef->n, ef->mx, ef->my - 1)
                + ef->vy_val (ef->n, ef->mx, ef->my))
              + 8 * ef->tau * cdiff_solver_mu_wave (ef->svr) / (ef->hy * ef->hy)),
           ef->vy_bot);

  fill_nz (ef->nz_values, ef->nz_cols, &ef->nnz,
           + 3 * ef->tau / (2 * ef->hx) * (
             + ef->vx_val (ef->n, ef->mx + 1, ef->my)
             + ef->vx_val (ef->n, ef->mx, ef->my))
           - 6 * ef->tau * cdiff_solver_mu_wave (ef->svr) / (ef->hx * ef->hx),
           ef->vy_right);

  fill_nz (ef->nz_values, ef->nz_cols, &ef->nnz,
           + ef->tau / (2 * ef->hy) * (
             + ef->vy_val (ef->n, ef->mx, ef->my + 1)
             + ef->vy_val (ef->n, ef->mx, ef->my))
           - 8 * ef->tau * cdiff_solver_mu_wave (ef->svr) / (ef->hy * ef->hy),
           ef->vy_top);

  fill_nz (ef->nz_values, ef->nz_cols, &ef->nnz,
           - 3 * ef->tau * ef->svr->p_drv (exp (ef->g_val (ef->n, ef->mx, ef->my))) / ef->hy,
           ef->g_bot);

  fill_nz (ef->nz_values, ef->nz_cols, &ef->nnz,
           3 * ef->tau * ef->svr->p_drv (exp (ef->g_val (ef->n, ef->mx, ef->my))) / ef->hy,
           ef->g_top);

  sparse_base_add_row (&ef->ws->matrix_base, ef->row, ef->nz_cols, ef->nz_values, ef->nnz);

  ef->ws->rhs_vector[ef->row] =
      + 6 * ef->vy_val (ef->n, ef->mx, ef->my)
      + 3 * ef->tau / (2 * ef->hx) * ef->vy_val (ef->n, ef->mx, ef->my) * (
        + ef->vx_val (ef->n, ef->mx + 1, ef->my)
        - ef->vy_val (ef->n, ef->mx - 1, ef->my))
      + 6 * ef->tau * (ef->svr->mu / exp (ef->g_val (ef->n, ef->mx, ef->my)) - cdiff_solver_mu_wave (ev->svr)) * 1 / (ef->hx * ef->hx ) * (
        + ef->vx_val (ef->n, ef->mx + 1, ef->my)
        - 2 * ef->vx_val (ef->n, ef->mx, ef->my)
        + ef->vx_val (ef->n, ef->mx - 1, ef->my))
      + 4. / (3 * ef->hy * ef->hy) * (
        + ef->vy_val (ef->n, ef->mx, ef->my + 1)
        - 2 * ef->vy_val (ef->n, ef->mx, ef->my)
        + ef->vy_val (ef->n, ef->mx, ef->my - 1))
      + ef->tau * ef->svr->mu / (2 * exp (ef->g_val (ef->n, ef->mx, ef->my)) * ef->hx * ef->hy) * (
        + ef->vx_val (ef->n, ef->mx + 1, ef->my + 1)
        - ef->vx_val (ef->n, ef->mx - 1, ef->my + 1)
        - ef->vx_val (ef->n, ef->mx + 1, ef->my - 1)
        + ef->vx_val (ef->n, ef->mx - 1, ef->my - 1))
      + 6 * ef->tau * ef->svr->f1 (ef->n, ef->mx, ef->my, ef->svr->mu, ef->svr->p_drv_type);

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
      double val = exp (-solver->ws.g[layer_index]);

      max = (max > val) ? max : val;

      layer_index++;
    }

  return max * solver->mu;
}

void cdiff_solver_fill_matrix_w_rhs (central_diff_solver *solver)
{

  eq_filler_t eqfill_obj;
  eq_filler_t *ef = &eqfill_obj;
  int mx_max = BOT_ROW_SQUARES_COUNT * ef->ws->MX;

  eq_filler_init (ef, solver);



  ef->my = 0;

  ef->mx = 0;
  /*body*/
  eq_filler_inc_layer_index (ef);

  for (ef->mx = 1; ef->mx < mx_max; ef->mx++)
    {
      /*body*/
      eq_filler_inc_layer_index (ef);
    }

  ef->mx = mx_max;
  /*body*/
  eq_filler_inc_layer_index (ef);

  for (ef->my = 1; ef->my < ef->ws->MY; ef->my++)
    {
      ef->mx = 0;
      /*body*/
      eq_filler_inc_layer_index (ef);

      for (ef->mx = 1; ef->mx < mx_max; ef->mx++)
        {
          /*body*/
          eq_filler_inc_layer_index (ef);
        }

      ef->mx = mx_max;
      /*body*/
      eq_filler_inc_layer_index (ef);
    }

  ef->my = ef->ws->MY;
  ef->mx = 0;
  /*body*/
  eq_filler_inc_layer_index (ef);

  for (ef->mx = 1; ef->mx < ef->ws->MX; ef->mx++)
    {
      /*body*/
      eq_filler_inc_layer_index (ef);
    }

  ef->mx = ef->ws->MX;
  /*body*/
  eq_filler_inc_layer_index (ef);

  for (ef->mx = ef->ws->MX + 1; ef->mx < 2 * ef->ws->MX; ef->mx++)
    {
      /*body*/
      eq_filler_inc_layer_index (ef);
    }

  ef->mx = 2 * ef->ws->MX;
  /*body*/
  eq_filler_inc_layer_index (ef);

  for (ef->mx = 2 * ef->ws->MX + 1; ef->mx < mx_max; ef->mx+)
    {
      /*body*/
      eq_filler_inc_layer_index (ef);
    }

  ef->mx = mx_max;
  /*body*/
  eq_filler_inc_layer_index (ef);

  for (ef->my = ef->ws->MY + 1; ef->my < 2 * ef->ws->MY; ef->my++)
    {
      ef->mx = ef->ws->MX;
      /*body*/
      eq_filler_inc_layer_index (ef);

      for (ef->mx = ef->ws->MX + 1; ef->mx < 2 * ef->ws->MX; ef->mx++)
        {
          /*body*/
          eq_filler_inc_layer_index (ef);
        }

      ef->mx = 2 * ef->ws->MX;
      /*body*/
      eq_filler_inc_layer_index (ef);
    }

  ef->my = 2 * ef->ws->MY;

  ef->mx = ef->ws->MX;
  /*body*/
  eq_filler_inc_layer_index (ef);

  for (ef->mx = ef->ws->MX + 1; ef->mx < 2 * ef->ws->MX; ef->mx++)
    {
      /*body*/
      eq_filler_inc_layer_index (ef);
    }

  ef->mx = 2 * ef->ws->MX;
  /*body*/
  eq_filler_inc_layer_index (ef);


  ef->mx = ef->ws->MX;
  /*body*/
  eq_filler_inc_layer_index (ef);
}

void cdiff_solver_solve_system (central_diff_solver *solver)
{
  FIX_UNUSED (solver);
}

#include "central_diff_solver_private.h"

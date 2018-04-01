#include "sokolov_solver_private.h"
#include "input/rhs.h"
#include "common/debug_utils.h"
#include "common/math_utils.h"

int sokolov_solver_init (sokolov_solver *solver,
                         mesh_info_t mesh_info,
                         solver_mode_t mode,
                         pressure_func_t p_func,
                         double solver_prec,
                         int solver_max_iter,
                         preconditioner_t precond,
                         linear_solver_t linear_solver)
{
  solver->mesh_info = mesh_info;
  solver->mode = mode;

  solver->h = &solver->h_obj;
  solver->vx = &solver->vx_obj;
  solver->vy = &solver->vy_obj;
  solver->composer_h = &solver->composer_h_obj;
  solver->composer_v = &solver->composer_v_obj;

  solver->p_drv_type = p_func;

  switch (p_func)
    {
    case pressure_linear:
      solver->p_drv = p_drv_linear;
      break;
    case pressure_polynomial:
      solver->p_drv = p_drv_polynomial;
      break;
    }

  hn_values_init (solver->h, mesh_info.MX, mesh_info.MY, mesh_info.N);
  nodes_values_init (solver->vx, mesh_info.MX, mesh_info.MY, mesh_info.N);
  nodes_values_init (solver->vy, mesh_info.MX, mesh_info.MY, mesh_info.N);

  system_composer_init (solver->composer_h, linear_solver, solver->h->layer_size, 1, solver_prec,
                        solver_max_iter, precond);

  system_composer_init (solver->composer_v, linear_solver, solver->vx->layer_size, 2,
                        solver_prec, solver_max_iter, precond);

  return 1;
}

void sokolov_solver_destroy (sokolov_solver *solver)
{
  hn_values_destroy (solver->h);
  nodes_values_destroy (solver->vx);
  nodes_values_destroy (solver->vy);
  system_composer_destroy (solver->composer_h);
  system_composer_destroy (solver->composer_v);
}

int sokolov_solver_compute (sokolov_solver *solver,
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
  solver->test_solution_g = test_solution_g;
  solver->test_solution_vx = test_solution_vx;
  solver->test_solution_vy = test_solution_vy;

  solver->f0 = f0;
  solver->f1 = f1;
  solver->f2 = f2;

  solver->start_g = start_g;
  solver->start_vx = start_vx;
  solver->start_vy = start_vy;

  sokolov_solver_fill_zero_layer (solver);
  sokolov_solver_fill_borders (solver);

  for (solver->layer = 1; solver->layer <= solver->mesh_info.N; solver->layer++)
    {
      sokolov_solver_fill_h_matrix_w_rhs (solver);
      sokolov_solver_solve_h_system (solver);
      sokolov_solver_fill_h_layer (solver);

      sokolov_solver_fill_v_matrix_w_rhs (solver);
      sokolov_solver_solve_v_system (solver);
      sokolov_solver_fill_v_layer (solver);
    }

  return 0;
}

void sokolov_solver_fill_zero_layer (sokolov_solver *solver)
{
  int mx;
  int my;
  double x;
  double y;
  double border_omega = solver->mesh_info.border_omega;
  int i = 0;

  for (i = 0; i < solver->h->layer_size; i++)
    {
      hn_values_get_mx_my (solver->h, i, &mx, &my);

      x = solver->mesh_info.hx * ((double)mx + 0.5);
      y = solver->mesh_info.hy * ((double)my + 0.5);

      solver->h->vals[i] = solver->start_g (x, y, border_omega);
    }

  for (i = 0; i < solver->vx->layer_size; i++)
    {
      nodes_values_get_mx_my (solver->vx, i, &mx, &my);

      x = solver->mesh_info.hx * mx;
      y = solver->mesh_info.hy * my;

      solver->vx->vals[i] = solver->start_vx (x, y, border_omega);
      solver->vy->vals[i] = solver->start_vy (x, y, border_omega);
    }
}

void sokolov_solver_fill_borders (sokolov_solver *solver)
{
  int layer = 1;
  for (layer = 1; layer <= solver->mesh_info.N; layer++)
    sokolov_solver_fill_borders_on_layer (solver, layer);
}

void sokolov_solver_fill_h_matrix_w_rhs (sokolov_solver *solver)
{
  int lli = 0;
  double rhs = 0;
  double coef = 0;

  int n = solver->layer - 1;
  int mx;
  int my;

  double t = solver->mesh_info.tau * n;
  double x;
  double y;

  int h_cur;
  int h_left;
  int h_top;
  int h_right;
  int h_bot;

  nz_row_t nz_row_obj;
  nz_row_t *nz_row = &nz_row_obj;

  nz_row_init (nz_row, 5);

  for (lli = 0; lli < solver->h->layer_size; lli++)
    {
      nz_row->nnz = 0;
      nz_row->row = lli;

      h_cur = lli;

      hn_values_get_mx_my (solver->h, lli, &mx, &my);

      x = solver->mesh_info.hx * mx + solver->mesh_info.hx * 0.5;
      y = solver->mesh_info.hy * my + solver->mesh_info.hy * 0.5;

      if (hn_values_is_border (solver->h, lli))
        {
          sparse_base_fill_nz_s (nz_row, 1, h_cur);
          rhs = hn_values_mx_my_val (solver->h, solver->layer, mx, my);
        }
      else
        {
          double gx = solver->mesh_info.tau / solver->mesh_info.hx;
          double gy = solver->mesh_info.tau / solver->mesh_info.hy;

          double val1 = nodes_avg_fwd_y (solver->vx, n, mx + 1, my);

          double val2 = nodes_avg_fwd_y (solver->vx, n, mx, my);

          double val3 = nodes_avg_fwd_x (solver->vy, n, mx, my + 1);

          double val4 = nodes_avg_fwd_x (solver->vy, n, mx, my);

          h_left = hn_values_index (solver->h, 0, mx - 1, my);
          h_top = hn_values_index (solver->h, 0, mx, my + 1);
          h_right = hn_values_index (solver->h, 0, mx + 1, my);
          h_bot = hn_values_index (solver->h, 0, mx, my - 1);

          /* mx my */
          coef =
              + 1
              + (fabs (val1) + val1) / (2 * fabs (val1)) * gx * val1
              - (fabs (val2) - val2) / (2 * fabs (val2)) * gx * val2
              + (fabs (val3) + val3) / (2 * fabs (val3)) * gy * val3
              - (fabs (val4) - val4) / (2 * fabs (val4)) * gy * val4;

          sparse_base_fill_nz_s (nz_row, coef, h_cur);

          /* mx - 1 my */
          coef =
              - (fabs (val2) + val2) / (2 * fabs (val2)) * gx * val2;

          sparse_base_fill_nz_s (nz_row, coef, h_left);

          /* mx my + 1 */
          coef =
              + (fabs (val3) - val3) / (2 * fabs (val3)) * gy * val3;

          sparse_base_fill_nz_s (nz_row, coef, h_top);

          /* mx + 1 my */
          coef =
              + (fabs (val1) - val1) / (2 * fabs (val1)) * gx * val1;

          sparse_base_fill_nz_s (nz_row, coef, h_right);

          /* mx my - 1 */
          coef =
              - (fabs (val4) + val4) / (2 * fabs (val4)) * gy * val4;

          sparse_base_fill_nz_s (nz_row, coef, h_bot);

          rhs = hn_values_mx_my_val (solver->h, n, mx, my) + solver->mesh_info.tau * solver->f0 (t, x, y, solver->mesh_info.mu,
                                                                                                 solver->p_drv_type);
        }

      sparse_base_add_row_s (solver->composer_h->matrix_base, nz_row);
      system_composer_set_rhs_val (solver->composer_h, rhs, nz_row->row);
    }

  nz_row_destroy (nz_row);
}

void sokolov_solver_solve_h_system (sokolov_solver *solver)
{
  system_composer_solve (solver->composer_h);
}

void sokolov_solver_fill_h_layer (sokolov_solver *solver)
{
  vector_double_t out[1];
  out[0] = solver->h->vals + solver->h->layer_size * solver->layer;
  system_composer_fill_nodes_values (solver->composer_h, out);
}

void sokolov_solver_fill_v_matrix_w_rhs (sokolov_solver *solver)
{
  FIX_UNUSED (solver);
}

void sokolov_solver_solve_v_system (sokolov_solver *solver)
{
  FIX_UNUSED (solver);
}

void sokolov_solver_fill_v_layer (sokolov_solver *solver)
{
  FIX_UNUSED (solver);
}

void sokolov_solver_fill_borders_on_layer (sokolov_solver *solver, int layer)
{
  hn_border_iter hn_iter;
  nodes_border_iter n_iter;
  if (solver->mode == solve_mode)
    {
      /*
       * H is filled with a constant RHO_LEFTMOST
       */

      int mx;
      int my;
      int index;


      /*border_leftmost*/
      mx = 0;
      for (my = 0; my <= solver->mesh_info.MY; my++)
        {
          index = nodes_values_index (solver->vx, layer, mx, my);
          solver->vx->vals[index] = solver->mesh_info.border_omega;
        }

      hn_border_iter_init (&hn_iter, solver->h);

      while (!hn_iter.is_end)
        {
          int index = hn_values_index (solver->h, layer, hn_iter.mx, hn_iter.my);
          solver->h->vals[index] = RHO_LEFTMOST;
          hn_border_iter_next (&hn_iter);
        }

      /*
       * Other borders are
       * already zero
       */
    }
  else
    {
      nodes_border_iter_init (&n_iter, solver->vx);

      while (!n_iter.is_end)
        {
          sokolov_solver_fill_v_values_from_functions (solver, layer, n_iter.mx, n_iter.my);
          nodes_border_iter_next (&n_iter);
        }


      hn_border_iter_init (&hn_iter, solver->h);

      while (!hn_iter.is_end)
        {
          sokolov_solver_fill_h_values_from_functions (solver, layer, hn_iter.mx, hn_iter.my);
          hn_border_iter_next (&hn_iter);
        }
    }
}

void sokolov_solver_fill_v_values_from_functions (sokolov_solver *solver, int layer, int mx, int my)
{
  double x = solver->mesh_info.hx * mx;
  double y = solver->mesh_info.hy * my;
  double t = solver->mesh_info.tau * layer;
  int index = nodes_values_index (solver->vx, layer, mx, my);
  solver->vx->vals[index] = solver->test_solution_vx (t, x, y);
  solver->vy->vals[index] = solver->test_solution_vy (t, x, y);
}

void sokolov_solver_fill_h_values_from_functions (sokolov_solver *solver, int layer, int mx, int my)
{
  double x = solver->mesh_info.hx * mx + 0.5 * solver->mesh_info.hx;
  double y = solver->mesh_info.hy * my + 0.5 * solver->mesh_info.hy;
  double t = solver->mesh_info.tau * layer;
  int index = hn_values_index (solver->h, layer, mx, my);
  solver->h->vals[index] = solver->test_solution_g (t, x, y);
}

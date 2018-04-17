#include "sokolov_solver_private.h"
#include "input/rhs.h"
#include "common/debug_utils.h"
#include "common/math_utils.h"
#include "common/test_macro/tests.h"

/* TEMPCODE_BEGIN */
#include "input/test_solutions.h"
/* TEMPCODE_END */

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

  solver->gamma = 1.4;

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
                            time_layer_func_t test_solution_h,
                            time_layer_func_t test_solution_vx,
                            time_layer_func_t test_solution_vy,
                            rhs_func_t f0,
                            rhs_func_t f1,
                            rhs_func_t f2,
                            layer_func_t start_vx,
                            layer_func_t start_vy,
                            layer_func_t start_h)
{
  solver->test_solution_h = test_solution_h;
  solver->test_solution_vx = test_solution_vx;
  solver->test_solution_vy = test_solution_vy;

  solver->f0 = f0;
  solver->f1 = f1;
  solver->f2 = f2;

  solver->start_h = start_h;
  solver->start_vx = start_vx;
  solver->start_vy = start_vy;

  sokolov_solver_fill_zero_layer (solver);
  sokolov_solver_fill_borders (solver);

  for (solver->layer = 1; solver->layer <= solver->mesh_info.N; solver->layer++)
    {
      fprintf (stdout, "Layer: %d, MX = %d, MY = %d, N = %d\n", solver->layer, solver->mesh_info.MX, solver->mesh_info.MY, solver->mesh_info.N);
      fprintf (stdout, "Filling H matrix...\n");

      sokolov_solver_fill_h_matrix_w_rhs (solver);

      fprintf (stdout, "Solving H system...\n");

      sokolov_solver_solve_h_system (solver);
      sokolov_solver_fill_h_layer (solver);

      fprintf (stdout, "Filling V matrix...\n");
      sokolov_solver_fill_v_matrix_w_rhs (solver);

      fprintf (stdout, "Solving V system...\n");
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
      sokolov_solver_set_hn_x_y (solver, mx, my, &x, &y);

      solver->h->vals[i] = solver->start_h (x, y, border_omega);
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
  double coef_cur;
  double coef_left;
  double coef_top;
  double coef_right;
  double coef_bot;

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

      sokolov_solver_set_hn_x_y (solver, mx, my, &x, &y);

      if (hn_values_is_border (solver->h, lli))
        {
          sparse_base_fill_nz_s (nz_row, 1, h_cur);
          if (solver->mode == test_mode || mx == 0)
            rhs = hn_values_mx_my_val (solver->h, solver->layer, mx, my);
          else
            {
              rhs = 0;
            }
        }
      else
        {
          double gx = solver->mesh_info.tau / solver->mesh_info.hx;
          double gy = solver->mesh_info.tau / solver->mesh_info.hy;

          double vx_r = nodes_avg_fwd_y (solver->vx, n, mx + 1, my);

          double vx_l = nodes_avg_fwd_y (solver->vx, n, mx, my);

          double vy_t = nodes_avg_fwd_x (solver->vy, n, mx, my + 1);

          double vy_b = nodes_avg_fwd_x (solver->vy, n, mx, my);

          h_left  = hn_values_index (solver->h, 0, mx - 1, my);
          h_top   = hn_values_index (solver->h, 0, mx, my + 1);
          h_right = hn_values_index (solver->h, 0, mx + 1, my);
          h_bot   = hn_values_index (solver->h, 0, mx, my - 1);

          /* mx my */
          coef_cur =
              + 1
              + math_is_pos_scaled (vx_r) * gx
              - math_is_neg_scaled (vx_l) * gx
              + math_is_pos_scaled (vy_t) * gy
              - math_is_neg_scaled (vy_b) * gy;

          DEBUG_ASSERT (!math_is_null (coef_cur));
          sparse_base_fill_nz_s (nz_row, coef_cur, h_cur);

          /* mx - 1 my */
          coef_left = - math_is_pos_scaled (vx_l) * gx;

          sparse_base_fill_nz_s (nz_row, coef_left, h_left);

          /* mx my + 1 */
          coef_top = + math_is_neg_scaled (vy_t) * gy;

          sparse_base_fill_nz_s (nz_row, coef_top, h_top);

          /* mx + 1 my */
          coef_right = + math_is_neg_scaled (vx_r) * gx;

          sparse_base_fill_nz_s (nz_row, coef_right, h_right);

          /* mx my - 1 */
          coef_bot = - math_is_pos_scaled (vy_b) * gy;

          sparse_base_fill_nz_s (nz_row, coef_bot, h_bot);

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
#ifdef DEBUG
  if (solver->mode == test_mode)
    sokolov_solver_fill_x_init_w_real_values (solver);
#endif
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
  int mx;
  int my;
  int n = solver->layer - 1;
  double hx = solver->mesh_info.hx;
  double hy = solver->mesh_info.hy;
  double tau = solver->mesh_info.tau;
  double mu = solver->mesh_info.mu;
  double gamma = solver->gamma;
  double x;
  double y;
  double coef = 0;
  double rhs = 0;
  int lli;

  int vx_c, vx_l, vx_r, vx_t, vx_b,
      vy_c, vy_l, vy_r, vy_t, vy_b;

  double gx, gy;

  double vxv_c, vxv_l, vxv_r, vxv_t, vxv_b,
      vyv_c, vyv_l, vyv_r, vyv_t, vyv_b;

  grid_area_t area;

  for (lli = 0; lli < solver->vx->layer_size; lli++)
    {
      nodes_values_get_mx_my (solver->vx, lli, &mx, &my);
      x = mx * hx;
      y = my * hy;

      area = nodes_values_get_area (solver->vx, mx, my);

      if (area != area_internal && solver->mode == test_mode)
        {
          int cols[1];
          double vals[1];
          int nnz = 1;

          cols[0] = 2 * lli;
          vals[0] = 1;

          sparse_base_add_row (solver->composer_v->matrix_base, 2 * lli, cols, vals, nnz);
          system_composer_set_rhs_val (solver->composer_v, solver->test_solution_vx ((n + 1) * tau, x, y), 2 * lli);

          cols[0] = 2 * lli + 1;

          sparse_base_add_row (solver->composer_v->matrix_base, 2 * lli + 1, cols, vals, nnz);
          system_composer_set_rhs_val (solver->composer_v, solver->test_solution_vy ((n + 1) * tau, x, y), 2 * lli + 1);
          continue;
        }

      if (area == border_leftmost)
        {
          int cols[1];
          double vals[1];
          int nnz = 1;

          cols[0] = 2 * lli;
          vals[0] = 1;

          sparse_base_add_row (solver->composer_v->matrix_base, 2 * lli, cols, vals, nnz);
          system_composer_set_rhs_val (solver->composer_v, solver->mesh_info.border_omega, 2 * lli);

          cols[0] = 2 * lli + 1;

          sparse_base_add_row (solver->composer_v->matrix_base, 2 * lli + 1, cols, vals, nnz);
          system_composer_set_rhs_val (solver->composer_v, 0, 2 * lli + 1);
          continue;
        }
      if (area == border_rightmost)
        {
          int cols[2];
          double vals[2];
          int nnz = 2;

          int vx_c = 2 * nodes_values_index (solver->vx, 0, mx, my);
          int vx_l = 2 * nodes_values_index (solver->vx, 0, mx - 1, my);

          cols[0] = vx_l;
          cols[1] = vx_c;

          vals[0] = -1;
          vals[1] = 1;

          sparse_base_add_row (solver->composer_v->matrix_base, 2 * lli, cols, vals, nnz);
          system_composer_set_rhs_val (solver->composer_v, 0, 2 * lli);

          nnz = 1;
          cols[0] = 2 * lli + 1;
          vals[0] = 1;

          sparse_base_add_row (solver->composer_v->matrix_base, 2 * lli + 1, cols, vals, nnz);
          system_composer_set_rhs_val (solver->composer_v, 0, 2 * lli + 1);
          continue;
        }
      if (area == border_topmost)
        {
          int cols[2];
          double vals[2];
          int nnz = 2;

          int vy_c = 2 * nodes_values_index (solver->vy, 0, mx, my) + 1;
          int vy_b = 2 * nodes_values_index (solver->vy, 0, mx, my - 1) + 1;

          nnz = 1;
          cols[0] = 2 * lli;
          vals[0] = 1;

          sparse_base_add_row (solver->composer_v->matrix_base, 2 * lli, cols, vals, nnz);
          system_composer_set_rhs_val (solver->composer_v, 0, 2 * lli);

          nnz = 2;

          cols[0] = vy_b;
          cols[1] = vy_c;

          vals[0] = -1;
          vals[1] = 1;

          sparse_base_add_row (solver->composer_v->matrix_base, 2 * lli + 1, cols, vals, nnz);
          system_composer_set_rhs_val (solver->composer_v, 0, 2 * lli + 1);


          continue;
        }
      if (area != area_internal)
        {
          int cols[1];
          double vals[1];
          int nnz = 1;

          cols[0] = 2 * lli;
          vals[0] = 1;

          sparse_base_add_row (solver->composer_v->matrix_base, 2 * lli, cols, vals, nnz);
          system_composer_set_rhs_val (solver->composer_v, 0, 2 * lli);

          cols[0] = 2 * lli + 1;

          sparse_base_add_row (solver->composer_v->matrix_base, 2 * lli + 1, cols, vals, nnz);
          system_composer_set_rhs_val (solver->composer_v, 0, 2 * lli + 1);
          continue;
        }

      vx_c = 2 * lli;
      vx_l = 2 * nodes_values_index (solver->vx, 0, mx - 1, my);
      vx_r = 2 * nodes_values_index (solver->vx, 0, mx + 1, my);
      vx_t = 2 * nodes_values_index (solver->vx, 0, mx, my + 1);
      vx_b = 2 * nodes_values_index (solver->vx, 0, mx, my - 1);

      vy_c = vx_c + 1;
      vy_l = vx_l + 1;
      vy_r = vx_r + 1;
      vy_t = vx_t + 1;
      vy_b = vx_b + 1;

      gx = tau / hx;
      gy = tau / hy;

      vxv_c = nodes_values_mx_my_val (solver->vx, n, mx, my);
      vxv_l = nodes_values_mx_my_val (solver->vx, n, mx - 1, my);
      vxv_r = nodes_values_mx_my_val (solver->vx, n, mx + 1, my);
      vxv_b = nodes_values_mx_my_val (solver->vx, n, mx, my - 1);
      vxv_t = nodes_values_mx_my_val (solver->vx, n, mx, my + 1);

      vyv_c = nodes_values_mx_my_val (solver->vy, n, mx, my);
      vyv_l = nodes_values_mx_my_val (solver->vy, n, mx - 1, my);
      vyv_r = nodes_values_mx_my_val (solver->vy, n, mx + 1, my);
      vyv_b = nodes_values_mx_my_val (solver->vy, n, mx, my - 1);
      vyv_t = nodes_values_mx_my_val (solver->vy, n, mx, my + 1);

      {
        if (math_is_null (hn_values_approx_in_node (solver->h, n + 1, mx, my)))
          {
            int cols[] = {vx_c};
            double vals[] = {1};
            int nnz = 1;

            sparse_base_add_row (solver->composer_v->matrix_base, 2 * lli, cols, vals, nnz);
            system_composer_set_rhs_val (solver->composer_v, 0, 2 * lli);
          }
        else
          {
            nz_row_t nz_row_obj;
            nz_row_t *nz_row = &nz_row_obj;

            nz_row_init (nz_row, 8);
            nz_row->row = 2 * lli;

            coef = (gx / 2.)
                   * (
                     - hn_values_avg_bwd_y (solver->h, n + 1, mx - 1, my)
                     * math_is_neg_scaled (vxv_c)
                     - hn_values_avg_bwd_y (solver->h, n + 1, mx - 2, my)
                     * math_is_neg_scaled (vxv_l)
                     )
                   - (4 * mu * tau) / (3 * hx * hx);

            sparse_base_fill_nz_s (nz_row, coef, vx_l);

            coef = + hn_values_approx_in_node (solver->h, n + 1, mx, my)
                   + (gx / 2.)
                   * (
                     + hn_values_avg_bwd_y (solver->h, n + 1, mx, my)
                     * math_is_neg_scaled (vxv_r)
                     - hn_values_avg_bwd_y (solver->h, n + 1, mx, my)
                     * math_is_pos_scaled (vxv_c)
                     + hn_values_avg_bwd_y (solver->h, n + 1, mx - 1, my)
                     * math_is_neg_scaled (vxv_c)
                     - hn_values_avg_bwd_y (solver->h, n + 1, mx - 1, my)
                     * math_is_pos_scaled (vxv_l)
                     )
                   + (2 * mu * tau) / (hx * hx)
                   + (8 * mu) / (3 * hy * hy);

            sparse_base_fill_nz_s (nz_row, coef, vx_c);

            coef = (gx / 2.)
                   * (
                     + hn_values_avg_bwd_y (solver->h, n + 1, mx + 1, my)
                     * math_is_pos_scaled (vxv_r)
                     + hn_values_avg_bwd_y (solver->h, n + 1, mx, my)
                     * math_is_pos_scaled (vxv_c)
                     )
                   - (4 * mu * tau) / (3 * hx * hx);

            sparse_base_fill_nz_s (nz_row, coef, vx_r);

            coef = (gy / 2.)
                   * (
                     - hn_values_avg_bwd_x (solver->h, n + 1, mx, my - 1)
                     * math_is_neg_scaled (vxv_c)
                     - hn_values_avg_bwd_x (solver->h, n + 1, mx, my - 2)
                     * math_is_neg_scaled (vxv_b)
                     );

            sparse_base_fill_nz_s (nz_row, coef, vy_b);

            coef = (gy / 2.)
                   * (
                     + hn_values_avg_bwd_x (solver->h, n + 1, mx, my)
                     * math_is_neg_scaled (vxv_t)
                     - hn_values_avg_bwd_x (solver->h, n + 1, mx, my)
                     * math_is_pos_scaled (vxv_c)
                     + hn_values_avg_bwd_x (solver->h, n + 1, mx, my - 1)
                     * math_is_neg_scaled (vxv_c)
                     - hn_values_avg_bwd_x (solver->h, n + 1, mx, my - 1)
                     * math_is_pos_scaled (vxv_b)
                     );

            sparse_base_fill_nz_s (nz_row, coef, vy_c);

            coef = (gy / 2.)
                   * (
                     + hn_values_avg_bwd_x (solver->h, n + 1, mx, my + 1)
                     * math_is_pos_scaled (vxv_t)
                     + hn_values_avg_bwd_y (solver->h, n + 1, mx, my)
                     * math_is_pos_scaled (vxv_c)
                     );

            sparse_base_fill_nz_s (nz_row, coef, vy_t);

            coef = - (mu * tau) / (hy * hy);
            sparse_base_fill_nz_s (nz_row, coef, vx_b);
            sparse_base_fill_nz_s (nz_row, coef, vx_t);

            rhs =
                + hn_values_approx_in_node (solver->h, n, mx, my) * vxv_c
                + gx * (gamma / (1 - gamma)) * hn_values_approx_in_node (solver->h, n + 1, mx, my)
                * (
                  + pow (hn_values_avg_bwd_y (solver->h, n + 1, mx, my), gamma - 1)
                  - pow (hn_values_avg_bwd_y (solver->h, n + 1, mx - 1, my), gamma - 1)
                  )
                + ((mu * tau) / (12 * hx * hy))
                * (
                  + nodes_values_mx_my_val (solver->vy, n, mx - 1, my - 1)
                  - nodes_values_mx_my_val (solver->vy, n, mx - 1, my + 1)
                  - nodes_values_mx_my_val (solver->vy, n, mx + 1,my - 1)
                  + nodes_values_mx_my_val (solver->vy, n, mx + 1, my + 1)
                  )
                + tau * solver->f1 (tau * (n + 1), x, y, mu, solver->p_drv_type) * hn_values_approx_in_node (solver->h, n + 1, mx, my);

            sparse_base_add_row_s (solver->composer_v->matrix_base, nz_row);
            system_composer_set_rhs_val (solver->composer_v, rhs, nz_row->row);
            nz_row_destroy (nz_row);
          }
      }

      {
        if (math_is_null (hn_values_approx_in_node (solver->h, n + 1, mx, my)))
          {
            int cols[] = {vy_c};
            double vals[] = {1};
            int nnz = 1;

            sparse_base_add_row (solver->composer_v->matrix_base, 2 * lli + 1, cols, vals, nnz);
            system_composer_set_rhs_val (solver->composer_v, 0, 2 * lli + 1);
          }
        else
          {
            nz_row_t nz_row_obj;
            nz_row_t *nz_row = &nz_row_obj;

            nz_row_init (nz_row, 8);
            nz_row->row = 2 * lli + 1;

            coef = (gx / 2.)
                   * (
                     - hn_values_avg_bwd_y (solver->h, n + 1, mx - 1, my)
                     * math_is_neg_scaled (vyv_c)
                     - hn_values_avg_bwd_y (solver->h, n + 1, mx - 2, my)
                     * math_is_neg_scaled (vyv_l)
                     );

            sparse_base_fill_nz_s (nz_row, coef, vx_l);

            coef = (gx / 2.)
                   * (
                     + hn_values_avg_bwd_y (solver->h, n + 1, mx, my)
                     * math_is_neg_scaled (vyv_r)
                     - hn_values_avg_bwd_y (solver->h, n + 1, mx, my)
                     * math_is_pos_scaled (vyv_c)
                     + hn_values_avg_bwd_y (solver->h, n + 1, mx - 1, my)
                     * math_is_neg_scaled (vyv_c)
                     - hn_values_avg_bwd_y (solver->h, n + 1, mx - 1, my)
                     * math_is_pos_scaled (vyv_l)
                     );

            sparse_base_fill_nz_s (nz_row, coef, vx_c);

            coef = (gx / 2.)
                   * (
                     + hn_values_avg_bwd_y (solver->h, n+1, mx + 1, my)
                     * math_is_pos_scaled (vyv_r)
                     + hn_values_avg_bwd_y (solver->h, n + 1, mx, my)
                     * math_is_pos_scaled (vyv_c)
                     );

            sparse_base_fill_nz_s (nz_row, coef, vx_r);

            coef = (gy / 2.)
                   * (
                     - hn_values_avg_bwd_x (solver->h, n + 1, mx, my - 1)
                     * math_is_neg_scaled (vyv_c)
                     - hn_values_avg_bwd_x (solver->h, n + 1, mx, my - 2)
                     * math_is_neg_scaled (vyv_b)
                     )
                   - (tau * mu * 4) / (3 * hy * hy);

            sparse_base_fill_nz_s (nz_row, coef, vy_b);

            coef = (gy / 2.)
                   * (
                     + hn_values_avg_bwd_x (solver->h, n + 1, mx, my)
                     * math_is_neg_scaled (vyv_t)
                     - hn_values_avg_bwd_x (solver->h, n + 1, mx, my)
                     * math_is_pos_scaled (vyv_c)
                     + hn_values_avg_bwd_x (solver->h, n + 1, mx, my - 1)
                     * math_is_neg_scaled (vyv_c)
                     - hn_values_avg_bwd_x (solver->h, n + 1, mx, my - 1)
                     * math_is_pos_scaled (vyv_b)
                     )
                   + hn_values_approx_in_node (solver->h, n + 1, mx, my)
                   + (tau * mu * 2) / (hx * hx)
                   + (tau * mu * 8) / (3 * hy * hy);

            sparse_base_fill_nz_s (nz_row, coef, vy_c);

            coef = (gy / 2.)
                   * (
                     + hn_values_avg_bwd_x (solver->h, n + 1, mx, my + 1)
                     * math_is_pos_scaled (vyv_t)
                     + hn_values_avg_bwd_x (solver->h, n + 1, mx, my)
                     * math_is_pos_scaled (vyv_c)
                     )
                   - (tau * mu * 4) / (3 * hy * hy);

            sparse_base_fill_nz_s (nz_row, coef, vy_t);

            coef = - (tau * mu) / (hx * hx);
            sparse_base_fill_nz_s (nz_row, coef, vy_l);
            sparse_base_fill_nz_s (nz_row, coef, vy_r);

            rhs =
                + hn_values_approx_in_node (solver->h, n, mx, my) * vyv_c
                + gy * (gamma / (1 - gamma)) * hn_values_approx_in_node (solver->h, n + 1, mx, my)
                * (
                  + pow (hn_values_avg_bwd_x (solver->h, n + 1, mx, my), gamma - 1)
                  - pow (hn_values_avg_bwd_x (solver->h, n + 1, mx, my - 1), gamma - 1)
                  )
                + ((mu * tau) / (12 * hx * hy))
                * (
                  + nodes_values_mx_my_val (solver->vx, n, mx - 1, my - 1)
                  - nodes_values_mx_my_val (solver->vx, n, mx - 1, my + 1)
                  - nodes_values_mx_my_val (solver->vx, n, mx + 1, my - 1)
                  + nodes_values_mx_my_val (solver->vx, n, mx + 1, my + 1)
                  )
                + tau * solver->f2 (tau * (n + 1), x, y, mu, solver->p_drv_type) * hn_values_approx_in_node (solver->h, n + 1, mx, my);

            sparse_base_add_row_s (solver->composer_v->matrix_base, nz_row);
            system_composer_set_rhs_val (solver->composer_v, rhs, nz_row->row);
            nz_row_destroy (nz_row);
          }
      }
    }
}

void sokolov_solver_solve_v_system (sokolov_solver *solver)
{
  system_composer_solve (solver->composer_v);
}

void sokolov_solver_fill_v_layer (sokolov_solver *solver)
{
  vector_double_t out[2];
  out[0] = solver->vx->vals + solver->vx->layer_size * solver->layer;
  out[1] = solver->vy->vals + solver->vy->layer_size * solver->layer;

  system_composer_fill_nodes_values (solver->composer_v, out);
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
  double x;
  double y;
  double t = solver->mesh_info.tau * layer;
  int index = hn_values_index (solver->h, layer, mx, my);
  sokolov_solver_set_hn_x_y (solver, mx, my, &x, &y);
  solver->h->vals[index] = solver->test_solution_h (t, x, y);
}

void sokolov_solver_set_hn_x_y (const sokolov_solver *solver, int mx, int my, double *x_ptr, double *y_ptr)
{
  double x = mx * solver->mesh_info.hx + 0.5 * solver->mesh_info.hx;
  double y = my * solver->mesh_info.hy + 0.5 * solver->mesh_info.hy;

  if (x_ptr)
    *x_ptr = x;

  if (y_ptr)
    *y_ptr = y;
}

void sokolov_solver_fill_x_init_w_real_values (sokolov_solver *solver)
{
  ASSERT_RETURN_VOID (solver->mode == test_mode);
  int i = 0;
  for (i = 0; i < solver->h->layer_size; i++)
    {
      int mx, my;
      double x, y;
      double t = solver->layer * solver->mesh_info.tau;
      hn_values_get_mx_my (solver->h, i, &mx, &my);
      sokolov_solver_set_hn_x_y (solver, mx, my, &x, &y);

      solver->composer_h->vector_to_compute[i] = solver->test_solution_h (t, x, y);
    }
  for (i = 0; i < solver->vx->layer_size; i++)
    {
      int mx, my;
      double x, y;
      double t = solver->layer * solver->mesh_info.tau;
      nodes_values_get_mx_my (solver->vx, i, &mx, &my);

      x = solver->mesh_info.hx * mx;
      y = solver->mesh_info.hy * my;

      solver->composer_v->vector_to_compute[2 * i] = solver->test_solution_vx (t, x, y);
      solver->composer_v->vector_to_compute[2 * i + 1] = solver->test_solution_vy (t, x, y);
    }
}

#include "solver_tester.h"
#include "common/vectors.h"
#include "linear_ops/vector_ops.h"
#include "common/debug_utils.h"
#include "common/math_utils.h"

#include "io/table_io.h"

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
                         layer_func_t start_g)
{
  int i;

  tester->N_start = N_start;
  tester->MX_start = MX_start;
  tester->MY_start = MY_start;
  tester->N_mult = N_mult;
  tester->MXY_mult = MXY_mult;
  tester->N_mult_count = N_mult_count;
  tester->MXY_mult_count = MXY_mult_count;

  tester->border_omega = border_omega;
  tester->mu = mu;
  tester->X = X;
  tester->Y = Y;
  tester->T = T;

  tester->g_func = g_func;
  tester->vx_func = vx_func;
  tester->vy_func = vy_func;

  tester->f0 = f0;
  tester->f1 = f1;
  tester->f2 = f2;

  tester->start_g = start_g;
  tester->start_vx = start_vx;
  tester->start_vy = start_vy;

  tester->cases_count = (N_mult_count) * (MXY_mult_count);
  tester->vector_size = tester->cases_count * UNKNOWN_FUNCTIONS_COUNT;

  tester->c_norms = VECTOR_CREATE (double, tester->vector_size);
  tester->l2_norms = VECTOR_CREATE (double, tester->vector_size);
  tester->w21_norms = VECTOR_CREATE (double, tester->vector_size);

  tester->cg_norms_text = VECTOR_CREATE (string_t, tester->cases_count);
  tester->l2g_norms_text = VECTOR_CREATE (string_t, tester->cases_count);
  tester->w21g_norms_text = VECTOR_CREATE (string_t, tester->cases_count);
  tester->cvx_norms_text = VECTOR_CREATE (string_t, tester->cases_count);
  tester->l2vx_norms_text = VECTOR_CREATE (string_t, tester->cases_count);
  tester->w21vx_norms_text = VECTOR_CREATE (string_t, tester->cases_count);
  tester->cvy_norms_text = VECTOR_CREATE (string_t, tester->cases_count);
  tester->l2vy_norms_text = VECTOR_CREATE (string_t, tester->cases_count);
  tester->w21vy_norms_text = VECTOR_CREATE (string_t, tester->cases_count);

  for (i = 0; i < tester->cases_count; i++)
    {
      tester->cg_norms_text[i] = VECTOR_CREATE (char, MAX_DOUBLE_LENGTH);
      tester->l2g_norms_text[i] = VECTOR_CREATE (char, MAX_DOUBLE_LENGTH);
      tester->w21g_norms_text[i] = VECTOR_CREATE (char, MAX_DOUBLE_LENGTH);
      tester->cvx_norms_text[i] = VECTOR_CREATE (char, MAX_DOUBLE_LENGTH);
      tester->l2vx_norms_text[i] = VECTOR_CREATE (char, MAX_DOUBLE_LENGTH);
      tester->w21vx_norms_text[i] = VECTOR_CREATE (char, MAX_DOUBLE_LENGTH);
      tester->cvy_norms_text[i] = VECTOR_CREATE (char, MAX_DOUBLE_LENGTH);
      tester->l2vy_norms_text[i] = VECTOR_CREATE (char, MAX_DOUBLE_LENGTH);
      tester->w21vy_norms_text[i] = VECTOR_CREATE (char, MAX_DOUBLE_LENGTH);
    }
}

void solver_tester_destroy (solver_tester *tester)
{
  int i;

  VECTOR_DESTROY (tester->c_norms);
  VECTOR_DESTROY (tester->l2_norms);
  VECTOR_DESTROY (tester->w21_norms);

  for (i = 0; i < tester->cases_count; i++)
    {
      VECTOR_DESTROY (tester->cg_norms_text[i]);
      VECTOR_DESTROY (tester->l2g_norms_text[i]);
      VECTOR_DESTROY (tester->w21g_norms_text[i]);
      VECTOR_DESTROY (tester->cvx_norms_text[i]);
      VECTOR_DESTROY (tester->l2vx_norms_text[i]);
      VECTOR_DESTROY (tester->w21vx_norms_text[i]);
      VECTOR_DESTROY (tester->cvy_norms_text[i]);
      VECTOR_DESTROY (tester->l2vy_norms_text[i]);
      VECTOR_DESTROY (tester->w21vy_norms_text[i]);
    }

  VECTOR_DESTROY (tester->cg_norms_text);
  VECTOR_DESTROY (tester->l2g_norms_text);
  VECTOR_DESTROY (tester->w21g_norms_text);
  VECTOR_DESTROY (tester->cvx_norms_text);
  VECTOR_DESTROY (tester->l2vx_norms_text);
  VECTOR_DESTROY (tester->w21vx_norms_text);
  VECTOR_DESTROY (tester->cvy_norms_text);
  VECTOR_DESTROY (tester->l2vy_norms_text);
  VECTOR_DESTROY (tester->w21vy_norms_text);
}

void solver_tester_test (solver_tester *tester,
                         double solver_prec,
                         int solver_max_iter,
                         preconditioner_t precond)
{
  int n_mult;
  int mxy_mult;
  int ind = 0;
  for (n_mult = 0; n_mult < tester->N_mult_count; n_mult++)
    {
      for (mxy_mult = 0; mxy_mult < tester->MXY_mult_count; mxy_mult++)
        {
          int N = tester->N_start;
          int MX = tester->MX_start;
          int MY = tester->MY_start;
          int i;

          for (i = 0; i < n_mult; i++)
            N *= tester->N_mult;

          for (i = 0; i < mxy_mult; i++)
            {
              MX *= tester->MXY_mult;
              MY *= tester->MXY_mult;
            }

          cdiff_solver_init (&tester->cds, test_mode, MX, MY, N,
                             tester->X, tester->Y, tester->T,
                             tester->border_omega, tester->mu,
                             solver_prec, solver_max_iter, precond,
                             custom_cgs);

          cdiff_solver_compute (&tester->cds, pressure_linear,
                                tester->g_func, tester->vx_func, tester->vy_func,
                                tester->f0, tester->f1, tester->f2,
                                tester->start_vx, tester->start_vy, tester->start_g);

          tester->c_norms[ind] = tester_grid_dif_c_norm (tester, grid_g);
          tester->c_norms[ind + 1] = tester_grid_dif_c_norm (tester, grid_vx);
          tester->c_norms[ind + 2] = tester_grid_dif_c_norm (tester, grid_vy);
          tester->l2_norms[ind] = tester_grid_dif_l2_norm (tester, grid_g);
          tester->l2_norms[ind + 1] = tester_grid_dif_l2_norm (tester, grid_vx);
          tester->l2_norms[ind + 2] = tester_grid_dif_l2_norm (tester, grid_vy);
          tester->w21_norms[ind] = tester_grid_dif_w21_norm (tester, grid_g);
          tester->w21_norms[ind + 1] = tester_grid_dif_w21_norm (tester, grid_vx);
          tester->w21_norms[ind + 2] = tester_grid_dif_w21_norm (tester, grid_vy);

          sprintf (tester->cg_norms_text[ind / 3], "%f", tester->c_norms[ind]);
          sprintf (tester->cvx_norms_text[ind / 3], "%f", tester->c_norms[ind + 1]);
          sprintf (tester->cvy_norms_text[ind / 3], "%f", tester->c_norms[ind + 2]);
          sprintf (tester->l2g_norms_text[ind / 3], "%f", tester->l2_norms[ind]);
          sprintf (tester->l2vx_norms_text[ind / 3], "%f", tester->l2_norms[ind + 1]);
          sprintf (tester->l2vy_norms_text[ind / 3], "%f", tester->l2_norms[ind + 2]);
          sprintf (tester->w21g_norms_text[ind / 3], "%f", tester->w21_norms[ind]);
          sprintf (tester->w21vx_norms_text[ind / 3], "%f", tester->w21_norms[ind + 1]);
          sprintf (tester->w21vy_norms_text[ind / 3], "%f", tester->w21_norms[ind + 2]);

          cdiff_solver_destroy (&tester->cds);
          ind += UNKNOWN_FUNCTIONS_COUNT;
        }
    }
}

double tester_grid_dif_c_norm (const solver_tester *tester, grid_func_t f)
{
  double max = 0;
  double val;
  int i;
  int n = tester->cds.ws.N;
  int mx;
  int my;
  for (i = 0; i < tester->cds.ws.layer_size; i++)
    {
      solver_workspace_get_mx_my (&tester->cds.ws, i, &mx, &my);
      val = fabs (solver_workspace_grid_val (&tester->cds.ws, n, mx, my, f) - tester_grid_true_val (tester, f, i));
      max = (max < val) ? val : max;
    }
  return max;
}

double tester_grid_true_val (const solver_tester *tester, grid_func_t f, int loc_layer_index)
{
  int mx;
  int my;
  double t;
  double x;
  double y;

  solver_workspace_get_mx_my (&tester->cds.ws, loc_layer_index, &mx, &my);

  t = tester->cds.ws.T;
  x = tester->cds.ws.hx * mx;
  y = tester->cds.ws.hy * my;

  switch (f)
    {
    case grid_g:
      return tester->g_func (t, x, y);
    case grid_vx:
      return tester->vx_func (t, x ,y);
    case grid_vy:
      return tester->vy_func (t, x, y);
    }

  ASSERT_RETURN (0, 0);
}

double tester_grid_dif_l2_norm (const solver_tester *tester, grid_func_t f)
{
  double coef = tester->cds.ws.hx * tester->cds.ws.hy;
  double sum = 0;
  int i;
  int n = tester->cds.ws.N;
  int mx;
  int my;

  for (i = 0; i < tester->cds.ws.layer_size; i++)
    {
      double val;

      solver_workspace_get_mx_my (&tester->cds.ws, i, &mx, &my);

      val = fabs (solver_workspace_grid_val (&tester->cds.ws, n, mx, my, f) - tester_grid_true_val (tester, f, i));

      if (solver_workspace_get_area (&tester->cds.ws, mx, my) == area_internal)
        sum += val * val;
      else
        sum += 0.5 * val * val;
    }

  return sqrt (coef * sum);
}

double tester_grid_dif_w21_norm (const solver_tester *tester, grid_func_t f)
{
  double coef = tester->cds.ws.hx * tester->cds.ws.hy;
  double sum = 0;
  int i;
  int n = tester->cds.ws.N;
  int mx;
  int my;
  grid_area_t area;

  for (i = 0; i < tester->cds.ws.layer_size; i++)
    {
      double val_right;
      double val_cur;
      double val_top;

      solver_workspace_get_mx_my (&tester->cds.ws, i, &mx, &my);

      val_cur = fabs (solver_workspace_grid_val (&tester->cds.ws, n, mx, my, f) - tester_grid_true_val (tester, f, i));

      area = solver_workspace_get_area (&tester->cds.ws, mx, my);

      if (area != border_right_top
          && area != border_rightmost
          && !(mx == 2 * tester->cds.ws.MX && my == 2 * tester->cds.ws.MY))
        {
          val_right = fabs (solver_workspace_grid_val (&tester->cds.ws, n, mx + 1, my, f) - tester_grid_true_val (tester, f, i + 1));
          sum += (val_right - val_cur) * (val_right - val_cur);
        }

      if (area != border_topmost
          && area != border_left_hor
          && area != border_right_hor
          && !(mx == 0 && my == tester->cds.ws.MY)
          && !(mx == BOT_ROW_SQUARES_COUNT * tester->cds.ws.MX && my == tester->cds.ws.MY))
        {
          int next_index = solver_workspace_final_index (&tester->cds.ws, 0, mx, my + 1);
          val_top = fabs (solver_workspace_grid_val (&tester->cds.ws, n, mx, my + 1, f) - tester_grid_true_val (tester, f, next_index));
          sum += (val_top - val_cur) * (val_top - val_cur);
        }

      if (solver_workspace_get_area (&tester->cds.ws, mx, my) == area_internal)
        sum += val_cur * val_cur;
      else
        sum += 0.5 * val_cur * val_cur;
    }

  return sqrt (coef * sum);
}

void solver_tester_print_results (const solver_tester *tester, FILE *fout)
{
  table_io table_obj;
  table_io *table = &table_obj;

  fprintf (fout, "CG norms:\n");
  table_io_init (table, tester->N_mult_count, tester->MXY_mult_count, tester->cg_norms_text, human_readable);
  fprintf (fout, table->table_text);
  fprintf (fout, "\n");
  table_io_destroy (table);

  table_io_init (table, tester->N_mult_count, tester->MXY_mult_count, tester->cg_norms_text, latex_format);
  fprintf (fout, table->table_text);
  fprintf (fout, "\n");
  table_io_destroy (table);
  fprintf (fout, "\n");

  fprintf (fout, "CVX norms:\n");
  table_io_init (table, tester->N_mult_count, tester->MXY_mult_count, tester->cvx_norms_text, human_readable);
  fprintf (fout, table->table_text);
  fprintf (fout, "\n");
  table_io_destroy (table);

  table_io_init (table, tester->N_mult_count, tester->MXY_mult_count, tester->cvx_norms_text, latex_format);
  fprintf (fout, table->table_text);
  fprintf (fout, "\n");
  table_io_destroy (table);
  fprintf (fout, "\n");

  fprintf (fout, "CVY norms:\n");
  table_io_init (table, tester->N_mult_count, tester->MXY_mult_count, tester->cvy_norms_text, human_readable);
  fprintf (fout, table->table_text);
  fprintf (fout, "\n");
  table_io_destroy (table);

  table_io_init (table, tester->N_mult_count, tester->MXY_mult_count, tester->cvy_norms_text, latex_format);
  fprintf (fout, table->table_text);
  fprintf (fout, "\n");
  table_io_destroy (table);
  fprintf (fout, "\n");


  fprintf (fout, "L2G norms:\n");
  table_io_init (table, tester->N_mult_count, tester->MXY_mult_count, tester->l2g_norms_text, human_readable);
  fprintf (fout, table->table_text);
  fprintf (fout, "\n");
  table_io_destroy (table);

  table_io_init (table, tester->N_mult_count, tester->MXY_mult_count, tester->l2g_norms_text, latex_format);
  fprintf (fout, table->table_text);
  fprintf (fout, "\n");
  table_io_destroy (table);
  fprintf (fout, "\n");

  fprintf (fout, "L2VX norms:\n");
  table_io_init (table, tester->N_mult_count, tester->MXY_mult_count, tester->l2vx_norms_text, human_readable);
  fprintf (fout, table->table_text);
  fprintf (fout, "\n");
  table_io_destroy (table);

  table_io_init (table, tester->N_mult_count, tester->MXY_mult_count, tester->l2vx_norms_text, latex_format);
  fprintf (fout, table->table_text);
  fprintf (fout, "\n");
  table_io_destroy (table);
  fprintf (fout, "\n");

  fprintf (fout, "L2VY norms:\n");
  table_io_init (table, tester->N_mult_count, tester->MXY_mult_count, tester->l2vy_norms_text, human_readable);
  fprintf (fout, table->table_text);
  fprintf (fout, "\n");
  table_io_destroy (table);

  table_io_init (table, tester->N_mult_count, tester->MXY_mult_count, tester->l2vy_norms_text, latex_format);
  fprintf (fout, table->table_text);
  fprintf (fout, "\n");
  table_io_destroy (table);
  fprintf (fout, "\n");

  fprintf (fout, "W21G norms:\n");
  table_io_init (table, tester->N_mult_count, tester->MXY_mult_count, tester->w21g_norms_text, human_readable);
  fprintf (fout, table->table_text);
  fprintf (fout, "\n");
  table_io_destroy (table);

  table_io_init (table, tester->N_mult_count, tester->MXY_mult_count, tester->w21g_norms_text, latex_format);
  fprintf (fout, table->table_text);
  fprintf (fout, "\n");
  table_io_destroy (table);
  fprintf (fout, "\n");

  fprintf (fout, "W21VX norms:\n");
  table_io_init (table, tester->N_mult_count, tester->MXY_mult_count, tester->w21vx_norms_text, human_readable);
  fprintf (fout, table->table_text);
  fprintf (fout, "\n");
  table_io_destroy (table);

  table_io_init (table, tester->N_mult_count, tester->MXY_mult_count, tester->w21vx_norms_text, latex_format);
  fprintf (fout, table->table_text);
  fprintf (fout, "\n");
  table_io_destroy (table);
  fprintf (fout, "\n");

  fprintf (fout, "W21VY norms:\n");
  table_io_init (table, tester->N_mult_count, tester->MXY_mult_count, tester->w21vy_norms_text, human_readable);
  fprintf (fout, table->table_text);
  fprintf (fout, "\n");
  table_io_destroy (table);

  table_io_init (table, tester->N_mult_count, tester->MXY_mult_count, tester->w21vy_norms_text, latex_format);
  fprintf (fout, table->table_text);
  fprintf (fout, "\n");
  table_io_destroy (table);
  fprintf (fout, "\n");
}

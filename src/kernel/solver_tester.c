#include "solver_tester.h"
#include "common/vectors.h"
#include "linear_ops/vector_ops.h"
#include "common/debug_utils.h"
#include "common/math_utils.h"

void solver_tester_init (solver_tester *tester,
                         const solver_core_workspace *ws,
                         time_layer_func_t g_func,
                         time_layer_func_t vx_func,
                         time_layer_func_t vy_func)
{
  int loc_layer_index;
  int mx, my;
  int n = ws->N;
  vector_double_t dif_g;
  vector_double_t dif_vx;
  vector_double_t dif_vy;
  tester->ws = ws;
  tester->g_func = g_func;
  tester->vx_func = vx_func;
  tester->vy_func = vy_func;


  dif_g  = VECTOR_CREATE (double, ws->layer_size);
  dif_vx = VECTOR_CREATE (double, ws->layer_size);
  dif_vy = VECTOR_CREATE (double, ws->layer_size);


  for (loc_layer_index = 0; loc_layer_index < ws->layer_size; loc_layer_index++)
    {
      double t;
      double x;
      double y;

      solver_workspace_get_mx_my (ws, loc_layer_index, &mx, &my);

      t = n * ws->tau;
      x = mx * ws->hx;
      y = my * ws->hy;

      dif_g[loc_layer_index] = solver_workspace_grid_g (ws, n, mx, my) - g_func (t, x, y);
      dif_vx[loc_layer_index] = solver_workspace_grid_vx (ws, n, mx, my) - vx_func (t, x ,y);
      dif_vy[loc_layer_index] = solver_workspace_grid_vy (ws, n, mx, my) - vy_func (t, x, y);
    }

  tester->g_l2_norm = grid_l2_norm (tester, dif_g);
  tester->vx_l2_norm = grid_l2_norm (tester, dif_vx);
  tester->vy_l2_norm = grid_l2_norm (tester, dif_vy);
  tester->g_c_norm = grid_c_norm (tester, dif_g);
  tester->vx_c_norm = grid_c_norm (tester, dif_vx);
  tester->vy_c_norm = grid_c_norm (tester, dif_vy);
  tester->g_w21_norm = grid_w21_norm (tester, dif_g);
  tester->vx_w21_norm = grid_w21_norm (tester, dif_vx);
  tester->vy_w21_norm = grid_w21_norm (tester, dif_vy);

  VECTOR_DESTROY (dif_g);
  VECTOR_DESTROY (dif_vx);
  VECTOR_DESTROY (dif_vy);
}

double grid_l2_norm (const solver_tester *tester, vector_double_t vec)
{
  double coef = tester->ws->hx * tester->ws->hy;
  double sum = 0;
  int i;
  int mx;
  int my;

  for (i = 0; i < tester->ws->layer_size; i++)
    {
      solver_workspace_get_mx_my (tester->ws, i, &mx, &my);

      if (solver_workspace_get_area (tester->ws, mx, my) == area_internal)
        sum += vec[i] * vec[i];
      else
        sum += 0.5 * vec[i] * vec[i];
    }

  return sqrt (coef * sum);
}

double grid_c_norm (const solver_tester *tester, vector_double_t vec)
{
  double max = 0;
  double val;
  int i;
  for (i = 0; i < tester->ws->layer_size; i++)
    {
      val = fabs (vec[i]);
      max = (max < val) ? val : max;
    }
  return max;
}

double grid_w21_norm (const solver_tester *tester, vector_double_t vec)
{
  double coef = tester->ws->hx * tester->ws->hy;
  double sum = 0;
  int i;
  int mx;
  int my;
  grid_area_t area;

  for (i = 0; i < tester->ws->layer_size; i++)
    {

      solver_workspace_get_mx_my (tester->ws, i, &mx, &my);
      area = solver_workspace_get_area (tester->ws, mx, my);

      if (area != border_right_top
          && area != border_rightmost
          && !(mx == 2 * tester->ws->MX && my == 2 * tester->ws->MY))
        {
          sum += (vec[i + 1] - vec[i]) * (vec[i + 1] - vec[i]);
        }

      if (area != border_topmost
          && area != border_left_hor
          && area != border_right_hor
          && !(mx == 0 && my == tester->ws->MY)
          && !(mx == BOT_ROW_SQUARES_COUNT * tester->ws->MX && my == tester->ws->MY))
        {
          int next_index = solver_workspace_final_index (tester->ws, 0, mx, my + 1);
          sum += (vec[next_index] - vec[i]) * (vec[next_index] - vec[i]);
        }

      if (solver_workspace_get_area (tester->ws, mx, my) == area_internal)
        sum += vec[i] * vec[i];
      else
        sum += 0.5 * vec[i] * vec[i];
    }

  return sqrt (coef * sum);
}

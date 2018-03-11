#include "solver_core_workspace_private.h"

#include <stddef.h>
#include "common/vectors.h"
#include "common/debug_utils.h"
#include "sparse/sparse_base_format.h"
#include "common/math_utils.h"
#include <stdio.h>

int solver_workspace_data_init (solver_core_workspace *data,
                                solver_mode_t mode,
                                int M1,
                                int M2,
                                int N,
                                double X,
                                double Y,
                                double T,
                                double border_omega,
                                linear_solver_t linear_solver)
{
  int success = 1;

  if (!data)
    return 0;

  if (!solver_workspace_data_check_input (M1, M2, N, X, Y, T))
    return 1;

  data->mode = mode;
  data->linear_solver = linear_solver;

  data->MX = M1;
  data->MY = M2;
  data->N = N;
  data->X = X;
  data->Y = Y;
  data->T = T;
  data->hx = X / M1;
  data->hy = Y / M2;
  data->tau = T / N;

  data->border_omega = border_omega;


  data->layer_size = (BOT_ROW_SQUARES_COUNT * M1 + 1) * (M2 + 1) +
                     M2 * (M1 + 1);

  data->vectors_size = (N + 1) * data->layer_size;

  data->matrix_size = UNKNOWN_FUNCTIONS_COUNT * data->layer_size;

  data->log_file = fopen ("log.log", "w");

  success &= ((data->g                 = VECTOR_CREATE (double, data->vectors_size)) != NULL);
  success &= ((data->vx                = VECTOR_CREATE (double, data->vectors_size)) != NULL);
  success &= ((data->vy                = VECTOR_CREATE (double, data->vectors_size)) != NULL);
  success &= ((data->vector_to_compute = VECTOR_CREATE (double, data->matrix_size))  != NULL);
  success &= ((data->rhs_vector        = VECTOR_CREATE (double, data->matrix_size))  != NULL);
  success &= (msr_init_empty (&data->matrix) == 0);
  success &= (sparse_base_init (&data->matrix_base, data->matrix_size, MAX_ROW_NZ) == 0);
  success &= (cgs_solver_init (&data->cgs_linear_solver, data->matrix_size, 1000, 1e-5, precond_jacobi) == 0);

  if (data->linear_solver == laspack_cgs)
    {
      success &= laspack_vector_init (&data->vector_to_compute_l, data->matrix_size);
      success &= laspack_vector_init (&data->rhs_vector_l, data->matrix_size);
    }

  return !success;
}

void solver_workspace_data_destroy (solver_core_workspace *data)
{
  if (!data)
    return;

  VECTOR_DESTROY (data->g);
  VECTOR_DESTROY (data->vx);
  VECTOR_DESTROY (data->vy);
  VECTOR_DESTROY (data->vector_to_compute);
  VECTOR_DESTROY (data->rhs_vector);

  msr_destroy (&data->matrix);
  sparse_base_destroy (&data->matrix_base);
  cgs_solver_destroy (&data->cgs_linear_solver);

  if (data->linear_solver == laspack_cgs)
    {
      laspack_vector_destroy (&data->vector_to_compute_l);
      laspack_vector_destroy (&data->rhs_vector_l);
    }
}


int solver_workspace_data_check_input (int M1, int M2, int N, double X1, double X2, double T)
{
  return !(M1 <= 2 || M2 <= 2 || N <= 2
      || X1 <= 0 || X2 <= 0 || T <= 0);

}

int solver_workspace_layer_begin_index (const solver_core_workspace *data, int n)
{
  return n * data->layer_size;
}

int solver_workspace_top_square_begin_index (const solver_core_workspace *data, int n)
{
  return n * data->layer_size + (BOT_ROW_SQUARES_COUNT * data->MX + 1) * (data->MY + 1);
}

int solver_workspace_final_index (const solver_core_workspace *data, int n, int mx, int my)
{
  int retval;
  int begin = solver_workspace_layer_begin_index (data, n);
  if (my <= data->MY)
    retval = begin + my * (BOT_ROW_SQUARES_COUNT  * data->MX + 1) + mx;
  else
    {
      int loc_mx = mx - data->MX;
      int loc_my = my - data->MY - 1;

      ASSERT_RETURN (mx >= data->MX && mx <= 2 * data->MX, 0);

      retval = solver_workspace_top_square_begin_index (data, n) + loc_my * (data->MX + 1) + loc_mx;
    }

  ASSERT_RETURN (retval < data->vectors_size, 0);
  return retval;
}

double solver_workspace_grid_vx (const solver_core_workspace *data, int n, int mx, int my)
{
  return data->vx[solver_workspace_final_index (data, n, mx, my)];
}

double solver_workspace_grid_vy (const solver_core_workspace *data, int n, int mx, int my)
{
  return data->vy[solver_workspace_final_index (data, n, mx, my)];
}

double solver_workspace_grid_g (const solver_core_workspace *data, int n, int mx, int my)
{
  return data->g[solver_workspace_final_index (data, n, mx, my)];
}

void solver_workspace_fill_layer (solver_core_workspace *data, int n)
{
  int i;
  int index_on_layer = solver_workspace_layer_begin_index (data, n);

  if (data->linear_solver == custom_cgs)
    for (i = 0; i < data->layer_size; i++)
      {
        data->g[index_on_layer]  = data->vector_to_compute[UNKNOWN_FUNCTIONS_COUNT * i];
        data->vx[index_on_layer] = data->vector_to_compute[UNKNOWN_FUNCTIONS_COUNT * i + 1];
        data->vy[index_on_layer] = data->vector_to_compute[UNKNOWN_FUNCTIONS_COUNT * i + 2];
        index_on_layer++;
    }
  else
    for (i = 0; i < data->layer_size; i++)
      {
        data->g[index_on_layer]  = V_GetCmp (&data->vector_to_compute_l.raw, UNKNOWN_FUNCTIONS_COUNT * i + 1);
        data->vx[index_on_layer] = V_GetCmp (&data->vector_to_compute_l.raw, UNKNOWN_FUNCTIONS_COUNT * i + 2);
        data->vy[index_on_layer] = V_GetCmp (&data->vector_to_compute_l.raw, UNKNOWN_FUNCTIONS_COUNT * i + 3);
        index_on_layer++;
      }

/*  if (data->mode == solve_mode && n != 0)
      solver_workspace_check_layer (data, n);
      */
}



void solver_workspace_get_mx_my (const solver_core_workspace *data, int loc_layer_index, int *mx_ptr, int *my_ptr)
{
  int mx = 0;
  int my = 0;

  int points_in_bot_row = BOT_ROW_SQUARES_COUNT * data->MX + 1;
  int points_in_bot_squares = (points_in_bot_row) * (data->MY + 1);

  if (loc_layer_index < points_in_bot_squares)
    {

      my = loc_layer_index / points_in_bot_row;
      mx = loc_layer_index - points_in_bot_row * my;
    }
  else
    {
      int my_loc = (loc_layer_index - points_in_bot_squares) / (data->MX + 1);
      int mx_loc = (loc_layer_index - points_in_bot_squares) - my_loc * (data->MX + 1);

      mx = data->MX + mx_loc;
      my = data->MY + my_loc + 1;
    }

  if (mx_ptr)
    *mx_ptr = mx;

  if (my_ptr)
    *my_ptr = my;
}

grid_area_t solver_workspace_get_area (const solver_core_workspace *data, int mx, int my)
{
  if (mx == 0)
    return border_leftmost;

  if (mx == BOT_ROW_SQUARES_COUNT * data->MX)
    return border_rightmost;

  if (my == 0)
    return border_botmost;

  if (my == 2 * data->MY)
    return border_topmost;

  if (my == data->MY)
    {
      if (mx > 0 && mx <= data->MX)
        return border_left_hor;

      if (mx >= 2 * data->MX && mx < BOT_ROW_SQUARES_COUNT * data->MX)
        return border_right_hor;

      return area_internal;
    }

  if (my > data->MY && my < 2 * data->MY)
    {
      if (mx == data->MX)
        return border_left_top;

      if (mx == 2 * data->MX)
        return border_right_top;

      return area_internal;
    }

  return area_internal;
}

void solver_workspace_check_layer (const solver_core_workspace *data, int n)
{
  int layer_index = solver_workspace_layer_begin_index (data, n);
  int i = 0;

  for (i = 0; i < data->layer_size; i++)
    {
      int mx;
      int my;
      double must_be;
      double actual_g;
      double actual_vx;
      double actual_vy;
      grid_area_t area;

      solver_workspace_get_mx_my (data, i, &mx, &my);

      actual_g = solver_workspace_grid_g (data, n, mx, my);
      actual_vx = solver_workspace_grid_vx (data, n, mx, my);
      actual_vy = solver_workspace_grid_vy (data, n, mx, my);

      area = solver_workspace_get_area (data, mx, my);

      must_be = actual_g;
      DEBUG_ASSERT (solver_workspace_check_border_value (data, area, grid_g, actual_g, &must_be) != 0);
      data->g[layer_index] = must_be;

      must_be = actual_vx;
      DEBUG_ASSERT (solver_workspace_check_border_value (data, area, grid_vx, actual_vx, &must_be) != 0);
      data->vx[layer_index] = must_be;

      must_be = actual_vy;
      DEBUG_ASSERT (solver_workspace_check_border_value (data, area, grid_vy, actual_vy, &must_be) != 0);
      data->vy[layer_index] = must_be;

      layer_index++;
    }
}

int solver_workspace_check_border_value  (const solver_core_workspace *data, grid_area_t area, grid_func_t func, double actual, double *must_be)
{
  double must_be_val;
  switch (area)
    {
    case border_leftmost:
      switch (func)
        {
        case grid_g:
          must_be_val = log (RHO_LEFTMOST);
          break;
        case grid_vx:
          must_be_val = data->border_omega;
          break;
        case grid_vy:
          must_be_val = 0;
          break;
        }
      break;
    case border_botmost:
      switch (func)
        {
        case grid_g:
          return 1;
        case grid_vx:
          must_be_val = 0;
          break;
        case grid_vy:
          must_be_val = 0;
          break;
        }
      break;
    case border_rightmost:
      switch (func)
        {
        case grid_g:
          return 1;
        case grid_vx:
          return 1;
        case grid_vy:
          must_be_val = 0;
          break;
        }
      break;
    case border_topmost:
      switch (func)
        {
        case grid_g:
          return 1;
        case grid_vx:
          must_be_val = 0;
          break;
        case grid_vy:
          return 1;
        }
      break;
    case border_left_hor:
      switch (func)
        {
        case grid_g:
          return 1;
        case grid_vx:
          must_be_val = 0;
          break;
        case grid_vy:
          must_be_val = 0;
          break;
        }
      break;
    case border_left_top:
      switch (func)
        {
        case grid_g:
          return 1;
        case grid_vx:
          must_be_val = 0;
          break;
        case grid_vy:
          must_be_val = 0;
          break;
        }
      break;
    case border_right_hor:
      switch (func)
        {
        case grid_g:
          return 1;
        case grid_vx:
          must_be_val = 0;
          break;
        case grid_vy:
          must_be_val = 0;
          break;
        }
      break;
    case border_right_top:
      switch (func)
        {
        case grid_g:
          return 1;
        case grid_vx:
          return 1;
        case grid_vy:
          return 1;
        }
      break;
    case area_internal:
      return 1;
    }

  if (must_be)
    *must_be = must_be_val;

  return math_fuzzy_eq (must_be_val, actual);
}

double solver_workspace_grid_val (const solver_core_workspace *data, int n, int mx, int my, grid_func_t func)
{
  switch (func)
    {
    case grid_g:
      return solver_workspace_grid_g (data, n, mx, my);
    case grid_vx:
      return solver_workspace_grid_vx (data, n, mx, my);
    case grid_vy:
      return solver_workspace_grid_vy (data, n, mx, my);
    }
  ASSERT_RETURN (0, 0);
}

int solver_workspace_func_col (const solver_core_workspace *data, int loc_layer_index, grid_func_t func)
{
  ASSERT_RETURN (loc_layer_index < data->matrix_size, -1);

  switch (func)
    {
    case grid_g:
      return UNKNOWN_FUNCTIONS_COUNT * loc_layer_index;
    case grid_vx:
      return UNKNOWN_FUNCTIONS_COUNT * loc_layer_index + 1;
    case grid_vy:
      return UNKNOWN_FUNCTIONS_COUNT * loc_layer_index + 2;
    }

  ASSERT_RETURN (0, -1);
}

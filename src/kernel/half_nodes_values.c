#include "half_nodes_values.h"
#include "common/debug_utils.h"
#include "kernel_typedefs.h"
#include "common/vectors.h"

int hn_values_init (half_nodes_values *vs, int MX, int MY, int N)
{
  vs->MX = MX;
  vs->MY = MY;
  vs->N = N;

/*  vs->layer_size = SQUARES_COUNT * vs->MX * vs->MY +
                   BOT_ROW_SQUARES_COUNT * 2 * vs->MX +
                   4 * vs->MY - 2;*/

  vs->layer_size = SQUARES_COUNT * MX * MY;
  vs->layer_size += BOT_ROW_SQUARES_COUNT * vs->MX;
  vs->layer_size += 2 * vs->MY;
  vs->layer_size += 2 * vs->MX;
  vs->layer_size += 2 * (vs->MY - 1);
  vs->layer_size += vs->MX;

  vs->vector_size = (N + 1) * vs->layer_size;

  vs->vals = VECTOR_CREATE (double, vs->vector_size);

  return vs->vals != NULL;
}

void hn_values_destroy (half_nodes_values *vs)
{
  VECTOR_DESTROY (vs->vals);
}

int hn_values_index (const half_nodes_values *vs, int n, int mx, int my)
{
  int row;
  int col;

  int retval = vs->vector_size;
  int begin = n * vs->layer_size;

  int low_sq_begin = BOT_ROW_SQUARES_COUNT * vs->MX;
  int mid_line_begin = low_sq_begin + vs->MY * (BOT_ROW_SQUARES_COUNT * vs->MX + 2);
  int top_sq_begin = mid_line_begin + BOT_ROW_SQUARES_COUNT * vs->MX;
  int top_line_begin = top_sq_begin + (vs->MY - 1) * (vs->MX + 2);

  if (my == -1)
    {
      row = 0;
      col = mx;

      ASSERT_RETURN (col >= 0 && col < BOT_ROW_SQUARES_COUNT * vs->MX, 0);

      retval = begin + col;
    }
  if (my < vs->MY && my > -1)
    {
      row = my;
      col = mx + 1;


      ASSERT_RETURN (col >= 0 && col <= BOT_ROW_SQUARES_COUNT * vs->MX + 1, 0);

      retval = begin
               + low_sq_begin
               + row * (BOT_ROW_SQUARES_COUNT  * vs->MX + 2) + col;
    }
  if (my == vs->MY)
    {
      row = 0;
      col = mx;

      ASSERT_RETURN (col >= 0 && col < BOT_ROW_SQUARES_COUNT * vs->MX, 0);

      retval = begin
               + mid_line_begin
               + col;
    }
  if (my >= vs->MY + 1 && my < 2 * vs->MY)
    {
      row = my  - vs->MY - 1;
      col = mx - vs->MX + 1;

      ASSERT_RETURN (col >= 0 && col <= vs->MX + 1
                     && row  < vs->MY - 1, 0);

      retval = begin
               + top_sq_begin
               + row * (vs->MX + 2) + col;
    }
  if (my == 2 * vs->MY)
    {
      row = 0;
      col = mx - vs->MX;

      ASSERT_RETURN (col >= 0 && col < vs->MX, 0);

      retval = begin
               + top_line_begin
               + col;
    }

  ASSERT_RETURN (retval < vs->vector_size, 0);
  return retval;
}

void hn_values_get_mx_my (const half_nodes_values *vs, int loc_layer_index, int *mx_ptr, int *my_ptr)
{
  int mx = -2;
  int my = -2;

  int row;
  int col;

  int low_sq_begin = BOT_ROW_SQUARES_COUNT * vs->MX;
  int mid_line_begin = low_sq_begin + vs->MY * (BOT_ROW_SQUARES_COUNT * vs->MX + 2);
  int top_sq_begin = mid_line_begin + BOT_ROW_SQUARES_COUNT * vs->MX;
  int top_line_begin = top_sq_begin + (vs->MY - 1) * (vs->MX + 2);

  if (loc_layer_index < low_sq_begin)
    {
      mx = loc_layer_index;
      my = -1;
    }
  if (loc_layer_index < mid_line_begin && loc_layer_index >= low_sq_begin)
    {
      int ind = loc_layer_index - low_sq_begin;
      row = ind / (BOT_ROW_SQUARES_COUNT * vs->MX + 2);
      col = ind - (BOT_ROW_SQUARES_COUNT * vs->MX + 2) * row;

      mx = col - 1;
      my = row;
    }
  if (loc_layer_index < top_sq_begin && loc_layer_index >= mid_line_begin)
    {
      int ind = loc_layer_index - mid_line_begin;
      my = vs->MY;
      mx = ind;
    }
  if (loc_layer_index < top_line_begin && loc_layer_index >= top_sq_begin)
    {
      int ind = loc_layer_index - top_sq_begin;
      row = ind / (vs->MX + 2);
      col = ind - (vs->MX + 2) * row;

      mx = col + vs->MX - 1;
      my = row + vs->MY + 1;
    }
  if (loc_layer_index >= top_line_begin)
    {
      int ind = loc_layer_index - top_line_begin;

      ASSERT_RETURN (loc_layer_index < vs->layer_size, );

      my = 2 * vs->MY;
      mx = ind + vs->MX;
    }

  DEBUG_ASSERT (mx != -2 && my != -2);

  if (mx_ptr)
    *mx_ptr = mx;

  if (my_ptr)
    *my_ptr = my;
}

void hn_border_iter_init (hn_border_iter *iter, const half_nodes_values *vs)
{
  iter->mx = -1;
  iter->my = 0;
  iter->vs = vs;
  iter->is_end = 0;
  iter->area = border_leftmost;
}

void hn_border_iter_next (hn_border_iter *iter)
{
  int next_mx = iter->mx;
  int next_my = iter->my;
  grid_area_t next_area = iter->area;

  switch (iter->area)
    {
    case border_leftmost:
      if (iter->my == iter->vs->MY - 1)
        {
          next_mx = 0;
          next_my = -1;
          next_area = border_botmost;
        }
      else
        next_my++;
      break;
    case border_botmost:
      if (iter->mx == BOT_ROW_SQUARES_COUNT * iter->vs->MX - 1)
        {
          next_mx = BOT_ROW_SQUARES_COUNT * iter->vs->MX;
          next_my = 0;
          next_area = border_rightmost;
        }
      else
        next_mx++;
      break;
    case border_rightmost:
      if (iter->my == iter->vs->MY - 1)
        {
          next_mx = iter->vs->MX;
          next_my = 2 * iter->vs->MY;
          next_area = border_topmost;
        }
      else
        next_my++;
      break;
    case border_topmost:
      if (iter->mx == 2 * iter->vs->MX - 1)
        {
          next_mx = 0;
          next_my = iter->vs->MY;
          next_area = border_left_hor;
        }
      else
        next_mx++;
      break;
    case border_left_hor:
      if (iter->mx == iter->vs->MX - 1)
        {
          next_my++;
          next_area = border_left_top;
        }
      else
        next_mx++;
      break;
    case border_left_top:
      if (iter->my == 2 * iter->vs->MY - 1)
        {
          next_mx = 2 * iter->vs->MX;
          next_my = iter->vs->MY;
          next_area = border_right_hor;
        }
      else
        next_my++;
      break;
    case border_right_hor:
      if (iter->mx == BOT_ROW_SQUARES_COUNT * iter->vs->MX - 1)
        {
          next_mx = 2 * iter->vs->MX;
          next_my++;
          next_area = border_right_top;
        }
      else
        next_mx++;
      break;
    case border_right_top:
      if (iter->my == 2 * iter->vs->MY - 1)
        iter->is_end = 1;
      else
        next_my++;
      break;
    case area_internal:
      DEBUG_ASSERT (0);
      break;
    }

  iter->mx = next_mx;
  iter->my = next_my;
  iter->area = next_area;
}

int hn_values_is_border (const half_nodes_values *vs, int lli)
{
  int mx;
  int my;

  hn_values_get_mx_my (vs, lli, &mx, &my);

  if (mx == -1 || my == -1)
    return 1;

  if (mx == BOT_ROW_SQUARES_COUNT * vs->MX || my == 2 * vs->MY)
    return 1;

  if (my == vs->MY)
    {
      if (mx >= 0 && mx < vs->MX)
        return 1;

      if (mx >= 2 * vs->MX && mx < BOT_ROW_SQUARES_COUNT * vs->MX)
        return 1;
    }

  if (mx == vs->MX - 1 || mx == 2 * vs->MX)
    {
      if (my >= vs->MY + 1 && my < 2 * vs->MY)
        return 1;
    }

  return 0;
}

grid_area_t hn_values_get_area (const half_nodes_values *vs, int mx, int my)
{
  if (mx == -1)
    return border_leftmost;

  if (mx == BOT_ROW_SQUARES_COUNT * vs->MX)
    return border_rightmost;

  if (my == 2 * vs->MY)
    return border_topmost;

  if (my == -1)
    return border_botmost;

  if (my == vs->MY
      && (0 <= mx && mx < vs->MX))
    return border_left_hor;

  if (my == vs->MY
      && (2 * vs->MX <= mx && mx < BOT_ROW_SQUARES_COUNT * vs->MX))
    return border_right_hor;

  if (mx == vs->MX - 1
      && vs->MY < my)
    return border_left_top;

  if (mx == 2 * vs->MX
      && vs->MY < my)
    return border_right_top;

  return area_internal;
}

double hn_values_val_by_index (const half_nodes_values *vs, int n, int lli)
{
  return vs->vals[n * vs->layer_size + lli];
}

double hn_values_mx_my_val (const half_nodes_values *vs, int n, int mx, int my)
{
  return vs->vals[hn_values_index (vs, n, mx, my)];
}

double hn_values_approx_in_node (const half_nodes_values *vs, int n, int mx, int my)
{
  double h1 = hn_values_mx_my_val (vs, n, mx, my);
  double h2 = hn_values_mx_my_val (vs, n, mx, my - 1);
  double h3 = hn_values_mx_my_val (vs, n, mx - 1, my - 1);
  double h4 = hn_values_mx_my_val (vs, n, mx - 1, my);

  return (h1 + h2 + h3 + h4) / 4;
}

double hn_values_avg_bwd_x (const half_nodes_values *vs, int n, int mx, int my)
{
  double h1 = hn_values_mx_my_val (vs, n, mx, my);
  double h2 = hn_values_mx_my_val (vs, n, mx - 1, my);

  return (h1 + h2) / 2;
}

double hn_values_avg_bwd_y (const half_nodes_values *vs, int n, int mx, int my)
{
  double h1 = hn_values_mx_my_val (vs, n, mx, my);
  double h2 = hn_values_mx_my_val (vs, n, mx, my - 1);

  return (h1 + h2) / 2;
}

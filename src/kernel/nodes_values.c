#include "nodes_values.h"
#include "kernel_typedefs.h"
#include "common/vectors.h"
#include "common/debug_utils.h"

int nodes_values_init (nodes_values *vs, int MX, int MY, int N)
{
  vs->MX = MX;
  vs->MY = MY;
  vs->N = N;

  vs->layer_size = (BOT_ROW_SQUARES_COUNT * MX + 1) * (MY + 1) +
                   MY * (MX + 1);

  vs->vector_size = (N + 1) * vs->layer_size;

  vs->vals = VECTOR_CREATE (double, vs->vector_size);

  return vs->vals != NULL;
}

void nodes_values_destroy (nodes_values *vs)
{
  VECTOR_DESTROY (vs->vals);
}

int nodes_values_index (const nodes_values *vs, int n, int mx, int my)
{
  int retval;
  int begin = n * vs->layer_size;
  if (my <= vs->MY)
    retval = begin + my * (BOT_ROW_SQUARES_COUNT  * vs->MX + 1) + mx;
  else
    {
      int loc_mx = mx - vs->MX;
      int loc_my = my - vs->MY - 1;

      ASSERT_RETURN (mx >= vs->MX && mx <= 2 * vs->MX, 0);

      retval = begin + (vs->MY + 1) * (BOT_ROW_SQUARES_COUNT * vs->MX + 1) + loc_my * (vs->MX + 1) + loc_mx;
    }

  ASSERT_RETURN (retval < vs->vector_size, 0);
  return retval;
}

void nodes_values_get_mx_my (const nodes_values *vs, int loc_layer_index, int *mx_ptr, int *my_ptr)
{
  int mx = 0;
  int my = 0;

  int points_in_bot_row = BOT_ROW_SQUARES_COUNT * vs->MX + 1;
  int points_in_bot_squares = (points_in_bot_row) * (vs->MY + 1);

  if (loc_layer_index < points_in_bot_squares)
    {

      my = loc_layer_index / points_in_bot_row;
      mx = loc_layer_index - points_in_bot_row * my;
    }
  else
    {
      int my_loc = (loc_layer_index - points_in_bot_squares) / (vs->MX + 1);
      int mx_loc = (loc_layer_index - points_in_bot_squares) - my_loc * (vs->MX + 1);

      mx = vs->MX + mx_loc;
      my = vs->MY + my_loc + 1;
    }

  if (mx_ptr)
    *mx_ptr = mx;

  if (my_ptr)
    *my_ptr = my;
}

double nodes_values_val (const nodes_values *vs, int n, int lli)
{
  return vs->vals[vs->layer_size * n + lli];
}

double nodes_values_mx_my_val (const nodes_values *vs, int n, int mx, int my)
{
  return vs->vals[nodes_values_index (vs, n, mx, my)];
}

void nodes_border_iter_init (nodes_border_iter *iter, const nodes_values *vs)
{
  iter->mx = 0;
  iter->my = 0;
  iter->is_end = 0;
  iter->area = border_leftmost;
  iter->vs = vs;
}

void nodes_border_iter_next (nodes_border_iter *iter)
{
  int next_mx = iter->mx;
  int next_my = iter->my;
  grid_area_t next_area = iter->area;

  switch (iter->area)
    {
    case border_leftmost:
      if (iter->my == iter->vs->MY)
        {
          next_mx = 1;
          next_my = 0;
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
      if (iter->my == iter->vs->MY)
        {
          next_mx = iter->vs->MX;
          next_my = 2 * iter->vs->MY;
          next_area = border_topmost;
        }
      else
        next_my++;
      break;
    case border_topmost:
      if (iter->mx == 2 * iter->vs->MX)
        {
          next_mx = 1;
          next_my = iter->vs->MY;
          next_area = border_left_hor;
        }
      else
        next_mx++;
      break;
    case border_left_hor:
      if (iter->mx == iter->vs->MX)
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

double nodes_avg_fwd_x (const nodes_values *vs, int n, int mx, int my)
{
  int i1 = nodes_values_index (vs, n, mx, my);
  int i2 = nodes_values_index (vs, n, mx + 1, my);

  return (vs->vals[i1] + vs->vals[i2]) / 2;
}

double nodes_avg_fwd_y (const nodes_values *vs, int n, int mx, int my)
{
  int i1 = nodes_values_index (vs, n, mx, my);
  int i2 = nodes_values_index (vs, n, mx, my + 1);

  return (vs->vals[i1] + vs->vals[i2]) / 2;
}

grid_area_t nodes_values_get_area (const nodes_values *vs, int mx, int my)
{
  if (mx == 0)
    return border_leftmost;

  if (mx == BOT_ROW_SQUARES_COUNT * vs->MX)
    return border_rightmost;

  if (my == 2 * vs->MY)
    return border_topmost;

  if (my == 0)
    return border_botmost;

  if (my == vs->MY &&
      0 < mx && mx <= vs->MX)
    return border_left_hor;

  if (my == vs->MY &&
      2 * vs->MX <= mx && mx < BOT_ROW_SQUARES_COUNT * vs->MX)
    return border_right_hor;

  if (mx == vs->MX &&
      vs->MY < my && my < 2 * vs->MY)
    return border_left_top;

  if (mx == 2 * vs->MX &&
      vs->MY < my && my < 2 * vs->MY)
    return border_right_top;

  return area_internal;
}

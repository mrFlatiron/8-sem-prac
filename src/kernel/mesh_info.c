#include "mesh_info.h"

int mesh_info_init (mesh_info_t *info, double X, double Y, double T, int MX, int MY, int N, double border_omega, double mu)
{
  info->X = X;
  info->Y = Y;
  info->T = T;
  info->MX = MX;
  info->MY = MY;
  info->N = N;
  info->border_omega = border_omega;
  info->mu = mu;

  if (X <= 0
      || Y <= 0
      || T <= 0
      || MX <= 2
      || MY <= 2
      || N <= 0)
    return 0;

  return 1;
}

grid_area_t mesh_info_get_area (const mesh_info_t *info, int mx, int my)
{
  if (mx == 0)
    return border_leftmost;

  if (mx == BOT_ROW_SQUARES_COUNT * info->MX)
    return border_rightmost;

  if (my == 0)
    return border_botmost;

  if (my == 2 * info->MY)
    return border_topmost;

  if (my == info->MY)
    {
      if (mx > 0 && mx <= info->MX)
        return border_left_hor;

      if (mx >= 2 * info->MX && mx < BOT_ROW_SQUARES_COUNT * info->MX)
        return border_right_hor;

      return area_internal;
    }

  if (my > info->MY && my < 2 * info->MY)
    {
      if (mx == info->MX)
        return border_left_top;

      if (mx == 2 * info->MX)
        return border_right_top;

      return area_internal;
    }

  return area_internal;
}

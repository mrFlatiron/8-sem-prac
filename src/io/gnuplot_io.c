#include "io/gnuplot_io.h"
#include "common/vectors.h"
#include <string.h>

void gnuplot_io_init (gnuplot_io *handle, const solver_core_workspace *ws, int layer)
{
  int size = ws->layer_size;
  int i;

  handle->full_size_g = size * 3;
  handle->full_size_v = size * 4;

  handle->coords_g = VECTOR_CREATE (string_t, handle->full_size_g);

  for (i = 0; i < handle->full_size_g; i++)
    handle->coords_g[i] = VECTOR_CREATE (char, 50);

  handle->coords_v = VECTOR_CREATE (string_t, handle->full_size_v);

  for (i = 0; i < handle->full_size_v; i++)
    handle->coords_v[i] = VECTOR_CREATE (char, 50);

  for (i = 0; i < size; i++)
    {
      int mx;
      int my;
      double x;
      double y;

      solver_workspace_get_mx_my (ws, i, &mx, &my);
      x = mx * ws->hx;
      y = my * ws->hy;

      sprintf (handle->coords_g[3 * i], "%.8f", x);
      sprintf (handle->coords_g[3 * i + 1], "%.8f", y);
      sprintf (handle->coords_g[3 * i + 2], "%.8f", solver_workspace_grid_g (ws, layer, mx, my));

      sprintf (handle->coords_v[4 * i], "%.8f", x);
      sprintf (handle->coords_v[4 * i + 1], "%.8f", y);
      sprintf (handle->coords_v[4 * i + 2], "%.8f", solver_workspace_grid_vx (ws, layer, mx, my));
      sprintf (handle->coords_v[4 * i + 3], "%.8f", solver_workspace_grid_vy (ws, layer, mx, my));
    }
}

void gnuplot_io_destroy (gnuplot_io *handle)
{
  int i;

  for (i = 0; i < handle->full_size_g; i++)
    VECTOR_DESTROY (handle->coords_g[i]);

  VECTOR_DESTROY (handle->coords_g);

  for (i = 0; i < handle->full_size_v; i++)
    VECTOR_DESTROY (handle->coords_v[i]);

  VECTOR_DESTROY (handle->coords_v);
}

#ifndef GNUPLOT_IO_H
#define GNUPLOT_IO_H

#include <stdio.h>
#include "common/vectors_fwd.h"
#include "kernel/solver_core_workspace.h"
#include "kernel/sokolov_solver.h"

typedef struct
{
  string_t *coords_g;
  string_t *coords_v;

  int full_size_g;
  int full_size_v;
} gnuplot_io;

void gnuplot_io_init_by_ws (gnuplot_io *handle, const solver_core_workspace *ws, int layer);
void gnuplot_io_init_by_sokolov (gnuplot_io *handle, const sokolov_solver *sok, int layer);
void gnuplot_io_destroy (gnuplot_io *handle);

#endif /* GNUPLOT_IO_H */

#include "central_diff_solver_private.h"
#include "input/rhs.h"
#include "cgs_solver.h"
#include "common/debug_utils.h"

int cdiff_solver_init (central_diff_solver *solver,
                       int M1,
                       int M2,
                       int N, double X,
                       double Y,
                       double T)
{
  solver_workspace_data_init (&solver->ws, M1, M2, N, X, Y, T);

  return 0;
}



void cdiff_solver_destroy (central_diff_solver *solver)
{
  solver_workspace_data_destroy (&solver->ws);
}



int cdiff_solver_compute (central_diff_solver *solver,
                          pressure_func_t p_func,
                          double mu,
                          area_func_t f0,
                          area_func_t f1, layer_func_t start_vx, layer_func_t start_vy, layer_func_t start_g)
{
  int i;

  solver->f0 = f0;
  solver->f1 = f1;
  solver->start_vx = start_vx;
  solver->start_vy = start_vy;
  solver->start_g  = start_g;

  switch (p_func)
    {
    case pressure_linear:
      solver->p_drv = &p_drv_linear;
      break;
    case pressure_polynomial:
      solver->p_drv = &p_drv_polynomial;
      break;
    }

  solver->mu = mu;

  cdiff_solver_init_first_layer (solver);

  for (i = 1; i <= solver->ws.N; i++)
    {
      cdiff_solver_fill_matrix (solver);
      cdiff_solver_fill_rhs (solver);
      cdiff_solver_solve_system (solver);
      solver_workspace_fill_layer (&solver->ws, i);
    }

  return 0;
}

void cdiff_solver_init_first_layer (central_diff_solver *solver)
{
  FIX_UNUSED (solver);
}

void cdiff_solver_fill_rhs (central_diff_solver *solver)
{
  FIX_UNUSED (solver);
}

void cdiff_solver_fill_matrix (central_diff_solver *solver)
{
  FIX_UNUSED (solver);
}

void cdiff_solver_solve_system (central_diff_solver *solver)
{
  FIX_UNUSED (solver);
}

#include "central_diff_solver_private.h"

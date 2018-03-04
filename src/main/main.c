#include <stdio.h>
#include <locale.h>

#include "common/vectors.h"
#include "common/debug_utils.h"
#include "kernel/solver_core_workspace.h"
#include "kernel/command_line_parser.h"
#include "kernel/cgs_solver.h"
#include "linear_ops/vector_ops.h"
#include "kernel/central_diff_solver.h"

#include "kernel/input/rhs.h"
#include "kernel/input/t0_functions.h"
#include "kernel/input/test_solutions.h"




int main (int argc, char *argv[])
{


  command_line_parser     parser_object;
  command_line_parser_ptr parser         = &parser_object;

  central_diff_solver      solver_object;
  central_diff_solver      *solver = &solver_object;

  int                     error_code;

  setlocale (LC_ALL, "en-GB.utf8");

  error_code = parse_command_line (parser, argc, argv);

  if (error_code)
    {
      fprintf (stderr, "%s\n", parser_info_str (parser, error_code));
      return 1;
    }


  error_code = cdiff_solver_init (solver,
                                  parser->solver_mode,
                                  parser->M1,
                                  parser->M2,
                                  parser->N,
                                  X_LEN,
                                  Y_LEN,
                                  parser->T,
                                  BORDER_OMEGA,
                                  test_g,
                                  test_vx,
                                  test_vy);

 if (error_code)
   {
     fprintf (stderr, "Could not initialize solver with given arguments\n");
     return error_code;
   }

 if (parser->solver_mode == test_mode)
   {
     cdiff_solver_compute (solver,
                           parser->p_func,
                           parser->mu,
                           rhs_test_f0,
                           rhs_test_f1,
                           rhs_test_f2,
                           t0_vx_test,
                           t0_vy_test,
                           t0_rho_test);
   }
 else
   {
     cdiff_solver_compute (solver,
                           parser->p_func,
                           parser->mu,
                           rhs_f0,
                           rhs_f1,
                           rhs_f2,
                           t0_vx,
                           t0_vy,
                           t0_rho);
   }

  cdiff_solver_destroy (solver);
#if 0
  double dense[] = {1, 1, 0, 0, 0,
                    0, 2, 0, 0, 0,
                    0, 2, 5, 0, 0,
                    0, 0, 4, 6, 8,
                    0, 0, 0, 0, 7};
  msr_matrix sparse_object;
  msr_matrix* sparse = &sparse_object;

  cgs_solver solver_object;
  cgs_solver *solver = &solver_object;

  double x[] = {1, 0, 1, 0, 1};

  double x_init[] = {0, 0, 0, 0, 0};

  double rhs[5];

  cgs_solver_error_t error;



  FIX_UNUSED (argc);
  FIX_UNUSED (argv);
  FIX_UNUSED (error);


  msr_init_from_vector (sparse, dense, 5);

  msr_dump (sparse, stdout);

  msr_mult_vector (sparse, x, rhs);

  cgs_solver_init (solver,
                   sparse->matrix_size,
                   300,
                   1e-6,
                   precond_jacobi);


  error = cgs_solver_solve (solver,
                    sparse,
                    rhs,
                    x_init,
                    x);

  cgs_solver_destroy (solver);
  msr_destroy (sparse);
#endif
  return 0;
}

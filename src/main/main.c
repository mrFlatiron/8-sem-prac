#include <stdio.h>

#include "common/vectors.h"
#include "common/debug_utils.h"
#include "kernel/solver_core_workspace.h"
#include "kernel/command_line_parser.h"
#include "kernel/cgs_solver.h"
#include "linear_ops/vector_ops.h"



int main (int argc, char *argv[])
{
#if 0
  command_line_parser     parser_object;
  command_line_parser_ptr parser         = &parser_object;

  general_solver_data     solver_data_object;
  general_solver_data     *solver_data         = &solver_data_object;

  int                     error_code;

  error_code = parse_command_line (parser, argc, argv);

  if (error_code)
    {
      fprintf (stderr, "%s\n", parser_info_str (parser, error_code));
      return 1;
    }


  error_code = solver_workspace_data_init (solver_data,
                                     parser->p_func,
                                     parser->M1,
                                     parser->M2,
                                     parser->N,
                                     parser->X,
                                     parser->Y,
                                     parser->T);

 if (error_code)
   {
     fprintf (stderr, "Could not initialize solver with given arguments\n");
     return error_code;
   }



  solver_workspace_data_destroy (solver_data);
#endif

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


  msr_init_frovector (sparse, dense, 5);

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

  return 0;
}

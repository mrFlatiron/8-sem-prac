#include <stdio.h>

#include "common/vectors.h"
#include "common/debug_utils.h"
#include "kernel/central_differences_solver.h"
#include "kernel/command_line_parser.h"

int main (int argc, char *argv[])
{
  command_line_parser     parser_object;
  command_line_parser_ptr parser         = &parser_object;
  cd_solver               solver_object;
  cd_solver_ptr           solver         = &solver_object;
  int                     error_code;

  error_code = parse_command_line (parser, argc, argv);

  if (error_code)
    {
      fprintf (stderr, "%s\n", parser_error_str (parser, error_code));
      return 1;
    }

  cdsolver_init (solver,
                 parser->p_func,
                 parser->M1,
                 parser->M2,
                 parser->N,
                 parser->X1,
                 parser->X2,
                 parser->T);


  error_code = cdsolver_compute (solver);

  if (error_code)
    {
      fprintf (stderr, "Failed to compute a solution\n");
      return error_code;
    }


  cdsolver_destroy (solver);

  return 0;
}

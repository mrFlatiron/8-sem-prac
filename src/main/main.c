#include <stdio.h>

#include "common/vectors.h"
#include "common/debug_utils.h"
#include "kernel/general_solver_data.h"
#include "kernel/command_line_parser.h"

#include <sys/sysinfo.h>

int main (int argc, char *argv[])
{
  command_line_parser     parser_object;
  command_line_parser_ptr parser         = &parser_object;
  general_solver_data     solver_data_object;
  general_solver_data     *solver_data         = &solver_data_object;
  int                     error_code;

  error_code = parse_command_line (parser, argc, argv);

  DEBUG_PAUSE ("Stop");

  if (error_code)
    {
      fprintf (stderr, "%s\n", parser_info_str (parser, error_code));
      return 1;
    }

  error_code = gen_solver_data_init (solver_data,
                                     parser->p_func,
                                     parser->M1,
                                     parser->M2,
                                     parser->N,
                                     parser->X1,
                                     parser->X2,
                                     parser->T,
                                     get_nprocs ());

 if (error_code)
   {
     fprintf (stderr, "Could not initialize solver with given arguments\n");
     return error_code;
   }



  gen_solver_data_destroy (solver_data);

  return 0;
}

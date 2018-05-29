#ifndef COMMAND_LINE_PARSER_H
#define COMMAND_LINE_PARSER_H

#include "kernel_typedefs.h"
#include "common/vectors_fwd.h"

#define OPTION_BUFFER_LEN 1024
#define VERSION_STR_BUFFER_LEN 100
#define HELP_STR_BUFFER_LEN 4096

typedef struct
{
  pressure_func_t       p_func;
  solver_t              solver;
  solver_mode_t         solver_mode;

  int                   N;
  int                   MX;
  int                   MY;

  double                T;
  double                border_omega;
  double                mu;

  int                   N_mult;
  int                   N_mult_count;
  int                   MXY_mult;
  int                   MXY_mult_count;

  preconditioner_t      precond;
  linear_solver_t       linear_solver;
  double                solver_precision;
  int                   solver_max_iter;

  char                  version_str[VERSION_STR_BUFFER_LEN];
  char                  help_str[HELP_STR_BUFFER_LEN];
} command_line_parser;
typedef command_line_parser* command_line_parser_ptr;

int parse_command_line (command_line_parser *parser, int argc, char *argv[]);
const char *parser_info_str (command_line_parser *parser, int error_code);
const char *parser_help_str (command_line_parser *parser);
const char *parser_version_str (command_line_parser *parser);

void encode_input_parameters (const command_line_parser *parser, string_t out);

#endif /* COMMAND_LINE_PARSER_H */

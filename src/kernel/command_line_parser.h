#ifndef COMMAND_LINE_PARSER_H
#define COMMAND_LINE_PARSER_H

#include "kernel_typedefs.h"

#define OPTION_BUFFER_LEN 1024

typedef struct
{
  pressure_func_t       p_func;
  solver_t              solver;
  table_output_format_t table_output_format;
  solver_mode_t         solver_mode;
  int                   N;
  int                   M1;
  int                   M2;
  double                T;
  double                X1;
  double                X2;

} command_line_parser;
typedef command_line_parser* command_line_parser_ptr;

int parse_command_line (command_line_parser *parser, int argc, char *argv[]);
const char *parser_error_str (command_line_parser *parser, int error_code);

//Below are private methods

const char *get_value (const char *option, int argc, char *argv[]);



#endif // COMMAND_LINE_PARSER_H

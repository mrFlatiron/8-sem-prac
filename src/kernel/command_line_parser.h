#ifndef COMMAND_LINE_PARSER_H
#define COMMAND_LINE_PARSER_H

#include "kernel_typedefs.h"

#define OPTION_BUFFER_LEN 1024
#define VERSION_STR_BUFFER_LEN 100
#define HELP_STR_BUFFER_LEN 1024

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
  char                  version_str[VERSION_STR_BUFFER_LEN];
  char                  help_str[HELP_STR_BUFFER_LEN];
} command_line_parser;
typedef command_line_parser* command_line_parser_ptr;

int parse_command_line (command_line_parser *parser, int argc, char *argv[]);
const char *parser_info_str (command_line_parser *parser, int error_code);
const char *parser_help_str (command_line_parser *parser);
const char *parser_version_str (command_line_parser *parser);

/*Below are private methods*/

const char *parser_get_value (command_line_parser *parser, const char *option, int argc, char *argv[]);
int parser_is_help_present (command_line_parser *parser, int argc, char *argv[]);
int parser_is_version_present (command_line_parser *parser, int argc, char *argv[]);


#endif /* COMMAND_LINE_PARSER_H */

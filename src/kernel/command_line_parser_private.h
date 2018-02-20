#ifndef COMMAND_LINE_PARSER_PRIVATE_H
#define COMMAND_LINE_PARSER_PRIVATE_H

#include "command_line_parser.h"

const char *parser_get_value (command_line_parser *parser, const char *option, int argc, char *argv[]);
int parser_is_help_present (command_line_parser *parser, int argc, char *argv[]);
int parser_is_version_present (command_line_parser *parser, int argc, char *argv[]);

#endif /* COMMAND_LINE_PARSER_PRIVATE_H */

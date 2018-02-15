#include "command_line_parser.h"
#include "common/debug_utils.h"

#include <string.h>
#include <stdlib.h>

int parse_command_line (command_line_parser *parser, int argc, char *argv[])
{
  const char *solver_value;
  const char *pressure_value;
  const char *table_format_value;
  const char *mode_value;
  const char *N_value;
  const char *M1_value;
  const char *M2_value;
  const char *T_value;
  const char *X1_value;
  const char *X2_value;


  if (!parser)
    return -1;

  if (argc == 1)
    return -2;

  solver_value       = parser_get_value (parser, "solver",       argc, argv);
  pressure_value     = parser_get_value (parser, "pressure",     argc, argv);
  table_format_value = parser_get_value (parser, "table-format", argc, argv);
  mode_value         = parser_get_value (parser, "solver-mode",  argc, argv);
  N_value            = parser_get_value (parser, "N",            argc, argv);
  M1_value           = parser_get_value (parser, "M1",           argc, argv);
  M2_value           = parser_get_value (parser, "M2",           argc, argv);
  T_value            = parser_get_value (parser, "T",            argc, argv);
  X1_value           = parser_get_value (parser, "X1",           argc, argv);
  X2_value           = parser_get_value (parser, "X2",           argc, argv);

  if (parser_is_help_present (parser, argc, argv))
    return -2;

  if (parser_is_version_present (parser, argc, argv))
    return -3;

  if (!solver_value)
    return 1;

  if (!strcmp ("central", solver_value))
    parser->solver = central_differences;
  else
    {
      if (!strcmp ("sokolov", solver_value))
        parser->solver = sokolov;
      else
        return 1;
    }

  if (!pressure_value)
    parser->p_func = pressure_linear;
  else
    {
      if (!strcmp ("linear", pressure_value))
        parser->p_func = pressure_linear;
      else
        {
          if (!strcmp ("polynomial", pressure_value))
            parser->p_func = pressure_polynomial;
          else
            return 2;
        }
    }

  if (!table_format_value)
    parser->table_output_format = human_readable;
  else
    {
      if (!strcmp ("latex", table_format_value))
        parser->table_output_format = latex_format;
      else
        {
          if (!strcmp ("simple", table_format_value))
            parser->table_output_format = human_readable;
          else
            return 3;
        }
    }

  if (!mode_value)
    return 4;

  if (!strcmp ("test", mode_value))
    parser->solver_mode = test_mode;
  else
    {
      if (!strcmp ("solve", mode_value))
        parser->solver = solve_mode;
      else
        return 4;
    }

  if (!N_value)
    return 5;

  parser->N = atoi (N_value);

  if (!M1_value)
    return 6;

  parser->M1 = atoi (M1_value);

  if (!M2_value)
    return 7;

  parser->M2 = atoi (M2_value);

  if (!T_value)
    return 8;

  parser->T = atof (T_value);

  if (!X1_value)
    return 9;

  parser->X1 = atof (X1_value);

  if (!X2_value)
    return 10;

  parser->X2 = atof (X2_value);

  return 0;
}

const char *parser_get_value (command_line_parser *parser, const char *option, int argc, char *argv[])
{
  char str_to_find[OPTION_BUFFER_LEN];
  int  str_to_find_len;
  const char *option_value = NULL;
  int i;

  FIX_UNUSED (parser);

  str_to_find[0] = 0;
  str_to_find_len = 0;

  strcat (str_to_find, "--");
  strcat (str_to_find, option);
  strcat (str_to_find, "=");
  str_to_find_len = strlen (str_to_find);


  for (i = 1; i < argc; i++)
    {
      option_value = strstr (argv[i], str_to_find);

      if (option_value)
        break;
    }

  if (!option_value || strlen (option_value) == 0)
    return NULL;

  option_value += str_to_find_len;

  return option_value;
}

const char *parser_info_str (command_line_parser *parser, int error_code)
{
  FIX_UNUSED (parser);
  switch (error_code)
    {
    case -3: return parser_version_str (parser);
    case -2: return parser_help_str (parser);
    case -1: return "Internal Error";
    case 0: return "No error";
    case 1:  return "Wrong --solver";
    case 2:  return "Wrong --pressure";
    case 3:  return "Wrong --table_format";
    case 4:  return "Wrong --solver-mode";
    case 5:  return "Wrong --N";
    case 6:  return "Wrong --M1";
    case 7:  return "Wrong --M2";
    case 8:  return "Wrong --T";
    case 9:  return "Wrong --X1";
    case 10:  return "Wrong --X2";
    default: return "Unknown Error";
    }

  return "Unknown Error";
}

int parser_is_help_present (command_line_parser *parser, int argc, char *argv[])
{
  int i;
  FIX_UNUSED (parser);


  for (i = 1; i < argc; i++)
    {
      if (!strcmp ("--help", argv[i]))
        return 1;
    }

  return 0;
}

int parser_is_version_present (command_line_parser *parser, int argc, char *argv[])
{
  int i;
  FIX_UNUSED (parser);


  for (i = 1; i < argc; i++)
    {
      if (!strcmp ("--version", argv[i]))
        return 1;
    }

  return 0;
}

const char *parser_help_str (command_line_parser *parser)
{
  /* dealing with -Woverlength-strings */
  if (!parser)
    return NULL;

  parser->help_str[0] = 0;

  strcat (parser->help_str,"--solver=[central, sokolov]     type=enum,   mandatory\n"
                           "--pressure=[linear, polynomial] type=enum,   optional, default=linear\n"
                           "--table-format=[simple, latex]  type=enum,   optional, default=simple\n"
                           "--solver-mode=[test, solve]     type=enum,   mandatory\n"
                           "--N=[3, 4, ...]                 type=int,    mandatory\n"
                           "--M1=[3, 4, ...]                type=int,    mandatory\n");

  strcat (parser->help_str, "--M2=[3, 4, ...]                type=int,    mandatory\n"
                            "--T=[double > 0]                type=double, mandatory\n"
                            "--X1=[any double]               type=double, mandatory\n"
                            "--X2=[any double]               type=double, mandatory");

  return parser->help_str;
}

const char *parser_version_str (command_line_parser *parser)
{
  sprintf (parser->version_str, "%s v%d.00", PROGRAM_NAME, PROGRAM_VERSION_NUM);
  return parser->version_str;
}

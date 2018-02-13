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

  solver_value       = get_value ("solver",       argc, argv);
  pressure_value     = get_value ("pressure",     argc, argv);
  table_format_value = get_value ("table-format", argc, argv);
  mode_value         = get_value ("solver-mode",  argc, argv);
  N_value            = get_value ("N",            argc, argv);
  M1_value           = get_value ("M1",           argc, argv);
  M2_value           = get_value ("M2",           argc, argv);
  T_value            = get_value ("T",            argc, argv);
  X1_value           = get_value ("X1",           argc, argv);
  X2_value           = get_value ("X2",           argc, argv);

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
    parser->p_func = linear;
  else
    {
      if (!strcmp ("linear", pressure_value))
        parser->p_func = linear;
      else
        {
          if (!strcmp ("polynomial", pressure_value))
            parser->p_func = polynomial;
          else
            return 2;
        }
    }

  if (!table_format_value)
    parser->table_output_format = human_readable;
  else
    {
      if (!strcmp ("latex", table_format_value))
        parser->table_output_format = latex;
      else
        {
          if (!strcmp ("simple", table_format_value))
            parser->table_output_format = human_readable;
          else
            return 3;
        }
    }

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

const char *get_value (const char *option, int argc, char *argv[])
{
  char str_to_find[OPTION_BUFFER_LEN];
  int  str_to_find_len;
  const char *option_value = NULL;

  str_to_find[0] = 0;
  str_to_find_len = 0;

  strcat (str_to_find, "--");
  strcat (str_to_find, option);
  strcat (str_to_find, "=");
  str_to_find_len = strlen (str_to_find);

  for (int i = 1; i < argc; i++)
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

const char *parser_error_str (command_line_parser *parser, int error_code)
{
  FIX_UNUSED (parser);
  switch (error_code)
    {
    case -1: return "Internal Error";
    case 0: return "No error";
    case 1:  return "Wrong --solver";
    case 2:  return "Wrong --pressure";
    case 3:  return "Wrong --table_format";
    case 4:  return "Wrong --N";
    case 5:  return "Wrong --M1";
    case 6:  return "Wrong --M2";
    case 7:  return "Wrong --T";
    case 8:  return "Wrong --X1";
    case 9:  return "Wrong --X2";
    default: return "Unknown Error";
    }

  return "Unknown Error";
}

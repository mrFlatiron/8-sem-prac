#include "command_line_parser_private.h"
#include "common/debug_utils.h"

#include <string.h>
#include <stdlib.h>
#include <sys/sysinfo.h>


int parse_command_line (command_line_parser *parser, int argc, char *argv[])
{
  const char *solver_value;
  const char *pressure_value;
  const char *mode_value;
  const char *N_value;
  const char *M1_value;
  const char *M2_value;
  const char *T_value;
  const char *omega_value;
  const char *mu_value;
  const char *linear_solver_value;


  if (!parser)
    return -1;

  if (argc == 1)
    return -2;

  solver_value       = parser_get_value (parser, "solver",       argc, argv);
  pressure_value     = parser_get_value (parser, "pressure",     argc, argv);
  mode_value         = parser_get_value (parser, "solver-mode",  argc, argv);
  N_value            = parser_get_value (parser, "N",            argc, argv);
  M1_value           = parser_get_value (parser, "MX",           argc, argv);
  M2_value           = parser_get_value (parser, "MY",           argc, argv);
  T_value            = parser_get_value (parser, "T",            argc, argv);
  omega_value        = parser_get_value (parser, "border-omega", argc, argv);
  mu_value           = parser_get_value (parser, "mu",           argc, argv);
  linear_solver_value= parser_get_value (parser, "linear-solver", argc, argv);


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


  if (!mode_value)
    return 4;

  if (!strcmp ("test", mode_value))
    parser->solver_mode = test_mode;
  else
    {
      if (!strcmp ("solve", mode_value))
        parser->solver_mode = solve_mode;
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

  if (!omega_value)
    return 11;

  parser->border_omega = atof (omega_value);

  if (!mu_value)
    return 12;

  parser->mu = atof (mu_value);

  if (!linear_solver_value)
    parser->linear_solver = custom_cgs;
  else
    {
      if (!strcmp ("laspack-cgs", linear_solver_value))
        parser->linear_solver = laspack_cgs;
      else
        {
          if (!strcmp ("custom-cgs", linear_solver_value))
            parser->linear_solver = custom_cgs;
          else
            return 13;
        }
}
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
    case 4:  return "Wrong --solver-mode";
    case 5:  return "Wrong --N";
    case 6:  return "Wrong --MX";
    case 7:  return "Wrong --MY";
    case 8:  return "Wrong --T";
    case 11: return "Wrong --border-omega";
    case 12: return "Wrong --mu";
    case 13: return "Wrong --linear-solver";
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

  char temp_str[HELP_STR_BUFFER_LEN];

  temp_str[0] = 0;

  if (!parser)
    return NULL;

  parser->help_str[0] = 0;

  strcat (parser->help_str,"--solver=[central, sokolov]     type=enum,   mandatory\n"
                           "--pressure=[linear, polynomial] type=enum,   optional, default=linear\n"
                           "--table-format=[simple, latex]  type=enum,   optional, default=simple\n"
                           "--solver-mode=[test, solve]     type=enum,   mandatory\n"
                           "--linear-solver=[laspack-cgs, custom-cgs] type=enum, optional, default=custom-cgs\n"
                           "--N=[3, 4, ...]                 type=int,    mandatory\n"
                           "--M1=[3, 4, ...]                type=int,    mandatory\n");

  strcat (temp_str,        "--M2=[3, 4, ...]                type=int,    mandatory\n"
                           "--T=[double > 0]                type=double, mandatory\n"
                           "--border-omega=[double > 0]     type=double, mandatory\n"
                           "--mu=[double > 0]               type=double, mandatory");

  strcat (parser->help_str, temp_str);

  return parser->help_str;
}

const char *parser_version_str (command_line_parser *parser)
{
  sprintf (parser->version_str, "%s v%d.00", PROGRAM_NAME, PROGRAM_VERSION_NUM);
  return parser->version_str;
}

#ifndef KERNEL_TYPEDEFS_H
#define KERNEL_TYPEDEFS_H

#include <stdio.h>

#define PROGRAM_VERSION_NUM 1
#define PROGRAM_NAME "8-sem-prac"

typedef enum
{
  pressure_linear,
  pressure_polynomial
} pressure_func_t;

typedef enum
{
  central_differences,
  sokolov
} solver_t;

typedef enum
{
  latex_format,
  human_readable
} table_output_format_t; /*for test mode*/

typedef enum
{
  test_mode, /*compare with known smooth solution*/
  solve_mode /*solve actual problem with zero rhs*/
} solver_mode_t;


#endif /* KERNEL_TYPEDEFS_H */

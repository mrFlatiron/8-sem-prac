#ifndef KERNEL_TYPEDEFS_H
#define KERNEL_TYPEDEFS_H

#include <stdio.h>

#define PROGRAM_VERSION_NUM 2
#define PROGRAM_NAME "8-sem-prac"

#define MY_PI                   3.14159265358979323846

#define RHO_LEFTMOST            1.0
#define BORDER_OMEGA            1.
#define X_LEN                   MY_PI
#define Y_LEN                   MY_PI

#define TOP_ROW_SQUARES_COUNT   1
#define BOT_ROW_SQUARES_COUNT   3
#define SQUARES_COUNT           (BOT_ROW_SQUARES_COUNT + TOP_ROW_SQUARES_COUNT)
#define DIMENSIONS              2
#define DIMENSIONS_W_TIME       DIMENSIONS + 1
#define UNKNOWN_FUNCTIONS_COUNT 3
#define MAX_ROW_NZ              9

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
  human_readable,
  gnuplot_xyz
} table_output_format_t; /*for test mode*/

typedef enum
{
  test_mode, /*compare with known smooth solution*/
  solve_mode /*solve actual problem with zero rhs*/
} solver_mode_t;

typedef enum
{
  laspack_cgs,
  custom_cgs
} linear_solver_t;

typedef enum
{
  precond_none,
  precond_jacobi
} preconditioner_t;

typedef enum
{
  border_leftmost,
  border_botmost,
  border_rightmost,
  border_topmost,
  border_left_hor,
  border_left_top,
  border_right_hor,
  border_right_top,
  area_internal

} grid_area_t;

typedef enum
{
  grid_g,
  grid_vx,
  grid_vy
} grid_func_t;

typedef double (*time_layer_func_t)  (double/*t*/, double/*x*/, double/*y*/);
typedef double (*layer_func_t) (double/*x*/, double/*y*/, double /*border_omega*/);
typedef double (*rhs_func_t) (double/*t*/, double/*x*/, double /*y*/, double /*mu*/, pressure_func_t);

#endif /* KERNEL_TYPEDEFS_H */

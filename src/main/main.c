#include <stdio.h>
#include <locale.h>

#include "common/vectors.h"
#include "common/debug_utils.h"
#include "kernel/solver_core_workspace.h"
#include "kernel/command_line_parser.h"
#include "kernel/cgs_solver.h"
#include "linear_ops/vector_ops.h"
#include "kernel/central_diff_solver.h"
#include "sparse/laspack_matrix.h"
#include "sparse/laspack_vector.h"
#include "io/table_io.h"
#include "kernel/solver_tester.h"
#include "3rd_party/laspack/itersolv.h"
#include "io/gnuplot_io.h"
#include <string.h>
#include "sys/stat.h"

#include "kernel/input/rhs.h"
#include "kernel/input/t0_functions.h"
#include "kernel/input/test_solutions.h"


int main (int argc, char *argv[])
{
/*Real main*/
  command_line_parser     parser_object;
  command_line_parser_ptr parser         = &parser_object;

  central_diff_solver      solver_object;
  central_diff_solver      *solver = &solver_object;

  solver_tester tester_obj;
  solver_tester *tester = &tester_obj;

  table_io table_obj;
  table_io *table = &table_obj;

  int                     error_code;

  setlocale (LC_ALL, "en-GB.utf8");

  error_code = parse_command_line (parser, argc, argv);

  if (error_code)
    {
      fprintf (stderr, "%s\n", parser_info_str (parser, error_code));
      return 1;
    }

 if (parser->solver_mode == test_mode)
   {

     solver_tester_init (tester,
                         parser->solver,
                         parser->N,
                         parser->MX,
                         parser->MY,
                         parser->N_mult,
                         parser->MXY_mult,
                         parser->N_mult_count,
                         parser->MXY_mult_count,
                         parser->border_omega,
                         parser->mu,
                         X_LEN, Y_LEN, parser->T,
                         test_g, test_vx, test_vy,
                         rhs_test_f0, rhs_test_f1, rhs_test_f2,
                         t0_vx_test, t0_vy_test, t0_g_test);

     solver_tester_test (tester,
                         parser->solver_precision,
                         parser->solver_max_iter,
                         parser->precond);

     solver_tester_print_results (tester, stdout);

     solver_tester_destroy (tester);
   }
 else
   {
     gnuplot_io gp_handle_obj;
     gnuplot_io *gp_io = &gp_handle_obj;

     int n;



     error_code = cdiff_solver_init (solver,
                                     parser->solver_mode,
                                     parser->MX,
                                     parser->MY,
                                     parser->N,
                                     X_LEN,
                                     Y_LEN,
                                     parser->T,
                                     parser->border_omega,
                                     parser->mu,
                                     parser->solver_precision,
                                     parser->solver_max_iter,
                                     parser->precond,
                                     parser->linear_solver);

     if (error_code)
       {
         fprintf (stderr, "Could not initialize solver with given arguments\n");
         return error_code;
       }

     cdiff_solver_compute (solver,
                           parser->p_func,
                           NULL,
                           NULL,
                           NULL,
                           rhs_f0,
                           rhs_f1,
                           rhs_f2,
                           t0_vx,
                           t0_vy,
                           t0_g);

     for (n = 0; n <= parser->N; n++)
       {

         char path_g[4096];
         char path_v[4096];
         FILE *gnu_out_g;
         FILE *gnu_out_v;

         if (n != 0 && n != parser->N)
           continue;

         mkdir ("out", S_IRWXU | S_IRWXG);
         mkdir ("out/gnuplot", S_IRWXU | S_IRWXG);

         mkdir ("out/gnuplot/g", S_IRWXU | S_IRWXG);
         mkdir ("out/gnuplot/v", S_IRWXU | S_IRWXG);

         sprintf (path_g, "out/gnuplot/g/%d", n);
         sprintf (path_v, "out/gnuplot/v/%d", n);
         gnu_out_g = fopen (path_g, "w");
         gnu_out_v = fopen (path_v, "w");

         if (!gnu_out_g || !gnu_out_v)
           {
             fprintf (stderr, "Couldn't open out file");
             return -1;
           }

         gnuplot_io_init (gp_io, &solver->ws, n);

         table_io_init (table, solver->ws.layer_size, 3, gp_io->coords_g, gnuplot_xyz);
         fprintf (gnu_out_g, table->table_text);
         table_io_destroy (table);

         table_io_init (table, solver->ws.layer_size, 4, gp_io->coords_v, gnuplot_xyz);
         fprintf (gnu_out_v, table->table_text);
         table_io_destroy (table);

         fclose (gnu_out_g);
         fclose (gnu_out_v);
         gnuplot_io_destroy (gp_io);
       }
     cdiff_solver_destroy (solver);
   }
#if 0 /*msr matrix from vector and cgs solver test*/
  double dense[] = {1, 1, 1, 0, 0,
                    6, 2, 3, 0, 0,
                    0, 2, 5, 8, 0,
                    0, 0, 4, 6, 8,
                    0, 0, 0, -1, 7};
  msr_matrix sparse_object;
  msr_matrix* sparse = &sparse_object;

  cgs_solver solver_object;
  cgs_solver *solver = &solver_object;

  double x[] = {1, 0, 1, 0, 1};
  double found_x[5];

/*  double x_init[] = {0, 0, 0, 0, 0};*/

  double rhs[5];

  cgs_solver_error_t error;



  FIX_UNUSED (argc);
  FIX_UNUSED (argv);
  FIX_UNUSED (error);


  msr_init_from_vector (sparse, dense, 5);

  msr_dump (sparse, stdout);

  msr_mult_vector (sparse, x, rhs);

  cgs_solver_init (solver,
                   sparse->matrix_size,
                   300,
                   1e-6,
                   precond_jacobi);


  error = cgs_solver_solve (solver,
                    sparse,
                    rhs,
                    NULL,
                    found_x);

  cgs_solver_destroy (solver);
  msr_destroy (sparse);
#endif
#if 0 /*sparse_base_format test*/
  sparse_base_format sparse_base_obj;
  sparse_base_format *sparse_base = &sparse_base_obj;
  msr_matrix sparse_object;
  msr_matrix* sparse = &sparse_object;

  double vals[3];
  int cols[3];
  int nnz;

  FIX_UNUSED (argc);
  FIX_UNUSED (argv);

  sparse_base_init (sparse_base,  5, 3);
  msr_init_empty (sparse);

  vals[0] = 1; vals[1] = 1;
  cols[0] = 0; cols[1] = 1;
  nnz = 2;

  sparse_base_add_row (sparse_base, 0, cols, vals, nnz);

  vals[0] = 2;
  cols[0] = 1;
  nnz = 1;

  sparse_base_add_row (sparse_base, 1, cols, vals, nnz);

  vals[0] = 2; vals[1] = 5;
  cols[0] = 1; cols[1] = 2;
  nnz = 2;

  sparse_base_add_row (sparse_base, 2, cols, vals, nnz);

  vals[0] = 4; vals[1] = 6; vals[2] = 8;
  cols[0] = 2; cols[1] = 3; cols[2] = 4;
  nnz = 3;

  sparse_base_add_row (sparse_base, 3, cols, vals, nnz);

  vals[0] = 7;
  cols[0] = 4;
  nnz = 1;

  sparse_base_add_row (sparse_base, 4, cols, vals, nnz);

  msr_fill_from_sparse_base (sparse, sparse_base);

  msr_dump (sparse, stdout);

  sparse_base_destroy (sparse_base);
  msr_destroy (sparse);
#endif
#if 0/*laspack_matrix test*/
  sparse_base_format sparse_base_obj;
  sparse_base_format *sparse_base = &sparse_base_obj;
  laspack_matrix matobj;
  laspack_matrix *mat = &matobj;
  laspack_vector vecobj;
  laspack_vector *rhs_l = &vecobj;
  laspack_vector x_obj_l;
  laspack_vector *x_l = &x_obj_l;

  double rhs[] = {2, 9, 5, 12, 7};

  double vals[3];
  int cols[3];
  int nnz;

  FIX_UNUSED (argc);
  FIX_UNUSED (argv);

  sparse_base_init (sparse_base,  5, 3);

  vals[0] = 1; vals[1] = 1; vals[2] = 1;
  cols[0] = 0; cols[1] = 1; cols[2] = 2;
  nnz = 3;

  sparse_base_add_row (sparse_base, 0, cols, vals, nnz);

  vals[0] = 6; vals[1] = 2; vals[2] = 3;
  cols[0] = 0; cols[1] = 1; cols[2] = 2;
  nnz = 3;

  sparse_base_add_row (sparse_base, 1, cols, vals, nnz);

  vals[0] = 2; vals[1] = 5; vals[2] = 8;
  cols[0] = 1; cols[1] = 2; cols[2] = 3;
  nnz = 3;

  sparse_base_add_row (sparse_base, 2, cols, vals, nnz);

  vals[0] = 4; vals[1] = 6; vals[2] = 8;
  cols[0] = 2; cols[1] = 3; cols[2] = 4;
  nnz = 3;

  sparse_base_add_row (sparse_base, 3, cols, vals, nnz);

  vals[0] = -1; vals[1] = 7;
  cols[0] = 3; cols[1] = 4;
  nnz = 2;

  sparse_base_add_row (sparse_base, 4, cols, vals, nnz);

  laspack_matrix_init (mat, sparse_base);

  laspack_vector_init (rhs_l, 5);
  laspack_vector_fill (rhs_l, rhs);

  laspack_vector_init (x_l, 5);

  laspack_matrix_dump (mat, stdout);
  fprintf (stdout, "\n");

  CGSIter (&mat->raw, &x_l->raw, &rhs_l->raw, 2000, NULL, 1.2);

  laspack_vector_dump (rhs_l, stdout);
  fprintf (stdout, "\n");

  laspack_vector_dump (x_l, stdout);
  fprintf (stdout, "\n");

  laspack_vector_destroy (x_l);
  laspack_vector_destroy (rhs_l);
  laspack_matrix_destroy (mat);
  sparse_base_destroy (sparse_base);
#endif
#if 0/*mapping test*/

  solver_core_workspace ws_obj;
  solver_core_workspace *ws = &ws_obj;

  int mx, my, index;
  int test_mx, test_my;

  FIX_UNUSED (argc);
  FIX_UNUSED (argv);

  solver_workspace_data_init (ws, test_mode, parser->N_mult_count, parser->MXY_mult_count, 3, X_LEN, Y_LEN, 1.5, 1);

  my = 0;
  mx = 5;
  index = 5;
  DEBUG_ASSERT (solver_workspace_final_index (ws, 0, mx, my) == index);
  solver_workspace_get_mx_my (ws, index, &test_mx, &test_my);
  DEBUG_ASSERT (test_mx == mx);
  DEBUG_ASSERT (test_my == my);

  mx = 9;
  my = 1;
  index = 19;
  DEBUG_ASSERT (solver_workspace_final_index (ws, 0, mx, my) == index);
  solver_workspace_get_mx_my (ws, index, &test_mx, &test_my);
  DEBUG_ASSERT (test_mx == mx);
  DEBUG_ASSERT (test_my == my);

  mx = 0;
  my = 1;
  index = 10;
  DEBUG_ASSERT (solver_workspace_final_index (ws, 0, mx, my) == index);
  solver_workspace_get_mx_my (ws, index, &test_mx, &test_my);
  DEBUG_ASSERT (test_mx == mx);
  DEBUG_ASSERT (test_my == my);

  mx = 0;
  my = 3;
  index = 30;
  DEBUG_ASSERT (solver_workspace_final_index (ws, 0, mx, my) == index);
  solver_workspace_get_mx_my (ws, index, &test_mx, &test_my);
  DEBUG_ASSERT (test_mx == mx);
  DEBUG_ASSERT (test_my == my);

  mx = 3;
  my = 3;
  index = 33;
  DEBUG_ASSERT (solver_workspace_final_index (ws, 0, mx, my) == index);
  solver_workspace_get_mx_my (ws, index, &test_mx, &test_my);
  DEBUG_ASSERT (test_mx == mx);
  DEBUG_ASSERT (test_my == my);

  mx = 6;
  my = 3;
  index = 36;
  DEBUG_ASSERT (solver_workspace_final_index (ws, 0, mx, my) == index);
  solver_workspace_get_mx_my (ws, index, &test_mx, &test_my);
  DEBUG_ASSERT (test_mx == mx);
  DEBUG_ASSERT (test_my == my);

  mx = 3;
  my = 4;
  index = 40;
  DEBUG_ASSERT (solver_workspace_final_index (ws, 0, mx, my) == index);
  solver_workspace_get_mx_my (ws, index, &test_mx, &test_my);
  DEBUG_ASSERT (test_mx == mx);
  DEBUG_ASSERT (test_my == my);

  mx = 6;
  my = 4;
  index = 43;
  DEBUG_ASSERT (solver_workspace_final_index (ws, 0, mx, my) == index);
  solver_workspace_get_mx_my (ws, index, &test_mx, &test_my);
  DEBUG_ASSERT (test_mx == mx);
  DEBUG_ASSERT (test_my == my);

  mx = 3;
  my = 5;
  index = 44;
  DEBUG_ASSERT (solver_workspace_final_index (ws, 0, mx, my) == index);
  solver_workspace_get_mx_my (ws, index, &test_mx, &test_my);
  DEBUG_ASSERT (test_mx == mx);
  DEBUG_ASSERT (test_my == my);

  mx = 4;
  my = 5;
  index = 45;
  DEBUG_ASSERT (solver_workspace_final_index (ws, 0, mx, my) == index);
  solver_workspace_get_mx_my (ws, index, &test_mx, &test_my);
  DEBUG_ASSERT (test_mx == mx);
  DEBUG_ASSERT (test_my == my);

  mx = 3;
  my = 6;
  index = 48;
  DEBUG_ASSERT (solver_workspace_final_index (ws, 0, mx, my) == index);
  solver_workspace_get_mx_my (ws, index, &test_mx, &test_my);
  DEBUG_ASSERT (test_mx == mx);
  DEBUG_ASSERT (test_my == my);

  mx = 6;
  my = 6;
  index = 51;
  DEBUG_ASSERT (solver_workspace_final_index (ws, 0, mx, my) == index);
  solver_workspace_get_mx_my (ws, index, &test_mx, &test_my);
  DEBUG_ASSERT (test_mx == mx);
  DEBUG_ASSERT (test_my == my);


  solver_workspace_data_destroy (ws);
#endif
#if 0/*table test*/
  const char *entries[] = {"one", "super two", "three", "four"};
  table_io table_obj;
  table_io *table = &table_obj;
  int must_be_len_human = strlen ("|  one  | super two |\n"
                                  "| three |    four   |");
  int must_be_len_latex = strlen ("begin{tabular}{|c|c|}\n"
                                  "\\hline\n"
                                     "one  & super two \\\\\n"
                                  "\\hline\n"
                                    "three &   four    \\\\\n"
                                  "\\hline\n"
                                  "end{tabular}");

  FIX_UNUSED (argc);
  FIX_UNUSED (argv);

  table_io_init (table, 2, 2, entries, human_readable);
  DEBUG_ASSERT (must_be_len_human <= table->length);
  fprintf (stdout, table->table_text);
  fprintf (stdout, "\n");
  table_io_destroy (table);

  table_io_init (table, 2, 2, entries, latex_format);
  fprintf (stdout, table->table_text);
  fprintf (stdout, "\n");
  DEBUG_ASSERT (must_be_len_latex <= table->length);
  table_io_destroy (table);
#endif
  return 0;
}

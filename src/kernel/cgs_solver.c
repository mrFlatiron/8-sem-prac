#include "cgs_solver.h"
#include "sparse/msr_matrix.h"
#include "common/math_utils.h"
#include "common/vectors.h"
#include "linear_ops/vector_ops.h"

#include <stdlib.h>



int cgs_solver_init (cgs_solver *solver,
                     int N,
                     int max_iter,
                     double exit_precision,
                     preconditioner_t precond)
{
  solver->iter = 0;
  solver->N = N;
  solver->max_iter = max_iter;
  solver->precond = precond;
  solver->preconditioned_matrix = (msr_matrix *) malloc (sizeof (msr_matrix));
  solver->precision = exit_precision;

  solver->preconditioned_rhs = VECTOR_CREATE (double, N);
  solver->residual_vec       = VECTOR_CREATE (double, N);
  solver->x_vec              = VECTOR_CREATE (double, N);
  solver->r_star_vec         = VECTOR_CREATE (double, N);
  solver->p_vec              = VECTOR_CREATE (double, N);
  solver->u_vec              = VECTOR_CREATE (double, N);
  solver->q_vec              = VECTOR_CREATE (double, N);
  solver->temp1_vec          = VECTOR_CREATE (double, N);
  solver->temp2_vec          = VECTOR_CREATE (double, N);

  return 0;
}

void cgs_solver_destroy (cgs_solver *solver)
{
  VECTOR_DESTROY (solver->preconditioned_rhs);
  VECTOR_DESTROY (solver->residual_vec);
  VECTOR_DESTROY (solver->x_vec);
  VECTOR_DESTROY (solver->r_star_vec);
  VECTOR_DESTROY (solver->p_vec);
  VECTOR_DESTROY (solver->u_vec);
  VECTOR_DESTROY (solver->q_vec);
  VECTOR_DESTROY (solver->temp1_vec);
  VECTOR_DESTROY (solver->temp2_vec);

  if (solver->preconditioned_matrix)
    {
      msr_destroy (solver->preconditioned_matrix);
      free (solver->preconditioned_matrix);
    }

}

cgs_solver_error_t cgs_solver_solve (cgs_solver *solver,
                                     const msr_matrix *matrix,
                                     const vector_double_t rhs,
                                     const vector_double_t init_x,
                                     vector_double_t out_x)
{
  cgs_solver_do_initialization (solver, matrix, rhs, init_x);

  if (solver->error_code)
    return solver->error_code;

  for (solver->iter = 0; solver->iter < solver->max_iter; solver->iter++)
    {
      cgs_solver_do_iter (solver);

      if (solver->error_code != cgs_error_ok)
        break;

      if (cgs_solver_check_converged (solver))
        break;
    }

  if (solver->iter == solver->max_iter)
    solver->error_code = cgs_error_max_iter_exceeded;

  if (solver->error_code == cgs_error_ok)
    VECTOR_COPY (double, solver->x_vec, out_x, solver->N);

  return solver->error_code;
}

void cgs_solver_apply_preconditioner (cgs_solver *solver)
{
  int i;
  msr_matrix *out_mat = solver->preconditioned_matrix;
  vector_double_t out_vec = solver->preconditioned_rhs;
  int N = solver->N;

  switch (solver->precond)
    {
    case precond_none:
      return;
    case precond_jacobi:
      {
        for (i = 0; i < N; i++)
          {
            int row_begin;
            int row_end;
            int ja_iter;
            double diag_val = out_mat->AA[i];

            if (math_is_null (diag_val))
              {
                solver->error_code = cgs_error_zero_on_diag;
                return;
              }

            out_mat->AA[i] = 1;

            msr_get_ja_row_bounds (out_mat, i, &row_begin, &row_end);

            for (ja_iter = row_begin; ja_iter < row_end; ja_iter++)
              out_mat->AA[ja_iter] /= diag_val;

            out_vec[i] /= diag_val;
          }
      }
    }
}

void cgs_solver_do_initialization (cgs_solver *solver,
                                      const msr_matrix *matrix,
                                      const vector_double_t rhs,
                                      const vector_double_t init_x)
{
  solver->error_code = cgs_error_ok;

  msr_init (solver->preconditioned_matrix, solver->N, matrix->array_size - 1 /*msr storage specific*/);

  VECTOR_COPY (double, matrix->AA, solver->preconditioned_matrix->AA, matrix->array_size);
  VECTOR_COPY (int,    matrix->JA, solver->preconditioned_matrix->JA, matrix->array_size);

  VECTOR_COPY (double, rhs, solver->preconditioned_rhs, solver->N);
  VECTOR_COPY (double, init_x, solver->x_vec, solver->N);

  cgs_solver_apply_preconditioner (solver);

  if (solver->error_code)
    return;

  VECTOR_SET (double, solver->r_star_vec, 1. / (double) solver->N, solver->N);

  msr_mult_vector (solver->preconditioned_matrix, solver->x_vec, solver->temp1_vec);

  linear_combination_1 (solver->preconditioned_rhs,
                        -1, solver->temp1_vec,
                        solver->residual_vec,
                        solver->N);

  VECTOR_COPY (double, solver->residual_vec, solver->p_vec, solver->N);
  VECTOR_COPY (double, solver->residual_vec, solver->u_vec, solver->N);
}

void cgs_solver_do_iter (cgs_solver *solver)
{
  double enumerator;
  double denominator;

  msr_mult_vector (solver->preconditioned_matrix, solver->p_vec, solver->temp1_vec);

  enumerator = dot_product (solver->residual_vec, solver->r_star_vec, solver->N);
  denominator = dot_product (solver->temp1_vec, solver->r_star_vec, solver->N);

  if (math_is_null (denominator))
    {
      solver->error_code = cgs_error_unknown;
      return;
    }

  solver->alpha = enumerator / denominator;

  linear_combination_1 (solver->u_vec,
                        -solver->alpha, solver->temp1_vec,
                        solver->q_vec,
                        solver->N);

  linear_combination_1 (solver->u_vec,
                        1, solver->q_vec,
                        solver->temp1_vec,
                        solver->N);

  linear_combination_w_override_1 (solver->x_vec,
                                   solver->alpha, solver->temp1_vec,
                                   solver->N);

  msr_mult_vector (solver->preconditioned_matrix, solver->temp1_vec, solver->temp2_vec);

  denominator = enumerator;

  if (math_is_null (denominator))
    {
      solver->error_code = cgs_error_unknown;
      return;
    }

  linear_combination_w_override_1 (solver->residual_vec,
                                   -solver->alpha, solver->temp2_vec,
                                   solver->N);

  enumerator = dot_product (solver->residual_vec, solver->r_star_vec, solver->N);

  solver->beta = enumerator / denominator;

  linear_combination_1 (solver->residual_vec,
                        solver->beta, solver->q_vec,
                        solver->u_vec,
                        solver->N);

  linear_combination_1 (solver->q_vec,
                        solver->beta, solver->p_vec,
                        solver->temp1_vec,
                        solver->N);

  linear_combination_1 (solver->u_vec,
                        solver->beta, solver->temp1_vec,
                        solver->p_vec,
                        solver->N);
}

int cgs_solver_check_converged (cgs_solver *solver)
{
  return l2_norm (solver->residual_vec, solver->N) <= solver->precision;
}

#ifndef CGS_SOLVER_H
#define CGS_SOLVER_H

#include "common/vectors_fwd.h"
#include "sparse/msr_matrix.h"

typedef enum
{
  precond_none,
  precond_jacobi
} preconditioner_t;

typedef enum
{
  cgs_error_ok,
  cgs_error_max_iter_exceeded,
  cgs_error_zero_on_diag,
  cgs_error_unknown
} cgs_solver_error_t;

typedef struct
{
  int              N;
  int              iter;
  int              max_iter;
  double           precision;

  msr_matrix*      preconditioned_matrix;
  vector_double_t  preconditioned_rhs;

  vector_double_t  r_star_vec;
  vector_double_t  p_vec;
  vector_double_t  u_vec;
  vector_double_t  q_vec;
  vector_double_t  temp1_vec;
  vector_double_t  temp2_vec;
  double alpha;
  double beta;

  vector_double_t  residual_vec;
  vector_double_t  x_vec;



  preconditioner_t precond;

  cgs_solver_error_t error_code;
} cgs_solver;

int cgs_solver_init (cgs_solver *solver,
                     int N,
                     int max_iter,
                     double exit_precision,
                     preconditioner_t precond
                     );

void cgs_solver_destroy (cgs_solver *solver);

cgs_solver_error_t cgs_solver_solve (cgs_solver *solver,
                                     const msr_matrix *matrix,
                                     const vector_double_t rhs,
                                     const vector_double_t init_x,
                                     vector_double_t out_x);

/* Below are private methods */

void cgs_solver_do_initialization (cgs_solver *solver,
                                      const msr_matrix *matrix,
                                      const vector_double_t rhs,
                                      const vector_double_t init_x);

void cgs_solver_apply_preconditioner (cgs_solver *solver);

void cgs_solver_do_iter (cgs_solver *solver);

int cgs_solver_check_converged (cgs_solver *solver);


#endif /* CGS_SOLVER_H */
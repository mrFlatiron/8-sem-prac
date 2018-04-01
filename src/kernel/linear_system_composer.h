#ifndef LINEAR_SYSTEM_COMPOSER_H
#define LINEAR_SYSTEM_COMPOSER_H

#include "common/vectors_fwd.h"
#include "kernel_typedefs.h"
#include "sparse/msr_matrix.h"
#include "sparse/sparse_base_format.h"
#include "cgs_solver.h"
#include "sparse/laspack_matrix.h"
#include "sparse/laspack_vector.h"

typedef struct
{
  linear_solver_t linear_solver;

  int nodes_count;
  int functions_count;

  int matrix_size;

  double precision;
  int max_iter;
  preconditioner_t precond;

  vector_double_t vector_to_compute;
  vector_double_t rhs_vector;

  laspack_vector rhs_vector_l_obj;
  laspack_vector  *rhs_vector_l;

  laspack_vector vector_to_compute_l_obj;
  laspack_vector  *vector_to_compute_l;

  sparse_base_format matrix_base_obj;
  sparse_base_format *matrix_base;

  msr_matrix matrix_obj;
  msr_matrix *matrix;

  laspack_matrix matrix_l_obj;
  laspack_matrix *matrix_l;

  cgs_solver cgs_linear_solver_obj;
  cgs_solver *cgs_linear_solver;

} linear_system_composer;

void system_composer_init (linear_system_composer *comp,
                           linear_solver_t solver,
                           int nodes_count,
                           int functions_count,
                           double prec,
                           int max_iter,
                           preconditioner_t precond);

void system_composer_destroy (linear_system_composer *comp);

void system_composer_solve (linear_system_composer *comp);

void system_composer_set_rhs_val (linear_system_composer *comp, double rhs, int row);

void system_composer_fill_nodes_values (const linear_system_composer *comp, vector_double_t *functions_as_vectors);

#endif /* LINEAR_SYSTEM_COMPOSER_H */

#include "linear_system_composer.h"
#include "common/vectors.h"
#include "3rd_party/laspack/itersolv.h"
#include "common/debug_utils.h"

void system_composer_init (linear_system_composer *comp,
                           linear_solver_t solver,
                           int nodes_count,
                           int functions_count,
                           double prec,
                           int max_iter,
                           preconditioner_t precond)
{
  comp->linear_solver = solver;
  comp->nodes_count = nodes_count;
  comp->functions_count = functions_count;
  comp->precision = prec;
  comp->max_iter = max_iter;
  comp->precond = precond;

  comp->rhs_vector_l = &comp->rhs_vector_l_obj;
  comp->vector_to_compute_l = &comp->vector_to_compute_l_obj;
  comp->matrix = &comp->matrix_obj;
  comp->matrix_base = &comp->matrix_base_obj;
  comp->matrix_l = &comp->matrix_l_obj;
  comp->cgs_linear_solver = &comp->cgs_linear_solver_obj;

  comp->matrix_size = functions_count * nodes_count;

  if (comp->linear_solver == laspack_cgs)
    {
      laspack_vector_init (comp->rhs_vector_l, comp->matrix_size);
      laspack_vector_init (comp->vector_to_compute_l, comp->matrix_size);
    }
  else
    {
      comp->rhs_vector = VECTOR_CREATE (double , comp->matrix_size);
      comp->vector_to_compute = VECTOR_CREATE (double, comp->matrix_size);
      msr_init_empty (comp->matrix);
      cgs_solver_init (comp->cgs_linear_solver, comp->matrix_size, comp->max_iter, comp->precision, comp->precond);
    }

  sparse_base_init (comp->matrix_base, comp->matrix_size, MAX_ROW_NZ);
}

void system_composer_destroy (linear_system_composer *comp)
{
  if (comp->linear_solver == laspack_cgs)
    {
      laspack_vector_destroy (comp->rhs_vector_l);
      laspack_vector_destroy (comp->vector_to_compute_l);
    }
  else
    {
      VECTOR_DESTROY (comp->rhs_vector);
      VECTOR_DESTROY (comp->vector_to_compute);
      msr_destroy (comp->matrix);
      sparse_base_destroy (comp->matrix_base);
      cgs_solver_destroy (comp->cgs_linear_solver);
    }
}

void system_composer_fill_nodes_values (const linear_system_composer *comp, vector_double_t *functions_as_vectors)
{
  int i;
  int j;
  for (i = 0; i < comp->nodes_count; i++)
    for (j = 0; j < comp->functions_count; j++)
      functions_as_vectors[j][i] = comp->vector_to_compute[comp->functions_count * i + j];
}

void system_composer_set_rhs_val (linear_system_composer *comp, double rhs, int row)
{
  if (comp->linear_solver == laspack_cgs)
    laspack_vector_set_el (comp->rhs_vector_l, rhs, row);
  else
    comp->rhs_vector[row] = rhs;
}

void system_composer_solve (linear_system_composer *comp)
{
  if (comp->linear_solver == laspack_cgs)
    {
      laspack_matrix_init (comp->matrix_l, comp->matrix_base);
      CGSIter (&comp->matrix_l->raw, &comp->vector_to_compute_l->raw, &comp->rhs_vector_l->raw,
               comp->max_iter, JacobiPrecond, 1);
      laspack_matrix_destroy (comp->matrix_l);
    }
  else
    {
      cgs_solver_error_t error;
      msr_fill_from_sparse_base (comp->matrix, comp->matrix_base);
      error = cgs_solver_solve (comp->cgs_linear_solver, comp->matrix,
                        comp->rhs_vector, comp->vector_to_compute,
                        comp->vector_to_compute);

      DEBUG_ASSERT (error == cgs_error_ok);
    }
}

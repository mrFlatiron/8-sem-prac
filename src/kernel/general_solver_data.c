#include "general_solver_data.h"

#include <stddef.h>
#include "common/vectors.h"
#include "common/debug_utils.h"

int gen_solver_data_init (general_solver_data *solver,
                          pressure_func_t p_func,
                          int M1,
                          int M2,
                          int N,
                          double X1,
                          double X2,
                          double T)
{
  int success = 1;
  int vector_size;

  if (!solver)
    return 0;

  if (!gen_solver_data_check_input (M1, M2, N, X1, X2, T))
    return 1;

  solver->m_p_func = p_func;
  solver->m_M1 = M1;
  solver->m_M2 = M2;
  solver->m_N = N;
  solver->m_X1 = X1;
  solver->m_X2 = X2;
  solver->m_T = T;

  solver->m_vectors_size = (N + 1) * (M1 + 1) * (M2 + 1);
  vector_size = solver->m_vectors_size;

  success &= ((solver->m_g1 = VECTOR_CREATE (double, vector_size)) != NULL);
  success &= ((solver->m_g2 = VECTOR_CREATE (double, vector_size)) != NULL);
  success &= ((solver->m_v1 = VECTOR_CREATE (double, vector_size)) != NULL);
  success &= ((solver->m_v2 = VECTOR_CREATE (double, vector_size)) != NULL);

  return !success;
}

void gen_solver_data_destroy (general_solver_data *solver)
{
  if (!solver)
    return;

  VECTOR_DESTROY (solver->m_g1);
  VECTOR_DESTROY (solver->m_g2);
  VECTOR_DESTROY (solver->m_v1);
  VECTOR_DESTROY (solver->m_v2);
}


int gen_solver_data_check_input (int M1, int M2, int N, double X1, double X2, double T)
{
  return !(M1 <= 2 || M2 <= 2 || N <= 2
      || X1 <= 0 || X2 <= 0 || T <= 0);

}

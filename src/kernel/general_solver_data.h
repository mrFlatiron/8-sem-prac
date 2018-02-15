#ifndef GENERAL_SOLVER_DATA_H
#define GENERAL_SOLVER_DATA_H

#include "common/vectors_fwd.h"
#include "kernel_typedefs.h"

typedef struct
{
  pressure_func_t m_p_func;

  int m_M1;
  int m_M2;
  int m_N;
  double m_X1;
  double m_X2;
  double m_T;
  double m_h1;
  double m_h2;
  double m_tau;

  int m_vectors_size;

  vector_double_t m_v1;
  vector_double_t m_v2;
  vector_double_t m_g1;
  vector_double_t m_g2;

  int m_cpus_available;
} general_solver_data;


int gen_solver_data_init (general_solver_data *solver,
                          pressure_func_t p_func,
                          int M1,
                          int M2,
                          int N,
                          double X1,
                          double X2,
                          double T,
                          int cpus_available);

void gen_solver_data_destroy (general_solver_data *solver);

/*Below are private methods*/

int gen_solver_data_check_input (int M1,
                           int M2,
                           int N,
                           double X1,
                           double X2,
                           double T);

#endif /* GENERAL_SOLVER_DATA_H */

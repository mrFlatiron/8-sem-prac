#ifndef CENTRAL_DIFFERENCES_SOLVER_H
#define CENTRAL_DIFFERENCES_SOLVER_H

#include "common/vectors_fwd.h"
#include "kernel_typedefs.h"

struct central_differences_solver
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
};

typedef struct central_differences_solver cd_solver;
typedef cd_solver* cd_solver_ptr;



int cdsolver_init (cd_solver *solver,
                     pressure_func_t p_func,
                     int M1,
                     int M2,
                     int N,
                     double X1,
                     double X2,
                     double T);

void cdsolver_destroy (cd_solver *solver);

int cdsolver_compute (cd_solver *solver);


//Below are private methods

int cdsolver_check_input (int M1,
                           int M2,
                           int N,
                           double X1,
                           double X2,
                           double T);

#endif // CENTRAL_DIFFERENCES_SOLVER_H

#ifndef CGS_SOLVER_PRIVATE_H
#define CGS_SOLVER_PRIVATE_H

#include "cgs_solver.h"

void cgs_solver_do_initialization (cgs_solver *solver,
                                      const msr_matrix *matrix,
                                      const vector_double_t rhs,
                                      const vector_double_t init_x);

void cgs_solver_apply_preconditioner (cgs_solver *solver);

void cgs_solver_do_iter (cgs_solver *solver);

int cgs_solver_check_converged (cgs_solver *solver);

#endif /* CGS_SOLVER_PRIVATE_H */

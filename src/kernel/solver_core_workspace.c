#include "solver_core_workspace_private.h"

#include <stddef.h>
#include "common/vectors.h"
#include "common/debug_utils.h"

int solver_workspace_data_init (solver_core_workspace *data,
                          int M1,
                          int M2,
                          int N,
                          double X,
                          double Y,
                          double T)
{
  int success = 1;

  if (!data)
    return 0;

  if (!solver_workspace_data_check_input (M1, M2, N, X, Y, T))
    return 1;

  data->M1 = M1;
  data->M2 = M2;
  data->N = N;
  data->X = X;
  data->Y = Y;
  data->T = T;

  data->layer_size = SQUARES_COUNT * (M1 + 1) * (M2 + 1);

  data->vectors_size = (N + 1) * data->layer_size;

  data->matrix_size = UNKNOWN_FUNCTIONS_COUNT * data->layer_size;


  success &= ((data->g                 = VECTOR_CREATE (double, data->vectors_size)) != NULL);
  success &= ((data->vx                = VECTOR_CREATE (double, data->vectors_size)) != NULL);
  success &= ((data->vy                = VECTOR_CREATE (double, data->vectors_size)) != NULL);
  success &= ((data->vector_to_compute = VECTOR_CREATE (double, data->matrix_size)) != NULL);
  success &= ((data->rhs_vector        = VECTOR_CREATE (double, data->matrix_size)) != NULL);

  return !success;
}

void solver_workspace_data_destroy (solver_core_workspace *data)
{
  if (!data)
    return;

  VECTOR_DESTROY (data->g);
  VECTOR_DESTROY (data->vx);
  VECTOR_DESTROY (data->vy);
  VECTOR_DESTROY (data->vector_to_compute);
  VECTOR_DESTROY (data->rhs_vector);
}


int solver_workspace_data_check_input (int M1, int M2, int N, double X1, double X2, double T)
{
  return !(M1 <= 2 || M2 <= 2 || N <= 2
      || X1 <= 0 || X2 <= 0 || T <= 0);

}

int solver_workspace_layer_begin_index (const solver_core_workspace *data, int n)
{
  return n * data->layer_size;
}

int solver_workspace_final_index (const solver_core_workspace *data, int n, int mx, int my)
{
  return solver_workspace_layer_begin_index (data, n) + my * (data->X + 1) + mx;
}

double solver_workspace_grid_vx (const solver_core_workspace *data, int n, int mx, int my)
{
  return data->vx[solver_workspace_final_index (data, n, mx, my)];
}

double solver_workspace_grid_vy (const solver_core_workspace *data, int n, int mx, int my)
{
  return data->vy[solver_workspace_final_index (data, n, mx, my)];
}

double solver_workspace_grid_g (const solver_core_workspace *data, int n, int mx, int my)
{
  return data->g[solver_workspace_final_index (data, n, mx, my)];
}

void solver_workspace_fill_layer (solver_core_workspace *data, int n)
{
  int i = 0;
  int ig = 0;
  int ivx = 0;
  int ivy = 0;
  int size = data->layer_size;

  for (i = 0; i < size; i++)
    {
      int index = solver_workspace_layer_begin_index (data, n) + i;
      data->g[index]  = data->vector_to_compute[ig];
      data->vx[index] = data->vector_to_compute[ivx];
      data->vy[index] = data->vector_to_compute[ivy];
      ig++;
      ivx++;
      ivy++;
    }
}

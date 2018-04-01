#ifndef MESH_INFO_H
#define MESH_INFO_H

#include "kernel_typedefs.h"

typedef struct
{
  int MX;
  int MY;
  int N;
  double X;
  double Y;
  double T;
  double hx;
  double hy;
  double tau;

  double border_omega;
  double mu;

} mesh_info_t;

int mesh_info_init (mesh_info_t *info,
                     double X,
                     double Y,
                     double T,
                     int MX,
                     int MY,
                     int N,
                     double border_omega,
                     double mu);

grid_area_t mesh_info_get_area (const mesh_info_t *info,
                                int mx, int my);

#endif /* MESH_INFO_H */

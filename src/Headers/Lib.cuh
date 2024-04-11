#ifndef __LIB_H__
#define __LIB_H__

/* Include Headers Here*/
#include "SPH.cuh"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>

/* Extern Function Here*/
extern __global__ void sph_mesh_cuda(double *x,double *y,double *accx,double *accy,double *drho,int *type,int *mesh,int *count,int ptc_num);
extern __global__ void sph_nnps_cuda(int *,double *,double *,int *,int *,int *,int *);
extern __global__ void sph_kernel_cuda(double *,double *,double *,double *,double *,double *,int *,int *,int *);
extern __global__ void sph_dummy_cuda(double *,double *,double *,double *,double *,double *,double *,double *,int *,int *,int *,double *,int *);

#endif
#ifndef __LIB_H__
#define __LIB_H__

/* Include Headers Here*/
#include "SPH.cuh"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>

/* Extern Function Here*/
extern __global__ void sph_mesh_cuda(SPH_CUDA *,SPH_ARG *,SPH_RIGID *);
/*extern __global__ void sph_nnps_cuda(int *,double *,double *,int *,int *,int *,int *);
extern __global__ void sph_kernel_cuda(double *,double *,double *,double *,double *,double *,int *,int *,int *);
extern __global__ void sph_dummy_cuda(double *,double *,double *,double *,double *,double *,double *,double *,int *,int *,int *,double *,int *);*/
extern __device__ void sph_cuda_lock(SPH_ARG *);
extern __device__ void sph_cuda_unlock(SPH_ARG *);

#endif
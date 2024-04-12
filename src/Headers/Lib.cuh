#ifndef __LIB_H__
#define __LIB_H__

/* Include Headers Here*/
#include "DataStructure.cuh"

/* Extern Function Here*/
extern __global__ void sph_mesh_cuda(SPH_CUDA *,SPH_ARG *);
extern __global__ void sph_nnps_cuda(SPH_CUDA *,SPH_ARG *,SPH_RIGID *);
extern void sph_nnps_cpu(SPH *);
/*extern __global__ void sph_nnps_cuda(int *,double *,double *,int *,int *,int *,int *);
extern __global__ void sph_kernel_cuda(double *,double *,double *,double *,double *,double *,int *,int *,int *);
extern __global__ void sph_dummy_c*/
extern __device__ void sph_lock_cuda(SPH_ARG *);
extern __device__ void sph_unlock_cuda(SPH_ARG *);

#endif
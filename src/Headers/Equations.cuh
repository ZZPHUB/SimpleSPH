#ifndef __EQUATIONS_H__
#define __EQUATIONS_H__

/* Include Headers Here */
#include "SPH.cuh"

/* Extern Functions Here */
extern __global__ void sph_governing_cuda(double *,double *,double *,double *,\
double *,double *,int *,int *,int *,double *,\
double *,double *,double *,double *,double *,int* ,int );


#endif
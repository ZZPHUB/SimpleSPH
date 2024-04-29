#ifndef __EQUATIONS_H__
#define __EQUATIONS_H__

#include "SPH.cuh"

extern __global__ void sph_governing_cuda(SPH_CUDA *,SPH_ARG *,SPH_RIGID *);
extern void sph_delta_cuda(SPH *);

#endif
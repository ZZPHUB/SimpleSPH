#include "Lib.cuh"

__device__ void sph_lock_cuda(SPH_ARG *arg)
{
    while (!atomicCAS(&(arg->lock),1,0))
    {
        continue;
    }
    
}

__device__ void sph_unlock_cuda(SPH_ARG *arg)
{
    atomicCAS(&(arg->lock),0,1);
}
#ifndef __SPH_H__
#define __SPH_H__

#include "DataStructure.cuh"


//#define PARA 0
#define PARA (0x01|0x02|0x04|0x08)
//0x01----------->density
//0x02----------->pressure
//0x04----------->velosity
//0x08----------->acceleration


#include "Lib.cuh"
#include "IO.cuh"
#include "Equations.cuh"


/* Extern Functions Here*/
void sph_init(SPH *);
void sph_free(SPH *);
__global__ void sph_predict_cuda(SPH_CUDA *,SPH_ARG *,SPH_RIGID *);
__global__ void sph_correct_cuda(SPH_CUDA *,SPH_ARG *,SPH_RIGID *);
__global__ void sph_rigid_cuda(SPH_CUDA *,SPH_ARG *,SPH_RIGID *);


#endif

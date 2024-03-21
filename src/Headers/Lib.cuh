#ifndef __LIB_H__
#define __LIB_H__

/* Include Headers Here*/
#include "SPH.cuh"
#include "sys/resource.h"
#include "cmath"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>

/* Extern Function Here*/
extern void ptc_nnps_direct(SPH *);
extern void ptc_nnps_mesh(SPH *);
extern void ptc_nnps_check(SPH_PAIR *,SPH_PAIR *,unsigned int *);
extern void ptc_kernel_serial(SPH *);
extern void ptc_kernel_parallel(SPH *);
extern void ptc_mesh_process(SPH *);
extern void ptc_info(SPH *);
extern void ptc_density_correct(SPH *);
extern void sph_avg_time(SPH *);
extern __global__ void ptc_mesh_cuda(double *,double *,double *,int );

#include "SPH.cuh"
#define CUDA_CHECK(call)             __cudaCheck(call, __FILE__, __LINE__)
#define LAST_KERNEL_CHECK(call)      __kernelCheck(__FILE__, __LINE__)

static void __cudaCheck(cudaError_t err, const char* file, const int line) {
    if (err != cudaSuccess) {
        printf("ERROR: %s:%d, ", file, line);
        printf("CODE:%s, DETAIL:%s\n", cudaGetErrorName(err), cudaGetErrorString(err));
        exit(1);
    }
}

static void __kernelCheck(const char* file, const int line) {
    cudaError_t err = cudaPeekAtLastError();
    if (err != cudaSuccess) {
        printf("ERROR: %s:%d, ", file, line);
        printf("CODE:%s, DETAIL:%s\n", cudaGetErrorName(err), cudaGetErrorString(err));
        exit(1);
    }
}

#endif
#include "Lib.cuh"

__global__ void sph_mesh_cuda(double *x,double *y,int *mesh,int ptc_num)
{
    //const int bid = blockIdx.x;
    //const int tid = threadIdx.x;
    const int id = threadIdx.x + blockIdx.x * blockDim.x;
    if(id >= ptc_num) return;

    int mid;

    if(y[id] < TOL_DOMAIN_DEEPTH && y[id] >= 0)
    {
        mid = __double2int_rz(y[id]/MESH_SPACING)*MESH_LENGTH_NUM_CUDA;
    }
    else
    {
        mid = (MESH_DEEPTH_NUM_CUDA - 1)*MESH_LENGTH_NUM_CUDA;
    }
    if(x[id] < TOL_DOMAIN_LENGTH && x[id] >= 0)
    {
        mid += __double2int_rz(x[id]/MESH_SPACING);
    }
    else
    {
        mid += MESH_LENGTH_NUM_CUDA - 1;
    }
    mid += MESH_DEEPTH_NUM_CUDA*MESH_LENGTH_NUM_CUDA*atomicAdd(&mesh[mid+(MESH_PTC_NUM-1)*MESH_DEEPTH_NUM_CUDA*MESH_LENGTH_NUM_CUDA],1);
    mesh[mid] = id;
}


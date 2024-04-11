#include "SPH.cuh"

__global__ void sph_fuck_you(SPH_CUDA *cuda,SPH_ARG *arg)
{
    //const int bid = blockIdx.x;
    //const int tid = threadIdx.x;
    const int id = threadIdx.x + blockIdx.x * blockDim.x;
    if(id >= arg->ptc_num) return;

    /*这里需要进行加速度和密度变化的初始化*/

    /*这里需要对pair_num进行初始化*/
    if(id == 0) arg->pair_num = 0;

    int mid = 0;

    if(cuda->y[id] < arg->domain_y && cuda->y[id] >= 0.0)
    {
        mid = __double2int_rz(cuda->y[id]/arg->mesh_dx)*arg->mesh_xnum;
    }
    else
    {
        mid = (arg->mesh_ynum - 1)*arg->mesh_xnum;
    }
    if(cuda->x[id] < arg->domain_x && cuda->x[id] >= 0.0)
    {
        mid += __double2int_rz(cuda->x[id]/arg->mesh_dx);
    }
    else
    {
        mid += arg->mesh_xnum - 1;
    }

    while(!atomicCAS(&arg->lock,1,0))
    {
        continue;
    } 
    cuda->mesh[mid + arg->mesh_num*cuda->mesh_count[mid]] = id;
    cuda->mesh_count[mid] += 1; 
    atomicCAS(&arg->lock,0,1);
}


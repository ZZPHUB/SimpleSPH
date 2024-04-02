#include "Lib.cuh"

__global__ void sph_mesh_cuda(double *x,double *y,double *accx,double *accy,double *drho,int *type,int *mesh,int ptc_num)
{
    //const int bid = blockIdx.x;
    //const int tid = threadIdx.x;
    const int id = threadIdx.x + blockIdx.x * blockDim.x;
    if(id >= ptc_num) return;
    accx[id] = 0.0;
    drho[id] = 0.0;
    if(type[id] == 0) accy[id] = -GRAVITY_ACC;
    else accy[id] = 0.0;

    int mid;

    if(y[id] < TOL_DOMAIN_DEEPTH && y[id] >= 0)
    {
        mid = __double2int_rz(y[id]/dev_mesh_spacing)*dev_mesh_lnum;
    }
    else
    {
        mid = (dev_mesh_dnum - 1)*dev_mesh_lnum;
    }
    if(x[id] < TOL_DOMAIN_LENGTH && x[id] >= 0)
    {
        mid += __double2int_rz(x[id]/dev_mesh_spacing);
    }
    else
    {
        mid += dev_mesh_lnum - 1;
    }
    mid += dev_mesh_tnum*atomicAdd(&mesh[mid+(MESH_PTC_NUM-1)*dev_mesh_tnum],1);
    mesh[mid] = id;
}


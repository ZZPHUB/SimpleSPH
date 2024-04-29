#include "Lib.cuh"

__global__ void sph_kernel_cuda(SPH_CUDA *cuda,SPH_ARG *arg,SPH_RIGID *rigid)
{
    double dx = 0.0;
    double dy = 0.0;
    double r = 0.0;
    double q = 0.0;
    int index_i = 0;
    int index_j = 0;
    int id = 0;
    const int mesh_id = blockIdx.x + blockIdx.y * gridDim.x;
    if( threadIdx.x < cuda->pair_count[mesh_id])
    {
        id = mesh_id * arg->pair_volume + threadIdx.x;
        index_i = cuda->pair_i[id];
        index_j = cuda->pair_j[id];

        dx = cuda->x[index_i] - cuda->x[index_j];
        dy = cuda->y[index_i] - cuda->y[index_j];
        r = sqrt(dx*dx+dy*dy);
        q = r/arg->h;

        if(q <= 2.0)
        {
            cuda->pair_w[id] = arg->alpha*((1.0+2.0*q)*pow((1.0-0.5*q),4));
            cuda->dwdx[id] = arg->alpha*(2.0*pow((1.0-0.5*q),4)-2.0*(1.0+2.0*q)*pow((1.0-0.5*q),3))*dx/(r*arg->h);
            cuda->dwdy[id] = arg->alpha*(2.0*pow((1.0-0.5*q),4)-2.0*(1.0+2.0*q)*pow((1.0-0.5*q),3))*dy/(r*arg->h);
        }
        else
        {
            cuda->pair_w[id] = 0;
            cuda->dwdx[id] = 0;
            cuda->dwdy[id] = 0;
        }
    }
    __syncthreads();
}

__global__ void sph_sum_w(SPH_CUDA *cuda,SPH_ARG *arg,SPH_RIGID *rigid)
{
    int index_i = 0;
    int index_j = 0;
    int id = 0;
    double ptc_w_i = 0.0;
    double ptc_w_j = 0.0;
    const int mesh_id = blockIdx.x + blockIdx.y * gridDim.x;
    if( threadIdx.x < cuda->pair_count[mesh_id])
    {
        id = mesh_id * arg->pair_volume + threadIdx.x;
        index_i = cuda->pair_i[id];
        index_j = cuda->pair_j[id];

        ptc_w_i = arg->m*cuda->pair_w[id]/cuda->rho[index_j];
        if(cuda->type[index_j] == 0) ptc_w_j = arg->m*cuda->pair_w[id]/cuda->rho[index_i];
        else ptc_w_j = cuda->pair_w[id];

        atomicAdd(&(cuda->ptc_w[index_i]),ptc_w_i);
        atomicAdd(&(cuda->ptc_w[index_j]),ptc_w_j);
    }
    __syncthreads();
}
/*if(q<=1.0)
    {
        cuda->pair_w[id] = arg->alpha*(2.0/3.0-q*q+0.5*q*q*q);
        cuda->dwdx[id] = arg->alpha*((-2.0+1.5*q)*dx/pow(arg->h,2));
        cuda->dwdy[id] = arg->alpha*((-2.0+1.5*q)*dy/pow(arg->h,2));
    }
    else if(1.0 <q && q < 2.0)
    {
        cuda->pair_w[id] = arg->alpha*((2.0-q)*(2.0-q)*(2.0-q)/6.0);
        cuda->dwdx[id] = -arg->alpha*0.5*((2.0-q)*(2.0-q)*dx/(arg->h*r));
        cuda->dwdy[id] = -arg->alpha*0.5*((2.0-q)*(2.0-q)*dy/(arg->h*r));
    }
    else
    {
        cuda->pair_w[id] = 0;
        cuda->dwdx[id] = 0;
        cuda->dwdy[id] = 0;
    }*/


//sum ptc w here,but get wrong result
/*atomicAdd(&(cuda->ptc_w[index_i]),cuda->pair_w[id]*arg->m/cuda->rho[index_j]);
    if(cuda->type[index_j] != 0)
    {
        atomicAdd(&(cuda->ptc_w[index_j]),cuda->pair_w[id]);
    }
    else
    {
        atomicAdd(&(cuda->ptc_w[index_j]),cuda->pair_w[id]*arg->m/cuda->rho[index_i]);
    }*/
#include "Lib.cuh"

__global__ void sph_dummy_cuda(SPH_CUDA *cuda,SPH_ARG *arg,SPH_RIGID *rigid)
{
    double dx = 0.0;
    double dy = 0.0;
    double tmp_accx = 0.0;
    double tmp_accy = 0.0;
    double tmp_prho = 0.0;
    double tmp_vx = 0.0;
    double tmp_vy = 0.0;
    int index_i = 0;
    int index_j = 0;
    int id = 0;
    const int mesh_id = blockIdx.x + blockIdx.y * gridDim.x;
    if( threadIdx.x < cuda->pair_count[mesh_id])
    {
        id = mesh_id * arg->pair_volume + threadIdx.x;
        index_i = cuda->pair_i[id];
        index_j = cuda->pair_j[id];

        if(cuda->type[index_j] != 0  && cuda->ptc_w[index_j] != 0.0)
        {
            if(cuda->type[index_j] == -1)
            {
                tmp_accx = 0.0;
                tmp_accy = 0.0;
            }
            else if(cuda->type[index_j] == 1)
            {
                tmp_accx = rigid->accx - pow(rigid->omega,2)*(cuda->x[index_j]-rigid->cogx)- \
                              rigid->alpha*(cuda->y[index_j]-rigid->cogy);
                tmp_accy = rigid->accy - pow(rigid->omega,2)*(cuda->y[index_j]-rigid->cogy)+ \
                              rigid->alpha*(cuda->x[index_j]-rigid->cogx);
            }
            dx = cuda->x[index_i]-cuda->y[index_j];
            dy = cuda->y[index_i]-cuda->y[index_j];

            tmp_prho = (cuda->p[index_i]+cuda->rho[index_i]*(tmp_accx*dx+(tmp_accy+arg->g)*dy))*cuda->pair_w[id]/cuda->ptc_w[index_j];
            atomicAdd(&(cuda->p[index_j]),tmp_prho);
            
            tmp_prho /= arg->c*arg->c;
            atomicAdd(&(cuda->rho[index_j]),tmp_prho);

            tmp_vx = cuda->vx[index_i]*cuda->pair_w[id]/cuda->ptc_w[index_j];
            atomicAdd(&(cuda->vx[index_j]),tmp_vx);

            tmp_vy = cuda->vy[index_i]*cuda->pair_w[id]/cuda->ptc_w[index_j];
            atomicAdd(&(cuda->vy[index_j]),tmp_vy);
        }   
    }
    __syncthreads();
}

__global__ void sph_check_rho(SPH_CUDA *cuda,SPH_ARG *arg,SPH_RIGID *rigid)
{
    const int id = threadIdx.x + blockIdx.x * blockDim.x;
    if(id < arg->ptc_num)
    {
        if(cuda->type[id] != 0)
        {
            if(cuda->rho[id] < arg->ref_rho) 
            {
                cuda->rho[id] = arg->ref_rho;
                cuda->p[id] = 0.0;
            }
        }
        /*if(cuda->rho[id] < arg->ref_rho)
        {
            cuda->rho[id] = arg->ref_rho;
            cuda->p[id] = 0.0;
        }
        else
        {
            cuda->p[id] = arg->c*arg->c*(cuda->rho[id] - arg->ref_rho);
        }*/
    }
}


__global__ void sph_filter_init(SPH_CUDA *cuda,SPH_ARG *arg,SPH_RIGID *rigid)
{
    const int id = blockIdx.x * blockDim.x + threadIdx.x;
    if(id < arg->ptc_num )
    {
        if(cuda->type[id] == 0 && cuda->ptc_w[id] != 0.0) cuda->rho[id] = arg->m*arg->alpha/cuda->ptc_w[id];
    }
}

__global__ void sph_filter_sum(SPH_CUDA *cuda,SPH_ARG *arg,SPH_RIGID *rigid)
{
    const int mesh_id = blockIdx.x + blockIdx.y * gridDim.x;
    int id = 0;
    int index_i = 0;
    int index_j = 0;
    double tmp_rho_i = 0.0;
    double tmp_rho_j = 0.0;
    if( threadIdx.x < cuda->pair_count[mesh_id])
    {
        id = mesh_id * arg->pair_volume + threadIdx.x;
        index_i = cuda->pair_i[id];
        index_j = cuda->pair_j[id];

        if(cuda->ptc_w[index_i] != 0.0) tmp_rho_i = arg->m*cuda->pair_w[id]/cuda->ptc_w[index_i];
        if(cuda->type[index_j] == 0)
        {
            if(cuda->ptc_w[index_j] != 0.0) tmp_rho_j = arg->m*cuda->pair_w[id]/cuda->ptc_w[index_j];
        }
        atomicAdd(&(cuda->rho[index_i]),tmp_rho_i);
        atomicAdd(&(cuda->rho[index_j]),tmp_rho_j);
        /*if(cuda->ptc_w[index_i] != 0.0) atomicAdd(&(cuda->rho[index_i]), arg->m*cuda->pair_w[id]/cuda->ptc_w[index_i]);
        if(cuda->type[index_j] == 0) 
        {
            if(cuda->ptc_w[index_j] != 0.0) atomicAdd(&(cuda->rho[index_j]),arg->m*cuda->pair_w[id]/cuda->ptc_w[index_j]);
        }*/
    }
}

void sph_rho_filter(SPH *sph)
{
    dim3 pair_block(sph->host_arg->pair_volume);
    dim3 pair_grid(sph->host_arg->mesh_xnum, sph->host_arg->mesh_ynum);
    dim3 ptc_block(256);
    dim3 ptc_grid((int)(sph->host_arg->ptc_num / 256) + 1);

    sph_filter_init<<<ptc_grid,ptc_block>>>(sph->cuda,sph->dev_arg,sph->dev_rigid);
    cudaDeviceSynchronize();
    sph_filter_sum<<<pair_grid,pair_block>>>(sph->cuda,sph->dev_arg,sph->dev_rigid);
    cudaDeviceSynchronize();
}
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
    if(threadIdx < cuda->pair_count[mesh_id])
    {
        id = mesh_id * arg->pair_volume + threadIdx.x;
        index_i = cuda->pair_i[id];
        index_j = cuda->pair_j[id];

        if(cuda->type[index_j] != 0 cuda->ptc_w[index_j] != 0.0)
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
            atomicadd(&(cuda->rho[index_j]),tmp_prho);

            tmp_vx = cuda->vx[index_i]*cuda->pair_w[id]/cuda->ptc_w[index_j];
            atomicAdd(&(cuda->vx[index_j]),tmp_vx);

            tmp_vy = cuda->vy[index_i]*cuda->pair_w[id]/cuda->ptc_w[index_j];
            atomicAdd(&(cuda->vy[index_j]),tmp_vy);
        }   
    }
    __syncthreads();
    if( threadIdx.x == 0)cuda->pair_count[mesh_id] = 0;
}

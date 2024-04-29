#include "Equations.cuh"

__global__ void sph_governing_cuda(SPH_CUDA *cuda,SPH_ARG *arg,SPH_RIGID *rigid)
{
    double accx = 0.0;
    double accy = 0.0;
    double drho_i = 0.0;
    double drho_j = 0.0;
    double tmp_acc_p = 0.0;
    double tmp_acc_v = 0.0;
    int index_i = 0.0;
    int index_j = 0.0;
    double dx = 0.0;
    double dy = 0.0;
    double dvx = 0.0;
    double dvy = 0.0;
    const int mesh_id = blockIdx.x + blockIdx.y * gridDim.x;
    int id = 0;
    if( threadIdx.x < cuda->pair_count[mesh_id]) 
    {
        id = mesh_id*arg->pair_volume + threadIdx.x;
        index_i = cuda->pair_i[id];
        index_j = cuda->pair_j[id];

        tmp_acc_p = cuda->p[index_i]/pow(cuda->rho[index_i],2) + cuda->p[index_j]/pow(cuda->rho[index_j],2);
        dx = cuda->x[index_i] - cuda->x[index_j];
        dy = cuda->y[index_i] - cuda->y[index_j];
        if(cuda->type[index_j] == 0)
        {
            dvx = cuda->vx[index_i] - cuda->vx[index_j];
            dvy = cuda->vy[index_i] - cuda->vy[index_j];
            drho_i = (cuda->vx[index_i]-cuda->vx[index_j])*cuda->dwdx[id]+(cuda->vy[index_i]-cuda->vy[index_j])*cuda->dwdy[id];
            drho_i *= arg->m;
        }
        else if(cuda->type[index_j] == -1)
        {
            dvx = cuda->vx[index_i] - (0.0 - cuda->vx[index_j]);
            dvy = cuda->vy[index_i] - (0.0 - cuda->vy[index_j]);
            drho_i = cuda->vx[index_i]*cuda->dwdx[id]+cuda->vy[index_i]*cuda->dwdy[id];
            drho_i *= arg->m;
        }
        else if(cuda->type[index_j] == 1)
        {
            dvx = cuda->vx[index_i] - (2.0*(rigid->vx - rigid->omega*(cuda->y[index_j]-rigid->cogy)) - cuda->vx[index_j]);
            dvy = cuda->vy[index_i] - (2.0*(rigid->vy + rigid->omega*(cuda->x[index_j]-rigid->cogx)) - cuda->vy[index_j]);
            drho_i = (cuda->vx[index_i] - (rigid->vx - rigid->omega*(cuda->y[index_j]-rigid->cogy)))*cuda->dwdx[id]+\
                      (cuda->vy[index_i] - (rigid->vy + rigid->omega*(cuda->x[index_j]-rigid->cogx)))*cuda->dwdy[id];
            drho_i *= arg->m;
        }

        tmp_acc_v = dx*dvx+dy*dvy;
        if(tmp_acc_v > 0.0) tmp_acc_v = 0.0;
        tmp_acc_v = (tmp_acc_v*0.01*arg->h*arg->c)/((dx*dx+dy*dy+0.01*arg->h)*0.5*(cuda->rho[index_i]+cuda->rho[index_j]));

        accx = arg->m * ( tmp_acc_v - tmp_acc_p) *cuda->dwdx[id];
        accy = arg->m * ( tmp_acc_v - tmp_acc_p) *cuda->dwdy[id];

        drho_j = drho_i;
        drho_i += 0.01*arg->h*arg->c*2*(cuda->rho[index_i]/cuda->rho[index_j]-1)*arg->m*(dx*cuda->dwdx[id]+dy*cuda->dwdy[id])/(dx*dx+dy*dy);
        drho_j += 0.01*arg->h*arg->c*2*(cuda->rho[index_j]/cuda->rho[index_i]-1)*arg->m*(dx*cuda->dwdx[id]+dy*cuda->dwdy[id])/(dx*dx+dy*dy); 
        
        atomicAdd(&(cuda->accx[index_i]),accx);
        atomicAdd(&(cuda->accx[index_j]),-accx);
        atomicAdd(&(cuda->accy[index_i]),accy);
        atomicAdd(&(cuda->accy[index_j]),-accy);
        atomicAdd(&(cuda->drho[index_i]),drho_i);
        atomicAdd(&(cuda->drho[index_j]),drho_j);
    }
    __syncthreads();
    if( threadIdx.x == 0)cuda->pair_count[mesh_id] = 0;
}
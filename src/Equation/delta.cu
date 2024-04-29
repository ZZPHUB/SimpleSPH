#include "Equations.cuh"

__global__ void sph_L_init(SPH_CUDA *cuda,SPH_ARG *arg,SPH_RIGID *rigid)
{
    const int id = threadIdx.x + blockIdx.x * blockDim.x;
    if(id < arg->ptc_num)
    {
        cuda->Lxx[id] = 0.0;
        cuda->Lxy[id] = 0.0;
        cuda->Lyx[id] = 0.0;
        cuda->Lyy[id] = 0.0;

        cuda->Lrho_x[id] = 0.0;
        cuda->Lrho_y[id] = 0.0;
    }
}

__global__ void sph_L_sum(SPH_CUDA *cuda,SPH_ARG *arg,SPH_RIGID *rigid)
{
    const int mesh_id = blockIdx.x + blockIdx.y * gridDim.x;
    int id = 0;
    int index_i = 0;
    int index_j = 0;
    double tmp_Lxx = 0.0;
    double tmp_Lxy = 0.0;
    double tmp_Lyx = 0.0;
    double tmp_Lyy = 0.0;

    if( threadIdx.x < cuda->pair_count[mesh_id])
    {
        id = mesh_id * arg->pair_volume + threadIdx.x;
        index_i = cuda->pair_i[id];
        index_j = cuda->pair_j[id];

        //dx = cuda->x[index_i] - cuda->x[index_j];
        //dy = cuda->y[index_i] - cuda->y[index_j];
        tmp_Lxx = (cuda->x[index_i] - cuda->x[index_j])*cuda->dwdx[id]*arg->m;
        tmp_Lxy = (cuda->x[index_i] - cuda->x[index_j])*cuda->dwdy[id]*arg->m;
        tmp_Lyx = (cuda->y[index_i] - cuda->y[index_j])*cuda->dwdx[id]*arg->m;
        tmp_Lyy = (cuda->y[index_i] - cuda->y[index_j])*cuda->dwdy[id]*arg->m;

        atomicAdd(&(cuda->Lxx[index_i]),tmp_Lxx/cuda->rho[index_j]);
        atomicAdd(&(cuda->Lxy[index_i]),tmp_Lxy/cuda->rho[index_j]);
        atomicAdd(&(cuda->Lyx[index_i]),tmp_Lyx/cuda->rho[index_j]);
        atomicAdd(&(cuda->Lyy[index_i]),tmp_Lyy/cuda->rho[index_j]);
        if(cuda->type[index_j] == 0)
        {
            atomicAdd(&(cuda->Lxx[index_j]),tmp_Lxx/cuda->rho[index_i]);
            atomicAdd(&(cuda->Lxy[index_j]),tmp_Lxy/cuda->rho[index_i]);
            atomicAdd(&(cuda->Lyx[index_j]),tmp_Lyx/cuda->rho[index_i]);
            atomicAdd(&(cuda->Lyy[index_j]),tmp_Lyy/cuda->rho[index_i]);
        }
    }
}

__global__ void sph_L_inver(SPH_CUDA *cuda,SPH_ARG *arg,SPH_RIGID *rigid)
{
    const int id = threadIdx.x + blockIdx.x * blockDim.x;
    double det = 0.0;
    double tmp_Lxx = 0.0;
    double tmp_Lyy = 0.0;

    if(id < arg->ptc_num)
    {
        if(cuda->type[id] == 0)
        {
            det = cuda->Lxx[id]*cuda->Lyy[id] - cuda->Lxy[id]*cuda->Lyx[id];
            if(det != 0.0)
            {
                tmp_Lxx = cuda->Lyy[id]/det;
                tmp_Lyy = cuda->Lxx[id]/det;

                cuda->Lxx[id] = tmp_Lxx;
                cuda->Lyy[id] = tmp_Lyy;
                cuda->Lxy[id] /= -det;
                cuda->Lyx[id] /= -det;
            }
        }
    }
}

__global__ void sph_L_rho(SPH_CUDA *cuda,SPH_ARG *arg,SPH_RIGID *rigid)
{
    const int mesh_id = blockIdx.x + blockIdx.y * gridDim.x;
    int id = 0;
    int index_i = 0;
    int index_j = 0;
    double tmp_Lrho_x_i = 0.0;
    double tmp_Lrho_y_i = 0.0;
    double tmp_Lrho_x_j = 0.0;
    double tmp_Lrho_y_j = 0.0;
    if( threadIdx.x < cuda->pair_count[mesh_id])
    {
        id = mesh_id * arg->pair_volume + threadIdx.x;
        index_i = cuda->pair_i[id];
        index_j = cuda->pair_j[id];

        tmp_Lrho_x_i = (cuda->rho[index_i]-cuda->rho[index_j])*(cuda->Lxx[index_i]*cuda->dwdx[id] + cuda->Lxy[index_i]*cuda->dwdy[id])*arg->m/cuda->rho[index_j];
        tmp_Lrho_y_i = (cuda->rho[index_i]-cuda->rho[index_j])*(cuda->Lyx[index_i]*cuda->dwdx[id] + cuda->Lyy[index_i]*cuda->dwdy[id])*arg->m/cuda->rho[index_j];
        atomicAdd(&(cuda->Lrho_x[index_i]),tmp_Lrho_x_i);
        atomicAdd(&(cuda->Lrho_y[index_i]),tmp_Lrho_y_i);

        tmp_Lrho_x_j = (cuda->rho[index_i]-cuda->rho[index_j])*(cuda->Lxx[index_j]*cuda->dwdx[id] + cuda->Lxy[index_j]*cuda->dwdy[id])*arg->m/cuda->rho[index_i];
        tmp_Lrho_y_j = (cuda->rho[index_i]-cuda->rho[index_j])*(cuda->Lyx[index_j]*cuda->dwdx[id] + cuda->Lyy[index_j]*cuda->dwdy[id])*arg->m/cuda->rho[index_i];
        atomicAdd(&(cuda->Lrho_x[index_j]),tmp_Lrho_x_j);
        atomicAdd(&(cuda->Lrho_y[index_j]),tmp_Lrho_y_j);
    }
}

__global__ void sph_delta_term(SPH_CUDA *cuda,SPH_ARG *arg,SPH_RIGID *rigid)
{
    const int mesh_id = blockIdx.x + blockIdx.y * gridDim.x;
    int id = 0;
    int index_i = 0;
    int index_j = 0;
    double drho = 0.0;
    double dx = 0.0;
    double dy = 0.0;
    if( threadIdx.x < cuda->pair_count[mesh_id])
    {
        id =mesh_id * arg->pair_volume + threadIdx.x;
        index_i = cuda->pair_i[id];
        index_j = cuda->pair_j[id];
        dx = cuda->x[index_i] - cuda->x[index_j];
        dy = cuda->y[index_i] - cuda->y[index_j];
        
        drho = 2.0*(cuda->rho[index_i] - cuda->rho[index_j]);
        drho += (cuda->Lrho_x[index_i] + cuda->Lrho_x[index_j])*dx + (cuda->Lrho_y[index_i]+cuda->Lrho_y[index_j])*dy;
        drho *= 0.01*arg->c*arg->h*(dx*cuda->dwdx[id] + dy*cuda->dwdy[id])*arg->m/(dx*dx+dy*dy);

        atomicAdd(&(cuda->drho[index_i]),-drho/cuda->rho[index_j]);
        atomicAdd(&(cuda->drho[index_j]),drho/cuda->rho[index_i]);
    }
    __syncthreads();
    if( threadIdx.x == 0) cuda->pair_count[mesh_id] = 0;
}

void sph_delta_cuda(SPH *sph)
{
    dim3 pair_block(sph->host_arg->pair_volume);
    dim3 pair_grid(sph->host_arg->mesh_xnum, sph->host_arg->mesh_ynum);
    dim3 ptc_block(256);
    dim3 ptc_grid((int)(sph->host_arg->ptc_num / 256) + 1);

    sph_L_init<<<ptc_grid,ptc_block>>>(sph->cuda,sph->dev_arg,sph->dev_rigid);
    cudaDeviceSynchronize();
    sph_L_sum<<<pair_grid,pair_block>>>(sph->cuda,sph->dev_arg,sph->dev_rigid);
    cudaDeviceSynchronize();
    sph_L_inver<<<ptc_grid,ptc_block>>>(sph->cuda,sph->dev_arg,sph->dev_rigid);
    cudaDeviceSynchronize();
    sph_L_rho<<<pair_grid,pair_block>>>(sph->cuda,sph->dev_arg,sph->dev_rigid);
    cudaDeviceSynchronize();
    sph_delta_term<<<pair_grid,pair_block>>>(sph->cuda,sph->dev_arg,sph->dev_rigid);
    cudaDeviceSynchronize();
}
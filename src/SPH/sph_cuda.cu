#include "SPH.cuh"
#include <assert.h>

int main(int argc,char *argv[])
{
    SPH_PARTICLE particle;
    SPH_KERNEL kernel;
    SPH_PAIR pair;
    SPH_RIGID wedge;
    SPH_MESH mesh;
    SPH_ARG arg;
    SPH_CUDA tmp_cuda;
    SPH sph;
    sph.particle = &particle;
    sph.kernel = &kernel;
    sph.pair = &pair;
    sph.host_rigid = &wedge;
    sph.host_arg = &arg;
    sph.mesh = &mesh;
    sph.tmp_cuda = &tmp_cuda;
    assert("argc == 2");
    arg.case_dir = argv[1];

    cudaSetDevice(0);
    sph_init(&sph);

    // define the seed for ptc data structure
    dim3 ptc_block(256);
    dim3 ptc_grid((int)(sph.host_arg->ptc_num / 256) + 1);
    // define the seed for mesh data structure
    dim3 mesh_block(32, 32);
    dim3 mesh_grid(sph.host_arg->mesh_xnum, sph.host_arg->mesh_ynum);
    // define the seed for pair data structre
    dim3 pair_block(128);
    dim3 pair_grid(sph.host_arg->mesh_xnum, sph.host_arg->mesh_ynum);

    // SPH_CUDA cuda;
    // SPH_ARG tmp_arg;
    // cudaMemcpy(&cuda,sph.cuda,sizeof(SPH_CUDA),cudaMemcpyDeviceToHost);
    // cudaDeviceSynchronize();
    // int *host_pair_count = (int *)calloc(sph.host_arg->mesh_num,sizeof(int));
    // int id = 0;
    // int *cpu_pair_i;
    // int *cpu_pair_j;
    // cudaMalloc(&cpu_pair_i,sizeof(int)*32*sph.host_arg->ptc_num);
    // cudaMalloc(&cpu_pair_j,sizeof(int)*32*sph.host_arg->ptc_num);

    for (sph.host_arg->init_step; sph.host_arg->init_step < sph.host_arg->total_step; sph.host_arg->init_step++)
    {
        printf("current step is:%d ", sph.host_arg->init_step);
        sph_mesh_cuda<<<ptc_grid, ptc_block>>>(sph.cuda, sph.dev_arg);
        if (sph.host_arg->init_step % sph.host_arg->print_step == 2 && sph.host_arg->init_impac_flag == 0)
        {
            sph_write_csv(&sph);
        }
        cudaDeviceSynchronize();
        sph_nnps_cuda<<<mesh_grid, mesh_block>>>(sph.cuda, sph.dev_arg, sph.dev_rigid);
        if (sph.host_arg->init_step % sph.host_arg->print_step == 1)
        {
            sph_save_single(&sph);
        }
        cudaDeviceSynchronize();
        sph_kernel_cuda<<<pair_grid, pair_block>>>(sph.cuda, sph.dev_arg, sph.dev_rigid);
        cudaDeviceSynchronize();
        sph_governing_cuda<<<pair_grid, pair_block>>>(sph.cuda, sph.dev_arg, sph.dev_rigid);
        cudaDeviceSynchronize();
        sph_predict_cuda<<<ptc_grid, ptc_block>>>(sph.cuda, sph.dev_arg, sph.dev_rigid);
        cudaDeviceSynchronize();
        sph_dummy_cuda<<<pair_grid, pair_block>>>(sph.cuda, sph.dev_arg, sph.dev_rigid);
        cudaDeviceSynchronize();

        sph_mesh_cuda<<<ptc_grid, ptc_block>>>(sph.cuda, sph.dev_arg);
        cudaDeviceSynchronize();
        sph_nnps_cuda<<<mesh_grid, mesh_block>>>(sph.cuda, sph.dev_arg, sph.dev_rigid);
        cudaDeviceSynchronize();
        sph_kernel_cuda<<<pair_grid, pair_block>>>(sph.cuda, sph.dev_arg, sph.dev_rigid);
        cudaDeviceSynchronize();
        sph_governing_cuda<<<pair_grid, pair_block>>>(sph.cuda, sph.dev_arg, sph.dev_rigid);
        cudaDeviceSynchronize();
        sph_correct_cuda<<<ptc_grid, ptc_block>>>(sph.cuda, sph.dev_arg, sph.dev_rigid);
        cudaDeviceSynchronize();
        sph_dummy_cuda<<<pair_grid, pair_block>>>(sph.cuda, sph.dev_arg, sph.dev_rigid);
        cudaDeviceSynchronize();
        if(sph.host_arg->init_impac_flag == 0)
        {
            sph_rigid_cuda<<<ptc_grid,ptc_block>>>(sph.cuda,sph.dev_arg,sph.dev_rigid);
            cudaDeviceSynchronize();
        }

        if (sph.host_arg->init_step % sph.host_arg->print_step == 0)
        {
            cudaMemcpy(sph.particle->x, sph.tmp_cuda->x, sizeof(double) * sph.host_arg->ptc_num, cudaMemcpyDeviceToHost);
            cudaDeviceSynchronize();
            cudaMemcpy(sph.particle->y, sph.tmp_cuda->y, sizeof(double) * sph.host_arg->ptc_num, cudaMemcpyDeviceToHost);
            cudaDeviceSynchronize();
            cudaMemcpy(sph.particle->vx, sph.tmp_cuda->vx, sizeof(double) * sph.host_arg->ptc_num, cudaMemcpyDeviceToHost);
            cudaDeviceSynchronize();
            cudaMemcpy(sph.particle->vy, sph.tmp_cuda->vy, sizeof(double) * sph.host_arg->ptc_num, cudaMemcpyDeviceToHost);
            cudaDeviceSynchronize();
            cudaMemcpy(sph.particle->accx, sph.tmp_cuda->accx, sizeof(double) * sph.host_arg->ptc_num, cudaMemcpyDeviceToHost);
            cudaDeviceSynchronize();
            cudaMemcpy(sph.particle->accy, sph.tmp_cuda->accy, sizeof(double) * sph.host_arg->ptc_num, cudaMemcpyDeviceToHost);
            cudaDeviceSynchronize();
            cudaMemcpy(sph.particle->density, sph.tmp_cuda->rho, sizeof(double) * sph.host_arg->ptc_num, cudaMemcpyDeviceToHost);
            cudaDeviceSynchronize();
            cudaMemcpy(sph.particle->pressure, sph.tmp_cuda->p, sizeof(double) * sph.host_arg->ptc_num, cudaMemcpyDeviceToHost);
            cudaDeviceSynchronize();
            if(sph.host_arg->init_impac_flag == 0)
            {
                cudaMemcpy(sph.host_rigid,sph.dev_rigid,sizeof(SPH_RIGID),cudaMemcpyDeviceToHost);
                cudaDeviceSynchronize();
            }
            
        }
        cudaError_t sph_error = cudaGetLastError();
        printf("%s\n", cudaGetErrorName(sph_error));
    }

    //save the last frame
    cudaMemcpy(sph.particle->x, sph.tmp_cuda->x, sizeof(double) * sph.host_arg->ptc_num, cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();
    cudaMemcpy(sph.particle->y, sph.tmp_cuda->y, sizeof(double) * sph.host_arg->ptc_num, cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();
    cudaMemcpy(sph.particle->vx, sph.tmp_cuda->vx, sizeof(double) * sph.host_arg->ptc_num, cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();
    cudaMemcpy(sph.particle->vy, sph.tmp_cuda->vy, sizeof(double) * sph.host_arg->ptc_num, cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();
    cudaMemcpy(sph.particle->accx, sph.tmp_cuda->accx, sizeof(double) * sph.host_arg->ptc_num, cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();
    cudaMemcpy(sph.particle->accy, sph.tmp_cuda->accy, sizeof(double) * sph.host_arg->ptc_num, cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();
    cudaMemcpy(sph.particle->density, sph.tmp_cuda->rho, sizeof(double) * sph.host_arg->ptc_num, cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();
    cudaMemcpy(sph.particle->pressure, sph.tmp_cuda->p, sizeof(double) * sph.host_arg->ptc_num, cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();
    cudaMemcpy(sph.particle->type,sph.tmp_cuda->type,sizeof(int)*sph.host_arg->ptc_num,cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();
    cudaMemcpy(sph.host_rigid,sph.dev_rigid,sizeof(SPH_RIGID),cudaMemcpyDeviceToHost);
    sph_save_last(&sph);

    sph_free(&sph);
    cudaDeviceReset();
    return 0;
}

__global__ void sph_predict_cuda(SPH_CUDA *cuda, SPH_ARG *arg, SPH_RIGID *rigid)
{
    const int id = threadIdx.x + blockIdx.x * blockDim.x;
    if (id < arg->ptc_num)
    {
        cuda->temp_x[id] = cuda->x[id];
        cuda->temp_y[id] = cuda->y[id];
        cuda->temp_vx[id] = cuda->vx[id];
        cuda->temp_vy[id] = cuda->vy[id];
        cuda->temp_rho[id] = cuda->rho[id];
        if (cuda->type[id] == 0)
        {
            cuda->x[id] += cuda->vx[id] * arg->dt * 0.5;
            cuda->y[id] += cuda->vy[id] * arg->dt * 0.5;
            cuda->vx[id] += cuda->accx[id] * arg->dt * 0.5;
            cuda->vy[id] += cuda->accy[id] * arg->dt * 0.5;
            cuda->rho[id] += cuda->drho[id] * arg->dt * 0.5;
            if (cuda->rho[id] < arg->ref_rho)
                cuda->rho[id] = arg->ref_rho;
            cuda->p[id] = arg->c * arg->c * (cuda->rho[id] - arg->ref_rho);
        }
        else
        {
            cuda->p[id] = 0.0;
            cuda->rho[id] = 0.0;
            cuda->vx[id] = 0.0;
            cuda->vy[id] = 0.0;
            cuda->rho[id] = arg->ref_rho;
        }
    }
}

__global__ void sph_correct_cuda(SPH_CUDA *cuda, SPH_ARG *arg, SPH_RIGID *rigid)
{
    const int id = threadIdx.x + blockIdx.x * blockDim.x;
    if (id < arg->ptc_num)
    {
        if (cuda->type[id] == 0)
        {
            cuda->x[id] = cuda->temp_x[id] + cuda->vx[id] * arg->dt;
            cuda->y[id] = cuda->temp_y[id] + cuda->vy[id] * arg->dt;
            cuda->vx[id] = cuda->temp_vx[id] + cuda->accx[id] * arg->dt;
            cuda->vy[id] = cuda->temp_vy[id] + cuda->accy[id] * arg->dt;
            cuda->rho[id] = cuda->temp_rho[id] + cuda->drho[id] * arg->dt;
            if (cuda->rho[id] < arg->ref_rho)
                cuda->rho[id] = arg->ref_rho;
            cuda->p[id] = arg->c * arg->c * (cuda->rho[id] - arg->ref_rho);
        }
        else if(cuda->type[id] == 1)
        {
            cuda->p[id] = 0.0;
            cuda->rho[id] = 0.0;
            cuda->vx[id] = 0.0;
            cuda->vy[id] = 0.0;
            cuda->rho[id] = arg->ref_rho;
            cuda->x[id] = cuda->temp_x[id] + arg->dt*(rigid->vx - rigid->omega*(cuda->y[id]-rigid->cogy));
            cuda->y[id] = cuda->temp_y[id] + arg->dt*(rigid->vy + rigid->omega*(cuda->x[id]-rigid->cogx));
        }
        else if(cuda->type[id] == -1)
        {
            cuda->p[id] = 0.0;
            cuda->rho[id] = 0.0;
            cuda->vx[id] = 0.0;
            cuda->vy[id] = 0.0;
            cuda->rho[id] = arg->ref_rho;
        }
    }
}

__global__ void sph_rigid_cuda(SPH_CUDA *cuda,SPH_ARG *arg,SPH_RIGID *rigid)
{
    __shared__ double accx;
    __shared__ double accy;
    __shared__ double alpha;
    const int id = threadIdx.x + blockIdx.x * blockDim.x;
    if(id == 0)
    {
        rigid->accx = 0.0;
        rigid->accy = -arg->g;
        rigid->alpha = 0.0;
        rigid->cogx = cuda->x[rigid->cog_ptc_id];
        rigid->cogy = cuda->y[rigid->cog_ptc_id];
    }
    if(threadIdx.x == 0)
    {
        accx = 0.0;
        accy = 0.0;
        alpha = 0.0;
    }
    __syncthreads();
    if(id < arg->ptc_num)
    {
        if(cuda->type[id] == 1)
        {
            accx += cuda->accx[id]*arg->m/rigid->mass;
            accy += cuda->accy[id]*arg->m/rigid->mass;
            alpha += (cuda->accy[id]*(cuda->x[id]-rigid->cogx)-cuda->accx[id]*(cuda->y[id]-rigid->cogy))*arg->m/rigid->moi;
        }   
    }
    __syncthreads();
    if( threadIdx.x == 0)
    {
        atomicAdd(&(rigid->accx),accx);
        atomicAdd(&(rigid->accy),accy);
        atomicAdd(&(rigid->alpha),alpha);
        atomicAdd(&(rigid->vx),accx*arg->dt);
        atomicAdd(&(rigid->vy),accy*arg->dt);
        atomicAdd(&(rigid->omega),alpha*arg->dt);
    }
    __syncthreads();
}
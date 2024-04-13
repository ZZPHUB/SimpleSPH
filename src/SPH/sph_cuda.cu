#include "SPH.cuh"

int main(void)
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

    cudaSetDevice(0);
    sph_init(&sph); 

    //define the seed for ptc data structure
    dim3 ptc_block(256);
    dim3 ptc_grid((int)(sph.particle->total/256)+1);
    //define the seed for mesh data structure
    dim3 mesh_block(32,32);
    dim3 mesh_grid(MESH_LENGTH_NUM,MESH_DEEPTH_NUM);
    //define the seed for pair data structre
    dim3 pair_block(128);
    dim3 pair_grid(MESH_LENGTH_NUM,MESH_DEEPTH_NUM);


    //SPH_CUDA cuda;
    //SPH_ARG tmp_arg;
    //cudaMemcpy(&cuda,sph.cuda,sizeof(SPH_CUDA),cudaMemcpyDeviceToHost);
    //cudaDeviceSynchronize();
    //int *host_pair_count = (int *)calloc(sph.host_arg->mesh_num,sizeof(int));
    //int id = 0;
    //int *cpu_pair_i;
    //int *cpu_pair_j;
    //cudaMalloc(&cpu_pair_i,sizeof(int)*32*sph.particle->total);
    //cudaMalloc(&cpu_pair_j,sizeof(int)*32*sph.particle->total);

    for(sph.host_arg->init_step;sph.host_arg->init_step<sph.host_arg->total_step;sph.host_arg->init_step++)
    {
        printf("current step is:%d ",sph.host_arg->init_step);
        
        sph_mesh_cuda<<<ptc_grid,ptc_block>>>(sph.cuda,sph.dev_arg);
        cudaDeviceSynchronize();
        sph_nnps_cuda<<<mesh_grid,mesh_block>>>(sph.cuda,sph.dev_arg,sph.dev_rigid);
        cudaDeviceSynchronize();
        sph_kernel_cuda<<<pair_grid,pair_block>>>(sph.cuda,sph.dev_arg,sph.dev_rigid);
        cudaDeviceSynchronize();
        sph_governing_cuda<<<pair_grid,pair_block>>>(sph.cuda,sph.dev_arg,sph.dev_rigid);
        cudaDeviceSynchronize();
        sph_predict_cuda<<<ptc_grid,ptc_block>>>(sph.cuda,sph.dev_arg,sph.dev_rigid);
        cudaDeviceSynchronize();
        sph_dummy_cuda<<<pair_grid,pair_block>>>(sph.cuda,sph.dev_arg,sph.dev_rigid);
        cudaDeviceSynchronize();
        
        sph_mesh_cuda<<<ptc_grid,ptc_block>>>(sph.cuda,sph.dev_arg);
        cudaDeviceSynchronize();
        sph_nnps_cuda<<<mesh_grid,mesh_block>>>(sph.cuda,sph.dev_arg,sph.dev_rigid);
        cudaDeviceSynchronize();
        sph_kernel_cuda<<<pair_grid,pair_block>>>(sph.cuda,sph.dev_arg,sph.dev_rigid);
        cudaDeviceSynchronize();
        sph_governing_cuda<<<pair_grid,pair_block>>>(sph.cuda,sph.dev_arg,sph.dev_rigid);
        cudaDeviceSynchronize();
        sph_correct_cuda<<<ptc_grid,ptc_block>>>(sph.cuda,sph.dev_arg,sph.dev_rigid);
        cudaDeviceSynchronize();
        sph_dummy_cuda<<<pair_grid,pair_block>>>(sph.cuda,sph.dev_arg,sph.dev_rigid);
        cudaDeviceSynchronize();

        if(sph.host_arg->init_step%400 == 0)
        {
            cudaMemcpy(sph.particle->x,sph.tmp_cuda->x,sizeof(double)*sph.particle->total,cudaMemcpyDeviceToHost);
            cudaDeviceSynchronize();
            cudaMemcpy(sph.particle->y,sph.tmp_cuda->y,sizeof(double)*sph.particle->total,cudaMemcpyDeviceToHost);
            cudaDeviceSynchronize();
            cudaMemcpy(sph.particle->vx,sph.tmp_cuda->vx,sizeof(double)*sph.particle->total,cudaMemcpyDeviceToHost);
            cudaDeviceSynchronize();
            cudaMemcpy(sph.particle->vy,sph.tmp_cuda->vy,sizeof(double)*sph.particle->total,cudaMemcpyDeviceToHost);
            cudaDeviceSynchronize();
            cudaMemcpy(sph.particle->accx,sph.tmp_cuda->accx,sizeof(double)*sph.particle->total,cudaMemcpyDeviceToHost);
            cudaDeviceSynchronize();
            cudaMemcpy(sph.particle->accy,sph.tmp_cuda->accy,sizeof(double)*sph.particle->total,cudaMemcpyDeviceToHost);
            cudaDeviceSynchronize();
            cudaMemcpy(sph.particle->density,sph.tmp_cuda->rho,sizeof(double)*sph.particle->total,cudaMemcpyDeviceToHost);
            cudaDeviceSynchronize();
            cudaMemcpy(sph.particle->pressure,sph.tmp_cuda->p,sizeof(double)*sph.particle->total,cudaMemcpyDeviceToHost);
            cudaDeviceSynchronize();

            sph_save_single(&sph);
        }

    
        cudaError_t sph_error = cudaGetLastError();
        printf("%s\n",cudaGetErrorName(sph_error));
    }

    sph_free(&sph);
    cudaDeviceReset();
    return 0;
}

__global__ void sph_predict_cuda(SPH_CUDA *cuda,SPH_ARG *arg,SPH_RIGID *rigid)
{
    const int id = threadIdx.x + blockIdx.x * blockDim.x;
    if(id < arg->ptc_num)
    {
        cuda->temp_x[id] = cuda->x[id];
        cuda->temp_y[id] = cuda->y[id];
        cuda->temp_vx[id] = cuda->vx[id];
        cuda->temp_vy[id] = cuda->vy[id];
        cuda->temp_rho[id] = cuda->rho[id];
        if(cuda->type[id] == 0)
        {
            cuda->x[id] += cuda->vx[id]*arg->dt*0.5;
            cuda->y[id] += cuda->vy[id]*arg->dt*0.5;
            cuda->vx[id] += cuda->accx[id]*arg->dt*0.5;
            cuda->vy[id] += cuda->accy[id]*arg->dt*0.5;
            cuda->rho[id] += cuda->drho[id]*arg->dt*0.5;
            if(cuda->rho[id] < arg->ref_rho) cuda->rho[id] = arg->ref_rho;
            cuda->p[id] = arg->c*arg->c*(cuda->rho[id] - arg->ref_rho);
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

__global__ void sph_correct_cuda(SPH_CUDA *cuda,SPH_ARG *arg,SPH_RIGID *rigid)
{
    const int id = threadIdx.x + blockIdx.x * blockDim.x;
    if(id < arg->ptc_num)
    {
        if(cuda->type[id] == 0)
        {
            cuda->x[id] = cuda->temp_x[id] + cuda->vx[id]*arg->dt;
            cuda->y[id] = cuda->temp_y[id] + cuda->vy[id]*arg->dt;
            cuda->vx[id] = cuda->temp_vx[id] + cuda->accx[id]*arg->dt;
            cuda->vy[id] = cuda->temp_vy[id] + cuda->accy[id]*arg->dt;
            cuda->rho[id] = cuda->temp_rho[id] + cuda->drho[id]*arg->dt;
            if(cuda->rho[id] < arg->ref_rho) cuda->rho[id] = arg->ref_rho;
            cuda->p[id] = arg->c*arg->c*(cuda->rho[id] - arg->ref_rho);
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

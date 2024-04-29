#include "Lib.cuh"

__global__ void sph_mesh_cuda(SPH_CUDA *cuda,SPH_ARG *arg,SPH_RIGID *rigid)
{
    //const int bid = blockIdx.x;
    //const int tid = threadIdx.x;
    const int id = threadIdx.x + blockIdx.x * blockDim.x;
    if(id >= arg->ptc_num) return;

    /*这里需要进行加速度和密度变化的初始化*/
    
    if(cuda->type[id] == 0)
    {
        cuda->accy[id] = -arg->g;
        cuda->ptc_w[id] = arg->alpha*arg->m/(1.0*cuda->rho[id]);
    }
    else 
    {
        cuda->ptc_w[id] = 0;
        cuda->accy[id] = 0.0;
        cuda->p[id] = 0.0;
        cuda->rho[id] = arg->ref_rho;
        cuda->vx[id] = 0.0;
        cuda->vy[id] = 0.0;
    }
    cuda->accx[id] = 0.0;
    cuda->drho[id] = 0.0;
    //cuda->ptc_w[id] = 0.0;

    /*这里需要对pair_num进行初始化*/
    if(id == 0) 
    {
        //printf("the arg tmp is:%d the pair num is:%d \n",arg->tmp,arg->pair_num);
        printf("the pair num is:%d \n",arg->pair_num);
        arg->tmp = 0;
        arg->pair_num = 0;
        rigid->cogx = cuda->x[rigid->cog_ptc_id];
        rigid->cogy = cuda->y[rigid->cog_ptc_id];
    }

    int mid = 0;
    int mesh_index = 0;

    if(cuda->y[id] < arg->mesh_y && cuda->y[id] >= 0.0)
    {
        mid = __double2int_rz(cuda->y[id]/arg->mesh_dx)*arg->mesh_xnum;
    }
    else
    {
        mid = (arg->mesh_ynum - 1)*arg->mesh_xnum;
    }
    if(cuda->x[id] < arg->mesh_x && cuda->x[id] >= 0.0)
    {
        mid += __double2int_rz(cuda->x[id]/arg->mesh_dx);
    }
    else
    {
        mid += arg->mesh_xnum - 1;
    }
    //printf("mid is:%lf\n",mid);
    //printf("xnum is:%d,ynum is:%d\n",__double2int_rz(cuda->x[id]/arg->mesh_dx),__double2int_rz(cuda->y[id]/arg->mesh_dx));
    //printf("x is:%lf,y is:%lf\n",cuda->x[id],cuda->y[id]);
    mesh_index = atomicAdd(&cuda->mesh_count[mid],1);
    mesh_index = mesh_index*arg->mesh_num + mid;
    cuda->mesh[mesh_index] = id;
}


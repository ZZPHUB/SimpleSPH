#include "SPH.cuh"
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <time.h>

using namespace std;

void sph_init(SPH *sph)
{
    SPH_PARTICLE *particle;
    SPH_PAIR *pair;
    SPH_KERNEL *kernel;
    SPH_MESH *mesh;
    SPH_CUDA *temp_cuda;
    SPH_ARG *arg;
    particle = sph->particle;
    pair = sph->pair;
    kernel = sph->kernel;
    mesh = sph->mesh;
    temp_cuda = sph->tmp_cuda;
    arg = sph->host_arg;

    sph_read_info(sph);

    /************stack is too small,so init data in heap***************/
    //particle data init
    particle->x = (double *)(calloc(arg->ptc_num,sizeof(double)));
    particle->y = (double *)(calloc(arg->ptc_num,sizeof(double)));
    particle->vx = (double *)(calloc(arg->ptc_num,sizeof(double)));
    particle->vy = (double *)(calloc(arg->ptc_num,sizeof(double)));
    particle->accx = (double *)(calloc(arg->ptc_num,sizeof(double)));
    particle->accy = (double *)(calloc(arg->ptc_num,sizeof(double)));
    particle->dif_density = (double *)(calloc(arg->ptc_num,sizeof(double)));
    particle->density = (double *)(calloc(arg->ptc_num,sizeof(double)));
    particle->temp_x = (double *)(calloc(arg->ptc_num,sizeof(double)));
    particle->temp_y = (double *)(calloc(arg->ptc_num,sizeof(double)));
    particle->temp_vx = (double *)(calloc(arg->ptc_num,sizeof(double)));
    particle->temp_vy = (double *)(calloc(arg->ptc_num,sizeof(double)));
    particle->temp_density = (double *)(calloc(arg->ptc_num,sizeof(double))); 
    //particle->mass = (double *)(calloc(arg->ptc_num,sizeof(double))); 
    particle->w = (double *)(calloc(arg->ptc_num,sizeof(double)));
    particle->pressure = (double *)(calloc(arg->ptc_num,sizeof(double)));
    particle->type = (int *)(calloc(arg->ptc_num,sizeof(int)));  

    //kernel data init
    kernel->w = (double *)(calloc(arg->pair_list_num,sizeof(double)));  //this code donnot use kernel value
    kernel->dwdx = (double *)(calloc(arg->pair_list_num,sizeof(double)));
    kernel->dwdy = (double *)(calloc(arg->pair_list_num,sizeof(double)));
   
    //pair data init
    pair->total = 0; 
    pair->i = (unsigned int *)(calloc(arg->pair_list_num,sizeof(unsigned int)));
    pair->j = (unsigned int *)(calloc(arg->pair_list_num,sizeof(unsigned int)));

    //mesh data init
    mesh->ptc = (int *)calloc(sph->host_arg->mesh_num*sph->host_arg->mesh_volume,sizeof(int));
    mesh->count = (int *)calloc(sph->host_arg->mesh_num,sizeof(int));
    sph->mesh = mesh;

    sph_read_vtk(sph);
    if(arg->init_impac_flag == 0)
    {
        for(int i=0;i<arg->ptc_num;i++)
        {
            if(particle->y[i] <= arg->fluid_x)
            {
                particle->pressure[i] = arg->ref_rho*arg->g*(arg->fluid_x - particle->y[i]);
                particle->density[i] = arg->ref_rho + particle->pressure/(arg->c * arg->c);
            }
        }
    }

    cudaMalloc(&(sph->dev_arg),sizeof(SPH_ARG));
    cudaMemcpy(sph->dev_arg,sph->host_arg,sizeof(SPH_ARG),cudaMemcpyHostToDevice);

    cudaMalloc(&(sph->dev_rigid),sizeof(SPH_RIGID));
    cudaMemcpy(sph->dev_rigid,sph->host_rigid,sizeof(SPH_RIGID),cudaMemcpyHostToDevice);

    /*cuda mem alloc*/
    cudaMalloc(&(sph->cuda),sizeof(SPH_CUDA));
    cudaMalloc(&(temp_cuda->x),arg->ptc_num*sizeof(double));
    cudaMalloc(&(temp_cuda->y),arg->ptc_num*sizeof(double));
    cudaMalloc(&(temp_cuda->temp_x),arg->ptc_num*sizeof(double));
    cudaMalloc(&(temp_cuda->temp_y),arg->ptc_num*sizeof(double));
    cudaMalloc(&(temp_cuda->vx),arg->ptc_num*sizeof(double));
    cudaMalloc(&(temp_cuda->vy),arg->ptc_num*sizeof(double));
    cudaMalloc(&(temp_cuda->temp_vx),arg->ptc_num*sizeof(double));
    cudaMalloc(&(temp_cuda->temp_vy),arg->ptc_num*sizeof(double));
    cudaMalloc(&(temp_cuda->accx),arg->ptc_num*sizeof(double));
    cudaMalloc(&(temp_cuda->accy),arg->ptc_num*sizeof(double));
    cudaMalloc(&(temp_cuda->rho),arg->ptc_num*sizeof(double));
    cudaMalloc(&(temp_cuda->drho),arg->ptc_num*sizeof(double));
    cudaMalloc(&(temp_cuda->temp_rho),arg->ptc_num*sizeof(double));
    cudaMalloc(&(temp_cuda->p),arg->ptc_num*sizeof(double));
    cudaMalloc(&(temp_cuda->type),arg->ptc_num*sizeof(int));
    cudaMalloc(&(temp_cuda->ptc_w),arg->ptc_num*sizeof(double));
    cudaMalloc(&(temp_cuda->Lxx),arg->ptc_num*sizeof(double));
    cudaMalloc(&(temp_cuda->Lxy),arg->ptc_num*sizeof(double));
    cudaMalloc(&(temp_cuda->Lyx),arg->ptc_num*sizeof(double));
    cudaMalloc(&(temp_cuda->Lyy),arg->ptc_num*sizeof(double));
    cudaMalloc(&(temp_cuda->Lrho_x),arg->ptc_num*sizeof(double));
    cudaMalloc(&(temp_cuda->Lrho_y),arg->ptc_num*sizeof(double));

    cudaMalloc(&(temp_cuda->pair_w),arg->pair_list_num*sizeof(double));
    cudaMalloc(&(temp_cuda->dwdx),arg->pair_list_num*sizeof(double));
    cudaMalloc(&(temp_cuda->dwdy),arg->pair_list_num*sizeof(double));
    cudaMalloc(&(temp_cuda->pair_i),arg->pair_list_num*sizeof(int));
    cudaMalloc(&(temp_cuda->pair_j),arg->pair_list_num*sizeof(int));
    cudaMalloc(&(temp_cuda->pair_count),sph->host_arg->pair_mesh_num*sizeof(int));
    cudaMalloc(&(temp_cuda->mesh),sph->host_arg->mesh_num*sph->host_arg->mesh_volume*sizeof(int));
    cudaMalloc(&(temp_cuda->mesh_count),sph->host_arg->mesh_num*sizeof(int));

    cudaMemcpy(temp_cuda->x, particle->x, arg->ptc_num*sizeof(double), cudaMemcpyHostToDevice); 
    cudaMemcpy(temp_cuda->y, particle->y, arg->ptc_num*sizeof(double), cudaMemcpyHostToDevice); 
    cudaMemcpy(temp_cuda->vx, particle->vx, arg->ptc_num*sizeof(double), cudaMemcpyHostToDevice); 
    cudaMemcpy(temp_cuda->vy, particle->vy, arg->ptc_num*sizeof(double), cudaMemcpyHostToDevice); 
    cudaMemcpy(temp_cuda->type, particle->type, arg->ptc_num*sizeof(int), cudaMemcpyHostToDevice); 
    cudaMemcpy(temp_cuda->rho, particle->density, arg->ptc_num*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(temp_cuda->accx, particle->accx, arg->ptc_num*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(temp_cuda->accy, particle->accx, arg->ptc_num*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(temp_cuda->drho, particle->dif_density, arg->ptc_num*sizeof(double), cudaMemcpyHostToDevice);
    
    cudaMemset(temp_cuda->p,0,arg->ptc_num*sizeof(double));
    cudaMemset(temp_cuda->temp_x,0,arg->ptc_num*sizeof(double));
    cudaMemset(temp_cuda->temp_y,0,arg->ptc_num*sizeof(double));
    cudaMemset(temp_cuda->temp_vx,0,arg->ptc_num*sizeof(double));
    cudaMemset(temp_cuda->temp_vy,0,arg->ptc_num*sizeof(double));
    cudaMemset(temp_cuda->temp_rho,0,arg->ptc_num*sizeof(double));
    cudaMemset(temp_cuda->ptc_w,0,arg->ptc_num*sizeof(double));
    cudaMemset(temp_cuda->Lxx,0,arg->ptc_num*sizeof(double));
    cudaMemset(temp_cuda->Lxy,0,arg->ptc_num*sizeof(double));
    cudaMemset(temp_cuda->Lyx,0,arg->ptc_num*sizeof(double));
    cudaMemset(temp_cuda->Lyy,0,arg->ptc_num*sizeof(double));
    cudaMemset(temp_cuda->Lrho_x,0,arg->ptc_num*sizeof(double));
    cudaMemset(temp_cuda->Lrho_y,0,arg->ptc_num*sizeof(double));

    cudaMemset(temp_cuda->pair_w,0,arg->pair_list_num*sizeof(double));
    cudaMemset(temp_cuda->dwdx,0,arg->pair_list_num*sizeof(double));
    cudaMemset(temp_cuda->dwdy,0,arg->pair_list_num*sizeof(double));
    cudaMemset(temp_cuda->pair_i,0,arg->pair_list_num*sizeof(int));
    cudaMemset(temp_cuda->pair_j,0,arg->pair_list_num*sizeof(int));
    cudaMemset(temp_cuda->pair_count,0,sph->host_arg->pair_mesh_num*sizeof(int));
    cudaMemset(temp_cuda->mesh,0,sph->host_arg->mesh_num*sph->host_arg->mesh_volume*sizeof(int));
    cudaMemset(temp_cuda->mesh_count,0,sph->host_arg->mesh_num*sizeof(int));

    cudaMemcpy(sph->cuda,temp_cuda,sizeof(SPH_CUDA),cudaMemcpyHostToDevice);
}

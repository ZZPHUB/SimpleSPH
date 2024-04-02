#include "SPH.cuh"
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <time.h>
using namespace std;


int main(void)
{
    SPH_PARTICLE particle;
    SPH_KERNEL kernel;
    SPH_PAIR pair;
    SPH_RIGID wedge;
    SPH_MESH mesh = NULL;
    SPH sph;
    sph.particle = &particle;
    sph.kernel = &kernel;
    sph.pair = &pair;
    sph.rigid = &wedge;
    sph.mesh = mesh;

    sph_init(&sph); 
    cudaSetDevice(0);

    double *dev_x = NULL;
    double *dev_y = NULL;
    double *dev_vx = NULL;
    double *dev_vy = NULL;
    double *dev_rho = NULL;
    double *dev_temp_x = NULL;
    double *dev_temp_y = NULL;
    double *dev_temp_vx = NULL;
    double *dev_temp_vy = NULL;
    double *dev_temp_rho = NULL;
    double *dev_p = NULL;
    double *dev_accx = NULL;
    double *dev_accy = NULL;
    double *dev_drho = NULL;
    double *dev_w = NULL;
    int *dev_type = NULL;

    int *dev_pair_i = NULL;
    int *dev_pair_j = NULL;

    double *dev_kernel_w = NULL;
    double *dev_kernel_dwdx = NULL;
    double *dev_kernel_dwdy = NULL;
    
    int *dev_mesh = NULL;

    double *dev_rigid = NULL;
    double host_rigid[10];

    host_rigid[VX] = sph.rigid->vx;
    host_rigid[VY] = sph.rigid->vy;
    host_rigid[ACCX] = sph.rigid->accx;
    host_rigid[ACCY] = sph.rigid->accy;
    host_rigid[OMEGA] = sph.rigid->omega;
    host_rigid[R_ALPHA] = sph.rigid->alpha;
    host_rigid[MASS] = sph.rigid->mass;
    host_rigid[MOI] = sph.rigid->moi;
    host_rigid[COGX] = sph.rigid->cogx;
    host_rigid[COGY] = sph.rigid->cogy;

    int host_count;
    int *dev_count;


    CUDA_CHECK(cudaMalloc((double**)&dev_x,sph.particle->total*sizeof(double)));
    CUDA_CHECK(cudaMalloc((double**)&dev_y,sph.particle->total*sizeof(double)));
    CUDA_CHECK(cudaMalloc((double**)&dev_vx,sph.particle->total*sizeof(double)));
    CUDA_CHECK(cudaMalloc((double**)&dev_vy,sph.particle->total*sizeof(double)));
    CUDA_CHECK(cudaMalloc((double**)&dev_rho,sph.particle->total*sizeof(double)));
    CUDA_CHECK(cudaMalloc((double**)&dev_temp_x,sph.particle->total*sizeof(double)));
    CUDA_CHECK(cudaMalloc((double**)&dev_temp_y,sph.particle->total*sizeof(double)));
    CUDA_CHECK(cudaMalloc((double**)&dev_temp_vx,sph.particle->total*sizeof(double)));
    CUDA_CHECK(cudaMalloc((double**)&dev_temp_vy,sph.particle->total*sizeof(double)));
    CUDA_CHECK(cudaMalloc((double**)&dev_temp_rho,sph.particle->total*sizeof(double)));
    CUDA_CHECK(cudaMalloc((double**)&dev_p,sph.particle->total*sizeof(double)));
    CUDA_CHECK(cudaMalloc((double**)&dev_accx,sph.particle->total*sizeof(double)));
    CUDA_CHECK(cudaMalloc((double**)&dev_accy,sph.particle->total*sizeof(double)));
    CUDA_CHECK(cudaMalloc((double**)&dev_drho,sph.particle->total*sizeof(double)));
    CUDA_CHECK(cudaMalloc((double**)&dev_w,sph.particle->total*sizeof(double)));
    CUDA_CHECK(cudaMalloc((int**)&dev_type,sph.particle->total*sizeof(double)));

    CUDA_CHECK(cudaMalloc((int**)&dev_pair_i,32*sph.particle->total*sizeof(double)));
    CUDA_CHECK(cudaMalloc((int**)&dev_pair_j,32*sph.particle->total*sizeof(double)));

    CUDA_CHECK(cudaMalloc((double**)&dev_kernel_w,32*sph.particle->total*sizeof(double)));
    CUDA_CHECK(cudaMalloc((double**)&dev_kernel_dwdx,32*sph.particle->total*sizeof(double)));
    CUDA_CHECK(cudaMalloc((double**)&dev_kernel_dwdy,32*sph.particle->total*sizeof(double)));
    
    CUDA_CHECK(cudaMalloc((int**)&dev_mesh,MESH_DEEPTH_NUM*MESH_LENGTH_NUM*MESH_PTC_NUM*sizeof(int)));

    CUDA_CHECK(cudaMalloc((int**)&dev_count,sizeof(int)));
    CUDA_CHECK(cudaMemset(dev_count,0,sizeof(int)));

    CUDA_CHECK(cudaMalloc((double**)&dev_rigid,sizeof(double)*10));

    /*----------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
    CUDA_CHECK(cudaMemcpy(dev_x, particle.x, sph.particle->total*sizeof(double), cudaMemcpyHostToDevice)); 
    CUDA_CHECK(cudaMemcpy(dev_y, particle.y, sph.particle->total*sizeof(double), cudaMemcpyHostToDevice)); 
    CUDA_CHECK(cudaMemcpy(dev_vx, particle.vx, sph.particle->total*sizeof(double), cudaMemcpyHostToDevice)); 
    CUDA_CHECK(cudaMemcpy(dev_vy, particle.vy, sph.particle->total*sizeof(double), cudaMemcpyHostToDevice)); 
    CUDA_CHECK(cudaMemcpy(dev_type, particle.type, sph.particle->total*sizeof(double), cudaMemcpyHostToDevice)); 
    CUDA_CHECK(cudaMemcpy(dev_rho, particle.density, sph.particle->total*sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(dev_accx, particle.accx, sph.particle->total*sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(dev_accy, particle.accx, sph.particle->total*sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(dev_drho, particle.accx, sph.particle->total*sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(dev_rigid,host_rigid,10*sizeof(double),cudaMemcpyHostToDevice));

    //define the seed for ptc data structure
    dim3 ptc_block(256);
    dim3 ptc_grid((int)(sph.particle->total/256)+1);
    //define the seed for mesh data structure
    dim3 mesh_block(32,32);
    dim3 mesh_grid(MESH_LENGTH_NUM,MESH_DEEPTH_NUM);
    //define the seed for pair data structre
    dim3 pair_block(512);
    dim3 pair_grid((int)(sph.particle->total/16)+1);

    // sph_avg_time(&sph);
    for(int i=sph.current_step;i<sph.total_step;i++)
    {    
    /*---------------------------------------Predict Step---------------------------------------Predict Step---------------------------------------Predict Step---------------------------------------Predict Step---------------------------------------Predict Step---------------------------------------Predict Step*/
        //CUDA_CHECK(cudaMemset(dev_mesh,0,MESH_DEEPTH_NUM*MESH_LENGTH_NUM*MESH_PTC_NUM*sizeof(int)));
        sph_mesh_cuda<<<ptc_grid,ptc_block>>>(dev_x,dev_y,dev_accx,dev_accy,dev_drho,dev_type,dev_mesh,sph.particle->total);
        CUDA_CHECK(cudaDeviceSynchronize());
        //__global__ void sph_mesh_cuda(double *x,double *y,double *accx,double *accy,double *drho,int *type,int *mesh,int ptc_num)

        sph_nnps_cuda<<<mesh_grid,mesh_block>>>(dev_mesh,dev_x,dev_y,dev_type,dev_pair_i,dev_pair_j,dev_count);
        //__global__ void sph_nnps_cuda(int *mesh,double *x,double *y,int *type,int *pair_i,int *pair_j)
        CUDA_CHECK(cudaDeviceSynchronize());

        sph_kernel_cuda<<<pair_grid,pair_block>>>(dev_x,dev_y,dev_kernel_w,dev_kernel_dwdx,dev_kernel_dwdy,dev_w,dev_pair_i,dev_pair_j,dev_count);
        CUDA_CHECK(cudaDeviceSynchronize());
        //sph_kernel_cuda(double *x,double *y,double *w,double *dwdx,double *dwdy,doubel *ptc_w,int *pair_i,int *pair_j,int pair_num)

        sph_governing_cuda<<<pair_grid,pair_block>>>(dev_x,dev_y,dev_vx,dev_vy,dev_rho,dev_p,dev_type,dev_pair_i,dev_pair_j,dev_kernel_dwdx,dev_kernel_dwdy,dev_accx,dev_accy,dev_drho,dev_rigid,dev_count,sph.particle->total);
        CUDA_CHECK(cudaDeviceSynchronize());
        //__global__ void sph_governing_cuda(double * x, double * y, double * vx, double * vy, double * rho, double * p, int * type, int * pair_i, int * pair_j, double * dwdx, double * dwdy, double * accx, double * accy, double * drho, double *rigid, int pair_num, int ptc_num)
        
        sph_predict_cuda<<<ptc_grid,ptc_block>>>(dev_x,dev_y,dev_temp_x,dev_temp_y,dev_vx,dev_vy,dev_temp_vx,dev_temp_vy,dev_accx,dev_accy,dev_rho,dev_temp_rho,dev_drho,dev_p,dev_type,sph.particle->total);
        CUDA_CHECK(cudaDeviceSynchronize());
        //__global__ void sph_predict_cuda(double *x,double *y,double *temp_x,double *temp_y,double *vx,double *vy,double *temp_vx,double *temp_vy,double *accx,double *accy,double *rho,double *temp_rho,double *drho,int ptc_num)
    
        CUDA_CHECK(cudaMemcpy(sph.particle->x,dev_x,sph.particle->total*sizeof(double),cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(sph.particle->y,dev_y,sph.particle->total*sizeof(double),cudaMemcpyDeviceToHost));
        sph_save_single(&sph);
    /*---------------------------------------Correct Step---------------------------------------Correct Step---------------------------------------Correct Step---------------------------------------Correct Step---------------------------------------Correct Step---------------------------------------Correct Step*/
        /*
        sph_mesh_cuda<<<ptc_grid,ptc_block>>>(dev_x,dev_y,dev_accx,dev_accy,dev_drho,dev_type,dev_mesh,sph.particle->total);
        CUDA_CHECK(cudaDeviceSynchronize());

        sph_nnps_cuda<<<mesh_grid,mesh_block>>>(dev_mesh,dev_x,dev_y,dev_type,dev_pair_i,dev_pair_j,dev_count);
        //__global__ void sph_nnps_cuda(int *mesh,double *x,double *y,int *type,int *pair_i,int *pair_j)
        CUDA_CHECK(cudaDeviceSynchronize());

        sph_kernel_cuda<<<pair_grid,pair_block>>>(dev_x,dev_y,dev_kernel_w,dev_kernel_dwdx,dev_kernel_dwdy,dev_w,dev_pair_i,dev_pair_j,dev_count);
        CUDA_CHECK(cudaDeviceSynchronize());
        //sph_kernel_cuda(double *x,double *y,double *w,double *dwdx,double *dwdy,double *w,int *pair_i,int *pair_j,int pair_num)

        sph_governing_cuda<<<pair_grid,pair_block>>>(dev_x,dev_y,dev_vx,dev_vy,dev_rho,dev_p,dev_type,dev_pair_i,dev_pair_j,dev_kernel_dwdx,dev_kernel_dwdy,dev_accx,dev_accy,dev_drho,dev_rigid,dev_count,sph.particle->total);
        CUDA_CHECK(cudaDeviceSynchronize());
        //__global__ void sph_governing_cuda(double * x, double * y, double * vx, double * vy, double * rho, double * p, int * type, int * pair_i, int * pair_j, double * dwdx, double * dwdy, double * accx, double * accy, double * drho, double *rigid, int pair_num, int ptc_num)
        
        sph_correct_cuda<<<ptc_grid,ptc_block>>>(dev_x,dev_y,dev_temp_x,dev_temp_y,dev_vx,dev_vy,dev_temp_vx,dev_temp_vy,dev_accx,dev_accy,dev_rho,dev_temp_rho,dev_drho,dev_p,dev_type,sph.particle->total);
        CUDA_CHECK(cudaDeviceSynchronize());
        //__global__ void sph_predict_cuda(double *x,double *y,double *temp_x,double *temp_y,double *vx,double *vy,double *temp_vx,double *temp_vy,double *accx,double *accy,double *rho,double *temp_rho,double *drho,int ptc_num)

        sph_dummy_cuda<<<pair_grid,pair_block>>>(dev_x,dev_y,dev_vx,dev_vy,dev_p,dev_rho,dev_w,dev_kernel_w,dev_pair_i,dev_pair_j,dev_type,dev_rigid,dev_count);
        CUDA_CHECK(cudaDeviceSynchronize());
        //__global__ void sph_dummy_cuda(double *vx,double *vy,double *p,double *rho,double *ptc_w,double *pair_w,int *pair_i,int *pair_j,int *type,double *rigid,int *pair_num)
        */
        printf("%d\n",sph.current_step);
    }
    sph_free(&sph);
    return 0;
}

__global__ void sph_predict_cuda(double *x,double *y,double *temp_x,double *temp_y,double *vx,double *vy,double *temp_vx,double *temp_vy,double *accx,double *accy,double *rho,double *temp_rho,double *drho,double *p,int *type,int ptc_num)
{
    const int id = threadIdx.x + blockIdx.x * blockDim.x;
    if(id >= ptc_num )return;

    if(type[id] == 0)
    {
        temp_x[id] = x[id];
        temp_y[id] = y[id];
        temp_vx[id] = vx[id];
        temp_vy[id] = vy[id];
        temp_rho[id] = rho[id];

        x[id] = vx[id]*dev_dt*0.5;
        y[id] = vy[id]*dev_dt*0.5;
        vx[id] = accx[id]*dev_dt*0.5;
        vy[id] = accy[id]*dev_dt*0.5;
        rho[id] = drho[id]*dev_dt*0.5;
        if(rho[id] < REF_DENSITY) rho[id]=REF_DENSITY;
    }
    /*
    else
    {
        vx[id] = 0.0;
        vy[id] = 0.0;
        p[id] = 0.0;
    }
    */
}


__global__ void sph_correct_cuda(double *x,double *y,double *temp_x,double *temp_y,double *vx,double *vy,double *temp_vx,double *temp_vy,double *accx,double *accy,double *rho,double *temp_rho,double *drho,double *p,int *type,int ptc_num)
{
    const int id = threadIdx.x + blockIdx.x * blockDim.x;
    if(id >= ptc_num )return; 

    if(type[id] == 0)
    {
        x[id] = temp_x[id] + vx[id]*dev_dt;
        y[id] = temp_y[id] + vy[id]*dev_dt;
        vx[id] = temp_vx[id] + accx[id]*dev_dt;
        vy[id] = temp_vy[id] + accy[id]*dev_dt;
        rho[id] = temp_rho[id] + drho[id]*dev_dt;
        if(rho[id] < REF_DENSITY) rho[id]=REF_DENSITY;
    }
    else
    {
        vx[id] = 0.0;
        vy[id] = 0.0;
        p[id] = 0.0;
    }
}
/*
    string filename = "../data/postprocess/vtk/sph"; 
    filename += to_string(sph.current_step/PRINT_TIME_STEP);
    filename += ".vtk";

    ofstream vtkfile;
    vtkfile.open(filename.c_str());

    vtkfile << "# vtk DataFile Version 3.0" << endl;
    vtkfile << "sph data" << endl;
    vtkfile << "ASCII" << endl;
    vtkfile << "DATASET UNSTRUCTURED_GRID" << endl;
    vtkfile << "POINTS " << sph.particle->total << " " << "double" << endl;

    for(unsigned int i=0;i<MESH_DEEPTH_NUM;i++)
    {
        for(unsigned int j=0;j<MESH_LENGTH_NUM;j++)
        {
            temp = sph.mesh[i*MESH_LENGTH_NUM+j+MESH_LENGTH_NUM*MESH_DEEPTH_NUM*(MESH_PTC_NUM-1)];
            for(unsigned int k=0;k<temp;k++)
            {
                temp_1 = sph.mesh[i*MESH_LENGTH_NUM+j+MESH_LENGTH_NUM*MESH_DEEPTH_NUM*k];
                vtkfile << setiosflags(ios::scientific) << sph.particle->x[temp_1] << " " \
                << sph.particle->y[temp_1] << " " << 0.0 << endl;
            }
        }
    }
    vtkfile.close();*/
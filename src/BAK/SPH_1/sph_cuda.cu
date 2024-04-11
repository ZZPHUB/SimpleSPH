#include "SPH.cuh"
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <time.h>
using namespace std;

__constant__  int c = ART_SOUND_VEL;
__constant__  int rho_0 = REF_DENSITY;
__constant__  int mesh_lnum = MESH_LENGTH_NUM;
__constant__  int mesh_dnum = MESH_DEEPTH_NUM;
__constant__  int mesh_pnum = MESH_PTC_NUM;
__constant__  int mesh_spacing = MESH_SPACING;

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

    double *dev_x;
    double *dev_y;
    double *dev_vx;
    double *dev_vy;
    double *dev_rho;
    double *dev_p;
    /*
    dev_pair_i,dev_pair_j,dev_pair_accx,dev_pair_accy,dev_pair_drho = NULL;
    */
    int *dev_mesh =NULL;

    int temp = 0;

    CUDA_CHECK(cudaMalloc((double**)&dev_x,particle.total*sizeof(double)));
    CUDA_CHECK(cudaMalloc((double**)&dev_y,particle.total*sizeof(double)));
    CUDA_CHECK(cudaMalloc((double**)&dev_vx,particle.total*sizeof(double)));
    CUDA_CHECK(cudaMalloc((double**)&dev_vy,particle.total*sizeof(double)));
    CUDA_CHECK(cudaMalloc((double**)&dev_rho,particle.total*sizeof(double)));
    CUDA_CHECK(cudaMalloc((double**)&dev_p,particle.total*sizeof(double)));
/*
    CUDA_CHECK(cudaMalloc((double**)&dev_pair_i,size*sizeof(double)));
    CUDA_CHECK(cudaMalloc((double**)&dev_pair_j,size*sizeof(double)));
    CUDA_CHECK(cudaMalloc((double**)&dev_pair_accx,size*sizeof(double)));
    CUDA_CHECK(cudaMalloc((double**)&dev_pair_accy,size*sizeof(double)));
    CUDA_CHECK(cudaMalloc((double**)&dev_pair_drho,size*sizeof(double)));
    */

    CUDA_CHECK(cudaMalloc((int**)&dev_mesh,MESH_DEEPTH_NUM*MESH_LENGTH_NUM*MESH_PTC_NUM*sizeof(int)));

       
    sph_avg_time(&sph);
    for(sph.current_step;sph.current_step<sph.total_step;sph.current_step++)
    {
        cudaMemcpy((void *)dev_x, (void *)particle.x, particle.total*sizeof(double), cudaMemcpyHostToDevice); 
        cudaMemcpy((void *)dev_y, (void *)particle.y, particle.total*sizeof(double), cudaMemcpyHostToDevice);
        cudaMemcpy((void *)dev_vx, (void *)particle.vx, particle.total*sizeof(double), cudaMemcpyHostToDevice);
        cudaMemcpy((void *)dev_vy, (void *)particle.vy, particle.total*sizeof(double), cudaMemcpyHostToDevice);
        cudaMemcpy((void *)dev_rho, (void *)particle.density, particle.total*sizeof(double), cudaMemcpyHostToDevice);
        cudaMemcpy((void *)dev_p, (void *)particle.pressure, particle.total*sizeof(double), cudaMemcpyHostToDevice);

        ptc_mesh_cuda<<<384,160>>>(dev_x,dev_y,dev_mesh,particle.total);
        cudaMemcpy(mesh, dev_mesh, MESH_DEEPTH_NUM*MESH_LENGTH_NUM*MESH_PTC_NUM,cudaMemcpyDeviceToHost);


    string filename = "../data/postprocess/vtk/sph"; 
    filename += to_string(sph.current_step/PRINT_TIME_STEP);
    filename += ".vtk";

    ofstream vtkfile;
    vtkfile.open(filename.c_str());

    vtkfile << "# vtk DataFile Version 3.0" << endl;
    vtkfile << "sph data" << endl;
    vtkfile << "ASCII" << endl;
    vtkfile << "DATASET UNSTRUCTURED_GRID" << endl;
    vtkfile << "POINTS " << particle.total << " " << "double" << endl;

    for(unsigned int i=0;i<MESH_DEEPTH_NUM;i++)
    {
        for(unsigned int j=0;j<MESH_LENGTH_NUM;j++)
        {
            for(unsigned int k=0;k<MESH_PTC_NUM;k++)
            temp = mesh[i*MESH_LENGTH_NUM+j+k];
            vtkfile << setiosflags(ios::scientific) << particle.x[temp] << " " \
            << particle.y[temp] << " " << 0.0 << endl;
        }
    }
    vtkfile.close();



        /*
        if(sph.current_step%PRINT_TIME_STEP == 0)
        {
            sph_save_single(&sph);
        }
        //calculate and integration
        sph_time_integral(&sph); 
        sph_save_rigid(&sph);
        ptc_info(&sph);
        sph_avg_time(&sph);
        */
    }
    sph_save_last(&sph);
    sph_free(&sph);
    return 0;
}
/*__global__ void sph_predict_cuda(double *x,double *y,double *temp_x,double *temp_y,double *vx,double *vy,double *temp_vx,double *temp_vy,double *accx,double *accy,double *rho,double *temp_rho,double *drho,double *p,int *type,int ptc_num)
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

        x[id] += vx[id]*dev_dt*0.5;
        y[id] += vy[id]*dev_dt*0.5;
        vx[id] += accx[id]*dev_dt*0.5;
        vy[id] += accy[id]*dev_dt*0.5;
        rho[id] += drho[id]*dev_dt*0.5;
        if(rho[id] < REF_DENSITY) rho[id]=REF_DENSITY;
    }
 
    else
    {
        vx[id] = 0.0;
        vy[id] = 0.0;
        p[id] = 0.0;
    }

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
}*/
/*
        CUDA_CHECK(cudaMemcpy(sph.mesh,dev_mesh,MESH_DEEPTH_NUM*MESH_LENGTH_NUM*MESH_PTC_NUM*sizeof(int),cudaMemcpyDeviceToHost));
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
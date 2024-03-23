#include "SPH.cuh"
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <time.h>
using namespace std;

__device__ int count;

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
    double *dev_pair_i;
    double *dev_pair_j;
    int *dev_type;
    /*
    dev_pair_i,dev_pair_j,dev_pair_accx,dev_pair_accy,dev_pair_drho = NULL;
    */
    int *dev_mesh =NULL;

    int temp = 0;
    int temp_1 = 0;
    int host_count = 0;


    CUDA_CHECK(cudaMalloc((double**)&dev_x,particle.total*sizeof(double)));
    CUDA_CHECK(cudaMalloc((double**)&dev_y,particle.total*sizeof(double)));
    CUDA_CHECK(cudaMalloc((double**)&dev_vx,particle.total*sizeof(double)));
    CUDA_CHECK(cudaMalloc((double**)&dev_vy,particle.total*sizeof(double)));
    CUDA_CHECK(cudaMalloc((double**)&dev_rho,particle.total*sizeof(double)));
    CUDA_CHECK(cudaMalloc((double**)&dev_p,particle.total*sizeof(double)));
    CUDA_CHECK(cudaMalloc((int**)&dev_type,particle.total*sizeof(double)));

    CUDA_CHECK(cudaMalloc((double**)&dev_pair_i,32*particle.total*sizeof(double)));
    CUDA_CHECK(cudaMalloc((double**)&dev_pair_j,32*particle.total*sizeof(double)));
    /*
    CUDA_CHECK(cudaMalloc((double**)&dev_pair_accx,size*sizeof(double)));
    CUDA_CHECK(cudaMalloc((double**)&dev_pair_accy,size*sizeof(double)));
    CUDA_CHECK(cudaMalloc((double**)&dev_pair_drho,size*sizeof(double)));
    */

    CUDA_CHECK(cudaMalloc((int**)&dev_mesh,MESH_DEEPTH_NUM*MESH_LENGTH_NUM*MESH_PTC_NUM*sizeof(int)));

       
   // sph_avg_time(&sph);
        CUDA_CHECK(cudaMencpy(&count,&host_count,sizeof(int),cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(dev_x, particle.x, particle.total*sizeof(double), cudaMemcpyHostToDevice)); 
        CUDA_CHECK(cudaMemcpy(dev_y, particle.y, particle.total*sizeof(double), cudaMemcpyHostToDevice)); 
        CUDA_CHECK(cudaMemcpy(dev_vx, particle.vx, particle.total*sizeof(double), cudaMemcpyHostToDevice)); 
        CUDA_CHECK(cudaMemcpy(dev_vy, particle.vy, particle.total*sizeof(double), cudaMemcpyHostToDevice)); 
        CUDA_CHECK(cudaMemcpy(dev_type, particle.type, particle.total*sizeof(double), cudaMemcpyHostToDevice)); 
        CUDA_CHECK(cudaMemcpy(dev_rho, particle.density, particle.total*sizeof(double), cudaMemcpyHostToDevice)); 
        dim3 block(64,64);
        dim3 grid(MESH_LENGTH_NUM,MESH_DEEPTH_NUM);
        sph_mesh_cuda<<<384,160>>>(dev_x,dev_y,dev_mesh,particle.total);
        cudaDeviceSynchronize();
        sph_nnps_cuda<<<grid,block>>>(dev_mesh,dev_x,dev_y,dev_type,dev_pair_i,dev_pair_j);
        cudaMemcpy(&host_count,&count,sizeof(int),cudaMemcpyDeviceToHost);
        
        //__global__ void sph_nnps_cuda(int *mesh,double *x,double *y,int *type,int *pair_i,int *pair_j)
        //CUDA_CHECK(cudaMemcpy(mesh, dev_mesh, sizeof(int)*MESH_DEEPTH_NUM*MESH_LENGTH_NUM*MESH_PTC_NUM,cudaMemcpyDeviceToHost));
        printf("gpu find :%d \n",host_count);
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
    vtkfile << "POINTS " << particle.total << " " << "double" << endl;

    for(unsigned int i=0;i<MESH_DEEPTH_NUM;i++)
    {
        for(unsigned int j=0;j<MESH_LENGTH_NUM;j++)
        {
            temp = sph.mesh[i*MESH_LENGTH_NUM+j+MESH_LENGTH_NUM*MESH_DEEPTH_NUM*(MESH_PTC_NUM-1)];
            for(unsigned int k=0;k<temp;k++)
            {
                temp_1 = sph.mesh[i*MESH_LENGTH_NUM+j+MESH_LENGTH_NUM*MESH_DEEPTH_NUM*k];
                vtkfile << setiosflags(ios::scientific) << particle.x[temp_1] << " " \
                << particle.y[temp_1] << " " << 0.0 << endl;
            }
        }
    }
    vtkfile.close();

*/
    sph_free(&sph);
    return 0;
}

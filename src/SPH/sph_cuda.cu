#include "SPH.cuh"


__global__ void check_pair(SPH_ARG *arg)
{
    printf("the pair num is:%d\n",arg->pair_num);
}

int main(void)
{
    SPH_PARTICLE particle;
    SPH_KERNEL kernel;
    SPH_PAIR pair;
    SPH_RIGID wedge;
    SPH_MESH mesh = NULL;
    SPH_ARG arg;
    SPH sph;
    sph.particle = &particle;
    sph.kernel = &kernel;
    sph.pair = &pair;
    sph.host_rigid = &wedge;
    sph.host_arg = &arg;
    sph.mesh = mesh;

    cudaSetDevice(0);
    sph_init(&sph); 

    //define the seed for ptc data structure
    dim3 ptc_block(256);
    dim3 ptc_grid((int)(sph.particle->total/256)+1);
    //define the seed for mesh data structure
    dim3 mesh_block(32);
    dim3 mesh_grid(MESH_LENGTH_NUM,MESH_DEEPTH_NUM);
    //define the seed for pair data structre
    dim3 pair_block(512);
    dim3 pair_grid((int)(sph.particle->total/16)+1);

    sph_fuck_you<<<ptc_grid,ptc_block>>>(sph.cuda,sph.dev_arg);
    cudaDeviceSynchronize();
    sph_nnps_cuda<<<mesh_grid,mesh_block>>>(sph.cuda,sph.dev_arg,sph.dev_rigid);
    cudaDeviceSynchronize();
    check_pair<<<1,1>>>(sph.dev_arg);
    cudaDeviceSynchronize();

    sph_free(&sph);
    cudaDeviceReset();
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
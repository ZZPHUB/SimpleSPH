#include "Lib.cuh"
using namespace std;

__global__ void ptc_mesh_cuda(double *x,double *y,int *mesh,int ptc_num)
{
    //const int bid = blockIdx.x;
    //const int tid = threadIdx.x;
    const int id = threadIdx.x + blockIdx.x * blockDim.x;
    if(id > ptc_num) return;

    int mid;

    if(y[id] < TOL_DOMAIN_DEEPTH && y[id] >= 0)
    {
        mid = __double2int_rz(y[id]/MESH_SPACING)*MESH_LENGTH_NUM;
    }
    else if(y[id] >= TOL_DOMAIN_DEEPTH)
    {
        mid = (MESH_DEEPTH_NUM - 1)*MESH_LENGTH_NUM;
    }
    if(x[id] < TOL_DOMAIN_LENGTH && x[id] >= 0)
    {
        mid += __double2int_rz(x[id]/MESH_SPACING);
    }
    else if(x[id] >= TOL_DOMAIN_LENGTH)
    {
        mid += MESH_LENGTH_NUM - 1;
    }
    mid += MESH_DEEPTH_NUM*MESH_LENGTH_NUM*(&mesh[mid+MESH_PTC_NUM],1);
    mesh[mid] = id;
    /*
    head = mesh[j][k][MESH_PTC_NUM-1];
    if(head<MESH_PTC_NUM-1)
    {
        mesh[j][k][head] = i;        
        mesh[j][k][MESH_PTC_NUM-1]++;
    }*/
    
    

}
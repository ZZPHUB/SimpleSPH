#include "Lib.cuh"
using namespace std;

void ptc_mesh_process(SPH *sph)
//void ptc_mesh_process(SPH_PARTICLE *particle,unsigned int ***mesh)
{
    SPH_PARTICLE *particle;
    SPH_MESH mesh;
    particle = sph->particle;
    mesh = sph->mesh;

    unsigned int mesh_ptc_tol=0;
    //init the head,which store the num of particle in grid
    for(unsigned int i=0;i<MESH_DEEPTH_NUM;i++)
    {
        for(unsigned int j=0;j<MESH_LENGTH_NUM;j++)
        {
            mesh[i][j][MESH_PTC_NUM-1] = 0;
        }
    }

    //mesh process
    for(unsigned int i=0;i<particle->total;i++)
    {
        unsigned int head;
        unsigned int j=0;
        unsigned int k=0;
        if(particle->y[i] < TOL_DOMAIN_DEEPTH && particle->y[i] >= 0)
        {
            j = (unsigned int)(particle->y[i]/MESH_SPACING);
        }
        else if(particle->y[i] >= TOL_DOMAIN_DEEPTH)
        {
            j = MESH_DEEPTH_NUM - 1;
        }
        else
        {
            j = 0;
        }
        if(particle->x[i] < TOL_DOMAIN_LENGTH && particle->x[i] >= 0)
        {
            k = (unsigned int)(particle->x[i]/MESH_SPACING);
        }
        else if(particle->x[i] >= TOL_DOMAIN_LENGTH)
        {
            k = MESH_LENGTH_NUM - 1;
        }
        else
        {
            k = 0;
        }
        head = mesh[j][k][MESH_PTC_NUM-1];
        if(head<MESH_PTC_NUM-1)
        {
            mesh[j][k][head] = i;        
            mesh[j][k][MESH_PTC_NUM-1]++;
        }
    }

    for(int i=0;i<MESH_DEEPTH_NUM;i++)
    {
        for(int j=0;j<MESH_LENGTH_NUM;j++)
        {
            mesh_ptc_tol = mesh_ptc_tol+mesh[i][j][MESH_PTC_NUM-1];
        }
    }
    if(mesh_ptc_tol == particle->total)
    {
        //cout << "num of particles in meshs equal to tolal particles " << endl;
    }
    else 
    {
        while(true)
        {
            cout << "some particles are not in the mesh" << endl;
            cout << "total particles is: " << particle->total << "mesh ptc tol: " << mesh_ptc_tol << endl;
        }
    }
}


__global__ void ptc_mesh_cuda(double *x,double *y,double *mesh,int ptc_num)
{
    //const int bid = blockIdx.x;
    //const int tid = threadIdx.x;
    const int id = threadIdx.x + blockIdx.x * blockDim.x;
    if(id > ptc_num) return;

    int mid;

    if(y[id] < TOL_DOMAIN_DEEPTH && y[id] >= 0)
    {
        mid = __double2in_rz(y[id]/mesh_spacing)*mesh_lnum;
    }
    else if(y[id] >= TOL_DOMAIN_DEEPTH)
    {
        mid = (mesh_dnum - 1)*mesh_lnum;
    }
    if(x[id] < TOL_DOMAIN_LENGTH && x[id] >= 0)
    {
        mid += __double2in_rz(x[id]/mesh_spacing);
    }
    else if(particle->x[id] >= TOL_DOMAIN_LENGTH)
    {
        mid += mesh_lnum - 1;
    }
    mesh[mid+atomicAdd(&mesh[mid][mesh_pnum], 1)] = id;
    /*
    head = mesh[j][k][MESH_PTC_NUM-1];
    if(head<MESH_PTC_NUM-1)
    {
        mesh[j][k][head] = i;        
        mesh[j][k][MESH_PTC_NUM-1]++;
    }*/
    
    

}
#include "SPH.H"
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
    SPH_PAIR pair_dircet;

    
    particle.total = FLUID_PTC_NUM+VIRTUAL_PTC_NUM+solid_ptc_num(); //get all of the particle number

    /************stack is too small,so init data in heap***************/
    //particle data init
    particle.x = (double *)(calloc(particle.total,sizeof(double)));
    particle.y = (double *)(calloc(particle.total,sizeof(double)));
    particle.vx = (double *)(calloc(particle.total,sizeof(double)));
    particle.vy = (double *)(calloc(particle.total,sizeof(double)));
    particle.accx = (double *)(calloc(particle.total,sizeof(double))); 
    particle.accy = (double *)(calloc(particle.total,sizeof(double)));
    particle.density = (double *)(calloc(particle.total,sizeof(double)));
    //particle.mass = (double *)(calloc(particle.total,sizeof(double))); //for particle's mass is constant
    particle.dif_density = (double *)(calloc(particle.total,sizeof(double)));
    particle.pressure = (double *)(calloc(particle.total,sizeof(double)));
    particle.visxx = (double *)(calloc(particle.total,sizeof(double)));
    particle.visyy = (double *)(calloc(particle.total,sizeof(double)));
    particle.visxy = (double *)(calloc(particle.total,sizeof(double)));
    particle.type = (char *)(calloc(particle.total,sizeof(char)));  

    //kernel data init
    kernel.w = (double *)(calloc(5*particle.total,sizeof(double)));  //this code donnot use kernel value
    kernel.dwdx = (double *)(calloc(5*particle.total,sizeof(double)));
    kernel.dwdy = (double *)(calloc(5*particle.total,sizeof(double)));
   
    //pair data init
    pair.total = 0; 
    pair.i = (unsigned int *)(calloc(5*particle.total,sizeof(unsigned int)));
    pair.j = (unsigned int *)(calloc(5*particle.total,sizeof(unsigned int)));

    //get time current time
    time_t current_time = 0;


    //mesh data init
    unsigned int ***mesh = (unsigned int ***)(calloc(MESH_DEEPTH_NUM,sizeof(unsigned int **)));
    for(int i=0;i<MESH_DEEPTH_NUM;i++)
    {
        mesh[i] = (unsigned int **)(calloc(MESH_LENGTH_NUM,sizeof(unsigned int *)));
        for(int j=0;j<MESH_LENGTH_NUM;j++)
        {
            mesh[i][j] = (unsigned int *)(calloc(MESH_PTC_NUM,sizeof(unsigned int)));
        }
    }

    ptc_generate(&particle);    //generate the fluid solid and dummy particles
    mesh_process(&particle,mesh);   //generate the mesh 

    //count the time of the nnps_mesh with parallel
    T_START
    nnps_mesh(&particle,&pair,mesh);    //use mesh to search interactive pairs
    T_END("nnps_mesh")

    cout << "total particle num is " << particle.total << endl;
    //cout << "nnds_direct find total pair is " << pair_dircet.total << endl;
    cout << "nnds_mesh find total pair is " << pair.total << endl;
    //cout << "the total same pair is " << total << endl;

    int temp_c = -1;
    int temp_tol = 0;
    for(int i=0;i<pair.total;i++)
    {
        if(temp_c != pair.i[i])
        {
            temp_tol++;
            temp_c = pair.i[i];
        }
    }

    cout << "fluid pair is " << temp_tol << endl;
    cout << "fluid particle is " << FLUID_PTC_NUM << endl;

    T_START
    ptc_kernel_parallel(&particle,&pair,&kernel);
    T_END("ptc_kernel_parallel")

    T_START
    ptc_init(&particle); //particles init values
    T_END("ptc_init")

    T_START
    ptc_vtk_mesh(&particle,mesh);
    T_END("ptc_vtk_mesh")


    free(particle.x);
    free(particle.y);
    free(particle.vx);
    free(particle.vy);
    free(particle.accx);
    free(particle.accy);
    free(particle.density);
    free(particle.dif_density);
    //free(particle.mass);
    free(particle.pressure);
    free(particle.type);

    free(kernel.w);
    free(kernel.dwdx);
    free(kernel.dwdy);
    
    free(pair.i);
    free(pair.j);
    free(mesh);
    return 0;
}



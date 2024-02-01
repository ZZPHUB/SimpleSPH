#include "SPH.H"
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <time.h>

using namespace std;

SPH_PARTICLE particle;
SPH_KERNEL kernel;
SPH_PAIR pair;
SPH_RIGID wall;
SPH_MESH mesh = NULL;
SPH sph;

SPH* sph_init(void)
{
    particle.fulid_ptc_num = FLUID_PTC_NUM;
    particle.wall_ptc_num = WALL_PTC_NUM;
    particle.total = particle.fulid_ptc_num+ particle.wall_ptc_num; //get all of the particle number

    /************stack is too small,so init data in heap***************/
    //particle data init
    particle.x = (double *)(calloc(particle.total,sizeof(double)));
    particle.y = (double *)(calloc(particle.total,sizeof(double)));
    particle.vx = (double *)(calloc(particle.total,sizeof(double)));
    particle.vy = (double *)(calloc(particle.total,sizeof(double)));
    particle.accx = (double *)(calloc(particle.total,sizeof(double))); 
    particle.accy = (double *)(calloc(particle.total,sizeof(double)));
    particle.density = (double *)(calloc(particle.total,sizeof(double)));
    particle.temp_x = (double *)(calloc(particle.total,sizeof(double)));
    particle.temp_y = (double *)(calloc(particle.total,sizeof(double)));
    particle.temp_vx = (double *)(calloc(particle.total,sizeof(double)));
    particle.temp_vy = (double *)(calloc(particle.total,sizeof(double)));
    particle.temp_density = (double *)(calloc(particle.total,sizeof(double))); 
    //particle.mass = (double *)(calloc(particle.total,sizeof(double))); //for particle's mass is constant
    particle.w = (double *)(calloc(particle.total,sizeof(double)));
    particle.dif_density = (double *)(calloc(particle.total,sizeof(double)));
    particle.pressure = (double *)(calloc(particle.total,sizeof(double)));
    #ifdef FLAG
    particle.visxx = (double *)(calloc(particle.total,sizeof(double)));
    particle.visyy = (double *)(calloc(particle.total,sizeof(double)));
    particle.visxy = (double *)(calloc(particle.total,sizeof(double)));
    #endif
    particle.type = (int *)(calloc(particle.total,sizeof(int)));  

    //kernel data init
    kernel.w = (double *)(calloc(30*particle.total,sizeof(double)));  //this code donnot use kernel value
    kernel.dwdx = (double *)(calloc(30*particle.total,sizeof(double)));
    kernel.dwdy = (double *)(calloc(30*particle.total,sizeof(double)));
   
    //pair data init
    pair.total = 0; 
    pair.i = (unsigned int *)(calloc(30*particle.total,sizeof(unsigned int)));
    pair.j = (unsigned int *)(calloc(30*particle.total,sizeof(unsigned int)));


    /* //mesh data init
    mesh = (SPH_MESH)(calloc(MESH_DEEPTH_NUM,sizeof(unsigned int **)));
    for(int i=0;i<MESH_DEEPTH_NUM;i++)
    {
        mesh[i] = (unsigned int **)(calloc(MESH_LENGTH_NUM,sizeof(unsigned int *)));
        for(int j=0;j<MESH_LENGTH_NUM;j++)
        {
            mesh[i][j] = (unsigned int *)(calloc(MESH_PTC_NUM,sizeof(unsigned int)));
        }
    }*/
    sph.particle = &particle;
    sph.kernel = &kernel;
    sph.pair = &pair;
    sph.rigid_0 = &wall;
    //sph.mesh = mesh;
    sph.d_time = DELTA_TIME;
    sph.c = ART_SOUND_VEL;
    sph.g = GRAVITY_ACC;
    sph.flag = 0.0;
    sph.current_step = 0;
    sph.total_step = INIT_TIME_STEP;

    ptc_generate(&sph);
    ptc_init(&sph);
    
    return &sph;
}
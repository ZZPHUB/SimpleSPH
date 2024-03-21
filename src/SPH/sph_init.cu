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
    SPH_MESH mesh;
    particle = sph->particle;
    pair = sph->pair;
    kernel = sph->kernel;

    particle->fulid_ptc_num = FLUID_PTC_NUM;  //fluid ptc num
    particle->wall_ptc_num = WALL_PTC_NUM;    //wall ptc num
    particle->rigid_ptc_num = ptc_rigid_num(); //rigid ptc num
    //get all of the particle number
    particle->total = particle->fulid_ptc_num+particle->wall_ptc_num+particle->rigid_ptc_num; 

    /************stack is too small,so init data in heap***************/
    //particle data init
    particle->x = (double *)(calloc(particle->total,sizeof(double)));
    particle->y = (double *)(calloc(particle->total,sizeof(double)));
    particle->vx = (double *)(calloc(particle->total,sizeof(double)));
    particle->vy = (double *)(calloc(particle->total,sizeof(double)));
    particle->density = (double *)(calloc(particle->total,sizeof(double)));
    particle->temp_x = (double *)(calloc(particle->total,sizeof(double)));
    particle->temp_y = (double *)(calloc(particle->total,sizeof(double)));
    particle->temp_vx = (double *)(calloc(particle->total,sizeof(double)));
    particle->temp_vy = (double *)(calloc(particle->total,sizeof(double)));
    particle->temp_density = (double *)(calloc(particle->total,sizeof(double))); 
    particle->mass = (double *)(calloc(particle->total,sizeof(double))); 
    particle->w = (double *)(calloc(particle->total,sizeof(double)));
    particle->pressure = (double *)(calloc(particle->total,sizeof(double)));
    particle->type = (int *)(calloc(particle->total,sizeof(int)));  

    //kernel data init
    kernel->w = (double *)(calloc(30*particle->total,sizeof(double)));  //this code donnot use kernel value
    kernel->dwdx = (double *)(calloc(30*particle->total,sizeof(double)));
    kernel->dwdy = (double *)(calloc(30*particle->total,sizeof(double)));
   
    //pair data init
    pair->total = 0; 
    pair->i = (unsigned int *)(calloc(30*particle->total,sizeof(unsigned int)));
    pair->j = (unsigned int *)(calloc(30*particle->total,sizeof(unsigned int)));

    //mesh data init
    mesh = (SPH_MESH)(calloc(MESH_DEEPTH_NUM*MESH_LENGTH_NUM*MESH_PTC_NUM,sizeof(int)));
    /*
    for(int i=0;i<MESH_DEEPTH_NUM;i++)
    {
        mesh[i] = (unsigned int **)(calloc(MESH_LENGTH_NUM,sizeof(unsigned int *)));
        for(int j=0;j<MESH_LENGTH_NUM;j++)
        {
            mesh[i][j] = (unsigned int *)(calloc(MESH_PTC_NUM,sizeof(unsigned int)));
        }
    }
    */
    
    sph->mesh = mesh;
    sph->d_time = DELTA_TIME;
    sph->c = ART_SOUND_VEL;
    sph->g = 0.0;
    sph->avg_time = 0.0;
    
    
    cout << "run a new case or an old case(press 1 for new,0 for old)" << endl;
    cin >> sph->new_case_flag;

    if(sph->new_case_flag == 1)
    {
        sph->current_step = 0;
        sph->total_step = INIT_TIME_STEP;
    }
    else if(sph->new_case_flag == 0)
    {
        cout << "the sph current time step is: " << endl;
        cin >> sph->current_step;
        cout << "the total sph time step is: " << endl;
        cin >> sph->total_step;
    }

    cout << "run a init case or a dynamic case(press 1 for init,0 for dynamic)" << endl;
    cin >> sph->init_impac_flag;
    
    cout << "save the last step or not(press 1 ta save,0 for not)" << endl;
    cin >> sph->save_last_flag;

    ptc_generate(sph);
    ptc_init(sph);
}
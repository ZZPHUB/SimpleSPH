#include "SPH.H"
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <time.h>

using namespace std;

void ptc_time_integral(SPH_PARTICLE *,SPH_PAIR *,SPH_KERNEL *,RIGID *,RIGID *,double);

int main(void)
{ 
    SPH_PARTICLE particle;
    SPH_KERNEL kernel;
    SPH_PAIR pair;
    SPH_PAIR pair_dircet;
    RIGID wall;
    RIGID wedge;

    
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
    particle.w = (double *)(calloc(particle.total,sizeof(double)));
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

    //rigid wall init
    wall.vx=wall.vy=wall.accx=wall.accy=wall.omega=wall.alpha=wall.cogx=wall.cogy=wall.mass=0;

    //rigid body init
    wedge.vx=wedge.vy=wedge.accx=wedge.accy=wedge.omega=wedge.alpha=0;
    wedge.cogx = TOL_DOMAIN_LENGTH/2;
    wedge.cogy = 1.024+4*PTC_SPACING;
    wedge.mass = 12.8;

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
    
    double delta_time = DELTA_TIME;  //the time step length
    unsigned int time_step = INIT_TIME_STEP; //time step num
    char filename[] = "../data/postprocess/sph001.vtk"; //filename 
    double scale[4] = {0,TOL_DOMAIN_LENGTH,0,TOL_DOMAIN_DEEPTH}; //output domain scale

    while (true)
    {
        for(unsigned int i=0;i<time_step;i++)
        {
            ptc_time_integral(&particle,&pair,&kernel,&wall,&wedge,delta_time);
            if(i%PRINT_TIME_STEP == 0)
            {
                ptc_info(&particle,&pair,&wedge,i);
                filename[23] = (i/PRINT_TIME_STEP)/100 + 48;
                filename[24] = ((i/PRINT_TIME_STEP)%100)/10 + 48;
                filename[25] = ((i/PRINT_TIME_STEP)/%10) + 48;
                ptc_vtk_direct(&particle,scale,filename);
            }
        }
        system("clear");
        cout << "press 0 to kill precess or num(>100) for more steps" << endl;
        cin >> time_step;
        if(time_step == 0) break;
    }

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

void ptc_time_integral(SPH_PARTICLE *particle,SPH_PAIR *pair,SPH_KERNEL *kernel,RIGID *wall,RIGID *wedge,double d_time)
{

}

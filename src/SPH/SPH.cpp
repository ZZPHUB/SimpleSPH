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

    

    //get time current time
    time_t current_time = 0;

    //if rigid body involve
    int rigid_flag = 0;


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
    ptc_init(&particle,&wall,&wedge);    //particle info init
    
    double delta_time = DELTA_TIME;  //the time step length
    unsigned int time_step = INIT_TIME_STEP; //time step num
    char filename[] = "../data/postprocess/sph001.vtk"; //filename 
    double scale[4] = {0,TOL_DOMAIN_LENGTH,0,TOL_DOMAIN_DEEPTH}; //output domain scale
    unsigned int step = 0;

    while (true)
    {
        for(step;step<time_step;step++)
        {
            //calculate and integration
            ptc_time_integral(&particle,&pair,&kernel,mesh,&wall,&wedge,delta_time,rigid_flag); 
            if(step%PRINT_TIME_STEP == 0)
            {
                ptc_info(&particle,&pair,&wedge,step);
                filename[23] = (step/PRINT_TIME_STEP)/100 + 48;
                filename[24] = ((step/PRINT_TIME_STEP)%100)/10 + 48;
                filename[25] = ((step/PRINT_TIME_STEP)%10) + 48;
                ptc_vtk_direct(&particle,scale,filename);
            }
        }
        system("clear");
        cout << "press 0 to kill precess or num(>100) for more steps" << endl;
        cin >> time_step;
        if(time_step == 0) break;
        time_step += INIT_TIME_STEP;
        rigid_flag = 1;
        scale[2] = FLUID_DOMAIN_DEEPTH/2;
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

void ptc_time_integral(SPH_PARTICLE *particle,SPH_PAIR *pair,SPH_KERNEL *kernel,unsigned int ***mesh,RIGID *wall,RIGID *wedge,double d_time,int flag)
{
    omp_lock_t lock;
    omp_init_lock(&lock);

    //ptc_mesh_process
    ptc_mesh_process(particle,mesh);

    //ptc_nnps_mesh
    ptc_nnps_mesh(particle,pair,mesh);

    //before the kernel generate,we need init the particle->w,for it donot involve the time integration
    

    //ptc_kernel_parallel
    ptc_kernel_parallel(particle,pair,kernel);

    //ptc_dif_density
    ptc_dif_density(particle,pair,kernel,wall,wedge);

    //ptc_viscous
    ptc_viscous(particle,pair,kernel,wall,wedge);

    //ptc_acceleration
    ptc_acc(particle,pair,kernel);
    
    //for virtual particles not in time integration,so the vx,vy,density,pressure need to be init every time step
    #pragma omp parallel for num_threads(TH_NUM)
    for(unsigned int i=0;i<particle->total;i++)
    {
        if(particle->type[i]!=0)
        {
            omp_set_lock(&lock);
            particle->vx[i] = 0;
            particle->vy[i] = 0;
            particle->accx[i] = 0;
            particle->accy[i] = 0;
            particle->density[i] = 0;
            particle->pressure[i] = 0;
            omp_unset_lock(&lock);
        }
    }

    //and alse,we need init the rigid body's acceleration and angular acceleration
    wedge->accx = 0;
    wedge->accy = 0;
    wedge->alpha = 0;

    //fluid_particles_info time integral and collect the rigid body's acceleration
    #pragma omp parallel for num_threads(TH_NUM)
    for(unsigned int i=0;i<particle->total;i++)
    {
        if(particle->type[i] == 0)
        {
            //fulid paritcles density,vx,vy,x,y time integral
            omp_set_lock(&lock);
            particle->density[i] += particle->dif_density[i]*d_time;
            particle->vx[i] += particle->accx[i]*d_time;
            particle->vy[i] += particle->accy[i]*d_time;
            particle->x[i] += particle->vx[i]*d_time;
            particle->y[i] += particle->vy[i]*d_time;
            omp_unset_lock(&lock);
        }
        else if(particle->type[i] == 1)
        {
            //rigid body acceleration and angular acceleration
            omp_set_lock(&lock);
            wedge->accx += particle->accx[i]*(PTC_MASS/wedge->mass);
            wedge->accy += particle->accy[i]*(PTC_MASS/wedge->mass);
            wedge->alpha += ((particle->x[i]-wedge->cogx)*particle->accy[i]-(particle->y[i]-wedge->cogy)*particle->accx[i])*(PTC_MASS/wedge->moi);
            omp_unset_lock(&lock);
        }
    }



    //rigid body velocity and angular velocity
    #pragma omp parallel for num_threads(TH_NUM)
    for(unsigned int i=0;i<particle->total;i++)
    {
        if(particle->type[i]==1)
        {
            omp_set_lock(&lock);
            wedge->vx += wedge->accx*d_time*flag;
            wedge->vy += (wedge->accy-GRAVITY_ACC)*d_time*flag;
            wedge->omega += wedge->alpha*d_time*flag;
            omp_unset_lock(&lock);
        }
    }

    //rigid body vx,vy,omega time integration
    #pragma omp parallel for num_threads(TH_NUM)
    for(unsigned int i=0;i<particle->total;i++)
    {
        if(particle->type[i]==1)
        {
            omp_set_lock(&lock);
            particle->x[i] += (wedge->vx - wedge->omega*(particle->y[i]-wedge->cogy))*d_time;
            particle->y[i] += (wedge->vy + wedge->omega*(particle->x[i]-wedge->cogx))*d_time;
            omp_unset_lock(&lock);
        }
    }

    /* -----------------------------------They are for next step calculation ------------------------------------------*/
    //fluid pressure
    fluid_ptc_pressure(particle);

    //rigid body particle pressure 
    #pragma omp parallel for num_threads(TH_NUM)
    for(unsigned int i=0;i<pair->total;i++)
    {
        if(particle->type[pair->j[i]]==1)
        {
            omp_set_lock(&lock);
            particle->pressure[pair->j[i]] += (particle->pressure[pair->i[i]]+particle->density[pair->i[i]]*\
            ((wedge->accx-wedge->alpha*(particle->y[pair->j[i]] - wedge->cogy)-wedge->omega*pow(particle->x[pair->j[i]]-wedge->cogx,2))*(particle->x[pair->i[i]]-particle->x[pair->j[i]])*kernel->w[i]+\
             (wedge->accy+wedge->alpha*(particle->x[pair->j[i]] - wedge->cogx)-wedge->omega*pow(particle->y[pair->j[i]]-wedge->cogy,2))*(particle->y[pair->i[i]]-particle->y[pair->j[i]])*kernel->w[i]))/particle->w[pair->j[i]];
            
            //correct virtual particles velocity for viscous calculation
            particle->vx[pair->j[i]] += particle->vx[pair->i[i]]/particle->w[pair->j[i]];
            particle->vy[pair->j[i]] += particle->vy[pair->i[i]]/particle->w[pair->j[i]];
            omp_unset_lock(&lock);
        }
        else if(particle->type[pair->j[i]]==-1)
        {
            //virtual particles pressure
            omp_set_lock(&lock);
            particle->pressure[pair->j[i]] += (particle->pressure[pair->i[i]]+particle->density[pair->i[i]]*\
            ((wall->accx-wall->alpha*(particle->y[pair->j[i]] - wall->cogy)-wall->omega*pow(particle->x[pair->j[i]]-wall->cogx,2))*(particle->x[pair->i[i]]-particle->x[pair->j[i]])*kernel->w[i]+\
             (wall->accy+wall->alpha*(particle->x[pair->j[i]] - wall->cogx)-wall->omega*pow(particle->y[pair->j[i]]-wall->cogy,2)+GRAVITY_ACC)*(particle->y[pair->i[i]]-particle->y[pair->j[i]])*kernel->w[i]))/particle->w[pair->j[i]];
            
            //correct virtual particle velocity for viscous calculation
            particle->vx[pair->j[i]] += particle->vx[pair->i[i]]/particle->w[pair->j[i]];
            particle->vy[pair->j[i]] += particle->vy[pair->i[i]]/particle->w[pair->j[i]];
            omp_unset_lock(&lock);
        }
    }

    //rigid body particles density
    #pragma omp parallel for num_threads(TH_NUM)
    for(unsigned int i=0;i<particle->total;i++)
    {
        if(particle->type[i]!=0)
        {
            omp_set_lock(&lock);
            particle->density[i] = particle->pressure[i]/pow(ART_SOUND_VEL,2) + REF_DENSITY;
            omp_unset_lock(&lock);
        }
    }
}

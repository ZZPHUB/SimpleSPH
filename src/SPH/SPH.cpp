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
    SPH_RIGID wall;
    SPH_RIGID wedge;
    SPH_MESH mesh = NULL;
    SPH sph;

    
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


    //mesh data init
    mesh = (SPH_MESH)(calloc(MESH_DEEPTH_NUM,sizeof(unsigned int **)));
    for(int i=0;i<MESH_DEEPTH_NUM;i++)
    {
        mesh[i] = (unsigned int **)(calloc(MESH_LENGTH_NUM,sizeof(unsigned int *)));
        for(int j=0;j<MESH_LENGTH_NUM;j++)
        {
            mesh[i][j] = (unsigned int *)(calloc(MESH_PTC_NUM,sizeof(unsigned int)));
        }
    }
    sph.particle = &particle;
    sph.kernel = &kernel;
    sph.pair = &pair;
    sph.rigid_0 = &wall;
    sph.rigid_1 = &wedge;
    sph.mesh = mesh;
    sph.d_time = DELTA_TIME;
    sph.c = ART_SOUND_VEL;
    sph.g = GRAVITY_ACC;
    sph.flag = 0.0;
    sph.current_step = 0;
    sph.total_step = INIT_TIME_STEP;



    ptc_generate(&sph);    //generate the fluid solid and dummy particles
    ptc_init(&sph);    //particle info init
    
    char filename[] = "../data/postprocess/sph000.vtk"; //filename 

    while (true)
    {
        for(sph.current_step;sph.current_step<sph.total_step;sph.current_step++)
        {
            if(sph.current_step%PRINT_TIME_STEP == 0)
            {
                filename[23] = (sph.current_step/PRINT_TIME_STEP)/100 + 48;
                filename[24] = ((sph.current_step/PRINT_TIME_STEP)%100)/10 + 48;
                filename[25] = ((sph.current_step/PRINT_TIME_STEP)%10) + 48;
                ptc_vtk_direct(&sph,filename);
            }
            //calculate and integration
            ptc_time_integral(&sph); 
            ptc_info(&sph);
        }
        system("clear");
        cout << "press 0 to kill precess or num(>100) for more steps" << endl;
        cin >> sph.total_step;
        if(sph.total_step == 0) break;
        sph.total_step += INIT_TIME_STEP;
        sph.flag = 1;
        wedge.vy = -5.0;
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
void ptc_density_correct(SPH *sph)
{
    SPH_PARTICLE *particle;
    SPH_PAIR *pair;
    SPH_KERNEL *kernel;
    particle = sph->particle;
    pair = sph->pair;
    kernel = sph->kernel;

    double m = PTC_MASS;    

    omp_lock_t lock;
    omp_init_lock(&lock);

    #pragma omp parallel for num_threads(TH_NUM)
    for(unsigned int i=0;i<particle->total;i++)
    {
        if(particle->type[i] == 0)
        {
            particle->density[i] = m*2.0*ALPHA/(3.0*particle->w[i]);
        }
    }

    #pragma omp parallel for num_threads(TH_NUM)
    for(unsigned int i=0;i<pair->total;i++)
    {
        particle->density[pair->i[i]] += m*kernel->w[i]/particle->w[pair->i[i]];
        if(particle->type[pair->j[i]] == 0)
        {
            particle->density[pair->j[i]] += m*kernel->w[i]/particle->w[pair->j[i]];
        }
    }
}

void ptc_dummy(SPH *sph)
{
    SPH_PARTICLE *particle;
    SPH_PAIR *pair;
    SPH_KERNEL *kernel;
    SPH_RIGID *wall;
    SPH_RIGID *wedge;
    particle = sph->particle;
    pair = sph->pair;
    kernel = sph->kernel;
    wall = sph->rigid_0;
    wedge = sph->rigid_1;

    omp_lock_t lock;
    omp_init_lock(&lock);
    //rigid body(wall & wedge)vx,vy,accx,accy,pressure init
    #pragma omp parallel for num_threads(TH_NUM)
    for(unsigned int i=0;i<particle->total;i++)
    {
        if(particle->type[i] != 0)
        {
            omp_set_lock(&lock);
            particle->vx[i] = 0;
            particle->vy[i] = 0;
            particle->pressure[i] = 0;
            particle->density[i] = 0;
            omp_unset_lock(&lock);
        }
    }
    //rigid body(wall & wedge) particle pressure 
    #pragma omp parallel for num_threads(TH_NUM)
    for(unsigned int i=0;i<pair->total;i++)
    {
        if(particle->type[pair->j[i]] != 0)
        {
            if(particle->w[pair->j[i]] != 0)
            {
                omp_set_lock(&lock);
                particle->pressure[pair->j[i]] += (particle->pressure[pair->i[i]]*kernel->w[i] +particle->density[pair->i[i]]*\
                ((wedge->accx-wedge->alpha*(particle->y[pair->j[i]] - wedge->cogy)-(particle->x[pair->j[i]]-wedge->cogx)*pow(wedge->omega,2))*(particle->x[pair->i[i]]-particle->x[pair->j[i]])*kernel->w[i]+\
                (wedge->accy+wedge->alpha*(particle->x[pair->j[i]] - wedge->cogx)-(particle->y[pair->j[i]]-wedge->cogy)*pow(wedge->omega,2))*(particle->y[pair->i[i]]-particle->y[pair->j[i]])*kernel->w[i]))/(particle->w[pair->j[i]]);
                //particle->pressure[pair->j[i]] += (particle->pressure[pair->i[i]]*kernel->w[i])/particle->w[pair->j[i]];
            
                //correct virtual particles velocity for viscous calculation
                particle->vx[pair->j[i]] += particle->vx[pair->i[i]]*kernel->w[i]/(particle->w[pair->j[i]]);
                particle->vy[pair->j[i]] += particle->vy[pair->i[i]]*kernel->w[i]/(particle->w[pair->j[i]]);
                omp_unset_lock(&lock);
            }
            else
            {
                //while(true) cout << pair->j[i] << " type is " <<particle->type[pair->j[i]] << " w = 0 " << endl;
            }
        }
    }
    //rigid body(wall & wedege) pressure 
    #pragma omp parallel for num_threads(TH_NUM)
    for(unsigned int i=0;i<particle->total;i++)
    {
        if(particle->type[i] != 0)
        {
            omp_set_lock(&lock);
            particle->density[i] = particle->pressure[i]/pow(sph->c,2)+REF_DENSITY;
            omp_unset_lock(&lock);
        }
    }
}

void ptc_time_integral(SPH *sph)
{
    SPH_PARTICLE *particle;
    SPH_PAIR *pair;
    SPH_KERNEL *kernel;
    SPH_RIGID *wall;
    SPH_RIGID *wedge;
    particle = sph->particle;
    pair = sph->pair;
    kernel = sph->kernel;
    wall = sph->rigid_0;
    wedge = sph->rigid_1;

    omp_lock_t lock;
    omp_init_lock(&lock);

    //generate mesh for nnps search
    ptc_mesh_process(sph);
    //search the interact particles pair
    ptc_nnps_mesh(sph);
    //generate the particle kernel value and differential kernel value
    ptc_kernel_parallel(sph);
    //get the particle pressure by eos
    fluid_ptc_pressure(sph);
    //get the ptc density change rate
    ptc_dif_density(sph);
    //get the acceleration of ptc 
    ptc_acc(sph);

    //PREDICT STEP
    #pragma omp parallel for num_threads(TH_NUM)
    for(unsigned int i=0;i<particle->total;i++)
    {
        omp_set_lock(&lock);
        particle->temp_x[i] = particle->x[i];
        particle->temp_y[i] = particle->y[i];
        particle->temp_vx[i] = particle->vx[i];
        particle->temp_vy[i] = particle->vy[i];
        particle->temp_density[i] = particle->density[i];
        omp_unset_lock(&lock);
    }
    #pragma omp parallel for num_threads(TH_NUM)
    for(unsigned int i=0;i<particle->total;i++)
    {
        if(particle->type[i] == 0)
        {
            omp_set_lock(&lock);
            particle->x[i] = particle->temp_x[i] + particle->temp_vx[i]*sph->d_time/2.0;
            particle->y[i] = particle->temp_y[i] + particle->temp_vy[i]*sph->d_time/2.0;
            particle->vx[i] = particle->temp_vx[i] + particle->accx[i]*sph->d_time/2.0;
            particle->vy[i] = particle->temp_vy[i] + particle->accy[i]*sph->d_time/2.0;
            particle->density[i] = particle->temp_density[i] + particle->dif_density[i]*sph->d_time/2.0;
            if(particle->density[i]<REF_DENSITY) particle->density[i] = REF_DENSITY;
            omp_unset_lock(&lock);
        }
    }
    //get rigid body's pressure and velosity
    ptc_dummy(sph);

    //generate mesh for nnps search
    ptc_mesh_process(sph);
    //search the interact particles pair
    ptc_nnps_mesh(sph);
    //generate the particle kernel value and differential kernel value
    ptc_kernel_parallel(sph);
    //get the particle pressure by eos
    fluid_ptc_pressure(sph);
    //get the ptc density change rate
    ptc_dif_density(sph);
    //get the acceleration of ptc 
    ptc_acc(sph);

    //CORRECT STEP
    #pragma omp parallel for num_threads(TH_NUM)
    for(unsigned int i=0;i<particle->total;i++)
    {
        if(particle->type[i] == 0)
        {
            omp_set_lock(&lock);
            particle->x[i] = particle->temp_x[i] + particle->vx[i]*sph->d_time;
            particle->y[i] = particle->temp_y[i] + particle->vy[i]*sph->d_time;
            particle->vx[i] = particle->temp_vx[i] + particle->accx[i]*sph->d_time;
            particle->vy[i] = particle->temp_vy[i] + particle->accy[i]*sph->d_time;
            particle->density[i] = particle->temp_density[i] + particle->dif_density[i]*sph->d_time;
            if(particle->density[i] < REF_DENSITY) particle->density[i] = REF_DENSITY; 
            omp_unset_lock(&lock);
        }
    }

    //get rigid body's pressure and velosity
    ptc_dummy(sph);
    
    //DENSITY CORRECT STEP
    if(sph->current_step%10 == 0)
    {
        ptc_density_correct(sph);
    }


    //and alse,we need init the rigid body's acceleration and angular acceleration
    wedge->accx = 0;
    wedge->accy = 0;
    wedge->alpha = 0;
    //collect the rigid body's acceleration
    #pragma omp parallel for num_threads(TH_NUM)
    for(unsigned int i=0;i<particle->total;i++)
    {
        if(particle->type[i] == 1)
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
            wedge->vx += wedge->accx*sph->d_time*sph->flag;
            wedge->vy += (wedge->accy-GRAVITY_ACC)*sph->d_time*sph->flag;
            wedge->omega += wedge->alpha*sph->d_time*sph->flag;
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
            particle->x[i] += (wedge->vx - wedge->omega*(particle->y[i]-wedge->cogy))*sph->d_time;
            particle->y[i] += (wedge->vy + wedge->omega*(particle->x[i]-wedge->cogx))*sph->d_time;
            omp_unset_lock(&lock);
        }
    }
}

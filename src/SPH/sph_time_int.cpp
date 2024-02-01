#include "SPH.H"
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <time.h>

void ptc_density_correct(SPH *sph)
{
    SPH_PARTICLE *particle;
    SPH_PAIR *pair;
    SPH_KERNEL *kernel;
    particle = sph->particle;
    pair = sph->pair;
    kernel = sph->kernel;

    double m = PTC_MASS;    
 
    for(unsigned int i=0;i<particle->total;i++)
    {
        if(particle->type[i] == 0)
        {
            particle->density[i] = m*2.0*ALPHA/(3.0*particle->w[i]);
        }
    }

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
 
    for(unsigned int i=0;i<particle->total;i++)
    {
        if(particle->type[i] != 0)
        {
            particle->vx[i] = 0;
            particle->vy[i] = 0;
            particle->pressure[i] = 0;
            particle->density[i] = 0;
        }
    }
    //rigid body(wall & wedge) particle pressure 
 
    for(unsigned int i=0;i<pair->total;i++)
    {
        if(particle->type[pair->j[i]] == -1)
        {
            if(particle->w[pair->j[i]] != 0)
            {
                particle->pressure[pair->j[i]] += (particle->pressure[pair->i[i]]*kernel->w[i] +particle->density[pair->i[i]]*\
                (GRAVITY_ACC)*(particle->y[pair->i[i]]-particle->y[pair->j[i]])*kernel->w[i])/(particle->w[pair->j[i]]);

                //correct virtual particles velocity for viscous calculation
                particle->vx[pair->j[i]] += particle->vx[pair->i[i]]*kernel->w[i]/(particle->w[pair->j[i]]);
                particle->vy[pair->j[i]] += particle->vy[pair->i[i]]*kernel->w[i]/(particle->w[pair->j[i]]);
            }
            else
            {
                //while(true) cout << pair->j[i] << " type is " <<particle->type[pair->j[i]] << " w = 0 " << endl;
            }
        }
    }

    //rigid body(wall & wedege) pressure 
    for(unsigned int i=0;i<particle->total;i++)
    {
        if(particle->type[i] != 0)
        {
            particle->density[i] = particle->pressure[i]/pow(sph->c,2)+REF_DENSITY;
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


    //search the interact particles pair
    ptc_nnps_direct(sph);
    //generate the particle kernel value and differential kernel value
    ptc_kernel_parallel(sph);
    //get the particle pressure by eos
    fluid_ptc_pressure(sph);
    //get the ptc density change rate
    ptc_dif_density(sph);
    //get the acceleration of ptc 
    ptc_acc(sph);

    //PREDICT STEP
    for(unsigned int i=0;i<particle->total;i++)
    {
        particle->temp_x[i] = particle->x[i];
        particle->temp_y[i] = particle->y[i];
        particle->temp_vx[i] = particle->vx[i];
        particle->temp_vy[i] = particle->vy[i];
        particle->temp_density[i] = particle->density[i];
    }
 
    for(unsigned int i=0;i<particle->total;i++)
    {
        if(particle->type[i] == 0)
        {
            particle->x[i] = particle->temp_x[i] + particle->temp_vx[i]*sph->d_time/2.0;
            particle->y[i] = particle->temp_y[i] + particle->temp_vy[i]*sph->d_time/2.0;
            particle->vx[i] = particle->temp_vx[i] + particle->accx[i]*sph->d_time/2.0;
            particle->vy[i] = particle->temp_vy[i] + particle->accy[i]*sph->d_time/2.0;
            particle->density[i] = particle->temp_density[i] + particle->dif_density[i]*sph->d_time/2.0;
            if(particle->density[i]<REF_DENSITY) particle->density[i] = REF_DENSITY;
        }
    }
    //get rigid body's pressure and velosity
    ptc_dummy(sph);


    //search the interact particles pair
    ptc_nnps_direct(sph);
    //generate the particle kernel value and differential kernel value
    ptc_kernel_parallel(sph);
    //get the particle pressure by eos
    fluid_ptc_pressure(sph);
    //get the ptc density change rate
    ptc_dif_density(sph);
    //get the acceleration of ptc 
    ptc_acc(sph);

    //CORRECT STEP
    for(unsigned int i=0;i<particle->total;i++)
    {
        if(particle->type[i] == 0)
        {
            particle->x[i] = particle->temp_x[i] + particle->vx[i]*sph->d_time;
            particle->y[i] = particle->temp_y[i] + particle->vy[i]*sph->d_time;
            particle->vx[i] = particle->temp_vx[i] + particle->accx[i]*sph->d_time;
            particle->vy[i] = particle->temp_vy[i] + particle->accy[i]*sph->d_time;
            particle->density[i] = particle->temp_density[i] + particle->dif_density[i]*sph->d_time;
            if(particle->density[i] < REF_DENSITY) particle->density[i] = REF_DENSITY; 
        }
    }

    //get rigid body's pressure and velosity
    ptc_dummy(sph);
    
    //DENSITY CORRECT STEP
    if(sph->current_step%10 == 0)
    {
        ptc_density_correct(sph);
    }
}

#include "SPH.H"
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <time.h>

void sph_time_integral(SPH *sph)
{
    SPH_PARTICLE *particle;
    SPH_PAIR *pair;
    SPH_KERNEL *kernel;
    particle = sph->particle;
    pair = sph->pair;
    kernel = sph->kernel;

    //generate the nnps mesh
    ptc_mesh_process(sph);
    //search the interact particles pair
    ptc_nnps_mesh(sph);
    //generate the particle kernel value and differential kernel value
    ptc_kernel_parallel(sph);
    //correct the ptc density 
    if(sph->current_step%20 == 0)
    {
        ptc_density_correct(sph);
    }
    //get the particle pressure by eos
    ptc_fluid_pressure(sph);
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

    //get the particle pressure by eos
    ptc_fluid_pressure(sph);
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
    sph_rigid_integral(sph);
}

void sph_rigid_integral(SPH *sph)
{
    SPH_PARTICLE *particle;
    SPH_PAIR *pair;
    SPH_RIGID *wedge;
    particle = sph->particle;
    pair = sph->pair;
    wedge = sph->rigid;

    double m = PTC_MASS;
    double temp_vx = 0.0;
    double temp_vy = 0.0;

    if(sph->init_impac_flag == 0)
    {
        //wedge acceleration and alpha init
        wedge->accx = wedge->alpha = 0.0;
        wedge->accy = -9.80;

        //calculate the wedge's acceleration and alpha
        for(unsigned int i=0;i<particle->total;i++)
        {
            if(particle->type[i] == 1)
            {
                wedge->accx += particle->accx[i]*m/wedge->mass;
                wedge->accy += particle->accy[i]*m/wedge->mass;
                wedge->alpha += (particle->accy[i]*(particle->x[i]-wedge->cogx)-\
                                 particle->accx[i]*(particle->y[i]-wedge->cogy))*m/wedge->mass;
            }
        }
        //rigid ptc time integral
        for(unsigned int i=0;i<particle->total;i++)
        {
            if(particle->type[i] == 1)
            {
                temp_vx = wedge->vx - wedge->omega*(particle->y[i]-wedge->cogy);
                temp_vy = wedge->vy + wedge->omega*(particle->x[i]-wedge->cogx);
                temp_vx += (wedge->accx-pow(wedge->omega,2)*(particle->x[i]-wedge->cogx)-wedge->alpha*(particle->y[i]-wedge->cogy))*sph->d_time;
                temp_vy += (wedge->accy-pow(wedge->omega,2)*(particle->y[i]-wedge->cogy)+wedge->alpha*(particle->x[i]-wedge->cogx))*sph->d_time;
                particle->x[i] += temp_vx*sph->d_time;
                particle->y[i] += temp_vy*sph->d_time;
            }
        }

        //calculate the center of gravity's moving
        wedge->vx += wedge->accx*sph->d_time;
        wedge->vy += wedge->accy*sph->d_time;
        wedge->omega += wedge->alpha*sph->d_time;
        wedge->cogx += wedge->vx*sph->d_time;
        wedge->cogy += wedge->vy*sph->d_time;
    }
}
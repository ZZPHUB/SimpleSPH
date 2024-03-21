#include "PreProcess.cuh"

void ptc_info_init(SPH *sph)
{
    SPH_PARTICLE *particle;
    SPH_PAIR *pair;
    particle = sph->particle;
    pair = sph->pair;

    double free_surf = PTC_SPACING*(FLUID_DEEPTH_NUM-1+4);

    if(sph->new_case_flag == 1)
    {
        for(int i=0;i<particle->total;i++)
        {
            particle->vx[i] = particle->vy[i] = 0;
            #ifndef ANALYSIS
                particle->pressure[i] = 0.0;
                particle->density[i] = REF_DENSITY;
                particle->mass[i] = PTC_MASS;
            #else
                if(particle->y[i] < free_surf)
                {
                    particle->pressure[i] = REF_DENSITY*GRAVITY_ACC*(free_surf-particle->y[i]);
                    particle->density[i] = particle->pressure[i]/pow(sph->c,2)+REF_DENSITY;
                    particle->mass[i] = particle->density[i]*pow(PTC_SPACING,2);
                }
                else
                {
                    particle->pressure[i] = 0.0;
                    particle->density[i] = REF_DENSITY;
                    particle->mass[i] = PTC_MASS;
                }
            #endif
        }
    }
}

void ptc_init(SPH *sph)
{
    ptc_info_init(sph);
    ptc_rigid_init(sph);
}
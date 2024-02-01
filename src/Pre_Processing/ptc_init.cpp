#include "PreProcess.H"

void ptc_init(SPH *sph)
{
    SPH_PARTICLE *particle;
    SPH_PAIR *pair;
    SPH_RIGID *wall;
    particle = sph->particle;
    pair = sph->pair;
    wall = sph->rigid_0;

    for(int i=0;i<particle->total;i++)
    {
        particle->vx[i] = particle->vy[i] = particle->accx[i] = \
        particle->accy[i] =particle->pressure[i] =particle->dif_density[i] = 0;
        particle->density[i] = REF_DENSITY;
    }

    //rigid wall init
    wall->vx=wall->vy=wall->accx=wall->accy=wall->omega=wall->alpha=wall->cogx=wall->cogy=wall->mass=0;

}

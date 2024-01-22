#include "PreProcess.H"

void ptc_init(SPH *sph)
//void ptc_init(SPH_PARTICLE *particle,RIGID *wall,RIGID *wedge)
{
    SPH_PARTICLE *particle;
    SPH_PAIR *pair;
    SPH_RIGID *wall;
    SPH_RIGID *wedge;
    particle = sph->particle;
    pair = sph->pair;
    wall = sph->rigid_0;
    wedge = sph->rigid_1;

    omp_lock_t lock;
    omp_init_lock(&lock);
    #pragma omp parallel for num_threads(TH_NUM)
    for(int i=0;i<particle->total;i++)
    {
        omp_set_lock(&lock);
        particle->vx[i] = particle->vy[i] = particle->accx[i] = \
        particle->accy[i] =particle->pressure[i] =particle->dif_density[i] = 0;
        //particle->mass[i] = REF_DENSITY*pow(PTC_SPACING,3);
        particle->density[i] = REF_DENSITY;
        omp_unset_lock(&lock);
    }

    //rigid body init
    wedge->vx=wedge->vy=wedge->accx=wedge->accy=wedge->omega=wedge->alpha=0;
    wedge->cogx = TOL_DOMAIN_LENGTH/2;
    wedge->cogy = 1.024+4*PTC_SPACING;
    wedge->mass = 12.8;
    wedge->moi = 0;

    //rigid wall init
   wall->vx=wall->vy=wall->accx=wall->accy=wall->omega=wall->alpha=wall->cogx=wall->cogy=wall->mass=0;

    unsigned int temp = 0;
    temp = solid_ptc_num();
    #pragma omp parallel for num_threads(TH_NUM)
    for(unsigned int i = FLUID_PTC_NUM+VIRTUAL_PTC_NUM;i<particle->total;i++)
    {
        if(particle->type[i]==1)
        {
            omp_set_lock(&lock);
            wedge->moi += (double)(wedge->mass/temp)*sqrt(pow(particle->x[i]-wedge->cogx,2)+pow(particle->y[i]-wedge->cogy,2));
            omp_unset_lock(&lock);
        }
    }
}

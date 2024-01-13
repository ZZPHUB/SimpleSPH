#include "Equations.H"

void fluid_ptc_pressure(SPH *sph)
{
    SPH_PARTICLE *particle;
    particle = sph->particle;
    
    omp_lock_t lock; //it seems that this donnot need a lock
    omp_init_lock(&lock);
    double c = ART_SOUND_VEL; //ariti_sound_velocity

    #pragma omp parallel for num_threads(TH_NUM)
    for(unsigned int i=0;i<particle->total;i++)
    {
        if(particle->type[i]==0)
        {
            /* the formula is $p=c \times (\rho-\rho_{ref})$ */
            omp_set_lock(&lock);
            particle->pressure[i] = c*c*(particle->density[i]-REF_DENSITY); 
            omp_unset_lock(&lock);
        }
    }
}
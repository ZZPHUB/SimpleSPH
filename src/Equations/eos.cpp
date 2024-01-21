#include "Equations.H"

void fluid_ptc_pressure(SPH *sph)
{
    SPH_PARTICLE *particle;
    particle = sph->particle;
    
    omp_lock_t lock; //it seems that this donnot need a lock
    omp_init_lock(&lock);
    double c = sph->c; //ariti_sound_velocity

    #ifdef LINEAR_EOS
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
    #else
    #pragma omp parallel for num_threads(TH_NUM)
    for(unsigned int i=0;i<particle->total;i++)
    {
        if(particle->type[i]==0)
        {
            omp_set_lock(&lock);
            particle->pressure[i] = (pow(ART_SOUND_VEL,2)*REF_DENSITY/7)*(pow(particle->density[i]/REF_DENSITY,7)-1);
            omp_unset_lock(&lock);
        }
    }
    #endif
}
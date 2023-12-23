#include "Equations.H"

void ptc_pressure(SPH_PARTICLE *particle)
{
    //omp_lock_t lock; //it seems that this donnot need a lock
    //omp_init_lock(&lock);
    double c = ART_SOUND_VEL; //ariti_sound_velocity

    #pragma omp parallel for num_threads(6)
    for(unsigned int i=0;i<PTC_TOL_NUM;i++)
    {
        if(particle->type[i]==0 || particle->type[i]==-1)
        {
            /* the formula is $p=c \times (\rho-\rho_{ref})$ */
            particle->pressure[i] = c*c*(particle->density[i]-REF_DENSITY); 
        }
    }
}
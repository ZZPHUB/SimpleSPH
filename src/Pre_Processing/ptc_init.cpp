#include "PreProcess.H"

void ptc_init(SPH_PARTICLE *particle)
{
    omp_lock_t lock;
    omp_init_lock(&lock);
    #pragma omp parallel for num_threads(6)
    for(int i=0;i<PTC_TOL_NUM;i++)
    {
        omp_set_lock(&lock);
        particle->vx[i] = particle->vy[i] = particle->accx[i] = \
        particle->accy[i] =particle->pressure[i] =particle->dif_density[i] = 0;
        particle->mass[i] = 1000.0*pow(PTC_SPACING,3);
        particle->density[i] = 1000;
        omp_unset_lock(&lock);
    }
}

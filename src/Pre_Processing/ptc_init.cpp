#include "PreProcess.H"

void ptc_init(SPH_PARTICLE *particle)
{
    #pragma omp parallel sections num_threads(3)
    {
        #pragma omp section
        {
            for(int i=0;i<PTC_TOL_NUM;i++)
            {
                particle->vx[i] = particle->vy[i] = particle->accx[i] = particle->accy[i] =particle->pressure[i] =particle->dif_density[i] = 0;
            }
        }
        #pragma omp section
        {
            for(int i=0;i<PTC_TOL_NUM;i++)
            {
                particle->mass[i] = 1000.0*pow(PTC_SPACING,3);
            }
        }
        #pragma omp section
        {
            for(int i=0;i<PTC_TOL_NUM;i++)
            {
                particle->density[i] = 1000;
            }
        }
    }
}

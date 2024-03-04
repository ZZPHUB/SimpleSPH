#include "PreProcess.H"

void ptc_info_init(SPH *sph)
{
    SPH_PARTICLE *particle;
    SPH_PAIR *pair;
    particle = sph->particle;
    pair = sph->pair;

    if(sph->new_case_flag == 1)
    {
        #pragma omp parallel for num_threads(TH_NUM)
        for(int i=0;i<particle->total;i++)
        {
            particle->vx[i] = particle->vy[i] = particle->accx[i] = \
            particle->accy[i]  = particle->dif_density[i] = 0;
            particle->pressure[i] = 0.0;
            particle->density[i] = REF_DENSITY;
        }
    }
}

void ptc_init(SPH *sph)
{
    ptc_info_init(sph);
    ptc_rigid_init(sph);
}
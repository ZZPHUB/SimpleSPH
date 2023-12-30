#include "Equations.H"

void ptc_dif_density(SPH_PARTICLE *particle,SPH_PAIR *pair,SPH_KERNEL *kernel)
{
    omp_lock_t lock;
    omp_init_lock(&lock);
    double m = PTC_MASS;
    double temp = 0; //particles velosity differentiation
    
    #pragma omp parallel for num_threads(TH_NUM)
    for(unsigned int i=0;i<pair->total;i++)
    {   
        if(particle->type[pair->i[i]]==0)    
        {
            omp_set_lock(&lock);
            temp = (particle->vx[pair->i[i]]-particle->vx[pair->j[i]])*kernel->dwdx[i] \
            + (particle->vy[pair->i[i]]-particle->vy[pair->j[i]])*kernel->dwdy[i];
            particle->dif_density[pair->i[i]] = particle->dif_density[pair->i[i]]+m*temp;
            particle->dif_density[pair->j[i]] = particle->dif_density[pair->j[i]]+m*temp;
            omp_unset_lock(&lock);
        }
    }
}
#include "Equations.H"

void ptc_dif_density(SPH *sph)
//void ptc_dif_density(SPH_PARTICLE *particle,SPH_PAIR *pair,SPH_KERNEL *kernel,RIGID *wall,RIGID *wedge)
{
    SPH_PARTICLE *particle;
    SPH_PAIR *pair;
    SPH_KERNEL *kernel;
    SPH_RIGID *wedge;
    particle = sph->particle;
    pair = sph->pair;
    kernel = sph->kernel;
    wedge = sph->rigid;

    //double m = PTC_MASS;
    omp_lock_t lock;
    omp_init_lock(&lock);

    #pragma omp parallel for num_threads(TH_NUM)
    for(unsigned int i=0;i<particle->total;i++)
    {
        particle->dif_density[i] = 0.0;
    }
    
    #pragma omp parallel for num_threads(TH_NUM)
    for(unsigned int i=0;i<pair->total;i++)
    {   
        double temp = 0.0; //particles velocity differentiation
        if(particle->type[pair->j[i]]==0)
        {
            temp = (particle->vx[pair->i[i]]-particle->vx[pair->j[i]])*kernel->dwdx[i] \
            + (particle->vy[pair->i[i]]-particle->vy[pair->j[i]])*kernel->dwdy[i];
        }
        else if(particle->type[pair->j[i]]==-1)
        {
            temp = (particle->vx[pair->i[i]])*kernel->dwdx[i] \
            + (particle->vy[pair->i[i]])*kernel->dwdy[i];
        }
        else if(particle->type[pair->j[i]]==1)
        {
            temp = (particle->vx[pair->i[i]]-(wedge->vx-wedge->omega*(particle->y[pair->j[i]]-wedge->cogy)))*kernel->dwdx[i] \
                    +(particle->vy[pair->i[i]]-(wedge->vy+wedge->omega*(particle->x[pair->j[i]]-wedge->cogx)))*kernel->dwdy[i];
        }
        omp_set_lock(&lock);
        particle->dif_density[pair->i[i]] += particle->mass[pair->j[i]]*temp;
        particle->dif_density[pair->j[i]] += particle->mass[pair->i[i]]*temp;
        omp_unset_lock(&lock);
    }
}
#include "Equations.H"

void ptc_dif_density(SPH *sph)
//void ptc_dif_density(SPH_PARTICLE *particle,SPH_PAIR *pair,SPH_KERNEL *kernel,RIGID *wall,RIGID *wedge)
{
    SPH_PARTICLE *particle;
    SPH_PAIR *pair;
    SPH_KERNEL *kernel;
    SPH_RIGID *wall;
    SPH_RIGID *wedge;
    particle = sph->particle;
    pair = sph->pair;
    kernel = sph->kernel;
    wall = sph->rigid_0;
    wedge = sph->rigid_1;

    omp_lock_t lock;
    omp_init_lock(&lock);
    double m = PTC_MASS;
    double temp = 0; //particles velosity differentiation

    #pragma omp parallel for num_threads(TH_NUM)
    for(unsigned int i=0;i<particle->total;i++)
    {
        omp_set_lock(&lock);
        particle->dif_density[i] = 0;
        omp_unset_lock(&lock);
    }
    
    #pragma omp parallel for num_threads(TH_NUM)
    for(unsigned int i=0;i<pair->total;i++)
    {   
        omp_set_lock(&lock);
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
        else
        {
            temp = (particle->vx[pair->i[i]] - (wedge->vx - wedge->omega*(particle->y[pair->j[i]] - wedge->cogy)))*kernel->dwdx[i] \
            + (particle->vy[pair->i[i]] - (wedge->vy + wedge->omega*(particle->x[pair->j[i]] - wedge->cogx)))*kernel->dwdy[i];
        }
        particle->dif_density[pair->i[i]] += m*temp;
        particle->dif_density[pair->j[i]] += m*temp;
        omp_unset_lock(&lock);
    }
}
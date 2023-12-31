#include "Equations.H"

void ptc_viscous(SPH_PARTICLE *particle,SPH_PAIR *pair,SPH_KERNEL *kernel,RIGID *wall,RIGID *wedge)
{
    omp_lock_t lock;
    omp_init_lock(&lock);
    double m = PTC_MASS;
    double div_vx = 0;
    double div_vy = 0;
    #pragma omp parallel for num_threads(TH_NUM)
    for(unsigned int i=0;i<pair->total;i++)
    {
        omp_set_lock(&lock);
        if(particle->type[pair->j[i]]==0)
        {
            div_vx = (particle->vx[pair->j[i]]-particle->vx[pair->i[i]])*kernel->dwdx[i];
            div_vy = (particle->vy[pair->j[i]]-particle->vy[pair->i[i]])*kernel->dwdy[i];
        }
        else if(particle->type[pair->j[i]]==-1)
        {
            div_vx = (2*(wall->vx - wall->omega*(particle->y[pair->j[i]] - wall->cogy)) - particle->vx[pair->j[i]] - particle->vx[pair->i[i]])*kernel->dwdx[i];
            div_vy = (2*(wall->vy + wall->omega*(particle->x[pair->j[i]] - wall->cogx)) - particle->vy[pair->j[i]] - particle->vy[pair->i[i]])*kernel->dwdy[i];
        }
        else if(particle->type[pair->j[i]]==1)
        {
            div_vx = (2*(wedge->vx - wedge->omega*(particle->y[pair->j[i]] - wedge->cogy)) - particle->vx[pair->j[i]] - particle->vx[pair->i[i]])*kernel->dwdx[i];
            div_vy = (2*(wedge->vy + wedge->omega*(particle->x[pair->j[i]] - wedge->cogx)) - particle->vy[pair->j[i]] - particle->vy[pair->i[i]])*kernel->dwdy[i];
        }

        particle->visxx[pair->i[i]] = particle->visxx[pair->i[i]] + \
        (4.0*m)/(3.0*particle->density[pair->j[i]])* div_vx - \
        (2.0*m)/(3.0*particle->density[pair->j[i]])* div_vy;

        particle->visxx[pair->j[i]] = particle->visxx[pair->j[i]] + \
        (4.0*m)/(3.0*particle->density[pair->i[i]])* div_vx - \
        (2.0*m)/(3.0*particle->density[pair->i[i]])* div_vy;

        particle->visyy[pair->i[i]] = particle->visyy[pair->i[i]] + \
        (4.0*m)/(3.0*particle->density[pair->j[i]])*div_vy - \
        (2.0*m)/(3.0*particle->density[pair->j[i]])*div_vx;

        particle->visyy[pair->j[i]] = particle->visyy[pair->j[i]] + \
        (4.0*m)/(3.0*particle->density[pair->i[i]])*div_vy - \
        (2.0*m)/(3.0*particle->density[pair->i[i]])*div_vx;

        particle->visxy[pair->i[i]] = particle->visxy[pair->i[i]] + \
        (m/particle->density[pair->j[i]])*(div_vx+div_vy);

        particle->visxy[pair->j[i]] = particle->visxy[pair->j[i]] + \
        (m/particle->density[pair->i[i]])*(div_vx+div_vy);
        omp_unset_lock(&lock);
    }
}
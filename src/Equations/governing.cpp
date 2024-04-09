#include "Equations.H"

void sph_governing(SPH *sph)
{
    SPH_PARTICLE *particle;
    SPH_PAIR *pair;
    SPH_KERNEL *kernel;
    SPH_RIGID *wedge;
    particle = sph->particle;
    pair = sph->pair;
    kernel = sph->kernel;
    wedge = sph->rigid;

    #pragma omp parallel for num_threads(TH_NUM)
    for(unsigned int i=0;i<particle->total;i++)
    {
        particle->dif_density[i] = 0.0;
        particle->accx[i] = 0.0;
        if(particle->type[i] == 0)
        {
            particle->accy[i] = -sph->g;
            particle->pressure[i] = sph->c*sph->c*(particle->density[i]-REF_DENSITY); 
        }
        else particle->accy[i] = 0.0;
    }


    //artificial viscous term
    #pragma omp parallel for num_threads(TH_NUM)
    for(unsigned int i=0;i<pair->total;i++)
    {
        double temp_acc = 0.0;
        double temp_rho = 0.0;
        double temp_p = 0.0;
        double dx = 0;
        double dy = 0;
        double dvx = 0;
        double dvy = 0;

        //omp_lock_t lock;
        //omp_init_lock(&lock);

        //to reduce the calculate tiems,define $p_i/\rho_i^2 + p_j/\rho_j^2$ to temp para
        temp_p = particle->pressure[pair->i[i]]/pow(particle->density[pair->i[i]],2) + particle->pressure[pair->j[i]]/pow(particle->density[pair->j[i]],2);
        
        dx = particle->x[pair->i[i]] - particle->x[pair->j[i]];
        dy = particle->y[pair->i[i]] - particle->y[pair->j[i]];
        if(particle->type[pair->j[i]] == 0)
        {
            dvx = particle->vx[pair->i[i]] - particle->vx[pair->j[i]];
            dvy = particle->vy[pair->i[i]] - particle->vy[pair->j[i]];
            temp_rho = (particle->vx[pair->i[i]]-particle->vx[pair->j[i]])*kernel->dwdx[i] \
            + (particle->vy[pair->i[i]]-particle->vy[pair->j[i]])*kernel->dwdy[i];
        }
        else if(particle->type[pair->j[i]] == -1) 
        {
            dvx = particle->vx[pair->i[i]] - (0.0 - particle->vx[pair->j[i]]);
            dvy = particle->vy[pair->i[i]] - (0.0 - particle->vy[pair->j[i]]);
            temp_rho = (particle->vx[pair->i[i]])*kernel->dwdx[i] \
            + (particle->vy[pair->i[i]])*kernel->dwdy[i];
        }
        else if(particle->type[pair->j[i]] == 1)
        {
            dvx = particle->vx[pair->i[i]] - (2.0*(wedge->vx-wedge->omega*(particle->y[pair->j[i]]-wedge->cogy)) - particle->vx[pair->j[i]]);
            dvy = particle->vy[pair->i[i]] - (2.0*(wedge->vy+wedge->omega*(particle->x[pair->j[i]]-wedge->cogx)) - particle->vy[pair->j[i]]);
            temp_rho = (particle->vx[pair->i[i]]-(wedge->vx-wedge->omega*(particle->y[pair->j[i]]-wedge->cogy)))*kernel->dwdx[i] \
                    +(particle->vy[pair->i[i]]-(wedge->vy+wedge->omega*(particle->x[pair->j[i]]-wedge->cogx)))*kernel->dwdy[i];
        }
        temp_acc = (dx*dvx+dy*dvy)/ \
               ((dx*dx+dy*dy+0.01*PTC_SML*PTC_SML)*\
               (particle->density[pair->i[i]]/2.0 + particle->density[pair->j[i]]/2.0));
        
        if(temp_acc < 0.0) temp_acc = 0.0;
        temp_acc *= 0.01*PTC_SML*sph->c;
        
        //omp_set_lock(&lock);
        #pragma omp atomic
        particle->accx[pair->i[i]] += particle->mass[pair->j[i]]*(temp_acc - temp_p)*kernel->dwdx[i];
        #pragma omp atomic
        particle->accx[pair->j[i]] -= particle->mass[pair->i[i]]*(temp_acc - temp_p)*kernel->dwdx[i];
        #pragma omp atomic
        particle->accy[pair->i[i]] += particle->mass[pair->j[i]]*(temp_acc - temp_p)*kernel->dwdy[i];
        #pragma omp atomic
        particle->accy[pair->j[i]] -= particle->mass[pair->i[i]]*(temp_acc - temp_p)*kernel->dwdy[i];
        #pragma omp atomic
        particle->dif_density[pair->i[i]] += particle->mass[pair->j[i]]*temp_rho;
        #pragma omp atomic
        particle->dif_density[pair->j[i]] += particle->mass[pair->i[i]]*temp_rho;      
        //omp_unset_lock(&lock);
    }
}
#include "Lib.H"

void ptc_density_correct(SPH *sph)
{
    SPH_PARTICLE *particle;
    SPH_PAIR *pair;
    SPH_KERNEL *kernel;
    particle = sph->particle;
    pair= sph->pair;
    kernel = sph->kernel;

    double a = ALPHA;
    //double m = PTC_MASS;

    omp_lock_t lock;
    omp_init_lock(&lock);

    #pragma omp parallel for num_threads(TH_NUM)
    for(unsigned int i=0;i<particle->total;i++)
    {
        if(particle->type[i] == 0)
        {
            particle->w[i] = (a*2.0*particle->mass[i])/(3.0*particle->density[i]);
        }
    }

    #pragma omp parallel for num_threads(TH_NUM)
    for(unsigned int i=0;i<pair->total;i++)
    {
        omp_set_lock(&lock);
        particle->w[pair->i[i]] += kernel->w[i]*particle->mass[pair->j[i]]/particle->density[pair->j[i]];
        if(particle->type[pair->j[i]]==0)
        {
            particle->w[pair->j[i]] += kernel->w[i]*particle->mass[pair->i[i]]/particle->density[pair->i[i]];
        }
        omp_unset_lock(&lock);
    }

    #pragma omp parallel for num_threads(TH_NUM)
    for(unsigned int i=0;i<particle->total;i++)
    {
        if(particle->type[i] == 0)
        {
            particle->density[i] = (a*2.0*particle->mass[i])/(3.0*particle->w[i]);
        }
    }

    #pragma omp parallel for num_threads(TH_NUM)
    for(unsigned int i=0;i<pair->total;i++)
    {
        omp_set_lock(&lock);
        particle->density[pair->i[i]] += particle->mass[pair->j[i]]*kernel->w[i]/particle->w[pair->i[i]];
        if(particle->type[pair->j[i]] == 0)
        {
            particle->density[pair->j[i]] += particle->mass[pair->i[i]]*kernel->w[i]/particle->w[pair->j[i]];
        }
        omp_unset_lock(&lock);
    }
}

void ptc_dummy(SPH *sph)
{
    SPH_PARTICLE *particle;
    SPH_PAIR *pair;
    SPH_KERNEL *kernel;
    SPH_RIGID *wedge;
    particle = sph->particle;
    pair = sph->pair;
    kernel = sph->kernel;
    wedge = sph->rigid;

    omp_lock_t lock;
    omp_init_lock(&lock);

    //rigid body(wall & wedge)vx,vy,accx,accy,pressure init
    #pragma omp parallel for num_threads(TH_NUM)
    for(unsigned int i=0;i<particle->total;i++)
    {
        if(particle->type[i] != 0)
        {
            particle->w[i] = 0;
            particle->vx[i] = 0;
            particle->vy[i] = 0;
            particle->pressure[i] = 0;
            particle->density[i] = 0;
        }
    }
    
    //the not fluid weight term 
    #pragma omp parallel for num_threads(TH_NUM)
    for(unsigned int i=0;i<pair->total;i++)
    {
        if(particle->type[pair->j[i]] != 0) 
        {
            omp_set_lock(&lock);
            particle->w[pair->j[i]] += kernel->w[i];
            omp_set_lock(&lock);
        }
    }
    
    //rigid body(wall & wedege) pressure and velocity
    #pragma omp parallel for num_threads(TH_NUM)
    for(unsigned int i=0;i<pair->total;i++)
    {
        double dx = 0.0;
        double dy = 0.0;
        double rigid_acc_x = 0.0;
        double rigid_acc_y = 0.0;
        if(particle->type[pair->j[i]] != 0 && particle->w[pair->j[i]] != 0.0)
        {
            if(particle->type[pair->j[i]] == -1)
            {
                rigid_acc_x = 0.0;
                rigid_acc_y = 0.0;
            }
            else if (particle->type[pair->j[i]] == 1)
            {
                rigid_acc_x = wedge->accx - pow(wedge->omega,2)*(particle->x[pair->j[i]]-wedge->cogx)- \
                              wedge->alpha*(particle->y[pair->j[i]]-wedge->cogy);
                rigid_acc_y = wedge->accy - pow(wedge->omega,2)*(particle->y[pair->j[i]]-wedge->cogy)+ \
                              wedge->alpha*(particle->x[pair->j[i]]-wedge->cogx);
            }
            dx = particle->x[pair->i[i]] - particle->x[pair->j[i]];
            dy = particle->y[pair->i[i]] - particle->y[pair->j[i]];
            omp_set_lock(&lock);
            particle->pressure[pair->j[i]] += (particle->pressure[pair->i[i]]+particle->density[pair->i[i]]*\
                        (rigid_acc_x*dx+(rigid_acc_y+GRAVITY_ACC)*dy))*kernel->w[i]/particle->w[pair->j[i]];
            particle->vx[pair->j[i]] += particle->vx[pair->i[i]]*kernel->w[i]/particle->w[pair->j[i]];
            particle->vy[pair->j[i]] += particle->vy[pair->i[i]]*kernel->w[i]/particle->w[pair->j[i]];
            omp_unset_lock(&lock);
        }
    }

    //rigid body(wall & wedege) densiy
    #pragma omp parallel for num_threads(TH_NUM)
    for(unsigned int i=0;i<particle->total;i++)
    {
        if(particle->type[i] != 0)
        {
            particle->density[i] = particle->pressure[i]/pow(sph->c,2)+REF_DENSITY;
        }
    }
}

#include "Equations.H"

void ptc_acc(SPH *sph)
//void ptc_acc(SPH_PARTICLE *particle,SPH_PAIR *pair,SPH_KERNEL *kernel)
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
    
    double m = PTC_MASS;
    double temp = 0.0;
    double temp_p = 0.0;
    double temp_rho_i = 0.0;
    double temp_rho_j = 0.0;
    double dx = 0;
    double dy = 0;
    double dvx = 0;
    double dvy = 0;

    
    for(unsigned int i=0;i<particle->total;i++)
    {
        particle->accx[i] = 0.0;
        if(particle->type[i] == 0)particle->accy[i] = -sph->g;
        else particle->accy[i] == 0.0;
    }

    for(unsigned int i=0;i<pair->total;i++)
    {
        //to reduce the calculate times,define $\rho^2$ to temp para
        temp_rho_i = pow(particle->density[pair->i[i]],2);
        temp_rho_j = pow(particle->density[pair->j[i]],2);

        //to reduce the calculate tiems,define $p_i/\rho_i^2 + p_j/\rho_j^2$ to temp para
        temp_p = particle->pressure[pair->i[i]]/temp_rho_i + particle->pressure[pair->j[i]]/temp_rho_j;
        particle->accx[pair->i[i]] -= m*temp_p*kernel->dwdx[i];
        particle->accx[pair->j[i]] += m*temp_p*kernel->dwdx[i];
        particle->accy[pair->i[i]] -= m*temp_p*kernel->dwdy[i];
        particle->accy[pair->j[i]] += m*temp_p*kernel->dwdy[i];
    }
    #ifdef FLAG
    #pragma omp parallel for num_threads(TH_NUM)
    for(unsigned int i=0;i<pair->total;i++)
    {
        omp_set_lock(&lock);

        //to reduce the calculate times,define $\rho^2$ to temp para
        temp_rho_i = pow(particle->density[pair->i[i]],2);
        temp_rho_j = pow(particle->density[pair->j[i]],2);

        //to reduce the calculate tiems,define $p_i/\rho_i^2 + p_j/\rho_j^2$ to temp para
        temp_p = particle->pressure[pair->i[i]]/temp_rho_i + particle->pressure[pair->j[i]]/temp_rho_j;

        particle->accx[pair->i[i]] = particle->accx[pair->i[i]]+ \
        MU*m*(particle->visxx[pair->i[i]]/temp_rho_i + particle->visxx[pair->j[i]]/temp_rho_j)*kernel->dwdx[i] + \
        MU*m*(particle->visxy[pair->i[i]]/temp_rho_i + particle->visxy[pair->j[i]]/temp_rho_j)*kernel->dwdy[i];

        particle->accx[pair->j[i]] = particle->accx[pair->j[i]] - \
        MU*m*(particle->visxx[pair->i[i]]/temp_rho_i + particle->visxx[pair->j[i]]/temp_rho_j)*kernel->dwdx[i] - \
        MU*m*(particle->visxy[pair->i[i]]/temp_rho_i + particle->visxy[pair->j[i]]/temp_rho_j)*kernel->dwdy[i];

        particle->accy[pair->i[i]] = particle->accy[pair->i[i]] + \
        MU*m*(particle->visyy[pair->i[i]]/temp_rho_i + particle->visyy[pair->j[i]]/temp_rho_j)*kernel->dwdy[i] + \
        MU*m*(particle->visxy[pair->i[i]]/temp_rho_i + particle->visxy[pair->j[i]]/temp_rho_j)*kernel->dwdx[i];

        particle->accy[pair->j[i]] = particle->accy[pair->j[i]] - \
        MU*m*(particle->visyy[pair->i[i]]/temp_rho_i + particle->visyy[pair->j[i]]/temp_rho_j)*kernel->dwdy[i] - \
        MU*m*(particle->visxy[pair->i[i]]/temp_rho_i + particle->visxy[pair->j[i]]/temp_rho_j)*kernel->dwdx[i];

        omp_unset_lock(&lock);
    }
    #else
    for(unsigned int i=0;i<pair->total;i++)
    {
        dx = particle->x[pair->i[i]] - particle->x[pair->j[i]];
        dy = particle->y[pair->i[i]] - particle->y[pair->j[i]];
        if(particle->type[pair->j[i]] == 0)
        {
            dvx = particle->vx[pair->i[i]] - particle->vx[pair->j[i]];
            dvy = particle->vy[pair->i[i]] - particle->vy[pair->j[i]];
        }
        else if(particle->type[pair->j[i]] == -1) 
        {
            dvx = particle->vx[pair->i[i]] - (0.0 - particle->vx[pair->j[i]]);
            dvy = particle->vy[pair->i[i]] - (0.0 - particle->vy[pair->j[i]]);
        }
        temp = (dx*dvx+dy*dvy)/ \
               ((dx*dx+dy*dy+0.01*PTC_SML*PTC_SML)*\
               (particle->density[pair->i[i]]/2.0 + particle->density[pair->j[i]]/2.0));
        
        particle->accx[pair->i[i]] += m*0.01*PTC_SML*sph->c*temp*kernel->dwdx[i];
        particle->accx[pair->j[i]] -= m*0.01*PTC_SML*sph->c*temp*kernel->dwdx[i];
        particle->accy[pair->i[i]] += m*0.01*PTC_SML*sph->c*temp*kernel->dwdy[i];
        particle->accy[pair->j[i]] -= m*0.01*PTC_SML*sph->c*temp*kernel->dwdy[i];
    }
    #endif
}
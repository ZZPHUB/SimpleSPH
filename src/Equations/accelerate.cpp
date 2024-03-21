#include "Equations.H"

void ptc_acc(SPH *sph)
//void ptc_acc(SPH_PARTICLE *particle,SPH_PAIR *pair,SPH_KERNEL *kernel)
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
    
    for(unsigned int i=0;i<particle->total;i++)
    {
        particle->accx[i] = 0.0;
        if(particle->type[i] == 0)particle->accy[i] = -sph->g;
        else particle->accy[i] = 0.0;
    }

    //pressure term
    for(unsigned int i=0;i<pair->total;i++)
    {
        //to reduce the calculate times,define $\rho^2$ to temp para
        double temp_p = 0.0;
        double temp_rho_i = 0.0;
        double temp_rho_j = 0.0;
        temp_rho_i = pow(particle->density[pair->i[i]],2);
        temp_rho_j = pow(particle->density[pair->j[i]],2);

        //to reduce the calculate tiems,define $p_i/\rho_i^2 + p_j/\rho_j^2$ to temp para
        temp_p = particle->pressure[pair->i[i]]/temp_rho_i + particle->pressure[pair->j[i]]/temp_rho_j;

        particle->accx[pair->i[i]] -= particle->mass[pair->j[i]]*temp_p*kernel->dwdx[i];
        particle->accx[pair->j[i]] += particle->mass[pair->i[i]]*temp_p*kernel->dwdx[i];
        particle->accy[pair->i[i]] -= particle->mass[pair->j[i]]*temp_p*kernel->dwdy[i];
        particle->accy[pair->j[i]] += particle->mass[pair->i[i]]*temp_p*kernel->dwdy[i];
    }

    //artificial viscous term
    for(unsigned int i=0;i<pair->total;i++)
    {
        double temp = 0.0;
        double dx = 0;
        double dy = 0;
        double dvx = 0;
        double dvy = 0;
        
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
        else if(particle->type[pair->j[i]] == 1)
        {
            dvx = particle->vx[pair->i[i]] - (2.0*(wedge->vx-wedge->omega*(particle->y[pair->j[i]]-wedge->cogy)) - particle->vx[pair->j[i]]);
            dvy = particle->vy[pair->i[i]] - (2.0*(wedge->vy+wedge->omega*(particle->x[pair->j[i]]-wedge->cogx)) - particle->vy[pair->j[i]]);
        }
        temp = (dx*dvx+dy*dvy)/ \
               ((dx*dx+dy*dy+0.01*PTC_SML*PTC_SML)*\
               (particle->density[pair->i[i]]/2.0 + particle->density[pair->j[i]]/2.0));
        
        if(temp < 0.0) temp = 0.0;
        
        particle->accx[pair->i[i]] += particle->mass[pair->j[i]]*0.01*PTC_SML*sph->c*temp*kernel->dwdx[i];
        particle->accx[pair->j[i]] -= particle->mass[pair->i[i]]*0.01*PTC_SML*sph->c*temp*kernel->dwdx[i];
        particle->accy[pair->i[i]] += particle->mass[pair->j[i]]*0.01*PTC_SML*sph->c*temp*kernel->dwdy[i];
        particle->accy[pair->j[i]] -= particle->mass[pair->i[i]]*0.01*PTC_SML*sph->c*temp*kernel->dwdy[i];      
    }
}
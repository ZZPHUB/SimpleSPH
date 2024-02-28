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
    double m = PTC_MASS;

    for(unsigned int i=0;i<particle->total;i++)
    {
        if(particle->type[i] == 0)
        {
            particle->w[i] = (a*2.0*m)/(3.0*particle->density[i]);
        }
    }

    for(unsigned int i=0;i<pair->total;i++)
    {
        particle->w[pair->i[i]] += kernel->w[i]*m/particle->density[pair->j[i]];
        if(particle->type[pair->j[i]]==0)
        {
            particle->w[pair->j[i]] += kernel->w[i]*m/particle->density[pair->i[i]];
        }
    }

    for(unsigned int i=0;i<particle->total;i++)
    {
        if(particle->type[i] == 0)
        {
            particle->density[i] = (a*2.0*m)/(3.0*particle->w[i]);
        }
    }

    for(unsigned int i=0;i<pair->total;i++)
    {
        particle->density[pair->i[i]] += m*kernel->w[i]/particle->w[pair->i[i]];
        if(particle->type[pair->j[i]] == 0)
        {
            particle->density[pair->j[i]] += m*kernel->w[i]/particle->w[pair->j[i]];
        }
    }
}


/*void ptc_density_correct(SPH *sph)
{
    SPH_PARTICLE *particle;
    SPH_PAIR *pair;
    SPH_KERNEL *kernel;
    particle = sph->particle;
    pair = sph->pair;
    kernel = sph->kernel;

    double m = PTC_MASS;    
 
    for(unsigned int i=0;i<particle->total;i++)
    {
        particle->w[i] = 0;
    }

    for(unsigned int i=0;i<pair->total;i++)
    {
        particle->w[pair->i[i]] += kernel->w[i]*m/particle->density[pair->j[i]];
        if(particle->type[pair->j[i]]==0) particle->w[pair->j[i]] += kernel->w[i]*m/particle->density[pair->i[i]];
    }

    for(unsigned int i=0;i<particle->total;i++)
    {
        if(particle->type[i] == 0)
        {
            particle->density[i] = m*2.0*ALPHA/(3.0*particle->w[i]);
        }
    }

    for(unsigned int i=0;i<pair->total;i++)
    {
        particle->density[pair->i[i]] += m*kernel->w[i]/particle->w[pair->i[i]];
        if(particle->type[pair->j[i]] == 0)
        {
            particle->density[pair->j[i]] += m*kernel->w[i]/particle->w[pair->j[i]];
        }
    }
}*/
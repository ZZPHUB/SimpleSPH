#include "PreProcess.H"


void fluid_ptc_generate(SPH_PARTICLE *particle)
{
    for(int i=0;i<(DOMAIN_LENGTH_NUM);i++)
    {
        for(int j=0;j<(DOMAIN_DEEPTH_NUM);j++)
        {
            particle->x[i*(DOMAIN_DEEPTH_NUM)+j] = i*PTC_SPACING;
            particle->y[i*(DOMAIN_DEEPTH_NUM)+j] = j*PTC_SPACING;
            particle->type[i*(DOMAIN_DEEPTH_NUM)+j] = 0;
        }
    }
}

void solid_ptc_generate(SPH_PARTICLE *particle)
{

}

void virtual_ptc_generate(SPH_PARTICLE *particle)
{

}

void ptc_generate(SPH_PARTICLE *particle)
{
    fluid_ptc_generate(particle);
    solid_ptc_generate(particle);
    virtual_ptc_generate(particle);
}
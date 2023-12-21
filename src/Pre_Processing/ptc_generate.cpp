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
    int k =0;
    for(int i=0;i<(DOMAIN_LENGTH_NUM+4);i++)
    {
        for(int j=0;j<(DOMAIN_DEEPTH_NUM+2);j++)
        {
            if(i<=1 || i>= DOMAIN_LENGTH_NUM+2)
            {
                particle->x[FLUID_PTC_NUM+k] = (i-2)*PTC_SPACING;
                particle->y[FLUID_PTC_NUM+k] = (j-2)*PTC_SPACING;
                particle->type[FLUID_PTC_NUM+k] = -1;
                k++;
            }
            else if (j <= 1)
            {
                particle->x[FLUID_PTC_NUM+k] = (i-2)*PTC_SPACING;
                particle->y[FLUID_PTC_NUM+k] = (j-2)*PTC_SPACING;
                particle->type[FLUID_PTC_NUM+k] = -1;
                k++;
            }
        }
    }
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
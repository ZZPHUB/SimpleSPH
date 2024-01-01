#include "PreProcess.H"


void fluid_ptc_generate(SPH_PARTICLE *particle)
{
    for(int i=0;i<(FLUID_LENGTH_NUM);i++)
    {
        for(int j=0;j<(FLUID_DEEPTH_NUM);j++)
        {
            particle->x[i*(FLUID_DEEPTH_NUM)+j] = (i+4)*PTC_SPACING;
            particle->y[i*(FLUID_DEEPTH_NUM)+j] = (j+2)*PTC_SPACING;
            particle->type[i*(FLUID_DEEPTH_NUM)+j] = 0;
        }
    }
}

void virtual_ptc_generate(SPH_PARTICLE *particle)
{
    int k =0;
    for(int i=0;i<(FLUID_LENGTH_NUM+8);i++)
    {
        for(int j=0;j<(FLUID_DEEPTH_NUM+8);j++)
        {
            if(i<2 || i> FLUID_LENGTH_NUM+5)
            {
                particle->x[FLUID_PTC_NUM+k] = i*PTC_SPACING;
                particle->y[FLUID_PTC_NUM+k] = j*PTC_SPACING;
                particle->type[FLUID_PTC_NUM+k] = -1;
                k++;
            }
            else if (j < 2)
            {
                particle->x[FLUID_PTC_NUM+k] = i*PTC_SPACING;
                particle->y[FLUID_PTC_NUM+k] = j*PTC_SPACING;
                particle->type[FLUID_PTC_NUM+k] = -1;
                k++;
            }
        }
    }
}

void ptc_generate(SPH_PARTICLE *particle)
{
    fluid_ptc_generate(particle);
    virtual_ptc_generate(particle);
    solid_ptc_generate(particle);
}
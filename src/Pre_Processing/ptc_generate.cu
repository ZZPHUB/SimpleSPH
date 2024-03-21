#include "PreProcess.cuh"


void ptc_fluid_generate(SPH *sph)
{
    SPH_PARTICLE *particle;
    particle = sph->particle;

    for(int i=0;i<(FLUID_LENGTH_NUM);i++)
    {
        for(int j=0;j<(FLUID_DEEPTH_NUM);j++)
        {
            particle->x[i*(FLUID_DEEPTH_NUM)+j] = (i+4)*PTC_SPACING;
            particle->y[i*(FLUID_DEEPTH_NUM)+j] = (j+4)*PTC_SPACING;
            particle->type[i*(FLUID_DEEPTH_NUM)+j] = 0;
        }
    }
}

void ptc_wall_generate(SPH *sph)
{
    SPH_PARTICLE *particle;
    particle = sph->particle;

    int k =0;

    for(int i=0;i<(FLUID_LENGTH_NUM+8);i++)
    {
        for(int j=0;j<(FLUID_DEEPTH_NUM+16);j++)
        {
            if(i<4 || i> FLUID_LENGTH_NUM+3)
            {
                particle->x[FLUID_PTC_NUM+k] = i*PTC_SPACING;
                particle->y[FLUID_PTC_NUM+k] = j*PTC_SPACING;
                particle->type[FLUID_PTC_NUM+k] = -1;
                k++;
            }
            else if (j < 4)
            {
                particle->x[FLUID_PTC_NUM+k] = i*PTC_SPACING;
                particle->y[FLUID_PTC_NUM+k] = j*PTC_SPACING;
                particle->type[FLUID_PTC_NUM+k] = -1;
                k++;
            }
        }
    }
}

void ptc_generate(SPH *sph)
{
    if(sph->new_case_flag == 1)
    {
        ptc_fluid_generate(sph);
        ptc_wall_generate(sph);
        ptc_rigid_generate(sph);
    }
   else if(sph->new_case_flag == 0)
   {
        ptc_read_vtk(sph);    
   } 
}
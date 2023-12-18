#include "PreProcess.H"



void fluid_ptc_generate(double *f_x,double *f_y)
{
    for(int i=0;i<(DOMAIN_LENGTH_NUM);i++)
    {
        for(int j=0;j<(DOMAIN_DEEPTH_NUM);j++)
        {
            f_x[i*(DOMAIN_DEEPTH_NUM)+j] = i*PTC_SPACING;
            f_y[i*(DOMAIN_DEEPTH_NUM)+j] = j*PTC_SPACING;
        }
    }
}

void solid_ptc_generate(double *s_x,double *s_y)
{

}

void virtual_ptc_generate(double *v_x,double *v_y)
{

}

void ptc_generate(double *p_x,double *p_y)
{
    fluid_ptc_generate(p_x,p_y);
    solid_ptc_generate(p_x,p_y);
    virtual_ptc_generate(p_x,p_y);
}
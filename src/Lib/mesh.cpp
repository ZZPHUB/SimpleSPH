#include "Lib.H"
using namespace std;

void mesh_process(SPH_PARTICLE *particle,unsigned ***mesh)
{
    unsigned int head;
    double x;
    double y;
    for(unsigned int i=0;i<PTC_TOL_NUM;i++)
    {
        if(particle->x[i]<0)
        {
            x = 0;
        }
        else if(particle->x[i]>=DOMAIN_LENGTH)
        {
            x = DOMAIN_LENGTH-PTC_SPACING;
        }
        else
        {
            x = particle->x[i];
        }
        if(particle->y[i]<0)
        {
            y = 0;
        }
        else if(particle->y[i]>=DOAMIN_DEEPTH)
        {
            y = DOMAIN_LENGTH-PTC_SPACING;
        }
        else
        {
            y = particle->y[i];
        }
        head = mesh[(unsigned int)(y/MESH_SPACING)][(unsigned int)(x/MESH_SPACING)][MESH_PTC_NUM-1];
        mesh[(unsigned int)(y/MESH_SPACING)][(unsigned int)(x/MESH_SPACING)][head] = i;        
        mesh[(unsigned int)(y/MESH_SPACING)][(unsigned int)(x/MESH_SPACING)][MESH_PTC_NUM-1]++;
    }

}
#include "Lib.H"
using namespace std;

void mesh_process(SPH_PARTICLE *particle,unsigned ***mesh)
{
    unsigned int head;
    unsigned int j=0;
    unsigned int k=0;
    for(unsigned int i=0;i<PTC_TOL_NUM;i++)
    {
        if(particle->y[i] < DOAMIN_DEEPTH && particle->y[i] >= 0)
        {
            j = (unsigned int)(particle->y[i]/MESH_SPACING);
        }
        else if(particle->y[i] >= DOAMIN_DEEPTH)
        {
            j = MESH_DEEPTH_NUM - 1;
        }
        else
        {
            j = 0;
        }
        if(particle->x[i] < DOMAIN_LENGTH && particle->x[i] >= 0)
        {
            k = (unsigned int)(particle->x[i]/MESH_SPACING);
        }
        else if(particle->x[i] >= DOMAIN_LENGTH)
        {
            k = MESH_LENGTH_NUM - 1;
        }
        else
        {
            k = 0;
        }
        head = mesh[j][k][MESH_PTC_NUM-1];
        mesh[j][k][head] = i;        
        mesh[j][k][MESH_PTC_NUM-1]++;
    }

}
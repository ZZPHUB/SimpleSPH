#include "Lib.H"
using namespace std;

void ptc_kernel(SPH_PARTICLE *particle,SPH_PAIR *pair,SPH_KERNEL *kernel)
{
    #pragma omp parallel sections num_threads(3)
    {
        #pragma omp section 
        {
            //temp particle distance
            double temp_dis = 0;
            //calculate the kernel value
            for(int i=0;i<pair->total;i++)
            {   temp_dis = PTC_DISTANCE(pair->i[i],pair->j[i]);
                if(0<temp_dis && temp_dis < PTC_SML)
                {
                    kernel->w[i] = ALPHA*(2.0/3.0-pow(temp_dis,2)+0.5*pow(temp_dis,3));
                }
                else if (PTC_SML <= temp_dis && temp_dis < PTC_REGION_RADIUS)
                {
                    kernel->w[i] = ALPHA*(pow(2.0-temp_dis,3)/6.0); 
                }
            }
        }
        #pragma omp section
        {
            //temp paritcle distance 
            double temp_dis;
            //calculate the x-direction differential kernel value
            for(int i=0;i<pair->total;i++)
            {
                temp_dis = PTC_DISTANCE(pair->i[i],pair->j[i]);
                if(0<temp_dis && temp_dis < PTC_SML)
                {
                    kernel->dwdx[i] = ALPHA*(-2.0+1.5*temp_dis/PTC_SML)*(particle->x[pair->i[i]]-particle->x[pair->j[i]])/pow(PTC_SML,2);
                }
                else if (PTC_SML <= temp_dis && temp_dis < PTC_REGION_RADIUS)
                {
                    kernel->dwdx[i] = ALPHA*0.5*pow(2-temp_dis/PTC_SML,2)*(particle->x[pair->i[i]]-particle->x[pair->j[i]])/(PTC_SML*temp_dis);
                }
            }
        }
        #pragma omp section
        {
            //temp paritcle distance 
            double temp_dis;
            //calculate the x-direction differential kernel value
            for(int i=0;i<pair->total;i++)
            {
                temp_dis = PTC_DISTANCE(pair->i[i],pair->j[i]);
                if(0<temp_dis && temp_dis < PTC_SML)
                {
                    kernel->dwdx[i] = ALPHA*(-2.0+1.5*temp_dis/PTC_SML)*(particle->y[pair->i[i]]-particle->y[pair->j[i]])/pow(PTC_SML,2);
                }
                else if (PTC_SML <= temp_dis && temp_dis < PTC_REGION_RADIUS)
                {
                    kernel->dwdx[i] = ALPHA*0.5*pow(2-temp_dis/PTC_SML,2)*(particle->y[pair->i[i]]-particle->y[pair->j[i]])/(PTC_SML*temp_dis);
                }
            }
        }
    }
}
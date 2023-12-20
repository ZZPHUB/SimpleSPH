#include "Lib.H"
using namespace std;

void ptc_kernel_serial(SPH_PARTICLE *particle,SPH_PAIR *pair,SPH_KERNEL *kernel)
{
    double temp_dis;
    for(int i=0;i<pair->total;i++)
    {   
        temp_dis = PTC_DISTANCE(pair->i[i],pair->j[i]);
        //temp_dis = sqrt(pow(particle->x[pair->i[i]]-particle->x[pair->j[i]],2)+pow(particle->y[pair->i[i]]-particle->y[pair->j[i]],2));
        if(0<temp_dis && temp_dis < PTC_SML)
        {
            kernel->w[i] = ALPHA*(2.0/3.0-pow(temp_dis,2)+0.5*pow(temp_dis,3));
            kernel->dwdx[i] = ALPHA*(-2.0+1.5*temp_dis/PTC_SML)*(particle->x[pair->i[i]]-particle->x[pair->j[i]])/pow(PTC_SML,2);
            kernel->dwdy[i] = ALPHA*(-2.0+1.5*temp_dis/PTC_SML)*(particle->y[pair->i[i]]-particle->y[pair->j[i]])/pow(PTC_SML,2);
        }
        else if (PTC_SML <= temp_dis && temp_dis < PTC_REGION_RADIUS)
        {
            kernel->w[i] = ALPHA*(pow(2.0-temp_dis,3)/6.0); 
            kernel->dwdx[i] = ALPHA*0.5*pow(2-temp_dis/PTC_SML,2)*(particle->x[pair->i[i]]-particle->x[pair->j[i]])/(PTC_SML*temp_dis);
            kernel->dwdy[i] = ALPHA*0.5*pow(2-temp_dis/PTC_SML,2)*(particle->y[pair->i[i]]-particle->y[pair->j[i]])/(PTC_SML*temp_dis);
        }
    }
}


void ptc_kernel_parallel(SPH_PARTICLE *particle,SPH_PAIR *pair,SPH_KERNEL *kernel)
{
    double temp_dis;
    #pragma omp parallel sections private(temp_dis) num_threads(3)
    {
        #pragma omp section 
        {
            //temp particle distance
            //calculate the kernel value
            for(int i=0;i<pair->total;i++)
            {   
                temp_dis = PTC_DISTANCE(pair->i[i],pair->j[i]);
                //temp_dis = sqrt(pow(particle->x[pair->i[i]]-particle->x[pair->j[i]],2)+pow(particle->y[pair->i[i]]-particle->y[pair->j[i]],2));
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
            //calculate the x-direction differential kernel value
            for(int i=0;i<pair->total;i++)
            {
                temp_dis = PTC_DISTANCE(pair->i[i],pair->j[i]);
                //temp_dis = sqrt(pow(particle->x[pair->i[i]]-particle->x[pair->j[i]],2)+pow(particle->y[pair->i[i]]-particle->y[pair->j[i]],2));
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
            //calculate the x-direction differential kernel value
            for(int i=0;i<pair->total;i++)
            {
                temp_dis = PTC_DISTANCE(pair->i[i],pair->j[i]);
                //temp_dis = sqrt(pow(particle->x[pair->i[i]]-particle->x[pair->j[i]],2)+pow(particle->y[pair->i[i]]-particle->y[pair->j[i]],2));
                if(0<temp_dis && temp_dis < PTC_SML)
                {
                    kernel->dwdy[i] = ALPHA*(-2.0+1.5*temp_dis/PTC_SML)*(particle->y[pair->i[i]]-particle->y[pair->j[i]])/pow(PTC_SML,2);
                }
                else if (PTC_SML <= temp_dis && temp_dis < PTC_REGION_RADIUS)
                {
                    kernel->dwdy[i] = ALPHA*0.5*pow(2-temp_dis/PTC_SML,2)*(particle->y[pair->i[i]]-particle->y[pair->j[i]])/(PTC_SML*temp_dis);
                }
            }
        }
    }
}

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
    omp_lock_t lock;
    omp_init_lock(&lock);
    
    #pragma omp parallel for num_threads(6)
    for(int i=0;i<pair->total;i++)
    {   
        temp_dis = PTC_DISTANCE(pair->i[i],pair->j[i]);
        if(0<temp_dis && temp_dis < PTC_SML)
        {
            omp_set_lock(&lock);
            kernel->w[i] = ALPHA*(2.0/3.0-pow(temp_dis,2)+0.5*pow(temp_dis,3));
            kernel->dwdx[i] = ALPHA*(-2.0+1.5*temp_dis/PTC_SML)*(particle->x[pair->i[i]]-particle->x[pair->j[i]])/pow(PTC_SML,2);
            kernel->dwdy[i] = ALPHA*(-2.0+1.5*temp_dis/PTC_SML)*(particle->y[pair->i[i]]-particle->y[pair->j[i]])/pow(PTC_SML,2);
            omp_unset_lock(&lock);
        }
        else if (PTC_SML <= temp_dis && temp_dis < PTC_REGION_RADIUS)
        {
            omp_set_lock(&lock);
            kernel->w[i] = ALPHA*(pow(2.0-temp_dis,3)/6.0); 
            kernel->dwdx[i] = ALPHA*0.5*pow(2-temp_dis/PTC_SML,2)*(particle->x[pair->i[i]]-particle->x[pair->j[i]])/(PTC_SML*temp_dis);
            kernel->dwdy[i] = ALPHA*0.5*pow(2-temp_dis/PTC_SML,2)*(particle->y[pair->i[i]]-particle->y[pair->j[i]])/(PTC_SML*temp_dis);
            omp_unset_lock(&lock);
        }
    }
}



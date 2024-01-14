#include "Lib.H"
using namespace std;

void ptc_kernel_parallel(SPH *sph)
//void ptc_kernel_parallel(SPH_PARTICLE *particle,SPH_PAIR *pair,SPH_KERNEL *kernel)
{
    SPH_PARTICLE *particle;
    SPH_PAIR *pair;
    SPH_KERNEL *kernel;
    particle = sph->particle;
    pair = sph->pair;
    kernel = sph->kernel;

    double temp_dis;
    omp_lock_t lock;
    omp_init_lock(&lock);
    
    //particle->w donot involve time integration
    #pragma omp parallel for num_threads(TH_NUM)
    for(unsigned int i=0;i<particle->total;i++)
    {
        if(particle->type[i]==0)
        {
            omp_set_lock(&lock);
            particle->w[i] = 2.0*ALPHA/3.0;
            omp_unset_lock(&lock);
        }
        else if(particle->type[i]==-1)
        {
            omp_set_lock(&lock);
            particle->w[i] = 0;
            omp_unset_lock(&lock);
        }
    }

    #pragma omp parallel for num_threads(TH_NUM)
    for(unsigned int i=0;i<pair->total;i++)
    {   
        temp_dis = PTC_DISTANCE(pair->i[i],pair->j[i]);
        if(0<temp_dis && temp_dis < PTC_SML)
        {
            omp_set_lock(&lock);
            //each pair's kernel value
            kernel->w[i] = ALPHA*(2.0/3.0-pow(temp_dis/PTC_SML,2)+0.5*pow(temp_dis/PTC_SML,3));

            //each pair's differential kernel value in x and y direction
            kernel->dwdx[i] = ALPHA*(-2.0+1.5*temp_dis/PTC_SML)*(particle->x[pair->i[i]]-particle->x[pair->j[i]])/pow(PTC_SML,2);
            kernel->dwdy[i] = ALPHA*(-2.0+1.5*temp_dis/PTC_SML)*(particle->y[pair->i[i]]-particle->y[pair->j[i]])/pow(PTC_SML,2);
            
            //each particles kernel value sum
            particle->w[pair->i[i]] += kernel->w[i];
            particle->w[pair->j[i]] += kernel->w[i];
            omp_unset_lock(&lock);
        }
        else if (PTC_SML <= temp_dis && temp_dis < PTC_REGION_RADIUS)
        {
            omp_set_lock(&lock);
            //each pair's kernel value
            kernel->w[i] = ALPHA*(pow(2.0-temp_dis/PTC_SML,3)/6.0); 

            //each pair's differential kernel value in x and y direction
            kernel->dwdx[i] = -1*ALPHA*0.5*pow(2-temp_dis/PTC_SML,2)*(particle->x[pair->i[i]]-particle->x[pair->j[i]])/(PTC_SML*temp_dis);
            kernel->dwdy[i] = -1*ALPHA*0.5*pow(2-temp_dis/PTC_SML,2)*(particle->y[pair->i[i]]-particle->y[pair->j[i]])/(PTC_SML*temp_dis);
            
            //each particles kernel value sum
            particle->w[pair->i[i]] += kernel->w[i];
            particle->w[pair->j[i]] += kernel->w[i];
            omp_unset_lock(&lock);
        }
    }
}



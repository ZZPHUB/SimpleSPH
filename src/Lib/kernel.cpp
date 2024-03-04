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

    double m = PTC_MASS;
    double r = 0;
    double q = 0;
    double a = ALPHA;
    double dx = 0;
    double dy = 0;

    omp_lock_t lock;
    omp_init_lock(&lock);

    #pragma omp parallel for num_threads(TH_NUM)
    for(unsigned int i=0;i<pair->total;i++)
    {   
        omp_set_lock(&lock);
        dx = particle->x[pair->i[i]] - particle->x[pair->j[i]];
        dy = particle->y[pair->i[i]] - particle->y[pair->j[i]];
        r = sqrt(dx*dx+dy*dy);
        q = r/PTC_SML;
        omp_unset_lock(&lock);

        if(0.0 <= q && q < 1.0)
        {
            //each pair's kernel value
            kernel->w[i] = a*(2.0/3.0-q*q+0.5*q*q*q);

            //each pair's differential kernel value in x and y direction
            kernel->dwdx[i] = a*((-2.0+1.5*q)*dx/pow(PTC_SML,2));
            kernel->dwdy[i] = a*((-2.0+1.5*q)*dy/pow(PTC_SML,2));
        }
        else if (1.0 <= q && q < 2.0)
        {
            //each pair's kernel value
            kernel->w[i] = a*((2.0-q)*(2.0-q)*(2.0-q)/6.0); 

            //each pair's differential kernel value in x and y direction
            kernel->dwdx[i] = -a*0.5*((2.0-q)*(2.0-q)*dx/(PTC_SML*r));
            kernel->dwdy[i] = -a*0.5*((2.0-q)*(2.0-q)*dy/(PTC_SML*r));
        }
        else
        {
            kernel->w[i] = kernel->dwdx[i] = kernel->dwdy[i] = 0.0;
        }
    }
}
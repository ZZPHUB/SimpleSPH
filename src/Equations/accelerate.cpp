#include "Equations.H"

void ptc_acc(SPH_PARTICLE *particle,SPH_PAIR *pair,SPH_KERNEL *kernel)
{
    omp_lock_t lock;
    omp_init_lock(&lock);
    double m = PTC_MASS;
    double temp_p = 0;
    double temp_rho_i = 0;
    double temp_rho_j = 0;

    #pragma omp parallel for num_threads(TH_NUM)
    for(unsigned int i=0;i<particle->total;i++)
    {
        omp_set_lock(&lock);
        particle->accx[i] = 0;
        particle->accy[i] = -GRAVITY_ACC;
        omp_unset_lock(&lock);
    }

    #pragma omp parallel for num_threads(TH_NUM)
    for(unsigned int i=0;i<pair->total;i++)
    {
        omp_set_lock(&lock);

        //to reduce the calculate times,define $\rho^2$ to temp para
        temp_rho_i = pow(particle->density[pair->i[i]],2);
        temp_rho_j = pow(particle->density[pair->j[i]],2);

        //to reduce the calculate tiems,define $p_i/\rho_i^2 + p_j/\rho_j^2$ to temp para
        temp_p = particle->pressure[pair->i[i]]/temp_rho_i + particle->pressure[pair->j[i]]/temp_rho_j;

        particle->accx[pair->i[i]] = particle->accx[pair->i[i]] - temp_p*kernel->dwdx[i] + \
        MU*(particle->visxx[pair->i[i]]/temp_rho_i + particle->visxx[pair->j[i]]/temp_rho_j)*kernel->dwdx[i] + \
        MU*(particle->visxy[pair->i[i]]/temp_rho_i + particle->visxy[pair->j[i]]/temp_rho_j)*kernel->dwdy[i];

        particle->accx[pair->j[i]] = particle->accx[pair->j[i]] + temp_p*kernel->dwdx[i] - \
        MU*(particle->visxx[pair->i[i]]/temp_rho_i + particle->visxx[pair->j[i]]/temp_rho_j)*kernel->dwdx[i] - \
        MU*(particle->visxy[pair->i[i]]/temp_rho_i + particle->visxy[pair->j[i]]/temp_rho_j)*kernel->dwdy[i];

        particle->accy[pair->i[i]] = particle->accy[pair->i[i]] - temp_p*kernel->dwdy[i] + \
        MU*(particle->visyy[pair->i[i]]/temp_rho_i + particle->visyy[pair->j[i]]/temp_rho_j)*kernel->dwdy[i] + \
        MU*(particle->visxy[pair->i[i]]/temp_rho_i + particle->visxy[pair->j[i]]/temp_rho_j)*kernel->dwdx[i];

        particle->accy[pair->j[i]] = particle->accy[pair->j[i]] + temp_p*kernel->dwdy[i] - \
        MU*(particle->visyy[pair->i[i]]/temp_rho_i + particle->visyy[pair->j[i]]/temp_rho_j)*kernel->dwdy[i] - \
        MU*(particle->visxy[pair->i[i]]/temp_rho_i + particle->visxy[pair->j[i]]/temp_rho_j)*kernel->dwdx[i];

        omp_unset_lock(&lock);
    }

    /*
    #pragma omp parallel for num_threads(TH_NUM)
    //add gravity acceleration in y-direction
    for(unsigned int i=0;i<particle->total;i++)
    {
        if(particle->type[i]==0)
        {
            omp_set_lock(&lock);
            particle->accy[i] -= GRAVITY_ACC;
            omp_unset_lock(&lock);
        }
    }
    */
    
}
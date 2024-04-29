#include "Lib.cuh"

__global__ void sph_sponge_cuda(SPH_CUDA *cuda,SPH_ARG *arg,SPH_RIGID *rigid)
{
    //do not impliment
    const int id = threadIdx.x + blockIdx.x * blockDim.x;
    double r = 0.0;
    double tmp_0 = 0.0;
    double tmp_1 = 0.0;
    double tmp_2 = 0.0;
    if(id < arg->ptc_num)
    {
        if(cuda->type[id] == 0)
        {
            r = sqrt(pow((cuda->x[id]-arg->fluid_x/2.0),2)+pow((cuda->y[id]-arg->fluid_y),2)) - arg->fluid_x/2.0;
            if(r > 0.0)
            {
                tmp_0 = 50.0*r/arg->sponge_dx;
                tmp_1 = pow(0.3,tmp_0);
                tmp_2 = 1.0 - pow(100,-tmp_1);

                cuda->accx[id] *= tmp_2;
                cuda->accy[id] *= tmp_2;
                cuda->drho[id] *= tmp_2;
            }
        }
    }
}
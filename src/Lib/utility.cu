#include "Lib.cuh"

/*__global__ void sph_dummy_cuda(double *x,double *y,double *vx,double *vy,double *p,double *rho,double *ptc_w,double *pair_w,int *pair_i,int *pair_j,int *type,double *rigid,int *pair_num)
{
    double rigid_accx = 0.0;
    double rigid_accy = 0.0;
    double dx = 0.0;
    double dy = 0.0;
    double temp_p = 0.0;
    double temp_vx = 0.0;
    double temp_vy = 0.0;

    const int id = threadIdx.x + blockIdx.x * blockDim.x;
    if(id >= pair_num[0]) return;

    if(type[pair_j[id]] != 0 && ptc_w[pair_j[id]] != 0.0)
    {
        rigid_accx = 0.0;
        rigid_accy = 0.0;
        dx = x[pair_i[id]]-x[pair_j[id]];
        dy = y[pair_i[id]]-y[pair_j[id]];
        if(type[pair_j[id]] == -1)
        {
            rigid_accx = 0.0;
            rigid_accy = 0.0;
        }
        else if (type[pair_j[id]] == 1)
        {
            rigid_accx = rigid[ACCX] - pow(rigid[OMEGA],2)*(x[pair_j[id]]-rigid[COGX])- \
                              rigid[R_ALPHA]*(y[pair_j[id]]-rigid[COGY]);
            rigid_accy = rigid[ACCY] - pow(rigid[OMEGA],2)*(x[pair_j[id]]-rigid[COGY])+ \
                              rigid[R_ALPHA]*(x[pair_j[id]]-rigid[COGX]);
        }
        temp_p = (p[pair_i[id]]+rho[pair_i[id]]*(rigid_accx*dx+(rigid_accy+GRAVITY_ACC)*dy))*pair_w[id]/ptc_w[pair_j[id]];
        temp_vx = vx[pair_i[id]]*pair_w[id]/ptc_w[pair_j[id]];
        temp_vy = vy[pair_i[id]]*pair_w[id]/ptc_w[pair_j[id]];

        atomicAdd(&p[pair_j[id]],temp_p);
        atomicAdd(&vx[pair_j[id]],temp_vx);
        atomicAdd(&vy[pair_j[id]],temp_vy);
    }
}*/


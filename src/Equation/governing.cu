#include "Equations.cuh"

__global__ void sph_governing_cuda(double *x,double *y,double *vx,double *vy,\
double *rho,double *p,int *type,int *pair_i,int *pair_j,double *dwdx,\
double *dwdy,double *accx,double *accy,double *drho,double *rigid,int* pair_num,int ptc_num)
{
    const int id = threadIdx.x + blockIdx.x * blockDim.x;
    if(id >= pair_num[0] )return;

    double dx;
    double dy;
    double dvx;
    double dvy;
    double rho_temp=0.0;
    double accx_temp=0.0;
    double accy_temp=0.0;
    double temp=0.0;

    dx = x[pair_i[id]]-x[pair_j[id]];
    dy = x[pair_i[id]]-y[pair_j[id]];

    accx_temp = (-dev_m*p[pair_i[id]]*dwdx[id])/pow(rho[pair_i[id]],2)+(-dev_m*p[pair_j[id]]*dwdx[id])/pow(rho[pair_j[id]],2);
    accy_temp = (-dev_m*p[pair_i[id]]*dwdy[id])/pow(rho[pair_i[id]],2)+(-dev_m*p[pair_j[id]]*dwdy[id])/pow(rho[pair_j[id]],2);
    //accx_temp = -dev_m*(p[pair_i[id]]/(rho[pair_i[id]]*rho[pair_i[id]])+p[pair_j[id]]/(rho[pair_j[id]]*rho[pair_j[id]]));

    //accx[id] = acc_temp*dwdx[id];
    //accy[id] = acc_temp*dwdy[id];

    if(type[pair_j[id]] == 0)
    {
        dvx = vx[pair_i[id]]-vx[pair_j[id]];
        dvy = vy[pair_i[id]]-vy[pair_j[id]];
        rho_temp = dvx*dwdx[id]+dvy*dwdy[id];
        rho_temp *= dev_m;
    }
    else if(type[pair_j[id]] == 1)
    {
        dvx = vx[pair_i[id]] - (2.0*(rigid[VX] - rigid[OMEGA]*(y[pair_j[id]]-rigid[COGY])) - vx[pair_j[id]]);
        dvy = vy[pair_i[id]] - (2.0*(rigid[VY] + rigid[OMEGA]*(x[pair_j[id]]-rigid[COGX])) - vy[pair_j[id]]);
        rho_temp = (vx[pair_i[id]]-(rigid[VX] - rigid[OMEGA]*(y[pair_j[id]]-rigid[COGY])))*dwdx[id]+\
                   (vy[pair_i[id]]-(rigid[VY] + rigid[OMEGA]*(x[pair_j[id]]-rigid[COGX])))*dwdy[id];
        rho_temp *= dev_m;
    }
    else if(type[pair_j[id]] == -1)
    {
        dvx = vx[pair_i[id]] - (0.0 - vx[pair_j[id]]); 
        dvy = vy[pair_i[id]] - (0.0 - vy[pair_j[id]]);
        rho_temp = vx[pair_i[id]]*dwdx[id]+vy[pair_j[id]]*dwdy[id];
        rho_temp *= dev_m;
    }
    /*
    accy_temp = dx*dvx+dy*dvy;
    if(accy_temp < 0.0) accy_temp = 0.0;
    
    accx_temp += accy_temp*dev_m*0.01*dev_h*dev_c/((dx*dx+dy*dy+0.01*dev_h*dev_h)*0.5*(rho[pair_i[id]]+rho[pair_j[id]]));
    accy_temp = accx_temp*dwdy[id];
    accx_temp *= dwdx[id];
    */

    atomicAdd(&accx[pair_i[id]], accx_temp);
    atomicAdd(&accx[pair_j[id]],-accx_temp);
    atomicAdd(&accy[pair_i[id]], accy_temp);
    atomicAdd(&accy[pair_j[id]],-accy_temp);
    atomicAdd(&drho[pair_i[id]],rho_temp);
    atomicAdd(&drho[pair_j[id]],rho_temp);
}
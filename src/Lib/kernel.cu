#include "Lib.cuh"

/*
__global__ void sph_kernel_cuda(double *x,double *y,double *w,double *dwdx,double *dwdy,double *ptc_w,int *pair_i,int *pair_j,int* pair_num)
{
    double dx,dy,q;
    const int id = threadIdx.x + blockIdx.x* blockDim.x;
    if(id >= pair_num[0]) return;

    dx = x[pair_i[id]]-x[pair_j[id]];
    dy = y[pair_i[id]]-y[pair_j[id]];
    q = sqrt(dx*dx+dy*dy)/dev_h;

    if(q<1.0)
    {
        w[id] = dev_a*(2.0/3.0-q*q+0.5*q*q*q);
        dwdx[id] = dev_a*((-2.0+1.5*q)*dx)/pow(dev_h,2);
        dwdy[id] = dev_a*((-2.0+1.5*q)*dy)/pow(dev_h,2);
        atomicAdd(&ptc_w[pair_i[id]],w[id]);
        atomicAdd(&ptc_w[pair_j[id]],w[id]);
        /*
            kernel->w[i] = a*(2.0/3.0-q*q+0.5*q*q*q);
            //each pair's differential kernel value in x and y direction
            kernel->dwdx[i] = a*((-2.0+1.5*q)*dx/pow(dev_h,2));
            kernel->dwdy[i] = a*((-2.0+1.5*q)*dy/pow(dev_h,2));
        */
    }
    else if(q<2.0)
    {
        w[id] = dev_a*((2.0-q)*(2.0-q)*(2.0-q))/6.0;
        dwdx[id] = -dev_a*0.5*((2.0-q)*(2.0-q)*dx)/(dev_h*dev_h*q);
        dwdy[id] = -dev_a*0.5*((2.0-q)*(2.0-q)*dy)/(dev_h*dev_h*q);
        atomicAdd(&ptc_w[pair_i[id]],w[id]);
        atomicAdd(&ptc_w[pair_j[id]],w[id]);
        /*
            //each pair's kernel value
            kernel->w[i] = a*((2.0-q)*(2.0-q)*(2.0-q)/6.0); 
            //each pair's differential kernel value in x and y direction
            kernel->dwdx[i] = -a*0.5*((2.0-q)*(2.0-q)*dx/(dev_h*r));
            kernel->dwdy[i] = -a*0.5*((2.0-q)*(2.0-q)*dy/(dev_h*r));
        */
    }
    else
    {
        w[id] = dwdx[id] = dwdy[id] = 0.0;
    }
}*/
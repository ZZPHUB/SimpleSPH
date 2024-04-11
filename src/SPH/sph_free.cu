#include "SPH.cuh"

void sph_free(SPH *sph)
{
    SPH_PARTICLE *particle;
    SPH_PAIR *pair;
    SPH_KERNEL *kernel;
    SPH_RIGID *wedge;
    particle = sph->particle;
    pair = sph->pair;
    kernel = sph->kernel;
    wedge = sph->host_rigid;

    cudaFree(sph->cuda);
    free(particle->x);
    free(particle->y);
    free(particle->vx);
    free(particle->vy);
    free(particle->accx);
    free(particle->accy);
    free(particle->density);
    free(particle->dif_density);
    free(particle->pressure);
    free(particle->type);
    free(particle->mass);

    free(kernel->w);
    free(kernel->dwdx);
    free(kernel->dwdy);
    
    free(pair->i);
    free(pair->j);
    //free(mesh);
}
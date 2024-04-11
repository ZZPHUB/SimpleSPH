#include "SPH.cuh"

void sph_free(SPH *sph)
{
    SPH_PARTICLE *particle;
    SPH_PAIR *pair;
    SPH_KERNEL *kernel;
    SPH_RIGID *wedge;
    SPH_CUDA cuda;
    particle = sph->particle;
    pair = sph->pair;
    kernel = sph->kernel;
    wedge = sph->host_rigid;
    cudaMemcpy(&cuda,sph->cuda,sizeof(SPH_CUDA),cudaMemcpyDeviceToHost);

    cudaFree(cuda->x);
    cudaFree(cuda->y);
    cudaFree(cuda->vx);
    cudaFree(cuda->vy);
    cudaFree(cuda->temp_x);
    cudaFree(cuda->temp_y);
    cudaFree(cuda->temp_vx);
    cudaFree(cuda->temp_vy);
    cudaFree(cuda->rho);
    cudaFree(cuda->drho);
    cudaFree(cuda->temp_rho);
    cudaFree(cuda->accx);
    cudaFree(cuda->accy);
    cudaFree(cuda->p);
    cudaFree(cuda->type);
    cudaFree(cuda->ptc_w);
    cudaFree(cuda->pair_i);
    cudaFree(cuda->pair_j);
    cudaFree(cuda->pair_w);
    cudaFree(cuda->dwdx);
    cudaFree(cuda->dwdy);
    cudaFree(cuda->mesh);
    cudaFree(cuda->mesh_count);
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
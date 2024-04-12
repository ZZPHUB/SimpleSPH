#include "SPH.cuh"

void sph_nnps_cpu(SPH *sph)
{
    SPH_MESH *mesh;
    SPH_PARTICLE *particle;
    SPH_ARG *arg;
    //SPH_PAIR *pair;
    mesh = sph->mesh;
    particle = sph->particle;
    arg = sph->host_arg;
    //pair = sph->pair;

    double q=0.0;
    double dx=0.0;
    double dy=0.0;
    int mesh_id=0;
    int index_i=0;
    int index_j=0;
    arg->pair_num = 0;
    for(int i=0;i<arg->mesh_xnum;i++)
    {
        for(int j=0;j<arg->mesh_ynum;j++)
        {
            mesh_id = i+j*arg->mesh_xnum;
            for(int m=0;m<mesh->count[mesh_id];m++)
            {
                //[x,y]->[x,y]
                index_i = mesh_id + m*arg->mesh_num;
                for(int n=m+1;n<mesh->count[mesh_id];n++)
                {
                    index_j = mesh_id + n*arg->mesh_num;
                    dx = particle->x[mesh->ptc[index_i]] -particle->x[mesh->ptc[index_j]];
                    dy = particle->y[mesh->ptc[index_i]] -particle->y[mesh->ptc[index_j]];
                    q = sqrt(dx*dx+dy*dy)/arg->h;
                    if(q<2.0)
                    {
                        if(particle->type[mesh->ptc[index_i]]==0 || particle->type[mesh->ptc[index_j]]==0)
                        {
                            arg->pair_num ++;
                        }
                    }
                }
                //[x,y]->[x+1,y]
                if(i<arg->mesh_xnum-1)
                {
                    for(int n=0;n<mesh->count[mesh_id+1];n++)
                    {
                        index_j = mesh_id+1+n*arg->mesh_num;
                        dx = particle->x[mesh->ptc[index_i]] -particle->x[mesh->ptc[index_j]];
                        dy = particle->y[mesh->ptc[index_i]] -particle->y[mesh->ptc[index_j]];
                        q = sqrt(dx*dx+dy*dy)/arg->h;
                        if(q<2.0)
                        {
                            if(particle->type[mesh->ptc[index_i]]==0 || particle->type[mesh->ptc[index_j]]==0)
                            {
                                arg->pair_num ++;
                            }
                        } 
                    }
                }
                //[x,y]->[x,y+1]
                if(j<arg->mesh_ynum-1)
                {
                    for(int n=0;n<mesh->count[mesh_id+arg->mesh_xnum];n++)
                    {
                        index_j = mesh_id+arg->mesh_xnum+n*arg->mesh_num;
                        dx = particle->x[mesh->ptc[index_i]] -particle->x[mesh->ptc[index_j]];
                        dy = particle->y[mesh->ptc[index_i]] -particle->y[mesh->ptc[index_j]];
                        q = sqrt(dx*dx+dy*dy)/arg->h;
                        if(q<2.0)
                        {
                            if(particle->type[mesh->ptc[index_i]]==0 || particle->type[mesh->ptc[index_j]]==0)
                            {
                                arg->pair_num ++;
                            }
                        } 
                    }
                }
                //[x,y]->[x+1,y+1]
                if(i<(arg->mesh_xnum-1) && j<(arg->mesh_ynum-1))
                {
                    for(int n=0;n<mesh->count[mesh_id+1+arg->mesh_xnum];n++)
                    {
                        index_j = mesh_id+1+arg->mesh_xnum+n*arg->mesh_num;
                        dx = particle->x[mesh->ptc[index_i]] -particle->x[mesh->ptc[index_j]];
                        dy = particle->y[mesh->ptc[index_i]] -particle->y[mesh->ptc[index_j]];
                        q = sqrt(dx*dx+dy*dy)/arg->h;
                        if(q<2.0)
                        {
                            if(particle->type[mesh->ptc[index_i]]==0 || particle->type[mesh->ptc[index_j]]==0)
                            {
                                arg->pair_num ++;
                            }
                        } 
                    }
                }
                //[x,y]->[x-1,y+1]
                if(i>0 && j<(arg->mesh_ynum-1))
                {
                    for(int n=0;n<mesh->count[mesh_id-1+arg->mesh_xnum];n++)
                    {
                        index_j = mesh_id-1+arg->mesh_xnum+n*arg->mesh_num;
                        dx = particle->x[mesh->ptc[index_i]] -particle->x[mesh->ptc[index_j]];
                        dy = particle->y[mesh->ptc[index_i]] -particle->y[mesh->ptc[index_j]];
                        q = sqrt(dx*dx+dy*dy)/arg->h;
                        if(q<2.0)
                        {
                            if(particle->type[mesh->ptc[index_i]]==0 || particle->type[mesh->ptc[index_j]]==0)
                            {
                                arg->pair_num ++;
                            }
                        }
                    }
                }
            }
        }
    }
    printf("the pair num in cpu is :%d\n",arg->pair_num);
}


__global__ void sph_nnps_cuda(SPH_CUDA *cuda,SPH_ARG *arg,SPH_RIGID *rigid)
{
    //blockIdx.x -> mesh x direction
    //blockIdx.y -> mesh y direction
    //threadIdx.x -> local mesh index
    //threadIdx.y -> near mesh index
    if( gridDim.x != arg->mesh_xnum || gridDim.y != arg->mesh_ynum) return;
    const int mesh_id = blockIdx.x + blockIdx.y * gridDim.x;
    int index_i = 0;
    int index_j = 0;
    double dx = 0.0;
    double dy = 0.0;
    double q = 0.0;
    //int tmp_count;
    __shared__ int count;
    if(threadIdx.x == 0 && threadIdx.y == 0) count=0;
    __syncthreads();
    
    index_i = cuda->mesh[mesh_id + threadIdx.x*arg->mesh_num];
    //(x,y)->(x,y)
    if( threadIdx.y> threadIdx.x && threadIdx.y<cuda->mesh_count[mesh_id])
    {
        index_j = cuda->mesh[mesh_id + threadIdx.y*arg->mesh_num];
        dx = cuda->x[index_i] - cuda->x[index_j];
        dy = cuda->y[index_i] - cuda->y[index_j];
        q = sqrt(dx*dx+dy*dy)/arg->h;
        if(q<2.0)
        {
            if(cuda->type[index_i] == 0 || cuda->type[index_j] == 0)
            {
                count+=1;
            }
        }
    }
    //(x,y)->(x+1,y)
    if( blockIdx.x < ( gridDim.x-1))
    {
        if( threadIdx.y < cuda->mesh_count[mesh_id+1])
        {
            index_j = cuda->mesh[mesh_id + 1 + threadIdx.y*arg->mesh_num];
            dx = cuda->x[index_i] - cuda->x[index_j];
            dy = cuda->y[index_i] - cuda->y[index_j];
            q = sqrt(dx*dx+dy*dy)/arg->h;
            if(q<2.0)
            {
                if(cuda->type[index_i] == 0 || cuda->type[index_j] == 0)
                {
                    count+=1;
                }
            }
        }
    }
    //(x,y)->(x,y+1)
    if( blockIdx.y < ( gridDim.y -1))
    {
        if( threadIdx.y < cuda->mesh_count[mesh_id+ gridDim.x])
        {
            index_j = cuda->mesh[mesh_id + gridDim.x + threadIdx.y*arg->mesh_num];
            dx = cuda->x[index_i] - cuda->x[index_j];
            dy = cuda->y[index_i] - cuda->y[index_j];
            q = sqrt(dx*dx+dy*dy)/arg->h;
            if(q<2.0)
            {
                if(cuda->type[index_i] == 0 || cuda->type[index_j] == 0)
                {
                    count+=1;
                }
            }
        }
    }
    //(x,y)->(x+1,y+1)
    if( blockIdx.x<( gridDim.x-1) && blockIdx.y<( gridDim.y-1))
    {
        if(threadIdx.y < cuda->mesh_count[mesh_id+ 1+ gridDim.x])
        {
           index_j = cuda->mesh[mesh_id + 1+  gridDim.x + threadIdx.y*arg->mesh_num];
            dx = cuda->x[index_i] - cuda->x[index_j];
            dy = cuda->y[index_i] - cuda->y[index_j];
            q = sqrt(dx*dx+dy*dy)/arg->h;
            if(q<2.0)
            {
                if(cuda->type[index_i] == 0 || cuda->type[index_j] == 0)
                {
                    count+=1;
                }
            } 
        }
    }
    //(x,y)->(x-1,y+1)
    if( blockIdx.x>0 && blockIdx.y<( gridDim.y-1))
    {
        if(threadIdx.y < cuda->mesh_count[mesh_id- 1+ gridDim.x])
        {
           index_j = cuda->mesh[mesh_id - 1+  gridDim.x + threadIdx.y*arg->mesh_num];
            dx = cuda->x[index_i] - cuda->x[index_j];
            dy = cuda->y[index_i] - cuda->y[index_j];
            q = sqrt(dx*dx+dy*dy)/arg->h;
            if(q<2.0)
            {
                if(cuda->type[index_i] == 0 || cuda->type[index_j] == 0)
                {
                    count+=1;
                }
            } 
        }
    }
    __syncthreads();
    if( threadIdx.x == 0 && threadIdx.y == 0)
    {
        atomicAdd(&(arg->pair_num),count);
        cuda->mesh_count[mesh_id]=0;
    }
    __syncthreads();
}
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

    //threadIdx.x---> (x)length direction
    //blockIdx.x ---> (y)depth direction

    double dx = 0.0;
    double dy = 0.0;
    double q = 0.0;
    int index_i = 0;
    int index_j = 0;
    int count_temp=0;
    if( blockIdx.x >= arg->mesh_ynum) return;
    if( threadIdx.x >= arg->mesh_xnum) return;
    const int mesh_id = threadIdx.x + blockIdx.x * blockDim.x;
    
    for(int i=0;i<cuda->mesh_count[mesh_id];i++)
    {
        index_i = mesh_id + i*arg->mesh_num;
        //(x,y)->(x,y)
        for(int j=i+1;j<cuda->mesh_count[mesh_id];j++)
        {
            index_j = mesh_id + j*arg->mesh_num;
            dx = cuda->x[cuda->mesh[index_i]] - cuda->x[cuda->mesh[index_j]];
            dy = cuda->y[cuda->mesh[index_i]] - cuda->y[cuda->mesh[index_j]];
            q = sqrt(dx*dx+dy*dy)/arg->h;
            if(q<2.0)
            {
                if(cuda->type[cuda->mesh[index_i]] == 0)
                {
                    count_temp = atomicAdd(&(arg->pair_num),1);
                    cuda->pair_i[count_temp] = cuda->mesh[index_i];
                    cuda->pair_j[count_temp] = cuda->mesh[index_j];
                }
                else if(cuda->type[cuda->mesh[index_j]] == 0)
                {
                    count_temp = atomicAdd(&(arg->pair_num),1);
                    cuda->pair_i[count_temp] = cuda->mesh[index_j];
                    cuda->pair_j[count_temp] = cuda->mesh[index_i];
                }
            }
        }
        //(x,y)->(x+1,y)
        if( threadIdx.x<(arg->mesh_xnum-1))
        {
            for(int j=0;j<cuda->mesh_count[mesh+1];j++)
            {
                index_j = mesh_id + 1 + j*arg->mesh_num;
                dx = cuda->x[cuda->mesh[index_i]] - cuda->x[cuda->mesh[index_j]];
                dy = cuda->y[cuda->mesh[index_i]] - cuda->y[cuda->mesh[index_j]];
                q = sqrt(dx*dx+dy*dy)/arg->h;
                if(q<2.0)
                {
                    if(cuda->type[cuda->mesh[index_i]] == 0)
                    {
                        count_temp = atomicAdd(&(arg->pair_num),1);
                        cuda->pair_i[count_temp] = cuda->mesh[index_i];
                        cuda->pair_j[count_temp] = cuda->mesh[index_j];
                    }
                    else if(cuda->type[cuda->mesh[index_j]] == 0)
                    {
                        count_temp = atomicAdd(&(arg->pair_num),1);
                        cuda->pair_i[count_temp] = cuda->mesh[index_j];
                        cuda->pair_j[count_temp] = cuda->mesh[index_i];
                    }
                }
            }   
        }

        //(x,y)->(x,y+1)
        if( blockIdx.x<(arg->mesh_ynum-1))
        {
            for(int j=0;j<cuda->mesh_count[mesh_id+arg->mesh_xnum];j++)
            {
                index_j = mesh_id + arg->mesh_xnum + j*arg->mesh_num;
                dx = cuda->x[cuda->mesh[index_i]] - cuda->x[cuda->mesh[index_j]];
                dy = cuda->y[cuda->mesh[index_i]] - cuda->y[cuda->mesh[index_j]];
                q = sqrt(dx*dx+dy*dy)/arg->h;
                if(q<2.0)
                {
                    if(cuda->type[cuda->mesh[index_i]] == 0)
                    {
                        count_temp = atomicAdd(&(arg->pair_num),1);
                        cuda->pair_i[count_temp] = cuda->mesh[index_i];
                        cuda->pair_j[count_temp] = cuda->mesh[index_j];
                    }
                    else if(cuda->type[cuda->mesh[index_j]] == 0)
                    {
                        count_temp = atomicAdd(&(arg->pair_num),1);
                        cuda->pair_i[count_temp] = cuda->mesh[index_j];
                        cuda->pair_j[count_temp] = cuda->mesh[index_i];
                    }
                }
            }
        }

        //(x,y)->(x+1,y+1)
        if( threadIdx.x<(arg->mesh_xnum-1) && blockIdx.x<(arg->mesh_ynum-1))
        {
            for(int j=0;j<cuda->mesh_count[mesh_id+1+arg->mesh_xnum];j++)
            {
                index_j = mesh_id + 1 + arg->mesh_xnum + j*arg->mesh_num;
                dx = cuda->x[cuda->mesh[index_i]] - cuda->x[cuda->mesh[index_j]];
                dy = cuda->y[cuda->mesh[index_i]] - cuda->y[cuda->mesh[index_j]];
                q = sqrt(dx*dx+dy*dy)/arg->h;
                if(q<2.0)
                {
                    if(cuda->type[cuda->mesh[index_i]] == 0)
                    {
                        count_temp = atomicAdd(&(arg->pair_num),1);
                        cuda->pair_i[count_temp] = cuda->mesh[index_i];
                        cuda->pair_j[count_temp] = cuda->mesh[index_j];
                    }
                    else if(cuda->type[cuda->mesh[index_j]] == 0)
                    {
                        count_temp = atomicAdd(&(arg->pair_num),1);
                        cuda->pair_i[count_temp] = cuda->mesh[index_j];
                        cuda->pair_j[count_temp] = cuda->mesh[index_i];
                    }
                }
            }
        }

        //(x,y)->(x-1,y+1)
        if( threadIdx.x>0 && blockIdx.x<(arg->mesh_ynum-1))
        {
            for(int j=0;j<cuda->mesh_count[mesh_id-1+arg->mesh_xnum];j++)
            {
                index_j = mesh_id - 1 + arg->mesh_xnum + j*arg->mesh_num;
                dx = cuda->x[cuda->mesh[index_i]] - cuda->x[cuda->mesh[index_j]];
                dy = cuda->y[cuda->mesh[index_i]] - cuda->y[cuda->mesh[index_j]];
                q = sqrt(dx*dx+dy*dy)/arg->h;
                if(q<2.0)
                {
                    if(cuda->type[cuda->mesh[index_i]] == 0)
                    {
                        count_temp = atomicAdd(&(arg->pair_num),1);
                        cuda->pair_i[count_temp] = cuda->mesh[index_i];
                        cuda->pair_j[count_temp] = cuda->mesh[index_j];
                    }
                    else if(cuda->type[cuda->mesh[index_j]] == 0)
                    {
                        count_temp = atomicAdd(&(arg->pair_num),1);
                        cuda->pair_i[count_temp] = cuda->mesh[index_j];
                        cuda->pair_j[count_temp] = cuda->mesh[index_i];
                    }
                }
            }
        }
    }
    cuda->mesh_count[mesh_id] = 0;
}

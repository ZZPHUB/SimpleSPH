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

    //blockIdx.x ---> length direction
    //blockIdx.y ---> deepth direction
    //threadIdx.x ---> search the mesh
    //threadIdx.y---> search the near mesh


    double q;
    int i,j;
    //int count_temp=0;
    const int mesh_id = blockIdx.x + blockIdx.y * gridDim.x;
    if( threadIdx.x >= cuda->mesh_count[mesh_id]) return;
    //if(threadIdx.y == 0)atomicAdd(&(arg->tmp),1);
    i = mesh_id + threadIdx.x*arg->mesh_num;
    
    //mesh[x,y]->mesh[x,y]
    if( threadIdx.y > threadIdx.x && threadIdx.y< cuda->mesh_count[mesh_id])
    {
        j = mesh_id + threadIdx.y*arg->mesh_num ;
        q = (cuda->x[cuda->mesh[i]]-cuda->x[cuda->mesh[j]])*(cuda->x[cuda->mesh[i]]-cuda->x[cuda->mesh[j]])+(cuda->y[cuda->mesh[i]]-cuda->y[cuda->mesh[j]])*(cuda->y[cuda->mesh[i]]-cuda->y[cuda->mesh[j]]);
        q = sqrt(q)/arg->h;
        if(q<2.0)
        {
            if(cuda->type[cuda->mesh[i]] == 0)
            {
                sph_lock_cuda(arg);
                cuda->pair_i[arg->pair_num] = cuda->mesh[i];
                cuda->pair_j[arg->pair_num] = cuda->mesh[j];
                arg->pair_num++;
                sph_unlock_cuda(arg);
                //count_temp = atomicAdd(&(arg->pair_num),1);
                //cuda->pair_i[count_temp] = cuda->mesh[i];
                //cuda->pair_j[count_temp] = cuda->mesh[j];
                if(cuda->mesh[i] == 0 || cuda->mesh[j] == 0)
                {
                    //printf("i:%d type:%d j:%d type:%d\n",i,cuda->type[cuda->mesh[i]],j,cuda->type[cuda->mesh[j]]);
                }
            }
            else if(cuda->type[cuda->mesh[j]] == 0)
            {
                sph_lock_cuda(arg);
                cuda->pair_i[arg->pair_num] = cuda->mesh[j];
                cuda->pair_j[arg->pair_num] = cuda->mesh[i];
                arg->pair_num++;
                sph_unlock_cuda(arg);
                //count_temp = atomicAdd(&(arg->pair_num),1);
                //cuda->pair_i[count_temp] = cuda->mesh[j];
                //cuda->pair_j[count_temp] = cuda->mesh[i];
                if(cuda->mesh[i] == 0 || cuda->mesh[j] == 0)
                {
                    //printf("i:%d type:%d j:%d type:%d\n",cuda->mesh[i],cuda->type[cuda->mesh[i]],cuda->mesh[j],cuda->type[cuda->mesh[j]]);
                    //printf("i:%d type:%d j:%d type:%d\n",i,cuda->type[cuda->mesh[i]],j,cuda->type[cuda->mesh[j]]);
                }
            }
        }
    }
    
    //mesh[x,y]->mesh[x+1,y]
    if( blockIdx.x < ( gridDim.x-1))
    {
        if( threadIdx.y< cuda->mesh_count[mesh_id+1] )
        {
            j = mesh_id + 1 + threadIdx.y*arg->mesh_num;
            q = (cuda->x[cuda->mesh[i]]-cuda->x[cuda->mesh[j]])*(cuda->x[cuda->mesh[i]]-cuda->x[cuda->mesh[j]])+(cuda->y[cuda->mesh[i]]-cuda->y[cuda->mesh[j]])*(cuda->y[cuda->mesh[i]]-cuda->y[cuda->mesh[j]]);
            q = sqrt(q)/arg->h;
            if(q<2.0)
            {
                if(cuda->type[cuda->mesh[i]] == 0)
                {
                    sph_lock_cuda(arg);
                    cuda->pair_i[arg->pair_num] = cuda->mesh[i];
                    cuda->pair_j[arg->pair_num] = cuda->mesh[j];
                    arg->pair_num++;
                    sph_unlock_cuda(arg);
                    //count_temp = atomicAdd(&(arg->pair_num),1);
                    //cuda->pair_i[count_temp] = cuda->mesh[i];
                    //cuda->pair_j[count_temp] = cuda->mesh[j];
                    if(cuda->mesh[i] == 0 || cuda->mesh[j] == 0)
                    {
                        //printf("i:%d type:%d j:%d type:%d\n",cuda->mesh[i],cuda->type[cuda->mesh[i]],cuda->mesh[j],cuda->type[cuda->mesh[j]]);
                        //printf("i:%d type:%d j:%d type:%d\n",i,cuda->type[cuda->mesh[i]],j,cuda->type[cuda->mesh[j]]);
                    }
                }
                else if(cuda->type[cuda->mesh[j]] == 0)
                {
                    sph_lock_cuda(arg);
                    cuda->pair_i[arg->pair_num] = cuda->mesh[j];
                    cuda->pair_j[arg->pair_num] = cuda->mesh[i];
                    arg->pair_num++;
                    sph_unlock_cuda(arg);
                    //count_temp = atomicAdd(&(arg->pair_num),1);
                    //cuda->pair_i[count_temp] = cuda->mesh[j];
                    //cuda->pair_j[count_temp] = cuda->mesh[i];
                    if(cuda->mesh[i] == 0 || cuda->mesh[j] == 0)
                    {
                        //printf("i:%d type:%d j:%d type:%d\n",cuda->mesh[i],cuda->type[cuda->mesh[i]],cuda->mesh[j],cuda->type[cuda->mesh[j]]);
                        //printf("i:%d type:%d j:%d type:%d\n",i,cuda->type[cuda->mesh[i]],j,cuda->type[cuda->mesh[j]]);
                    }
                }
            }
        }
    }

    //mesh[x,y]->mesh[x,y+1]
    if( blockIdx.y < ( gridDim.y-1))
    {
        if( threadIdx.y< cuda->mesh[ mesh_id + gridDim.x] )
        {
            j = mesh_id + gridDim.x + threadIdx.y*arg->mesh_num;
            q = (cuda->x[cuda->mesh[i]]-cuda->x[cuda->mesh[j]])*(cuda->x[cuda->mesh[i]]-cuda->x[cuda->mesh[j]])+(cuda->y[cuda->mesh[i]]-cuda->y[cuda->mesh[j]])*(cuda->y[cuda->mesh[i]]-cuda->y[cuda->mesh[j]]);
            q = sqrt(q)/arg->h;
            if(q<2.0)
            {
                if(cuda->type[cuda->mesh[i]] == 0)
                {
                    sph_lock_cuda(arg);
                    cuda->pair_i[arg->pair_num] = cuda->mesh[i];
                    cuda->pair_j[arg->pair_num] = cuda->mesh[j];
                    arg->pair_num++;
                    sph_unlock_cuda(arg);
                    //count_temp = atomicAdd(&(arg->pair_num),1);
                    //cuda->pair_i[count_temp] = cuda->mesh[i];
                    //cuda->pair_j[count_temp] = cuda->mesh[j];
                    if(cuda->mesh[i] == 0 || cuda->mesh[j] == 0)
                    {
                        //printf("i:%d type:%d j:%d type:%d\n",cuda->mesh[i],cuda->type[cuda->mesh[i]],cuda->mesh[j],cuda->type[cuda->mesh[j]]);
                        //printf("i:%d type:%d j:%d type:%d\n",i,cuda->type[cuda->mesh[i]],j,cuda->type[cuda->mesh[j]]);
                    }
                }
                else if(cuda->type[cuda->mesh[j]] == 0)
                {
                    sph_lock_cuda(arg);
                    cuda->pair_i[arg->pair_num] = cuda->mesh[j];
                    cuda->pair_j[arg->pair_num] = cuda->mesh[i];
                    arg->pair_num++;
                    sph_unlock_cuda(arg);
                    //count_temp = atomicAdd(&(arg->pair_num),1);
                    //cuda->pair_i[count_temp] = cuda->mesh[j];
                    //cuda->pair_j[count_temp] = cuda->mesh[i];
                    if(cuda->mesh[i] == 0 || cuda->mesh[j] == 0)
                    {
                        //printf("i:%d type:%d j:%d type:%d\n",cuda->mesh[i],cuda->type[cuda->mesh[i]],cuda->mesh[j],cuda->type[cuda->mesh[j]]);
                        //printf("i:%d type:%d j:%d type:%d\n",i,cuda->type[cuda->mesh[i]],j,cuda->type[cuda->mesh[j]]);
                    }
                }
            }
        }
    }

    //mesh[x,y]->mesh[x+1,y+1]
    if( blockIdx.x < ( gridDim.x-1) && blockIdx.y < ( gridDim.y-1))
    {
        if( threadIdx.y< cuda->mesh_count[mesh_id + 1 + gridDim.x])
        {
            j = mesh_id + 1 + gridDim.x + threadIdx.y*arg->mesh_num;
            q = (cuda->x[cuda->mesh[i]]-cuda->x[cuda->mesh[j]])*(cuda->x[cuda->mesh[i]]-cuda->x[cuda->mesh[j]])+(cuda->y[cuda->mesh[i]]-cuda->y[cuda->mesh[j]])*(cuda->y[cuda->mesh[i]]-cuda->y[cuda->mesh[j]]);
            q = sqrt(q)/arg->h;
            if(q<2.0)
            {
                if(cuda->type[cuda->mesh[i]] == 0)
                {
                    sph_lock_cuda(arg);
                    cuda->pair_i[arg->pair_num] = cuda->mesh[i];
                    cuda->pair_j[arg->pair_num] = cuda->mesh[j];
                    arg->pair_num++;
                    sph_unlock_cuda(arg);
                    //count_temp = atomicAdd(&(arg->pair_num),1);
                    //cuda->pair_i[count_temp] = cuda->mesh[i];
                    //cuda->pair_j[count_temp] = cuda->mesh[j];
                    if(cuda->mesh[i] == 0 || cuda->mesh[j] == 0)
                    {
                       // printf("i:%d type:%d j:%d type:%d\n",cuda->mesh[i],cuda->type[cuda->mesh[i]],cuda->mesh[j],cuda->type[cuda->mesh[j]]);
                       // printf("i:%d type:%d j:%d type:%d\n",i,cuda->type[cuda->mesh[i]],j,cuda->type[cuda->mesh[j]]);
                    }
                }
                else if(cuda->type[cuda->mesh[j]] == 0)
                {
                    sph_lock_cuda(arg);
                    cuda->pair_i[arg->pair_num] = cuda->mesh[j];
                    cuda->pair_j[arg->pair_num] = cuda->mesh[i];
                    arg->pair_num++;
                    sph_unlock_cuda(arg);
                    //count_temp = atomicAdd(&(arg->pair_num),1);
                    //cuda->pair_i[count_temp] = cuda->mesh[j];
                    //cuda->pair_j[count_temp] = cuda->mesh[i];
                    if(cuda->mesh[i] == 0 || cuda->mesh[j] == 0)
                    {
                        //printf("i:%d type:%d j:%d type:%d\n",cuda->mesh[i],cuda->type[cuda->mesh[i]],cuda->mesh[j],cuda->type[cuda->mesh[j]]);
                        //printf("i:%d type:%d j:%d type:%d\n",i,cuda->type[cuda->mesh[i]],j,cuda->type[cuda->mesh[j]]);
                    }
                }
            }
        }
    }

    //mesh[x,y]->mesh[x+1,y-1]
    if( blockIdx.x < ( gridDim.x-1) && blockIdx.y > 0)
    {
        if( threadIdx.y< cuda->mesh[mesh_id + 1 - gridDim.x])
        {
            j = mesh_id + 1 - gridDim.x + threadIdx.y*arg->mesh_num;
            q = (cuda->x[cuda->mesh[i]]-cuda->x[cuda->mesh[j]])*(cuda->x[cuda->mesh[i]]-cuda->x[cuda->mesh[j]])+(cuda->y[cuda->mesh[i]]-cuda->y[cuda->mesh[j]])*(cuda->y[cuda->mesh[i]]-cuda->y[cuda->mesh[j]]);
            q = sqrt(q)/arg->h;
            if(q<2.0)
            {
                if(cuda->type[cuda->mesh[i]] == 0)
                {
                    sph_lock_cuda(arg);
                    cuda->pair_i[arg->pair_num] = cuda->mesh[i];
                    cuda->pair_j[arg->pair_num] = cuda->mesh[j];
                    arg->pair_num++;
                    sph_unlock_cuda(arg);
                    //count_temp = atomicAdd(&(arg->pair_num),1);
                    //cuda->pair_i[count_temp] = cuda->mesh[i];
                    //cuda->pair_j[count_temp] = cuda->mesh[j];
                    if(cuda->mesh[i] == 0 || cuda->mesh[j] == 0)
                    {
                        //printf("i:%d type:%d j:%d type:%d\n",cuda->mesh[i],cuda->type[cuda->mesh[i]],cuda->mesh[j],cuda->type[cuda->mesh[j]]);
                        //printf("i:%d type:%d j:%d type:%d\n",i,cuda->type[cuda->mesh[i]],j,cuda->type[cuda->mesh[j]]);
                    }
                }
                else if(cuda->type[cuda->mesh[j]] == 0)
                {
                    sph_lock_cuda(arg);
                    cuda->pair_i[arg->pair_num] = cuda->mesh[j];
                    cuda->pair_j[arg->pair_num] = cuda->mesh[i];
                    arg->pair_num++;
                    sph_unlock_cuda(arg);
                    //count_temp = atomicAdd(&(arg->pair_num),1);
                    //cuda->pair_i[count_temp] = cuda->mesh[j];
                    //cuda->pair_j[count_temp] = cuda->mesh[i];
                    if(cuda->mesh[i] == 0 || cuda->mesh[j] == 0)
                    {
                        //printf("i:%d type:%d j:%d type:%d\n",cuda->mesh[i],cuda->type[cuda->mesh[i]],cuda->mesh[j],cuda->type[cuda->mesh[j]]);
                        //printf("i:%d type:%d j:%d type:%d\n",i,cuda->type[cuda->mesh[i]],j,cuda->type[cuda->mesh[j]]);
                    }
                }
            }
        }
    }
    cuda->mesh_count[mesh_id] = 0;
}

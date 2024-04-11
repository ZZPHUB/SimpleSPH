#include "Lib.cuh"

__global__ void sph_nnps_cuda(SPH_CUDA *cuda,SPH_ARG *arg,SPH_RIGID *rigid)
{

    //blockIdx.x ---> length direction
    //blockIdx.y ---> deepth direction
    //threadIdx.x ---> search the mesh
    //threadIdx.y---> search the near mesh


    double q;
    int i,j;
    int count_temp=0;
    const int mesh_id = blockIdx.x + blockIdx.y * gridDim.x;
    if( threadIdx.x >= cuda->mesh_count[mesh_id]) return;
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
                count_temp = atomicAdd(&(arg->pair_num),1);
                cuda->pair_i[count_temp] = cuda->mesh[i];
                cuda->pair_j[count_temp] = cuda->mesh[j];
            }
            else if(cuda->type[cuda->mesh[j]] == 0)
            {
                count_temp = atomicAdd(&(arg->pair_num),1);
                cuda->pair_i[count_temp] = cuda->mesh[j];
                cuda->pair_j[count_temp] = cuda->mesh[i];
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
                    count_temp = atomicAdd(&(arg->pair_num),1);
                    cuda->pair_i[count_temp] = cuda->mesh[i];
                    cuda->pair_j[count_temp] = cuda->mesh[j];
                }
                else if(cuda->type[cuda->mesh[j]] == 0)
                {
                    count_temp = atomicAdd(&(arg->pair_num),1);
                    cuda->pair_i[count_temp] = cuda->mesh[j];
                    cuda->pair_j[count_temp] = cuda->mesh[i];
                }
            }
        }
    }

    //mesh[x,y]->mesh[x,y+1]
    if( blockIdx.y < ( gridDim.y-1))
    {
        mesh_near_num = mesh[ mesh_id + gridDim.x];
        if( threadIdx.y< mesh_near_num )
        {
            j = mesh_id + gridDim.x + threadIdx.y*arg->mesh_num;
            q = (cuda->x[cuda->mesh[i]]-cuda->x[cuda->mesh[j]])*(cuda->x[cuda->mesh[i]]-cuda->x[cuda->mesh[j]])+(cuda->y[cuda->mesh[i]]-cuda->y[cuda->mesh[j]])*(cuda->y[cuda->mesh[i]]-cuda->y[cuda->mesh[j]]);
            q = sqrt(q)/arg->h;
            if(q<2.0)
            {
                if(cuda->type[cuda->mesh[i]] == 0)
                {
                    count_temp = atomicAdd(&(arg->pair_num),1);
                    cuda->pair_i[count_temp] = cuda->mesh[i];
                    cuda->pair_j[count_temp] = cuda->mesh[j];
                }
                else if(cuda->type[cuda->mesh[j]] == 0)
                {
                    count_temp = atomicAdd(&(arg->pair_num),1);
                    cuda->pair_i[count_temp] = cuda->mesh[j];
                    cuda->pair_j[count_temp] = cuda->mesh[i];
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
                    count_temp = atomicAdd(&(arg->pair_num),1);
                    cuda->pair_i[count_temp] = cuda->mesh[i];
                    cuda->pair_j[count_temp] = cuda->mesh[j];
                }
                else if(cuda->type[cuda->mesh[j]] == 0)
                {
                    count_temp = atomicAdd(&(arg->pair_num),1);
                    cuda->pair_i[count_temp] = cuda->mesh[j];
                    cuda->pair_j[count_temp] = cuda->mesh[i];
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
                if(cuda->type[cuda->cuda->mesh[i]] == 0)
                {
                    count_temp = atomicAdd(&(arg->pair_num),1);
                    cuda->pair_i[count_temp] = cuda->mesh[i];
                    cuda->pair_j[count_temp] = cuda->mesh[j];
                }
                else if(cuda->type[cuda->mesh[j]] == 0)
                {
                    count_temp = atomicAdd(&(arg->pair_num),1);
                    cuda->pair_i[count_temp] = cuda->mesh[j];
                    cuda->pair_j[count_temp] = cuda->mesh[i];
                }
            }
        }
    }
    cuda->mesh[mesh_id] = 0;
}

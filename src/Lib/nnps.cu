#include "Lib.cuh"

__global__ void sph_nnps_cuda(int *mesh,double *x,double *y,int *type,int *pair_i,int *pair_j)
{
    /*
    blockIdx.x ---> length direction
    blockIdx.y ---> deepth direction
    threadIdx.x ---> search the mesh
    threadIdx.y ---> search the near mesh
    */
    double q;
    int i,j;
    int count_temp;
    int mesh_ptc_num;
    int mesh_near_ptc_num;
    if( blockIdx.x >= MESH_LENGTH_NUM_CUDA || blockIdx.y >= MESH_DEEPTH_NUM_CUDA) return;
    mesh_ptc_num = mesh[ blockIdx.x + blockIdx.y*MESH_LENGTH_NUM_CUDA + MESH_DEEPTH_NUM_CUDA*MESH_LENGTH_NUM_CUDA*(MESH_PTC_NUM-1)];
    if( threadIdx.x >= mesh_ptc_num)return;
    i = blockIdx.x + blockIdx.y*MESH_LENGTH_NUM_CUDA + threadIdx.x*MESH_DEEPTH_NUM_CUDA*MESH_LENGTH_NUM_CUDA;
    

    //mesh[i,j]->mesh[i,j]
    mesh_near_ptc_num = mesh_ptc_num;
    if( threadIdx.y > threadIdx.x && threadIdx.y < mesh_near_ptc_num)
    {
        j = blockIdx.x + blockIdx.y*MESH_LENGTH_NUM_CUDA + threadIdx.y*MESH_DEEPTH_NUM_CUDA*MESH_LENGTH_NUM_CUDA;
        q = (x[mesh[i]]-x[mesh[j]])*(x[mesh[i]]-x[mesh[j]])+(y[mesh[i]]-y[mesh[j]])*(y[mesh[i]]-y[mesh[j]]);
        q = sqrt(q);
        if(q<2.0)
        {
            if(type[mesh[i]] == 0)
            {
                count_temp = atomicAdd(&dev_count,1);
                pair_i[count_temp] = i;
                pair_j[count_temp] = j;
            }
            else if(type[mesh[j]] == 0)
            {
                count_temp = atomicAdd(&dev_count,1);
                pair_i[count_temp] = j;
                pair_j[count_temp] = i;
            }
        }
    }

    //mesh[i,j]->mesh[i,j+1]
    if( blockIdx.x < (MESH_LENGTH_NUM_CUDA-1))
    {
        mesh_near_ptc_num = mesh[ (blockIdx.x+1) + blockIdx.y*MESH_LENGTH_NUM_CUDA + MESH_DEEPTH_NUM_CUDA*MESH_LENGTH_NUM_CUDA*(MESH_PTC_NUM-1)];
        if( threadIdx.y < mesh_near_ptc_num )
        {
            j = (blockIdx.x+1) + blockIdx.y*MESH_LENGTH_NUM_CUDA + threadIdx.y*MESH_DEEPTH_NUM_CUDA*MESH_LENGTH_NUM_CUDA;
            q = (x[mesh[i]]-x[mesh[j]])*(x[mesh[i]]-x[mesh[j]])+(y[mesh[i]]-y[mesh[j]])*(y[mesh[i]]-y[mesh[j]]);
            q = sqrt(q);
            if(q<2.0)
            {
                if(type[mesh[i]] == 0)
                {
                    count_temp = atomicAdd(&dev_count,1);
                    pair_i[count_temp] = i;
                    pair_j[count_temp] = j;
                }
                else if(type[mesh[j]] == 0)
                {
                    count_temp = atomicAdd(&dev_count,1);
                    pair_i[count_temp] = j;
                    pair_j[count_temp] = i;
                }
            }
        }
    }

    //mesh[i,j]->mesh[i+1,j]
    if( blockIdx.y < (MESH_DEEPTH_NUM_CUDA-1))
    {
        mesh_near_ptc_num = mesh[ blockIdx.x + ( blockIdx.y+1)*MESH_LENGTH_NUM_CUDA + MESH_DEEPTH_NUM_CUDA*MESH_LENGTH_NUM_CUDA*(MESH_PTC_NUM-1)];
        if( threadIdx.y < mesh_near_ptc_num )
        {
            j = blockIdx.x +( blockIdx.y+1)*MESH_LENGTH_NUM_CUDA + threadIdx.y*MESH_DEEPTH_NUM_CUDA*MESH_LENGTH_NUM_CUDA;
            q = (x[mesh[i]]-x[mesh[j]])*(x[mesh[i]]-x[mesh[j]])+(y[mesh[i]]-y[mesh[j]])*(y[mesh[i]]-y[mesh[j]]);
            q = sqrt(q);
            if(q<2.0)
            {
                if(type[mesh[i]] == 0)
                {
                    count_temp = atomicAdd(&dev_count,1);
                    pair_i[count_temp] = i;
                    pair_j[count_temp] = j;
                }
                else if(type[mesh[j]] == 0)
                {
                    count_temp = atomicAdd(&dev_count,1);
                    pair_i[count_temp] = j;
                    pair_j[count_temp] = i;
                }
            }
        }
    }

    //mesh[i,j]->mesh[i+1,j+1]
    if( blockIdx.x < (MESH_LENGTH_NUM_CUDA-1) && blockIdx.y < (MESH_DEEPTH_NUM_CUDA-1))
    {
        mesh_near_ptc_num = mesh[( blockIdx.x+1) + ( blockIdx.y+1)*MESH_LENGTH_NUM_CUDA + MESH_DEEPTH_NUM_CUDA*MESH_LENGTH_NUM_CUDA*(MESH_PTC_NUM-1)];
        if( threadIdx.y < mesh_near_ptc_num)
        {
            j = ( blockIdx.x+1) +( blockIdx.y+1)*MESH_LENGTH_NUM_CUDA + threadIdx.y*MESH_DEEPTH_NUM_CUDA*MESH_LENGTH_NUM_CUDA;
            q = (x[mesh[i]]-x[mesh[j]])*(x[mesh[i]]-x[mesh[j]])+(y[mesh[i]]-y[mesh[j]])*(y[mesh[i]]-y[mesh[j]]);
            q = sqrt(q);
            if(q<2.0)
            {
                if(type[mesh[i]] == 0)
                {
                    count_temp = atomicAdd(&dev_count,1);
                    pair_i[count_temp] = i;
                    pair_j[count_temp] = j;
                }
                else if(type[mesh[j]] == 0)
                {
                    count_temp = atomicAdd(&dev_count,1);
                    pair_i[count_temp] = j;
                    pair_j[count_temp] = i;
                }
            }
        }
    }

    //mesh[i,j]->mesh[i-1,j+1]
    if( blockIdx.x < (MESH_LENGTH_NUM_CUDA-1) && blockIdx.y > 0)
    {
        mesh_near_ptc_num = mesh[( blockIdx.x+1) + ( blockIdx.y-1)*MESH_LENGTH_NUM_CUDA + MESH_DEEPTH_NUM_CUDA*MESH_LENGTH_NUM_CUDA*(MESH_PTC_NUM-1)];
        if( threadIdx.y < mesh_near_ptc_num)
        {
            j = ( blockIdx.x+1) +( blockIdx.y-1)*MESH_LENGTH_NUM_CUDA + threadIdx.y*MESH_DEEPTH_NUM_CUDA*MESH_LENGTH_NUM_CUDA;
            q = (x[mesh[i]]-x[mesh[j]])*(x[mesh[i]]-x[mesh[j]])+(y[mesh[i]]-y[mesh[j]])*(y[mesh[i]]-y[mesh[j]]);
            q = sqrt(q);
            if(q<2.0)
            {
                if(type[mesh[i]] == 0)
                {
                    count_temp = atomicAdd(&dev_count,1);
                    pair_i[count_temp] = i;
                    pair_j[count_temp] = j;
                }
                else if(type[mesh[j]] == 0)
                {
                    count_temp = atomicAdd(&dev_count,1);
                    pair_i[count_temp] = j;
                    pair_j[count_temp] = i;
                }
            }
        }
    }
}

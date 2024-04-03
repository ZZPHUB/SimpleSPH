#include "Lib.cuh"

__global__ void sph_nnps_cuda(int *mesh,double *x,double *y,int *type,int *pair_i,int *pair_j,int *count)
{
    /*
    blockIdx.x ---> length direction
    blockIdx.y ---> deepth direction
    threadIdx.x ---> search the mesh
    blockIdx.z---> search the near mesh
    */
    double q;
    int i,j;
    int count_temp;
    int mesh_ptc_num;
    int mesh_near_ptc_num;
    if( blockIdx.x >= dev_mesh_lnum || blockIdx.y >= dev_mesh_dnum ) return;
    mesh_ptc_num = mesh[ blockIdx.x + blockIdx.y*dev_mesh_lnum + dev_mesh_tnum*(MESH_PTC_NUM-2)];
    mesh[ blockIdx.x + blockIdx.y*dev_mesh_lnum + dev_mesh_tnum*(MESH_PTC_NUM-2)] = 0;
    if( threadIdx.x >= mesh_ptc_num)return;
    i = blockIdx.x + blockIdx.y*dev_mesh_lnum + threadIdx.x*dev_mesh_tnum;
    

    //mesh[i,j]->mesh[i,j]
    mesh_near_ptc_num = mesh_ptc_num;
    if( blockIdx.z> threadIdx.x && blockIdx.z< mesh_near_ptc_num)
    {
        j = blockIdx.x + blockIdx.y*dev_mesh_lnum + blockIdx.z*dev_mesh_tnum;
        q = (x[mesh[i]]-x[mesh[j]])*(x[mesh[i]]-x[mesh[j]])+(y[mesh[i]]-y[mesh[j]])*(y[mesh[i]]-y[mesh[j]]);
        q = sqrt(q)/PTC_SML;
        if(q<2.0)
        {
            if(type[mesh[i]] == 0)
            {
                count_temp = atomicAdd(count,1);
                pair_i[count_temp] = mesh[i];
                pair_j[count_temp] = mesh[j];
            }
            else if(type[mesh[j]] == 0)
            {
                count_temp = atomicAdd(count,1);
                pair_i[count_temp] = mesh[j];
                pair_j[count_temp] = mesh[i];
            }
        }
    }

    //mesh[i,j]->mesh[i,j+1]
    if( blockIdx.x < (dev_mesh_lnum-1))
    {
        mesh_near_ptc_num = mesh[ (blockIdx.x+1) + blockIdx.y*dev_mesh_lnum + dev_mesh_tnum*(MESH_PTC_NUM-2)];
        if( blockIdx.z< mesh_near_ptc_num )
        {
            j = ( blockIdx.x +1) + blockIdx.y*dev_mesh_lnum + blockIdx.z*dev_mesh_tnum;
            q = (x[mesh[i]]-x[mesh[j]])*(x[mesh[i]]-x[mesh[j]])+(y[mesh[i]]-y[mesh[j]])*(y[mesh[i]]-y[mesh[j]]);
            q = sqrt(q)/PTC_SML;
            if(q<2.0)
            {
                if(type[mesh[i]] == 0)
                {
                    count_temp = atomicAdd(count,1);
                    pair_i[count_temp] = mesh[i];
                    pair_j[count_temp] = mesh[j];
                }
                else if(type[mesh[j]] == 0)
                {
                    count_temp = atomicAdd(count,1);
                    pair_i[count_temp] = mesh[j];
                    pair_j[count_temp] = mesh[i];
                }
            }
        }
    }

    //mesh[i,j]->mesh[i+1,j]
    if( blockIdx.y < (dev_mesh_dnum-1))
    {
        mesh_near_ptc_num = mesh[ blockIdx.x + ( blockIdx.y+1)*dev_mesh_lnum + dev_mesh_tnum*(MESH_PTC_NUM-2)];
        if( blockIdx.z< mesh_near_ptc_num )
        {
            j = blockIdx.x +( blockIdx.y+1)*dev_mesh_lnum + blockIdx.z*dev_mesh_tnum;
            q = (x[mesh[i]]-x[mesh[j]])*(x[mesh[i]]-x[mesh[j]])+(y[mesh[i]]-y[mesh[j]])*(y[mesh[i]]-y[mesh[j]]);
            q = sqrt(q)/PTC_SML;
            if(q<2.0)
            {
                if(type[mesh[i]] == 0)
                {
                    count_temp = atomicAdd(count,1);
                    pair_i[count_temp] = mesh[i];
                    pair_j[count_temp] = mesh[j];
                }
                else if(type[mesh[j]] == 0)
                {
                    count_temp = atomicAdd(count,1);
                    pair_i[count_temp] = mesh[j];
                    pair_j[count_temp] = mesh[i];
                }
            }
        }
    }

    //mesh[i,j]->mesh[i+1,j+1]
    if( blockIdx.x < (dev_mesh_lnum-1) && blockIdx.y < (dev_mesh_dnum-1))
    {
        mesh_near_ptc_num = mesh[( blockIdx.x+1) + ( blockIdx.y+1)*dev_mesh_lnum + dev_mesh_tnum*(MESH_PTC_NUM-2)];
        if( blockIdx.z< mesh_near_ptc_num)
        {
            j = ( blockIdx.x+1) +( blockIdx.y+1)*dev_mesh_lnum + blockIdx.z*dev_mesh_tnum;
            q = (x[mesh[i]]-x[mesh[j]])*(x[mesh[i]]-x[mesh[j]])+(y[mesh[i]]-y[mesh[j]])*(y[mesh[i]]-y[mesh[j]]);
            q = sqrt(q)/PTC_SML;
            if(q<2.0)
            {
                if(type[mesh[i]] == 0)
                {
                    count_temp = atomicAdd(count,1);
                    pair_i[count_temp] = mesh[i];
                    pair_j[count_temp] = mesh[j];
                }
                else if(type[mesh[j]] == 0)
                {
                    count_temp = atomicAdd(count,1);
                    pair_i[count_temp] = mesh[j];
                    pair_j[count_temp] = mesh[i];
                }
            }
        }
    }

    //mesh[i,j]->mesh[i-1,j+1]
    if( blockIdx.x < (dev_mesh_lnum-1) && blockIdx.y > 0)
    {
        mesh_near_ptc_num = mesh[( blockIdx.x+1) + ( blockIdx.y-1)*dev_mesh_lnum + dev_mesh_tnum*(MESH_PTC_NUM-2)];
        if( blockIdx.z< mesh_near_ptc_num)
        {
            j = ( blockIdx.x+1) +( blockIdx.y-1)*dev_mesh_lnum + blockIdx.z*dev_mesh_tnum;
            q = (x[mesh[i]]-x[mesh[j]])*(x[mesh[i]]-x[mesh[j]])+(y[mesh[i]]-y[mesh[j]])*(y[mesh[i]]-y[mesh[j]]);
            q = sqrt(q)/PTC_SML;
            if(q<2.0)
            {
                if(type[mesh[i]] == 0)
                {
                    count_temp = atomicAdd(count,1);
                    pair_i[count_temp] = mesh[i];
                    pair_j[count_temp] = mesh[j];
                }
                else if(type[mesh[j]] == 0)
                {
                    count_temp = atomicAdd(count,1);
                    pair_i[count_temp] = mesh[j];
                    pair_j[count_temp] = mesh[i];
                }
            }
        }
    }
}

/*#include "Lib.cuh"

__global__ void sph_nnps_cuda(int *mesh,double *x,double *y,int *type,int *pair_i,int *pair_j,int *count)
{

    //blockIdx.x ---> length direction
    //blockIdx.y ---> deepth direction
    //threadIdx.x ---> search the mesh
    //threadIdx.y---> search the near mesh
    k ---> search the near mesh


    double q;
    int i,j;
    int count_temp=0;
    int mesh_local_num;
    int mesh_near_num;
    //if( blockIdx.x >= dev_mesh_lnum || blockIdx.y >= dev_mesh_dnum ) return;
    const int mesh_id = blockIdx.x + blockIdx.y * gridDim.x;


    mesh_local_num = mesh[ mesh_id + dev_mesh_tnum*(MESH_PTC_NUM-2)];
    mesh[ mesh_id + dev_mesh_tnum*(MESH_PTC_NUM-2)] = 0;
    if( threadIdx.x >= mesh_local_num)return;
    i = mesh_id + threadIdx.x*dev_mesh_tnum;
    
    //mesh[x,y]-->mesh[x,y]
    mesh_near_num = mesh_local_num;
    for(int k=blockIdx.x+1;k<mesh_near_num;k++)
    {
        j = mesh_id + k * dev_mesh_tnum;
        q = (x[mesh[i]]-x[mesh[j]])*(x[mesh[i]]-x[mesh[j]])+(y[mesh[i]]-y[mesh[j]])*(y[mesh[i]]-y[mesh[j]]);
        q = sqrt(q)/dev_h;
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

    //mesh[x,y]-->mesh[x+1,y]
    if( blockIdx.x < ( gridDim.x-1))
    {
        mesh_near_num = mesh[mesh_id + 1 + dev_mesh_tnum*(MESH_PTC_NUM-2)];
        for(int k=0;k<mesh_near_num;k++)
        {
            j = mesh_id + 1 +dev_mesh_tnum*k;
            q = (x[mesh[i]]-x[mesh[j]])*(x[mesh[i]]-x[mesh[j]])+(y[mesh[i]]-y[mesh[j]])*(y[mesh[i]]-y[mesh[j]]);
            q = sqrt(q)/dev_h;
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

    //mesh[x,y]-->mesh[x+1,y+1]
    if( blockIdx.x<( gridDim.x-1) && blockIdx.y<( gridDim.y-1))
    {
        mesh_near_num = mesh[mesh_id + 1 + gridDim.x + dev_mesh_tnum*(MESH_PTC_NUM-2)];
        for(int k=0;k<mesh_near_num;k++)
        {
            j = mesh_id + 1 + gridDim.x + dev_mesh_tnum*k;
            q = (x[mesh[i]]-x[mesh[j]])*(x[mesh[i]]-x[mesh[j]])+(y[mesh[i]]-y[mesh[j]])*(y[mesh[i]]-y[mesh[j]]);
            q = sqrt(q)/dev_h;
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

    //mesh[x,y]-->mesh[x,y+1]
    if( blockIdx.y<( gridDim.y-1))
    {
        mesh_near_num = mesh[mesh_id + gridDim.x + dev_mesh_tnum*(MESH_PTC_NUM-2)];
        for(int k=0;k<mesh_near_num;k++)
        {
            j = mesh_id + gridDim.x + dev_mesh_tnum*k;
            q = (x[mesh[i]]-x[mesh[j]])*(x[mesh[i]]-x[mesh[j]])+(y[mesh[i]]-y[mesh[j]])*(y[mesh[i]]-y[mesh[j]]);
            q = sqrt(q)/dev_h;
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

    //mesh[x,y]-->mesh[x-1,y+1]
    if( blockIdx.x>0 && blockIdx.y<( gridDim.y-1))
    {
        mesh_near_num = mesh[mesh_id - 1 + gridDim.x + dev_mesh_tnum*(MESH_PTC_NUM-2)];
        for(int k=0;k<mesh_near_num;k++)
        {
            j = mesh_id - 1 + gridDim.x + dev_mesh_tnum*k;
            q = (x[mesh[i]]-x[mesh[j]])*(x[mesh[i]]-x[mesh[j]])+(y[mesh[i]]-y[mesh[j]])*(y[mesh[i]]-y[mesh[j]]);
            q = sqrt(q)/dev_h;
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
*/


    /*
    //mesh[i,j]->mesh[i,j]
    mesh_near_num = mesh_local_num;
    if( threadIdx.y > threadIdx.x && threadIdx.y< mesh_near_num)
    {
        j = blockIdx.x + blockIdx.y*dev_mesh_lnum + threadIdx.y*dev_mesh_tnum;
        q = (x[mesh[i]]-x[mesh[j]])*(x[mesh[i]]-x[mesh[j]])+(y[mesh[i]]-y[mesh[j]])*(y[mesh[i]]-y[mesh[j]]);
        q = sqrt(q)/dev_h;
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
        mesh_near_num = mesh[ (blockIdx.x+1) + blockIdx.y*dev_mesh_lnum + dev_mesh_tnum*(MESH_PTC_NUM-2)];
        if( threadIdx.y< mesh_near_num )
        {
            j = ( blockIdx.x +1) + blockIdx.y*dev_mesh_lnum + threadIdx.y*dev_mesh_tnum;
            q = (x[mesh[i]]-x[mesh[j]])*(x[mesh[i]]-x[mesh[j]])+(y[mesh[i]]-y[mesh[j]])*(y[mesh[i]]-y[mesh[j]]);
            q = sqrt(q)/dev_h;
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
        mesh_near_num = mesh[ blockIdx.x + ( blockIdx.y+1)*dev_mesh_lnum + dev_mesh_tnum*(MESH_PTC_NUM-2)];
        if( threadIdx.y< mesh_near_num )
        {
            j = blockIdx.x +( blockIdx.y+1)*dev_mesh_lnum + threadIdx.y*dev_mesh_tnum;
            q = (x[mesh[i]]-x[mesh[j]])*(x[mesh[i]]-x[mesh[j]])+(y[mesh[i]]-y[mesh[j]])*(y[mesh[i]]-y[mesh[j]]);
            q = sqrt(q)/dev_h;
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
        mesh_near_num = mesh[( blockIdx.x+1) + ( blockIdx.y+1)*dev_mesh_lnum + dev_mesh_tnum*(MESH_PTC_NUM-2)];
        if( threadIdx.y< mesh_near_num)
        {
            j = ( blockIdx.x+1) +( blockIdx.y+1)*dev_mesh_lnum + threadIdx.y*dev_mesh_tnum;
            q = (x[mesh[i]]-x[mesh[j]])*(x[mesh[i]]-x[mesh[j]])+(y[mesh[i]]-y[mesh[j]])*(y[mesh[i]]-y[mesh[j]]);
            q = sqrt(q)/dev_h;
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
        mesh_near_num = mesh[( blockIdx.x+1) + ( blockIdx.y-1)*dev_mesh_lnum + dev_mesh_tnum*(MESH_PTC_NUM-2)];
        if( threadIdx.y< mesh_near_num)
        {
            j = ( blockIdx.x+1) +( blockIdx.y-1)*dev_mesh_lnum + threadIdx.y*dev_mesh_tnum;
            q = (x[mesh[i]]-x[mesh[j]])*(x[mesh[i]]-x[mesh[j]])+(y[mesh[i]]-y[mesh[j]])*(y[mesh[i]]-y[mesh[j]]);
            q = sqrt(q)/dev_h;
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

}*/

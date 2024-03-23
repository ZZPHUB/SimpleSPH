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
        q = (x[mesh[i]]-x[mesh[j]])**(x[mesh[i]]-x[mesh[j]])+(y[mesh[i]]-y[mesh[j]])**(y[mesh[i]]-y[mesh[j]]);
        q = sqrt(q);
        if(q<2.0)
        {
            if(type[mesh[i]] == 0)
            {
                count_temp = atomicAdd(&count,1);
                pair_i[count_temp] = i;
                pair_j[count_temp] = j;
            }
            else if(type[mesh[j]] == 0)
            {
                count_temp = atomicAdd(&count,1);
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
            q = (x[mesh[i]]-x[mesh[j]])**(x[mesh[i]]-x[mesh[j]])+(y[mesh[i]]-y[mesh[j]])**(y[mesh[i]]-y[mesh[j]]);
            q = sqrt(q);
            if(q<2.0)
            {
                if(type[mesh[i]] == 0)
                {
                    count_temp = atomicAdd(&count,1);
                    pair_i[count_temp] = i;
                    pair_j[count_temp] = j;
                }
                else if(type[mesh[j]] == 0)
                {
                    count_temp = atomicAdd(&count,1);
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
            q = (x[mesh[i]]-x[mesh[j]])**(x[mesh[i]]-x[mesh[j]])+(y[mesh[i]]-y[mesh[j]])**(y[mesh[i]]-y[mesh[j]]);
            q = sqrt(q);
            if(q<2.0)
            {
                if(type[mesh[i]] == 0)
                {
                    count_temp = atomicAdd(&count,1);
                    pair_i[count_temp] = i;
                    pair_j[count_temp] = j;
                }
                else if(type[mesh[j]] == 0)
                {
                    count_temp = atomicAdd(&count,1);
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
            q = (x[mesh[i]]-x[mesh[j]])**(x[mesh[i]]-x[mesh[j]])+(y[mesh[i]]-y[mesh[j]])**(y[mesh[i]]-y[mesh[j]]);
            q = sqrt(q);
            if(q<2.0)
            {
                if(type[mesh[i]] == 0)
                {
                    count_temp = atomicAdd(&count,1);
                    pair_i[count_temp] = i;
                    pair_j[count_temp] = j;
                }
                else if(type[mesh[j]] == 0)
                {
                    count_temp = atomicAdd(&count,1);
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
            q = (x[mesh[i]]-x[mesh[j]])**(x[mesh[i]]-x[mesh[j]])+(y[mesh[i]]-y[mesh[j]])**(y[mesh[i]]-y[mesh[j]]);
            q = sqrt(q);
            if(q<2.0)
            {
                if(type[mesh[i]] == 0)
                {
                    count_temp = atomicAdd(&count,1);
                    pair_i[count_temp] = i;
                    pair_j[count_temp] = j;
                }
                else if(type[mesh[j]] == 0)
                {
                    count_temp = atomicAdd(&count,1);
                    pair_i[count_temp] = j;
                    pair_j[count_temp] = i;
                }
            }
        }
    }
}

void ptc_nnps_mesh(SPH *sph)
//void ptc_nnps_mesh(SPH_PARTICLE *particle,SPH_PAIR *pair,unsigned int ***mesh)
{
    SPH_PARTICLE *particle;
    SPH_PAIR *pair;
    SPH_MESH mesh;
    particle = sph->particle;
    pair = sph->pair;
    mesh = sph->mesh;

    pair->total = 0;

    omp_lock_t lock;
    omp_init_lock(&lock);

    for(int j=0;j<MESH_LENGTH_NUM;j++)
    {
        for(int i=0;i<MESH_DEEPTH_NUM;i++)
        {
            for(unsigned int k=0;k<mesh[i][j][MESH_PTC_NUM-1];k++)
            {
                //mesh[i][j]-->mesh[i][j]
                for(unsigned int m=k+1;m<mesh[i][j][MESH_PTC_NUM-1];m++)
                {
                    if(PTC_DISTANCE(mesh[i][j][k],mesh[i][j][m])<=PTC_REGION_RADIUS)
                    {
                        if(particle->type[mesh[i][j][k]]==0)
			            {
                            pair->i[pair->total] = mesh[i][j][k];
                            pair->j[pair->total] = mesh[i][j][m];
                            pair->total++;
                        }
                        else if (particle->type[mesh[i][j][m]]==0)
                        {
                            pair->i[pair->total] = mesh[i][j][m];
                            pair->j[pair->total] = mesh[i][j][k];
                            pair->total++;
                        }
                    }
                }
                //mesh[i][j]-->mesh[i][j+1]
                if(j<(MESH_LENGTH_NUM-1))
                {
                    for(unsigned int m=0;m<mesh[i][j+1][MESH_PTC_NUM-1];m++)
                    {
                        if(PTC_DISTANCE(mesh[i][j][k],mesh[i][j+1][m])<=PTC_REGION_RADIUS)
                        {
                            if(particle->type[mesh[i][j][k]]==0)
			                {
                                pair->i[pair->total] = mesh[i][j][k];
                                pair->j[pair->total] = mesh[i][j+1][m];
                                pair->total++;
                            }
                            else if (particle->type[mesh[i][j+1][m]]==0)
                            {
                                pair->i[pair->total] = mesh[i][j+1][m];
                                pair->j[pair->total] = mesh[i][j][k];
                                pair->total++;
                            }
                        }
                    }
                }
                //mesh[i][j]-->mesh[i+1][j]
                if(i<(MESH_DEEPTH_NUM-1))
                {
                    for(unsigned int m=0;m<mesh[i+1][j][MESH_PTC_NUM-1];m++)
                    {
                        if(PTC_DISTANCE(mesh[i][j][k],mesh[i+1][j][m])<=PTC_REGION_RADIUS)
                        {
                            if(particle->type[mesh[i][j][k]]==0)
			                {
                                pair->i[pair->total] = mesh[i][j][k];
                                pair->j[pair->total] = mesh[i+1][j][m];
                                pair->total++;
                            }
                            else if (particle->type[mesh[i+1][j][m]]==0)
                            {
                                pair->i[pair->total] = mesh[i+1][j][m];
                                pair->j[pair->total] = mesh[i][j][k];
                                pair->total++;
                            }
                        }
                    }
                }
                //mesh[i][j]-->mesh[i+1][j+1]
                if(i < (MESH_DEEPTH_NUM-1) && j < (MESH_LENGTH_NUM-1))
                {
                    for(unsigned int m=0;m<mesh[i+1][j+1][MESH_PTC_NUM-1];m++)
                    {
                        if(PTC_DISTANCE(mesh[i][j][k],mesh[i+1][j+1][m])<=PTC_REGION_RADIUS)
                        {
                            if(particle->type[mesh[i][j][k]]==0)
			                {
                                pair->i[pair->total] = mesh[i][j][k];
                                pair->j[pair->total] = mesh[i+1][j+1][m];
                                pair->total++;
                            }
                            else if (particle->type[mesh[i+1][j+1][m]]==0)
                            {
                                pair->i[pair->total] = mesh[i+1][j+1][m];
                                pair->j[pair->total] = mesh[i][j][k];
                                pair->total++;
                            }
                        }
                    }

                }
                //mesh[i][j]-->mesh[i-1][j+1]
                if(i > 0 && j<(MESH_LENGTH_NUM-1))
                {
                    for(unsigned int m=0;m<mesh[i-1][j+1][MESH_PTC_NUM-1];m++)
                    {
                        if(PTC_DISTANCE(mesh[i][j][k],mesh[i-1][j+1][m])<=PTC_REGION_RADIUS)
                        {
                            if(particle->type[mesh[i][j][k]]==0)
			                {
                                pair->i[pair->total] = mesh[i][j][k];
                                pair->j[pair->total] = mesh[i-1][j+1][m];
                                pair->total++;
                            }
                            else if (particle->type[mesh[i-1][j+1][m]]==0)
                            {
                                pair->i[pair->total] = mesh[i-1][j+1][m];
                                pair->j[pair->total] = mesh[i][j][k];
                                pair->total++;
                            }
                        }
                    }
                }
            }
        }
    }
}
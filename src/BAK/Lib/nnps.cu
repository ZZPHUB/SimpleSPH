#include "Lib.cuh"
using namespace std;

void ptc_nnps_direct(SPH *sph)
//void ptc_nnps_direct(SPH_PARTICLE *particle,SPH_PAIR *pair)
{   
    SPH_PARTICLE *particle;
    SPH_PAIR *pair;
    particle = sph->particle;
    pair = sph->pair;

    pair->total = 0;
    for(int i=0;i<particle->total;i++)
    {   
        for(int j=i+1;j<particle->total;j++)
        {
        /*  if the distance between the particles i and j is less or equ to PTC_RADIUS,then they are a pair*/
            if(PTC_DISTANCE(i,j) <= PTC_REGION_RADIUS) 
            {
                if(particle->type[i] == 0)
                {
                    pair->i[pair->total] = i;
                    pair->j[pair->total] = j;
                    pair->total = pair->total+1;
                }
                else if (particle->type[j] == 0)
                {
                    pair->i[pair->total] = j;
                    pair->j[pair->total] = i;
                    pair->total = pair->total+1;
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

void ptc_nnps_check(SPH_PAIR *pair,SPH_PAIR *pair_direct,unsigned int *total)
{
    *total = 0;
    if(pair->total == pair_direct->total)
    {
        for(int i=0;i<pair->total;i++)
        {   
            for(int j=0;j<pair_direct->total;j++)
            {
                if(pair->i[i]==pair_direct->i[j] && pair->j[i] == pair_direct->j[j])
                {
                    (*total)++;
                }
                else if(pair->i[i] == pair_direct->j[j] && pair->j[i] == pair_direct->i[j])
                {
                    (*total)++;
                }
            }
        }
        }
}
__global__ void sph_nnps_cuda(SPH_CUDA *cuda,SPH_ARG *arg,SPH_RIGID *rigid)
{

    //n---> (x)length direction
    //m ---> (y)depth direction

    double dx = 0.0;
    double dy = 0.0;
    double q = 0.0;
    int index_i = 0;
    int index_j = 0;
    //int arg->pair_num=0;
    int mesh_id = 0;
    //if( blockIdx.x >= arg->mesh_ynum) return;
    //if( threadIdx.x >= arg->mesh_xnum) return;
    //const int mesh_id = threadIdx.x + blockIdx.x * blockDim.x;
    for(int m =0;m<arg->mesh_ynum;m++)
    {
        for(int n=0;n<arg->mesh_xnum;n++)
        {
            mesh_id = n+m*arg->mesh_xnum;   
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
                            //arg->pair_num = atomicAdd(&(arg->pair_num),1);
                            cuda->pair_i[arg->pair_num] = cuda->mesh[index_i];
                            cuda->pair_j[arg->pair_num] = cuda->mesh[index_j];
                            arg->pair_num +=1;
                        }
                        else if(cuda->type[cuda->mesh[index_j]] == 0)
                        {
                            //arg->pair_num = atomicAdd(&(arg->pair_num),1);
                            cuda->pair_i[arg->pair_num] = cuda->mesh[index_j];
                            cuda->pair_j[arg->pair_num] = cuda->mesh[index_i];
                            arg->pair_num +=1;
                        }
                    }
                }
                //(x,y)->(x+1,y)
                if( n<(arg->mesh_xnum-1))
                {
                    for(int j=0;j<cuda->mesh_count[mesh_id+1];j++)
                    {
                        index_j = mesh_id + 1 + j*arg->mesh_num;
                        dx = cuda->x[cuda->mesh[index_i]] - cuda->x[cuda->mesh[index_j]];
                        dy = cuda->y[cuda->mesh[index_i]] - cuda->y[cuda->mesh[index_j]];
                        q = sqrt(dx*dx+dy*dy)/arg->h;
                        if(q<2.0)
                        {
                            if(cuda->type[cuda->mesh[index_i]] == 0)
                            {
                                //arg->pair_num = atomicAdd(&(arg->pair_num),1);
                                cuda->pair_i[arg->pair_num] = cuda->mesh[index_i];
                                cuda->pair_j[arg->pair_num] = cuda->mesh[index_j];
                                arg->pair_num +=1;
                            }
                            else if(cuda->type[cuda->mesh[index_j]] == 0)
                            {
                                //arg->pair_num = atomicAdd(&(arg->pair_num),1);
                                cuda->pair_i[arg->pair_num] = cuda->mesh[index_j];
                                cuda->pair_j[arg->pair_num] = cuda->mesh[index_i];
                                arg->pair_num +=1;
                            }
                        }
                    }   
                }

                //(x,y)->(x,y+1)
                if( m<(arg->mesh_ynum-1))
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
                                //arg->pair_num = atomicAdd(&(arg->pair_num),1);
                                cuda->pair_i[arg->pair_num] = cuda->mesh[index_i];
                                cuda->pair_j[arg->pair_num] = cuda->mesh[index_j];
                                arg->pair_num +=1;
                            }
                            else if(cuda->type[cuda->mesh[index_j]] == 0)
                            {
                                //arg->pair_num = atomicAdd(&(arg->pair_num),1);
                                cuda->pair_i[arg->pair_num] = cuda->mesh[index_j];
                                cuda->pair_j[arg->pair_num] = cuda->mesh[index_i];
                                arg->pair_num +=1;
                            }
                        }
                    }
                }

                //(x,y)->(x+1,y+1)
                if( n<(arg->mesh_xnum-1) && m<(arg->mesh_ynum-1))
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
                                //arg->pair_num = atomicAdd(&(arg->pair_num),1);
                                cuda->pair_i[arg->pair_num] = cuda->mesh[index_i];
                                cuda->pair_j[arg->pair_num] = cuda->mesh[index_j];
                                arg->pair_num +=1;
                            }
                            else if(cuda->type[cuda->mesh[index_j]] == 0)
                            {
                                //arg->pair_num = atomicAdd(&(arg->pair_num),1);
                                cuda->pair_i[arg->pair_num] = cuda->mesh[index_j];
                                cuda->pair_j[arg->pair_num] = cuda->mesh[index_i];
                                arg->pair_num +=1;
                            }
                        }
                    }
                }

                //(x,y)->(x-1,y+1)
                if( n>0 && m<(arg->mesh_ynum-1))
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
                                //arg->pair_num = atomicAdd(&(arg->pair_num),1);
                                cuda->pair_i[arg->pair_num] = cuda->mesh[index_i];
                                cuda->pair_j[arg->pair_num] = cuda->mesh[index_j];
                                arg->pair_num +=1;
                            }
                            else if(cuda->type[cuda->mesh[index_j]] == 0)
                            {
                                //arg->pair_num = atomicAdd(&(arg->pair_num),1);
                                cuda->pair_i[arg->pair_num] = cuda->mesh[index_j];
                                cuda->pair_j[arg->pair_num] = cuda->mesh[index_i];
                                arg->pair_num +=1;
                            }
                        }
                    }
                }
            }
            cuda->mesh_count[mesh_id] = 0;
        }
    }
}

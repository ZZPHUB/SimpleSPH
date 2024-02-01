#include "Lib.H"
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
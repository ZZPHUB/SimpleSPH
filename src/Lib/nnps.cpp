#include "Lib.H"
using namespace std;

void nnps_direct(SPH_PARTICLE *particle,SPH_PAIR *pair)
{   cout << "here in nnps" << endl;
    pair->total = 0;
    #pragma omp parallel num_threads(7)
    {
        #pragma omp for schedule(dynamic,PTC_TOL_NUM/35)
        for(int i=0;i<PTC_TOL_NUM;i++)
        {   
            //if(omp_get_thread_num()==6) cout << "i is " << i << "thid is " << omp_get_thread_num() << endl;
            for(int j=i+1;j<PTC_TOL_NUM;j++)
            {
            /*  if the distance between the particles i and j is less or equ to PTC_RADIUS,then they are a pair*/
                if(PTC_DISTANCE(i,j) <= PTC_REGION_RADIUS) 
                {
                    pair->i[pair->total] = i;
                    pair->j[pair->total] = j;
                    pair->total = pair->total+1;

                    //if(omp_get_thread_num()==6) cout << "total is " << pair->total << "thid is " << omp_get_thread_num() << endl;

                }
            }
        }
    }
}

void nnps_mesh(SPH_PARTICLE *particle,SPH_PAIR *pair,unsigned int ***mesh)
{
    //#pragma omp parallel for num_threads(6)
    for(int i=0;i<MESH_DEEPTH_NUM;i++)
    {
        for(int j=0;j<MESH_LENGTH_NUM;j++)
        {
            for(unsigned int k=0;k<mesh[i][j][MESH_PTC_NUM-1];k++)
            {
                //mesh[i][j]-->mesh[i][j]
                for(unsigned int m=0;m<mesh[i][j][MESH_PTC_NUM-1];m++)
                {
                    if(PTC_DISTANCE(mesh[i][j][k],mesh[i][j][m])<=PTC_REGION_RADIUS)
                    {
                        pair->i[pair->total] = mesh[i][j][k];
                        pair->j[pair->total] = mesh[i][j][m];
                        pair->total++;
                    }
                }
                //mesh[i][j]-->mesh[i][j+1]
                if(j<(MESH_LENGTH_NUM-1))
                {
                    for(unsigned int m=0;m<mesh[i][j+1][MESH_PTC_NUM-1];m++)
                    {
                        if(PTC_DISTANCE(mesh[i][j][k],mesh[i][j+1][m])<=PTC_REGION_RADIUS)
                        {
                            pair->i[pair->total] = mesh[i][j][k];
                            pair->j[pair->total] = mesh[i][j+1][m];
                            pair->total++;
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
                            pair->i[pair->total] = mesh[i][j][k];
                            pair->j[pair->total] = mesh[i+1][j][m];
                            pair->total++;
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
                            pair->i[pair->total] = mesh[i][j][k];
                            pair->j[pair->total] = mesh[i+1][j+1][m];
                            pair->total++;
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
                            pair->i[pair->total] = mesh[i][j][k];
                            pair->j[pair->total] = mesh[i-1][j+1][m];
                            pair->total++;
                        }
                    }
                }
            }
        }
    }
}
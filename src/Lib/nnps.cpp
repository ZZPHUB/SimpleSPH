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

                    if(omp_get_thread_num()==6) cout << "total is " << pair->total << "thid is " << omp_get_thread_num() << endl;

                }
            }
        }
    }
}
#include "SPH.H"
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <time.h>
using namespace std;


int main(void)
{ 
    SPH *sph = sph_init();    
    
    char filename[] = "../data/postprocess/sph000.vtk"; //filename 

    while (true)
    {
        for(sph->current_step;sph->current_step<sph->total_step;sph->current_step++)
        {
            if(sph->current_step%PRINT_TIME_STEP == 0)
            {
                filename[23] = (sph->current_step/PRINT_TIME_STEP)/100 + 48;
                filename[24] = ((sph->current_step/PRINT_TIME_STEP)%100)/10 + 48;
                filename[25] = ((sph->current_step/PRINT_TIME_STEP)%10) + 48;
                ptc_vtk_direct(sph,filename);
            }
            //calculate and integration
            ptc_time_integral(sph); 
            ptc_info(sph);
        }
        system("clear");
        cout << "press 0 to kill precess or num(>100) for more steps" << endl;
        cin >> sph->total_step;
        if(sph->total_step == 0) break;
    }

    sph_free(sph);
    return 0;
}

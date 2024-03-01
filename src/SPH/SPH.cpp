#include "SPH.H"
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <time.h>
using namespace std;


int main(void)
{ 
    SPH_PARTICLE particle;
    SPH_KERNEL kernel;
    SPH_PAIR pair;
    //SPH_RIGID wall;
    SPH_MESH mesh = NULL;
    SPH sph;
    sph.particle = &particle;
    sph.kernel = &kernel;
    sph.pair = &pair;
    //sph.rigid_0 = &wall;
    sph.mesh = mesh;

    sph_init(&sph);    
    
    char filename[] = "../data/postprocess/vtk/sph000.vtk"; //filename 

    for(sph.current_step;sph.current_step<sph.total_step;sph.current_step++)
    {
        if(sph.current_step%PRINT_TIME_STEP == 0)
        {
            filename[27] = (sph.current_step/PRINT_TIME_STEP)/100 + 48;
            filename[28] = ((sph.current_step/PRINT_TIME_STEP)%100)/10 + 48;
            filename[29] = ((sph.current_step/PRINT_TIME_STEP)%10) + 48;
            sph_save_single(&sph,filename);
        }
        //calculate and integration
        sph_time_integral(&sph); 
        ptc_info(&sph);
    }
    sph_save_last(&sph);
    sph_free(&sph);
    return 0;
}

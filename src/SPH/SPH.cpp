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
    SPH_RIGID wedge;
    SPH_MESH mesh = NULL;
    SPH sph;
    sph.particle = &particle;
    sph.kernel = &kernel;
    sph.pair = &pair;
    sph.rigid = &wedge;
    sph.mesh = mesh;

    sph_init(&sph);    
    
    for(sph.current_step;sph.current_step<sph.total_step;sph.current_step++)
    {
        if(sph.current_step%PRINT_TIME_STEP == 0)
        {
            sph_save_single(&sph);
        }
        //calculate and integration
        sph_time_integral(&sph); 
        ptc_info(&sph);
    }
    sph_save_last(&sph);
    sph_free(&sph);
    return 0;
}

#include "SPH.H"
#include <fstream>
#include <iomanip>
using namespace std;

int main(void)
{ 
    //get_max_men();
    SPH_PARTICLE particle;
    SPH_KERNEL kernel;
    SPH_PARE pare;


    /************stack is too small***************
    double x[PTC_TOL_NUM]={0};
    double y[PTC_TOL_NUM]={0};
    double vx[PTC_TOL_NUM]={0};
    double vy[PTC_TOL_NUM]={0};
    double pressure[PTC_TOL_NUM]={0};
    double density[PTC_TOL_NUM]={0};
    char type[PTC_TOL_NUM]={0};
    ***********************************************/

   

    particle.x = (double *)(malloc(sizeof(double)*PTC_TOL_NUM));
    particle.y = (double *)(malloc(sizeof(double)*PTC_TOL_NUM));
    particle.vx = (double *)(malloc(sizeof(double)*PTC_TOL_NUM));
    particle.vy = (double *)(malloc(sizeof(double)*PTC_TOL_NUM));
    particle.pressure = (double *)(malloc(sizeof(double)*PTC_TOL_NUM));
    particle.type = (char *)(malloc(sizeof(char)*PTC_TOL_NUM));  
    ptc_generate(particle.x,particle.y);

    ofstream writefile;
    writefile.open("../data/init.vtk");

    writefile << "# vtk DataFile Version 3.0" << endl;
    writefile << "init data" << endl;
    writefile << "ASCII" << endl;
    writefile << "DATASET UNSTRUCTURED_GRID" << endl;
    writefile << "POINTS " << FLUID_PTC_NUM<< " " << "double" << endl;

    for(int i=0;i<FLUID_PTC_NUM;i++)
    {
        writefile << setiosflags(ios::scientific) << particle.x[i]<< " " << particle.y[i]<< " " << 0.0 << endl;
    }
    
    writefile.close();
    return 0;
}



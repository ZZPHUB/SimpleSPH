#include "SPH.H"
#include <fstream>
using namespace std;

int main(void)
{ 
    SPH_PARTICLE particle;
    SPH_KERNEL kernel;
    SPH_PARE pare;

    double x[PTC_TOL_NUM]={0};
    double y[PTC_TOL_NUM]={0};
    double vx[PTC_TOL_NUM]={0};
    double vy[PTC_TOL_NUM]={0};
    double pressure[PTC_TOL_NUM]={0};
    double density[PTC_TOL_NUM]={0};
    char type[PTC_TOL_NUM]={0};

    particle.x = x;
    particle.y = y;
    particle.vx = vx;
    particle.vy = vy;
    particle.pressure = pressure;
    particle.type = type;  
    ptc_generate(particle.x,particle.y);

    ofstream writefile;
    writefile.open("../data/init.vtk");

    writefile << "vtk DataFile Version 3.0" << endl;
    writefile << "init data" << endl;
    writefile << "ASCII" << endl;
    writefile << "DATASET UNSTRUCTURED_GRID" << endl;
    writefile << "POINTS " << FLUID_PTC_NUM << "double" << endl;

    for(int i=0;i<FLUID_PTC_NUM;i++)
    {
        writefile << particle.x[i]<< " " << particle.y[i] << 0.0 << endl;
    }
    
    writefile.close();
    return 0;
}



#include "SPH.H"
#include <fstream>
#include <iomanip>
#include <stdlib.h>
using namespace std;

int main(void)
{ 
    SPH_PARTICLE particle;
    SPH_KERNEL kernel;
    SPH_PAIR pair;

    /************stack is too small,so init data in heap***************/
    particle.x = (double *)(calloc(PTC_TOL_NUM,sizeof(double)));
    particle.y = (double *)(calloc(PTC_TOL_NUM,sizeof(double)));
    particle.vx = (double *)(calloc(PTC_TOL_NUM,sizeof(double)));
    particle.vy = (double *)(calloc(PTC_TOL_NUM,sizeof(double)));
    particle.accx = (double *)(calloc(PTC_TOL_NUM,sizeof(double))); 
    particle.accy = (double *)(calloc(PTC_TOL_NUM,sizeof(double)));
    particle.density = (double *)(calloc(PTC_TOL_NUM,sizeof(double)));
    particle.mass = (double *)(calloc(PTC_TOL_NUM,sizeof(double)));
    particle.dif_density = (double *)(calloc(PTC_TOL_NUM,sizeof(double)));
    particle.pressure = (double *)(calloc(PTC_TOL_NUM,sizeof(double)));
    particle.type = (char *)(calloc(PTC_TOL_NUM,sizeof(char)));  

    kernel.w = (double *)(calloc(4*PTC_TOL_NUM,sizeof(double)));
    kernel.dwdx = (double *)(calloc(4*PTC_TOL_NUM,sizeof(double)));
    kernel.dwdy = (double *)(calloc(4*PTC_TOL_NUM,sizeof(double)));
   
    pair.total = 0; 
    pair.i = (int *)(calloc(4*PTC_TOL_NUM,sizeof(int)));
    pair.j = (int *)(calloc(4*PTC_TOL_NUM,sizeof(int)));

    ptc_generate(&particle);//generate the fluid solid and dummy particles
    cout << particle.x[3] << " " << particle.y[3] << endl;
    nnps_direct(&particle,&pair);
    cout << pair.total << endl;
/*
    for(int i=0;i<pair.total;i=i+1){
    //cout << sqrt(pow(particle.x[pair.i[i]]-particle.x[pair.j[i]],2)+pow(particle.y[pair.i[i]]-particle.y[pair.j[i]],2)) << endl;
    cout << pair.i[i] << " " << pair.j[i] << endl;
    }
    //cout << PTC_DISTANCE(pair.i[2],pair.j[2]) << endl;
    ptc_kernel(&particle,&pair,&kernel);
    cout << "total interact particle num is " << pair.total << endl;

    ptc_init(&particle); //particles init values

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
    */

    free((void *)particle.x);
    free((void *)particle.y);
    free((void *)particle.vx);
    free((void *)particle.vy);
    free((void *)particle.accx);
    free((void *)particle.accy);
    free((void *)particle.density);
    free((void *)particle.dif_density);
    free((void *)particle.mass);
    free((void *)particle.pressure);
    free((void *)particle.type);

    free((void *)kernel.w);
    free((void *)kernel.dwdx);
    free((void *)kernel.dwdy);
    
    free((void *)pair.i);
    free((void *)pair.j);
    return 0;
}



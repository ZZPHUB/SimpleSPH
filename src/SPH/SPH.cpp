#include "SPH.H"



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

    return 0;
}



#include "Lib.H"

void ptc_info(SPH_PARTICLE *particle,SPH_PAIR *pair,RIGID *wedge,unsigned int time)
{
    system("clear");
    cout << "***************************SPH*************************" << endl;
    cout << "Total Particles Num: " << particle->total << endl;
    cout << "------------------------------------------------------"
    cout << "Fluid Particles Num: " << FLUID_PTC_NUM << endl;
    cout << "------------------------------------------------------"
    cout << "Virtual Particles Num: " << VIRTUAL_PTC_NUM << endl;
    cout << "------------------------------------------------------"
    cout << "Solid Particles Num: " << particle->total-FLUID_PTC_NUM-VIRTUAL_PTC_NUM << endl;
    cout << "------------------------------------------------------"
    cout << "Total Particles Pair Num: " << pair->total << endl;
    cout << "------------------------------------------------------"
    cout << "Current Time Step: " << time << endl;
    cout << "------------------------------------------------------"
    cout << "Rigid Body Velocity in X-direction: " << wedge->vx << endl;
    cout << "------------------------------------------------------"
    cout << "Rigid Body Velocity in Y-direction: " << wedge->vy << endl;
    cout << "------------------------------------------------------"
    cout << "Rigid Body Acceleration in X-direction: " << wedge->accx << endl;
    cout << "------------------------------------------------------"
    cout << "Rigid Body Acceleration in Y-direction: " << wedge->accy << endl;
    cout << "------------------------------------------------------"
    cout << "Rigid Body Angular Velocity: " << wedge->omega << endl;
    cout << "------------------------------------------------------"
    cout << "Rigid Body Angular Acceleration: "<< wedge->alpha << endl;
    cout << "******************************************************"

}
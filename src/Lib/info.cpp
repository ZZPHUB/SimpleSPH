#include "Lib.H"

void ptc_info(SPH_PARTICLE *particle,SPH_PAIR *pair,RIGID *wedge,unsigned int time)
{
    system("clear");
    cout << "***************************SPH*************************" << endl;
    cout << "Total Particles Num: " << particle->total << endl;
    cout << "------------------------------------------------------" << endl;
    cout << "Fluid Particles Num: " << FLUID_PTC_NUM << endl;
    cout << "------------------------------------------------------" << endl;
    cout << "Virtual Particles Num: " << VIRTUAL_PTC_NUM << endl;
    cout << "------------------------------------------------------" << endl;
    cout << "Solid Particles Num: " << particle->total-FLUID_PTC_NUM-VIRTUAL_PTC_NUM << endl;
    cout << "------------------------------------------------------" << endl;
    cout << "Total Particles Pair Num: " << pair->total << endl;
    cout << "------------------------------------------------------" << endl;
    cout << "Current Time Step: " << time << endl;
    cout << "------------------------------------------------------" << endl;
    cout << "Rigid Body Velocity in X-direction: " << wedge->vx << endl;
    cout << "------------------------------------------------------" << endl;
    cout << "Rigid Body Velocity in Y-direction: " << wedge->vy << endl;
    cout << "------------------------------------------------------" << endl;
    cout << "Rigid Body Acceleration in X-direction: " << wedge->accx << endl;
    cout << "------------------------------------------------------" << endl;
    cout << "Rigid Body Acceleration in Y-direction: " << wedge->accy << endl;
    cout << "------------------------------------------------------" << endl;
    cout << "Rigid Body Angular Velocity: " << wedge->omega << endl;
    cout << "------------------------------------------------------" << endl;
    cout << "Rigid Body Angular Acceleration: "<< wedge->alpha << endl;
    cout << "******************************************************" << endl;

}
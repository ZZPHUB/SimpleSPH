#include "Lib.H"

void ptc_info(SPH *sph)
//void ptc_info(SPH_PARTICLE *particle,SPH_PAIR *pair,RIGID *wedge,unsigned int time)
{
    SPH_PARTICLE *particle;
    SPH_PAIR *pair;
    SPH_RIGID *wedge;
    particle = sph->particle;
    pair = sph->pair;
    wedge = sph->rigid;

    system("clear");
    cout << "***************************SPH*************************" << endl;
    cout << "Total Particles Num: " << particle->total << endl;
    cout << "------------------------------------------------------" << endl;
    cout << "Fluid Particles Num: " << particle->fulid_ptc_num << endl;
    cout << "------------------------------------------------------" << endl;
    cout << "Wall Particles Num: " << particle->wall_ptc_num << endl;
    cout << "------------------------------------------------------" << endl;
    cout << "Solid Particles Num: " << particle->rigid_ptc_num << endl;
    cout << "------------------------------------------------------" << endl;
    cout << "Total Particles Pair Num: " << pair->total << endl;
    cout << "------------------------------------------------------" << endl;
    cout << "Current Time Step: " << sph->current_step << endl;
    //cout << "------------------------------------------------------" << endl;
    //cout << "One Step Use Time: " << (double)(sph->end_tiem-sph->start_time)/(double)sph->current_step << endl;
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
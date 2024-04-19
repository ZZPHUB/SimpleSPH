#include "IO.cuh"
using namespace std;
void sph_write_csv(SPH *sph)
{
    SPH_RIGID *rigid;
    rigid = sph->host_rigid;
    SPH_ARG *arg;
    arg = sph->host_arg;

    string filename = arg->case_dir;
    filename += "/csv/rigid.csv";

    ofstream rigid_info;
    rigid_info.open(filename.c_str(),ios::app);

    rigid_info << arg->init_step*arg->dt << "," << rigid->cogx << "," << rigid->cogy << "," << rigid->vx << "," <<
    rigid->vy << "," << rigid->omega << "," << rigid->accx << "," << rigid->accy << "," << rigid->alpha << endl;

}
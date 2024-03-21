#include "PostProcess.cuh"
using namespace std;

void sph_save_rigid(SPH *sph)
{
    SPH_RIGID *wedge;
    wedge = sph->rigid;

    static vector<double> rigid_vy;
    static vector<double> rigid_accy;

    if(sph->init_impac_flag == 0)
    {
        if(sph->current_step%PRINT_TIME_STEP == 0)
        {
            rigid_vy.push_back(wedge->vy);
            rigid_accy.push_back(wedge->accy);
        }
        if(sph->current_step==sph->total_step-1)
        {
            string filename = "../data/postprocess/csv/rigid.csv";
            ofstream rigidfile;
            rigidfile.open(filename.c_str());
            rigidfile << "time, vy, accy" << endl;
            for(unsigned int i=0;i<rigid_vy.size();i++)
            {
                rigidfile << setiosflags(ios::scientific)<< sph->d_time*((double)sph->current_step) \
                << ", " << rigid_vy.at(i) << ", " << rigid_accy.at(i) << endl;
            }
            rigidfile.close();
        }
    }
}
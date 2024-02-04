#include "PostProcess.H"
using namespace std;

void ptc_pressure_conv(SPH *sph)
{
    SPH_PARTICLE *particle;
    particle = sph->particle;
    char filename[] = "../data/postprocess/csv/p_conv.csv";

    ofstream writefile;
    writefile.open(filename,ios::out|ios::app);

    writefile << setiosflags(ios::scientific) <<    \
    ((double)sph->current_step)*sph->d_time <<", "  \
    <<particle->pressure[1994] << ", " <<9800*particle->y[2226]<< endl;

    writefile.close();
}
#include "PostProcess.H"
using namespace std;


void ptc_vtk_direct(SPH *sph,char *filename)
{
    SPH_PARTICLE *particle;
    particle = sph->particle;

    unsigned int ptc_num = 0;

    ptc_num = particle->total;

    ofstream writefile;
    writefile.open(filename);

    writefile << "# vtk DataFile Version 3.0" << endl;
    writefile << "sph data" << endl;
    writefile << "ASCII" << endl;
    writefile << "DATASET UNSTRUCTURED_GRID" << endl;
    writefile << "POINTS " << ptc_num << " " << "double" << endl;

    for(unsigned int i=0;i<particle->total;i++)
    {
        writefile << setiosflags(ios::scientific) << particle->x[i] << " " \
        << particle->y[i] << " " << 0.0 << endl;
    }

    writefile << "POINT_DATA" << " " << ptc_num << endl;

    //density
    if(PARA&0x01)
    {
        writefile << "SCALARS "<< "density double 1" << endl;
        writefile << "LOOKUP_TABLE DEFAULT" << endl;
        for(unsigned int i=0;i<particle->total;i++)
        {
            writefile << setiosflags(ios::scientific) << particle->density[i] << endl;
        }
    }
    //pressure
    if(PARA&0x02)
    {
        writefile << "SCALARS "<< "pressure double 1" << endl;
        writefile << "LOOKUP_TABLE DEFAULT" << endl;
        for(unsigned int i=0;i<particle->total;i++)
        {
            writefile << setiosflags(ios::scientific) << particle->pressure[i] << endl;
        }
    }
    //velocity
    if(PARA&0x04)
    {
        writefile << "VECTORS "<< "velocity double" << endl;
        for(unsigned int i=0;i<particle->total;i++)
        {
            writefile << setiosflags(ios::scientific) << particle->vx[i] <<" " << particle->vy[i] << " " \
            << 0.0 << endl;

        }
    }
    //acceleration
    if(PARA&0x08)
    {
        writefile << "VECTORS "<< "acceleration double" << endl;
        for(unsigned int i=0;i<particle->total;i++)
        {
            writefile << setiosflags(ios::scientific) << particle->accx[i] <<" " << particle->accy[i] << " " \
            << 0.0 << endl;
        }
    }

    writefile.close();

}
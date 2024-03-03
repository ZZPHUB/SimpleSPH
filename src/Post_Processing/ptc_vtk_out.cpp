#include "PostProcess.H"
using namespace std;


void sph_save_single(SPH *sph)
{
    SPH_PARTICLE *particle;
    particle = sph->particle;

    unsigned int ptc_num = 0;
    ptc_num = particle->total;


    string filename = "../data/postprocess/vtk/sph"; 
    filename += to_string(sph->current_step);
    filename += ".vtk";

    ofstream vtkfile;
    vtkfile.open(filename.c_str());

    vtkfile << "# vtk DataFile Version 3.0" << endl;
    vtkfile << "sph data" << endl;
    vtkfile << "ASCII" << endl;
    vtkfile << "DATASET UNSTRUCTURED_GRID" << endl;
    vtkfile << "POINTS " << ptc_num << " " << "double" << endl;

    for(unsigned int i=0;i<particle->total;i++)
    {
        vtkfile << setiosflags(ios::scientific) << particle->x[i] << " " \
        << particle->y[i] << " " << 0.0 << endl;
    }

    vtkfile << "POINT_DATA" << " " << ptc_num << endl;

    //density
    if(PARA&0x01)
    {
        vtkfile << "SCALARS "<< "density double 1" << endl;
        vtkfile << "LOOKUP_TABLE DEFAULT" << endl;
        for(unsigned int i=0;i<particle->total;i++)
        {
            vtkfile << setiosflags(ios::scientific) << particle->density[i] << endl;
        }
    }
    //pressure
    if(PARA&0x02)
    {
        vtkfile << "SCALARS "<< "pressure double 1" << endl;
        vtkfile << "LOOKUP_TABLE DEFAULT" << endl;
        for(unsigned int i=0;i<particle->total;i++)
        {
            vtkfile << setiosflags(ios::scientific) << particle->pressure[i] << endl;
        }
    }
    //velocity
    if(PARA&0x04)
    {
        vtkfile << "VECTORS "<< "velocity double" << endl;
        for(unsigned int i=0;i<particle->total;i++)
        {
            vtkfile << setiosflags(ios::scientific) << particle->vx[i] <<" " << particle->vy[i] << " " \
            << 0.0 << endl;

        }
    }
    //acceleration
    if(PARA&0x08)
    {
        vtkfile << "VECTORS "<< "acceleration double" << endl;
        for(unsigned int i=0;i<particle->total;i++)
        {
            vtkfile << setiosflags(ios::scientific) << particle->accx[i] <<" " << particle->accy[i] << " " \
            << 0.0 << endl;
        }
    }
    vtkfile.close();

}


void sph_save_last(SPH *sph)
{
    SPH_PARTICLE *particle;
    SPH_RIGID *wedge;
    particle = sph->particle;
    wedge = sph->rigid;

    unsigned int ptc_num = 0;

    ptc_num = particle->total;

    ofstream vtkfile;
    vtkfile.open("../data/postprocess/save.vtk");

    vtkfile << "# vtk DataFile Version 3.0" << endl;
    vtkfile << "sph data" << endl;
    vtkfile << "ASCII" << endl;
    vtkfile << "DATASET UNSTRUCTURED_GRID" << endl;
    vtkfile << "POINTS " << ptc_num << " " << "double" << endl;

    for(unsigned int i=0;i<particle->total;i++)
    {
        vtkfile << setiosflags(ios::scientific) << particle->x[i] << " " \
        << particle->y[i] << " " << 0.0 << endl;
    }

    vtkfile << "POINT_DATA" << " " << ptc_num << endl;

    //ptc type
    vtkfile << "SCALARS "<< "type int 1" << endl;
    vtkfile << "LOOKUP_TABLE DEFAULT" << endl;
    for(unsigned int i=0;i<particle->total;i++)
    {
        vtkfile  << particle->type[i] << endl;
    }
    //density
    vtkfile << "SCALARS "<< "density double 1" << endl;
    vtkfile << "LOOKUP_TABLE DEFAULT" << endl;
    for(unsigned int i=0;i<particle->total;i++)
    {
        vtkfile << setiosflags(ios::scientific) << particle->density[i] << endl;
    }
    //pressure
    vtkfile << "SCALARS "<< "pressure double 1" << endl;
    vtkfile << "LOOKUP_TABLE DEFAULT" << endl;
    for(unsigned int i=0;i<particle->total;i++)
    {
        vtkfile << setiosflags(ios::scientific) << particle->pressure[i] << endl;
    }
    //velocity
    vtkfile << "VECTORS "<< "velocity double" << endl;
    for(unsigned int i=0;i<particle->total;i++)
    {
        vtkfile << setiosflags(ios::scientific) << particle->vx[i] <<" " << particle->vy[i] << " " \
        << 0.0 << endl;

    }
    vtkfile.close();

    ofstream infofile;
    infofile.open("../data/postprocess/info.txt");

    //infofile << "#the wedge's velocity in x-direction" << end;
    infofile << setiosflags(ios::scientific) << wedge->vx << endl;
    //infofile << "#the wedge's velocity in y-direction" << endl;
    infofile << setiosflags(ios::scientific) << wedge->vy << endl;
    //infofile << "#the wedge's omega" << endl;
    infofile << setiosflags(ios::scientific) << wedge->omega << endl;
    //infofile << "#the wedge's center of gravity in x-direction" << endl;
    infofile << setiosflags(ios::scientific) << wedge->cogx << endl;
    //infofile << "#the wedge's center of gravity in y-direction" << endl;
    infofile << setiosflags(ios::scientific) << wedge->cogy << endl;
    //infofile << "#the wedge's moi" << endl;
    infofile << setiosflags(ios::scientific) << wedge->moi << endl;
}
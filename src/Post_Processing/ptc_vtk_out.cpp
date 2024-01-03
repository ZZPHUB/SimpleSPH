#include "PostProcess.H"
using namespace std;

void ptc_vtk_mesh(SPH_PARTICLE *particle,unsigned int ***mesh)
{
    int ptc_num = particle->total;
    
    ofstream writefile;
    writefile.open("../data/init.vtk");

    writefile << "# vtk DataFile Version 3.0" << endl;
    writefile << "init data" << endl;
    writefile << "ASCII" << endl;
    writefile << "DATASET UNSTRUCTURED_GRID" << endl;
    writefile << "POINTS " << ptc_num << " " << "double" << endl;
    for(int i=0;i<MESH_DEEPTH_NUM;i++)
    {
        for(int j=0;j<MESH_LENGTH_NUM;j++)
        {
            for(int m=0;m<mesh[i][j][MESH_PTC_NUM-1];m++)
            {
                writefile << setiosflags(ios::scientific) <<particle->x[mesh[i][j][m]] << " " \
                << particle->y[mesh[i][j][m]] << " " << 0.0 << endl;
            }
        }
    }
    /*
    for(unsigned int i=0;i<ptc_num;i++)
    {
        writefile << setiosflags(ios::scientific) << particle->x[i] << " " \
        << particle->y[i] << " " << 0.0 << endl;
    }
    */

    writefile << "POINT_DATA" << " " << ptc_num << endl;
    writefile << "SCALARS "<< "type double 1" << endl;
    writefile << "LOOKUP_TABLE DEFAULT" << endl;
    for(int i=0;i<MESH_DEEPTH_NUM;i++)
    {
        for(int j=0;j<MESH_LENGTH_NUM;j++)
        {
            for(int m=0;m<mesh[i][j][MESH_PTC_NUM-1];m++)
            {
                writefile << setiosflags(ios::scientific) << (double)(0.0001*particle->type[mesh[i][j][m]]) << endl;
            }
        }
    }
   /*
   for(unsigned int i=0;i<ptc_num;i++)
   {
    writefile << setiosflags(ios::scientific) << (double)(0.0000001*particle->type[i]) << endl;
   }
   */

    writefile.close();    
}

void ptc_vtk_direct(SPH_PARTICLE *particle,double *scale,char *filename)
{
    unsigned int ptc_num = 0;
    //get the number of the particles
    for(unsigned int i=0;i<particle->total;i++)
    {
        if(particle->x[i]>= scale[0] && particle->x[i] <= scale[1] && \
           particle->y[i]>= scale[2] && particle->y[i] <= scale[3])
           {
               ptc_num++;
           } 
    }

    ofstream writefile;
    writefile.open(filename);

    writefile << "# vtk DataFile Version 3.0" << endl;
    writefile << "sph data" << endl;
    writefile << "ASCII" << endl;
    writefile << "DATASET UNSTRUCTURED_GRID" << endl;
    writefile << "POINTS " << ptc_num << " " << "double" << endl;

    for(unsigned int i=0;i<particle->total;i++)
    {
        if(particle->x[i]>= scale[0] && particle->x[i] <= scale[1] && \
           particle->y[i]>= scale[2] && particle->y[i] <= scale[3])
           {
                writefile << setiosflags(ios::scientific) << particle->x[i] << " " \
                << particle->y[i] << " " << 0.0 << endl;
           } 
    }

    writefile << "POINT_DATA" << " " << ptc_num << endl;

    //density
    if(PARA&0x01)
    {
        writefile << "SCALARS "<< "density double 1" << endl;
        writefile << "LOOKUP_TABLE DEFAULT" << endl;
        for(unsigned int i=0;i<particle->total;i++)
        {
            if(particle->x[i]>= scale[0] && particle->x[i] <= scale[1] && \
                particle->y[i]>= scale[2] && particle->y[i] <= scale[3])
                {
                    writefile << setiosflags(ios::scientific) << particle->density[i] << endl;
                } 
        }
    }
    //pressure
    if(PARA&0x02)
    {
        writefile << "SCALARS "<< "pressure double 1" << endl;
        writefile << "LOOKUP_TABLE DEFAULT" << endl;
        for(unsigned int i=0;i<particle->total;i++)
        {
            if(particle->x[i]>= scale[0] && particle->x[i] <= scale[1] && \
                particle->y[i]>= scale[2] && particle->y[i] <= scale[3])
                {
                    writefile << setiosflags(ios::scientific) << particle->pressure[i] << endl;
                } 
        }
    }
    //velocity
    if(PARA&0x04)
    {
        writefile << "VECTORS "<< "velocity double" << endl;
        //writefile << "LOOKUP_TABLE DEFAULT" << endl;
        for(unsigned int i=0;i<particle->total;i++)
        {
            if(particle->x[i]>= scale[0] && particle->x[i] <= scale[1] && \
                particle->y[i]>= scale[2] && particle->y[i] <= scale[3])
                {
                    writefile << setiosflags(ios::scientific) << particle->vx[i] <<" " << particle->vy[i] << " " \
                    << 0.0 << endl;
                } 
        }
    }
    //acceleration
    if(PARA&0x08)
    {
        writefile << "VECTORS "<< "acceleration double" << endl;
        //writefile << "LOOKUP_TABLE DEFAULT" << endl;
        for(unsigned int i=0;i<particle->total;i++)
        {
            if(particle->x[i]>= scale[0] && particle->x[i] <= scale[1] && \
                particle->y[i]>= scale[2] && particle->y[i] <= scale[3])
                {
                    writefile << setiosflags(ios::scientific) << particle->accx[i] <<" " << particle->accy[i] << " " \
                    << 0.0 << endl;
                } 
        }
    }

    writefile.close();

}
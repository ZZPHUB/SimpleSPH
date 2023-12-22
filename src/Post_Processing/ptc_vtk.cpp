#include "PostProcess.H"
using namespace std;

void ptc_vtk_mesh(SPH_PARTICLE *particle,unsigned int ***mesh)
{
    int ptc_num = PTC_TOL_NUM;
    
    ofstream writefile;
    writefile.open("../data/init.vtk");

    writefile << "# vtk DataFile Version 3.0" << endl;
    writefile << "init data" << endl;
    writefile << "ASCII" << endl;
    writefile << "DATASET UNSTRUCTURED_GRID" << endl;
    writefile << "POINTS " << ptc_num << " " << "double" << endl;
    for(int i=0;i<MESH_DEEPTH_NUM;i=i++)
    {
        for(int j=0;j<MESH_LENGTH_NUM;j=j++)
        {
            for(int m=0;m<mesh[i][j][MESH_PTC_NUM-1];m++)
            {
                writefile << setiosflags(ios::scientific) <<particle->x[mesh[i][j][m]] << " " \
                << particle->y[mesh[i][j][m]] << " " << 0.0 << endl;
            }
        }
    }

    writefile << "POINT_DATA" << " " << ptc_num << endl;
    writefile << "SCALARS "<< "name double 1" << endl;
    writefile << "LOOKUP_TABLE DEFAULT" << endl;
    for(int i=0;i<MESH_DEEPTH_NUM;i=i++)
    {
        for(int j=0;j<MESH_LENGTH_NUM;j=j++)
        {
            for(int m=0;m<mesh[i][j][MESH_PTC_NUM-1];m++)
            {
                writefile << setiosflags(ios::scientific) << (double)(0.000001*particle->type[mesh[i][j][m]]) << endl;
            }
        }
    }

    writefile.close();    
}
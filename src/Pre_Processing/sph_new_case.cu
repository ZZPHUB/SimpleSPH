#include "SPH.cuh"
#include <iostream>
#include <iomanip>
#include <string>
#include <assert.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkSmartPointer.h>
#include <vtkType.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPointSet.h>
#include <vtkDataSetReader.h>
#include <vtkPointData.h>
#include <vtkDataArray.h>
#include <vtkDoubleArray.h>

using namespace std;

void get_input(SPH *);
void get_rigid_num(SPH *);
void fluid_ptc_generate(SPH *);
void rigid_ptc_generate(SPH *);
void write_vtk(SPH *);
void rigid_init(SPH *);

int main(int argc,char *argv[])
{
    SPH_ARG arg;
    SPH_RIGID rigid;
    SPH_PARTICLE particle;
    SPH sph;
    sph.host_arg = &arg;
    sph.host_rigid = &rigid;
    sph.particle = &particle;
    //if(argc != 2) printf("\033[0;32;31m Error in %s:%d\033[m\n",__FILE__,__LINE__);
    assert(argc == 2);
    arg.case_dir = argv[1];
    get_input(&sph);
    get_rigid_num(&sph);
    arg.pair_volume = (int)(64*arg.ptc_num/arg.mesh_num);
    particle.x = (double *)calloc(arg.ptc_num,sizeof(double));
    particle.y = (double *)calloc(arg.ptc_num,sizeof(double));
    particle.type = (int *)calloc(arg.ptc_num,sizeof(int));
    fluid_ptc_generate(&sph);
    rigid_ptc_generate(&sph);

    rigid_init(&sph);

    sph_write_info(&sph);
    write_vtk(&sph);

    return 0;
}

void get_input(SPH *sph)
{
    SPH_ARG *arg;
    arg = sph->host_arg;
    SPH_RIGID *rigid;
    rigid = sph->host_rigid;
    //SPH_PARTICLE *particle;
    //particle = sph->particle;

    cout << "fluid length (fluid_x) is:" << endl;
    cin >> arg->fluid_x;
    cout << "fluid depth (fluid_y) is:" << endl;
    cin >> arg->fluid_y;
    cout << "particle spacing (ptc_dx) is:" << endl;
    cin >> arg->ptc_dx;

    arg->h = 2.0*arg->ptc_dx;
    arg->r = 2.0*arg->h;
    arg->wall_layer = 4;
    arg->fluid_xnum = (int)(arg->fluid_x/arg->ptc_dx)+1-2*arg->wall_layer;
    arg->fluid_ynum = (int)(arg->fluid_y/arg->ptc_dx)+1-arg->wall_layer;

    arg->rigid_ptc_num = 0;
    arg->fluid_ptc_num = arg->fluid_xnum*arg->fluid_ynum;
    arg->wall_ptc_num = ((int)(arg->fluid_x/arg->ptc_dx)+1)*((int)(1.1*arg->fluid_y/arg->ptc_dx)+1)- \
                             ((int)(arg->fluid_x/arg->ptc_dx)+1-2*arg->wall_layer)*((int)(1.1*arg->fluid_y/arg->ptc_dx)+1-arg->wall_layer); 
    arg->ptc_num = arg->rigid_ptc_num + arg->fluid_ptc_num + arg->wall_ptc_num;

    arg->mesh_dx = arg->r;
    arg->mesh_x = arg->fluid_x;
    arg->mesh_y = arg->fluid_y * 1.5;
    arg->mesh_xnum = (int)(arg->mesh_x/arg->mesh_dx)+1;
    arg->mesh_ynum = (int)(arg->mesh_y/arg->mesh_dx)+1;
    arg->mesh_num = arg->mesh_xnum*arg->mesh_ynum;
    arg->mesh_volume = 33;

    arg->g = 9.8;
    arg->c = 10.0*sqrt(arg->g*arg->fluid_y);
    arg->ref_rho = 1000.0;
    arg->m = arg->ref_rho*arg->ptc_dx*arg->ptc_dx;
    arg->dt = 0.00002;
    arg->sst = 0.0;
    arg->alpha = 7.0/(4.0*3.14159265358*pow(arg->h,2));

    arg->init_step = 0;
    arg->total_step = 80000;
    arg->print_step = 400;

    arg->new_case_flag =1;
    arg->init_impac_flag = 1;
    arg->save_last_flag = 1;

    rigid->vx = 0.0;
    rigid->vy = 0.0;
    rigid->omega = 0.0;
    rigid->accx = 0.0;
    rigid->accy = 0.0;
    rigid->alpha = 0.0;
    rigid->mass = 12.8;
    rigid->moi = 0.0;
    rigid->offset_x = 0.0;
    rigid->offset_y = 0.0;
    rigid->offset_angl = 0.0;
    rigid->cog_ptc_id = 0.0;
    rigid->cogx = 0.0;
    rigid->cogy = 0.0; 
}

void get_rigid_num(SPH *sph)
{
    SPH_ARG *arg;
    //SPH_PARTICLE *particle;
    //SPH_RIGID *rigid;
    arg = sph->host_arg;
    //particle = sph->particle;
    //rigid = sph->host_rigid;

    std::string filename = arg->case_dir;
    filename += "/wedge.vtk";

    arg->rigid_ptc_num = 0;
    double x[3];

    vtkSmartPointer<vtkUnstructuredGridReader> reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
    reader->SetFileName(filename.c_str());
    reader->Update();
    
    vtkUnstructuredGrid *vtkdata;
    vtkdata = reader->GetOutput();
    for(vtkIdType i=0;i<vtkdata->GetNumberOfPoints();i++)
    {
        vtkdata->GetPoint(i,x);
        if(x[2]==0)
        {
            arg->rigid_ptc_num++;
        }
    }
    arg->ptc_num += arg->rigid_ptc_num;
}

void fluid_ptc_generate(SPH *sph)
{
    SPH_ARG *arg;
    SPH_PARTICLE *particle;
    //SPH_RIGID *rigid;
    arg = sph->host_arg;
    particle = sph->particle;
    //rigid = sph->host_rigid;
    int index = 0;
    for(int x=0;x<arg->fluid_xnum;x++)
    {
        for(int y=0;y<arg->fluid_ynum;y++)
        {
            particle->x[index] = (x+arg->wall_layer)*arg->ptc_dx;
            particle->y[index] = (y+arg->wall_layer)*arg->ptc_dx;
            particle->type[index] = 0;
            index++;
        }
    }
    for(int x=0;x<(arg->fluid_xnum+2*arg->wall_layer);x++)
    {
        for(int y=0;y<((int)(1.1*arg->fluid_y/arg->ptc_dx)+1);y++)
        {
            if(x < arg->wall_layer || x > (arg->fluid_xnum+arg->wall_layer-1))
            {
                particle->x[index] = x*arg->ptc_dx;
                particle->y[index] = y*arg->ptc_dx;
                particle->type[index] = -1;
                index++;
            }
            else if (y < arg->wall_layer)
            {
                particle->x[index] = x*arg->ptc_dx;
                particle->y[index] = y*arg->ptc_dx;
                particle->type[index] = -1;
                index++;
            }
        }
    }
    //if(index != (particle->fluid_ptc_num+particle->wall_ptc_num)) printf("\033[0;32;31m Error in %s:%d\033[m\n",__FILE__,__LINE__);
    assert(index == (arg->fluid_ptc_num + arg->wall_ptc_num));
}

void rigid_ptc_generate(SPH *sph)
{
    SPH_ARG *arg;
    SPH_PARTICLE *particle;
    //SPH_RIGID *rigid;
    arg = sph->host_arg;
    particle = sph->particle;
    //rigid = sph->host_rigid;

    std::string filename = arg->case_dir;
    filename += "/wedge.vtk"; 

    double x[3];
    int index = arg->fluid_ptc_num + arg->wall_ptc_num;

    vtkSmartPointer<vtkUnstructuredGridReader> reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
    reader->SetFileName(filename.c_str());
    reader->Update();
    
    vtkUnstructuredGrid *vtkdata;
    vtkdata = reader->GetOutput();
    for(vtkIdType i=0;i<vtkdata->GetNumberOfPoints();i++)
    {
        vtkdata->GetPoint(i,x);
        if(x[2]==0)
        {
            particle->x[index] = x[0] + arg->fluid_x/2.0;
            particle->y[index] = x[1] + arg->fluid_y*1.1;
            particle->type[index] = 1;
            index++;
        }
    }
    assert(index == arg->ptc_num);
}

void write_vtk(SPH *sph)
{
    SPH_ARG *arg;
    SPH_PARTICLE *particle;
    //SPH_RIGID *rigid;
    arg = sph->host_arg;
    particle = sph->particle;
    //rigid = sph->host_rigid;
    
    string filename = arg->case_dir; 
    filename += "/init.vtk";

    ofstream vtkfile;
    vtkfile.open(filename.c_str());

    vtkfile << "# vtk DataFile Version 3.0" << endl;
    vtkfile << "sph data" << endl;
    vtkfile << "ASCII" << endl;
    vtkfile << "DATASET UNSTRUCTURED_GRID" << endl;
    vtkfile << "POINTS " << arg->ptc_num << " " << "double" << endl;

    for(unsigned int i=0;i<arg->ptc_num;i++)
    {
        vtkfile << setiosflags(ios::scientific) << particle->x[i] << " " \
        << particle->y[i] << " " << 0.0 << endl;
    }

    vtkfile << "POINT_DATA" << " " << arg->ptc_num << endl;

    vtkfile << "SCALARS "<< "density double 1" << endl;
    vtkfile << "LOOKUP_TABLE DEFAULT" << endl;
    for(unsigned int i=0;i<arg->ptc_num;i++)
    {
        vtkfile << setiosflags(ios::scientific) << arg->ref_rho << endl;
    }
    vtkfile << "SCALARS "<< "pressure double 1" << endl;
    vtkfile << "LOOKUP_TABLE DEFAULT" << endl;
    for(unsigned int i=0;i<arg->ptc_num;i++)
    {
        vtkfile << setiosflags(ios::scientific) << 0.0 << endl;
    }
    vtkfile << "SCALARS " << "type int 1" << endl;
    vtkfile << "LOOKUP_TABLE DEFAULT" << endl;
    for(unsigned int i=0;i<arg->ptc_num;i++)
    {
        vtkfile << particle->type[i] << endl; 
    }
    vtkfile << "VECTORS "<< "velocity double" << endl;
    for(unsigned int i=0;i<arg->ptc_num;i++)
    {
        vtkfile << setiosflags(ios::scientific) << 0.0 <<" " << 0.0 << " " \
        << 0.0 << endl;

    }
    vtkfile << "VECTORS "<< "acceleration double" << endl;
    for(unsigned int i=0;i<arg->ptc_num;i++)
    {
        vtkfile << setiosflags(ios::scientific) << 0.0 <<" " << 0.0 << " " \
        << 0.0 << endl;
    }

    vtkfile.close();

}

void rigid_init(SPH *sph)
{
    SPH_ARG *arg;
    SPH_PARTICLE *particle;
    SPH_RIGID *rigid;
    arg = sph->host_arg;
    particle = sph->particle;
    rigid = sph->host_rigid;

    double tmp_cogx = 0.0;
    double tmp_cogy = 0.0;
    double tmp_r = 10000.0;

    for(int i=0;i<arg->ptc_num;i++)
    {
        if(particle->type[i] == 1)
        {
            tmp_cogx += particle->x[i];
            tmp_cogy += particle->y[i];
        }
    }
    tmp_cogx /= (double)arg->rigid_ptc_num;
    tmp_cogy /= (double)arg->rigid_ptc_num;

    for(int i=0;i<arg->ptc_num;i++)
    {
        if(particle->type[i] == 1)
        {
            if(tmp_r >= (pow((particle->x[i]-tmp_cogx),2)+pow((particle->y[i]-tmp_cogy),2)) )
            {
                tmp_r = pow((particle->x[i]-tmp_cogx),2)+pow((particle->y[i]-tmp_cogy),2);
                rigid->cog_ptc_id = i;
            }
        }
    }
    rigid->cogx = particle->x[rigid->cog_ptc_id];
    rigid->cogy = particle->y[rigid->cog_ptc_id];

    for(int i=0;i<arg->ptc_num;i++)
    {
        if(particle->type[i] == 1)
        {
            rigid->moi += (rigid->mass)*(pow((particle->x[i]-rigid->cogx),2)+pow((particle->y[i]-rigid->cogy),2));
        }
    }
    rigid->moi /= arg->rigid_ptc_num;
}
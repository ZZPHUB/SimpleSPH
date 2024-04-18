#include <vtkUnstructuredGridReader.h>
#include <vtkSmartPointer.h>
#include <vtkType.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPointSet.h>
#include <vtkDataSetReader.h>
#include <vtkPointData.h>
#include <vtkDataArray.h>
#include <vtkDoubleArray.h>
#include "IO.cuh"
using namespace std;

void sph_read_vtk(SPH *sph)
{
    SPH_PARTICLE *particle;
    SPH_ARG *arg;
    SPH_RIGID *rigid;
    particle = sph->particle;
    arg = sph->host_arg;
    rigid = sph->host_rigid;

    string filename = arg->case_dir;
    filename += "/init.vtk";

    vtkSmartPointer<vtkUnstructuredGridReader> reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
    reader->SetFileName(filename.c_str());
    reader->ReadAllScalarsOn();
    reader->ReadAllVectorsOn();
    reader->Update();

    vtkUnstructuredGrid *vtkdata;
    vtkdata = reader->GetOutput();

    if(vtkdata->GetNumberOfPoints() != particle->total)
    {
        while(true)
        {
            std::cout << "here is ptc_read_vtk,the num of vtk file not equal to the defined ptc num" << endl;
        }
    }

    vtkPointData* pointdata = vtkdata->GetPointData();  

    vtkDataArray* pressure_array = pointdata->GetScalars("pressure");
    vtkDataArray* density_array = pointdata->GetScalars("density");
    //vtkDataArray* mass_array = pointdata->GetScalars("mass");
    vtkDataArray* type_array = pointdata->GetScalars("type");
    vtkDataArray* velocity_array = pointdata->GetVectors("velocity");
    vtkDataArray* acc_array = pointdata->GetVectors("acceleration");

    vtkDoubleArray* pressure_data = vtkDoubleArray::SafeDownCast(pressure_array);
    vtkDoubleArray* density_data = vtkDoubleArray::SafeDownCast(density_array);
    //vtkDoubleArray* mass_data = vtkDoubleArray::SafeDownCast(mass_array);
    vtkIntArray* type_data = vtkIntArray::SafeDownCast(type_array);
    vtkDoubleArray* velocity_data = vtkDoubleArray::SafeDownCast(velocity_array);
    vtkDoubleArray* acc_data = vtkDoubleArray::SafeDownCast(acc_array);

    if(pressure_data != nullptr && density_data != nullptr \
        && type_data != nullptr && velocity_data != nullptr && acc_data != nullptr)
    {
	    for(vtkIdType i=0;i<vtkdata->GetNumberOfPoints();i++)
	     {
            double p = 0.0;
            double d = 0.0;
            //double m = 0.0;
            double v[3] = {0.0};
            double x[3] = {0.0};
            double a[3] = {0.0};
            int t = 0;
	    
            pressure_data->GetTuple(i,&p);
            density_data->GetTuple(i,&d);
            //mass_data->GetTuple(i,&m);
            t=type_data->GetValue(i);
            velocity_data->GetTuple(i,v);
            vtkdata->GetPoint(i,x);
            acc_data->GetTuple(i,a);
            
            particle->x[i] = x[0];
            particle->y[i] = x[1];
            particle->pressure[i] = p;
            particle->density[i] = d;
            //particle->mass[i] = m;
            particle->type[i] = t;
            particle->vx[i] = v[0];
            particle->vy[i] = v[1];
            particle->accx[i] = a[0];
            particle->accy[i] = a[1];
	    }
    }
    else
    {
        while (true)
        {
            cout << "some case are null" << endl;
        }
        
    }
}

void sph_save_single(SPH *sph)
{
    SPH_PARTICLE *particle;
    SPH_ARG *arg;
    SPH_RIGID *rigid;
    particle = sph->particle;
    arg = sph->host_arg;
    rigid = sph->host_rigid;

    unsigned int ptc_num = 0;
    ptc_num = particle->total;


    string filename = arg->case_dir; 
    filename += "/vtk/sph";
    filename += to_string(sph->host_arg->init_step/sph->host_arg->print_step);
    //filename += to_string(sph->current_step/PRINT_TIME_STEP);
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
    SPH_RIGID *rigid;
    SPH_ARG *arg;
    particle = sph->particle;
    rigid = sph->host_rigid;
    arg = sph->host_arg;

    string filename = arg->case_dir;
    filename += "/last.vtk";
    unsigned int ptc_num = 0;
    ptc_num = particle->total;

    if(sph->host_arg->save_last_flag == 1)
    {
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

        //ptc type
        vtkfile << "SCALARS "<< "type int 1" << endl;
        vtkfile << "LOOKUP_TABLE DEFAULT" << endl;
        for(unsigned int i=0;i<particle->total;i++)
        {
            vtkfile  << particle->type[i] << endl;
        }
        //mass
        /*vtkfile << "SCALARS "<< "mass double 1" << endl;
        vtkfile << "LOOKUP_TABLE DEFAULT" << endl;
        for(unsigned int i=0;i<particle->total;i++)
        {
            vtkfile << particle->mass[i] << endl;
        }*/
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
        //acceleration
        vtkfile << "VECTORS "<< "acceleration double" << endl;
        for(unsigned int i=0;i<particle->total;i++)
        {
            vtkfile << setiosflags(ios::scientific) << particle->accx[i] <<" " << particle->accy[i] << " " \
            << 0.0 << endl;
        }

        vtkfile.close();

        sph_write_info(sph);
    }
}
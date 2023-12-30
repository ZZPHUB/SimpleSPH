#include "PreProcess.H"

void solid_ptc_generate(SPH_PARTICLE *particle)
{
    std::string filename = "../data/preprocess/wedge.vtk";
    double x[3];
    unsigned int tol=0;

    vtkSmartPointer<vtkUnstructuredGridReader> reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
    reader->SetFileName(filename.c_str());
    reader->Update();
    
    vtkUnstructuredGrid *vtkdata;
    vtkdata = reader->GetOutput();
    for(vtkIdType i=0;i<vtkdata->GetNumberOfPoints();i++)
    {
        vtkdata->GetPoint(i,x);
        if(x[2]==0 &&(x[0]!=0 || x[1]!= 0))
        {
            particle->x[FLUID_PTC_NUM+VIRTUAL_PTC_NUM+tol] = (double)(x[0]/1000)+0.2;
            particle->y[FLUID_PTC_NUM+VIRTUAL_PTC_NUM+tol] = (double)(x[1]/1000)+1.0+4*PTC_SPACING;
            particle->type[FLUID_PTC_NUM+VIRTUAL_PTC_NUM+tol] = 1;
            tol++;
        }
    }
}

unsigned int solid_ptc_num(void)
{
    std::string filename = "../data/preprocess/wedge.vtk";
    double x[3];
    unsigned int tol=0;

    vtkSmartPointer<vtkUnstructuredGridReader> reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
    reader->SetFileName(filename.c_str());
    reader->Update();
    
    vtkUnstructuredGrid *vtkdata;
    vtkdata = reader->GetOutput();
    for(vtkIdType i=0;i<vtkdata->GetNumberOfPoints();i++)
    {
        vtkdata->GetPoint(i,x);
        if(x[2]==0 && (x[0]!=0 || x[1]!= 0))
        {
            tol++;
        }
    }
    return tol;
}
#include "PreProcess.H"

void solid_ptc_generate(SPH_PARTICLE *particle)
{
    std::string filename = "./wedge.vtk";
    double x[3];

    vtkSmartPointer<vtkUnstructuredGridReader> reader =      vtkSmartPointer<vtkUnstructuredGridReader>::New();
    reader->SetFileName(filename.c_str());
    reader->Update();
    
    vtkUnstructuredGrid *vtkdata;
    vtkdata = reader->GetOutput();
    //std::cout << "tol p is: " << vtkdata->GetNumberOfPoints() << std::endl;
    for(vtkIdType i=0;i<vtkdata->GetNumberOfPoints();i++)
    {
        vtkdata->GetPoint(i,x);
        if(x[2]==0 &&(x[0]!=0 || x[1]!= 0))
        {
            particle->x[FLUID_PTC_NUM+VIRTUAL_PTC_NUM+(int)i] = x[0];
            particle->y[FLUID_PTC_NUM+VIRTUAL_PTC_NUM+(int)i] = x[1];
        }
    }

    return 0;
}
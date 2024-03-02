#include "PreProcess.H"
using namespace std;

void ptc_rigid_generate(SPH *sph)
//void solid_ptc_generate(SPH_PARTICLE *particle)
{
    SPH_PARTICLE *particle;
    particle = sph->particle;
    
    #if(PTC_SPACING == 0.005)
    std::string filename = "../data/preprocess/wedge005.vtk";
    #elif(PTC_SPACING == 0.001)
    std::string filename = "../data/preprocess/wedge001.vtk";
    #endif

    double x[3] = {0};
    unsigned int tol=0;
    tol = particle->fulid_ptc_num+particle->wall_ptc_num;

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
            particle->x[tol] = (double)(x[0]/1000)+FLUID_DOMAIN_LENGTH/2.0;
            particle->y[tol] = (double)(x[1]/1000)+FLUID_DOMAIN_DEEPTH+4.0*PTC_SPACING;
            particle->type[tol] = 1;
            tol++;
        }
    }
}

unsigned int ptc_rigid_num(void)
{
    #if(PTC_SPACING == 0.005)
    std::string filename = "../data/preprocess/wedge005.vtk";
    #elif(PTC_SPACING == 0.001)
    std::string filename = "../data/preprocess/wedge001.vtk";
    #endif

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

void ptc_rigid_init(SPH *sph)
{
    SPH_PARTICLE *particle;
    SPH_RIGID *wedge;
    particle = sph->particle;
    wedge = sph->rigid;

    wedge->accx = wedge->accy = wedge->alpha = 0.0;
    wedge->mass = 12.8;

    if(sph->new_case_flag == 1)
    {
        wedge->vx = wedge->vy = wedge->omega = 0.0;
        wedge->cogx = FLUID_DOMAIN_LENGTH/2.0;
        wedge->cogy = FLUID_DOMAIN_DEEPTH+4.0*PTC_SPACING+0.032;

        //calculate the moi of the wedge
        for(unsigned int i=0;i<particle->total;i++)
        {
            if(particle->type[i] == 1)
            {
                wedge->moi = (wedge->mass/particle->rigid_ptc_num)*(pow((particle->x[i]-wedge->cogx),2)+pow((particle->y[i]-wedge->cogy),2))
            }
        }
    }
    else if(sph->new_case_flag == 0)
    {
        cout << "the wedge's velocity in x-direction: " << endl;
        cin >> wedge->vx;
        cout << "the wedge's velocity in y-direction: " << endl;
        cin >> wedge->vy;
        cout << "the wedge's omega is : " << endl;
        cin >> wedge->omega;
        cout << "the wedge's center of gravity in x-direction: " << endl;
        cin >> wedge->cogx;
        cout << "the wedge's center of gravity in y-direction: " << endl;
        cin >> wedge->cogy;
        cout << "the wedge's moi is : " << endl;
        cin >> wedge->moi;
    }   
}

void ptc_read_vtk(SPH *sph)
{
    SPH_PARTICLE *particle;
    particle = sph->particle;

    std::string filename = "../Preprocess/init.vtk";

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
    vktDataArray* density_array = pointdata->GetScalars("density");
    vktDataArray* type_array = pointdata->GetScalars("type");
    vktDataArray* velocity_array = pointdata->GetVectors("velocity");

    vtkDoubleArray* pressure_data = vtkDoubleArray::SafeDownCast(pressure_array);
    vtkDoubleArray* density_data = vtkDoubleArray::SafeDownCast(density_array);
    vtkIntArray* type_data = vtkIntArray::SafeDownCast(type_array);
    vtkDoubleArray* velocity_data = vtkDoubleArray::SafeDownCast(velocity_array);

    if(pressure_data != nullptr && density_data != nullptr && type_data != nullptr && velocity_data != nullptr)
    {
	    for(unsigned int i=0;i<vtk->GetNumberOfPoints();i++)
	     {
            double p = 0.0;
            double d = 0.0;
            double v[3] = {0.0};
            double x[3] = {0.0};
            int t = 0;
	        
            pressure_data->GetTuple(i,&p);
            density_data->GetTuple(i,&d);
            type_data->GetTuple(i,&t);
            velocity_data->GetTuple(i,v);
            vtkdata->GetPoint(i,x);
            
            particle->x[i] = x[0];
            particle->y[i] = x[1];
            particle->pressure[i] = p;
            particle->density[i] = d;
            particle->type[i] = t;
            particle->vx[i] = v[0];
            particle->vy[i] = v[1];
	    }
    }
}
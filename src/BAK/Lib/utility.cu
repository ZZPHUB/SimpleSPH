#include "Lib.cuh"

void ptc_density_correct(SPH *sph)
{
    SPH_PARTICLE *particle;
    SPH_PAIR *pair;
    SPH_KERNEL *kernel;
    particle = sph->particle;
    pair= sph->pair;
    kernel = sph->kernel;

    double a = ALPHA;
    //double m = PTC_MASS;


    for(unsigned int i=0;i<particle->total;i++)
    {
        if(particle->type[i] == 0)
        {
            particle->w[i] = (a*2.0*particle->mass[i])/(3.0*particle->density[i]);
        }
    }

    for(unsigned int i=0;i<pair->total;i++)
    {
        particle->w[pair->i[i]] += kernel->w[i]*particle->mass[pair->j[i]]/particle->density[pair->j[i]];
        if(particle->type[pair->j[i]]==0)
        {
            particle->w[pair->j[i]] += kernel->w[i]*particle->mass[pair->i[i]]/particle->density[pair->i[i]];
        }
    }

    for(unsigned int i=0;i<particle->total;i++)
    {
        if(particle->type[i] == 0)
        {
            particle->density[i] = (a*2.0*particle->mass[i])/(3.0*particle->w[i]);
        }
    }

    for(unsigned int i=0;i<pair->total;i++)
    {
        particle->density[pair->i[i]] += particle->mass[pair->j[i]]*kernel->w[i]/particle->w[pair->i[i]];
        if(particle->type[pair->j[i]] == 0)
        {
            particle->density[pair->j[i]] += particle->mass[pair->i[i]]*kernel->w[i]/particle->w[pair->j[i]];
        }
    }
}

void ptc_dummy(SPH *sph)
{
    SPH_PARTICLE *particle;
    SPH_PAIR *pair;
    SPH_KERNEL *kernel;
    SPH_RIGID *wedge;
    particle = sph->particle;
    pair = sph->pair;
    kernel = sph->kernel;
    wedge = sph->rigid;


    //rigid body(wall & wedge)vx,vy,accx,accy,pressure init
    for(unsigned int i=0;i<particle->total;i++)
    {
        if(particle->type[i] != 0)
        {
            particle->w[i] = 0;
            particle->vx[i] = 0;
            particle->vy[i] = 0;
            particle->pressure[i] = 0;
            particle->density[i] = 0;
        }
    }
    
    //the not fluid weight term 
    for(unsigned int i=0;i<pair->total;i++)
    {
        if(particle->type[pair->j[i]] != 0) 
        {
            particle->w[pair->j[i]] += kernel->w[i];
        }
    }
    
    //rigid body(wall & wedege) pressure and velocity
    for(unsigned int i=0;i<pair->total;i++)
    {
        double dx = 0.0;
        double dy = 0.0;
        double rigid_acc_x = 0.0;
        double rigid_acc_y = 0.0;
        if(particle->type[pair->j[i]] != 0 && particle->w[pair->j[i]] != 0.0)
        {
            if(particle->type[pair->j[i]] == -1)
            {
                rigid_acc_x = 0.0;
                rigid_acc_y = 0.0;
            }
            else if (particle->type[pair->j[i]] == 1)
            {
                rigid_acc_x = wedge->accx - pow(wedge->omega,2)*(particle->x[pair->j[i]]-wedge->cogx)- \
                              wedge->alpha*(particle->y[pair->j[i]]-wedge->cogy);
                rigid_acc_y = wedge->accy - pow(wedge->omega,2)*(particle->y[pair->j[i]]-wedge->cogy)+ \
                              wedge->alpha*(particle->x[pair->j[i]]-wedge->cogx);
            }
            dx = particle->x[pair->i[i]] - particle->x[pair->j[i]];
            dy = particle->y[pair->i[i]] - particle->y[pair->j[i]];
            particle->pressure[pair->j[i]] += (particle->pressure[pair->i[i]]+particle->density[pair->i[i]]*\
                        (rigid_acc_x*dx+(rigid_acc_y+GRAVITY_ACC)*dy))*kernel->w[i]/particle->w[pair->j[i]];
            particle->vx[pair->j[i]] += particle->vx[pair->i[i]]*kernel->w[i]/particle->w[pair->j[i]];
            particle->vy[pair->j[i]] += particle->vy[pair->i[i]]*kernel->w[i]/particle->w[pair->j[i]];
        }
    }

    //rigid body(wall & wedege) densiy
    for(unsigned int i=0;i<particle->total;i++)
    {
        if(particle->type[i] != 0)
        {
            particle->density[i] = particle->pressure[i]/pow(sph->c,2)+REF_DENSITY;
        }
    }
}

void sph_avg_time(SPH *sph)
{
    static unsigned int step = 0;
    static double start;
    static double end;
    if(step == 0)
    {
        start = (double)time(nullptr);
    }
    else
    {
        end = (double)time(nullptr);
        sph->avg_time = (end-start)/(double)step;
    }
    step++;
}
void sph_read_info(SPH *sph)
{
    string filename = sph->arg->case_dir;
    filename += "/case_info"
    ifstream case_info;
    case_info.open(filename.c_str());
    string line;
    //case_info << sph->arg->c << endl;
    getline(case_info,line);
    sph->arg->c = stod(line.c_str());
    //case_info << sph->arg->g <<endl;   //gravity acceleration
    getline(case_info,line);
    sph->arg->g = stod(line.c_str());
    //case_info << sph->arg->ref_rho << endl;     //reference density
    getline(case_info,line);
    sph->arg->ref_rho = stod(line.c_str());
    //case_info << sph->arg->ptc_dx << endl;      //ptc delta spacing
    getline(case_info,line);
    sph->arg->ptc_dx = stod(line.c_str());
    //case_info << sph-arg->r << endl;       //ptc radius
    getline(case_info,line);
    sph->arg->r = stod(line.c_str());
    //case_info << sph->arg->h << endl;   //smoothed length
    getline(case_info,line);
    sph->arg->h = stod(line.c_str());
    //case_info << sph->arg->m << endl;   //ptc mass
    getline(case_info,line);
    sph->arg->m = stod(line.c_str());
    //case_info << sph->arg->alpha << endl;   //kernel function's para
    getline(case_info,line);
    sph->arg->alpha = stod(line.c_str());
    //case_info << sph->arg->sst << endl;    //single step time
    getline(case_info,line);
    sph->arg->sst = stod(line.c_str());
    //case_info << sph->arg->dt << endl;  //delta t
    getline(case_info,line);
    sph->arg->dt = stod(line.c_str());

    //case_info << sph->arg->fluid_x << endl;     //fluid length
    getline(case_info,line);
    sph->arg->fluid_x = stod(line.c_str());
    //case_info << sph->arg->fluid_y << endl;     //fluid depth
    getline(case_info,line);
    sph->arg->fluid_y = stod(line.c_str());
    //case_info << sph->arg->fluid_xnum << endl;     //fluid length direction ptc num
    getline(case_info,line);
    sph->arg->fluid_xnum = stoi(line.c_str());
    //case_info << sph->arg->fluid_ynum << endl;     //fluid depth direction ptc num
    getline(case_info,line);
    sph->arg->fluid_ynum = stoi(line.c_str());
    
    //case_info << sph->arg->domain_x << endl;    //total domain length
    getline(case_info,line);
    sph->arg->domain_x = stod(line.c_str());
    //case_info << sph->arg->domain_y << endl;    //total domain depth
    getline(case_info,line);
    sph->arg->domain_y = stod(line.c_str());

    //case_info << sph->arg->mesh_dx << endl;     //mesh delta spacing
    getline(case_info,line);
    sph->arg->mesh_dx = stod(line.c_str());
    //case_info << sph->arg->mesh_xnum << endl;      //mesh length direction num
    getline(case_info,line);
    sph->arg->mesh_xnum = stoi(line.c_str());
    //case_info << sph->arg->mesh_ynum << endl;      //mesh depth direction num
    getline(case_info,line);
    sph->arg->mesh_ynum = stoi(line.c_str());
    //case_info << sph->arg->mesh_num << endl;       //total mesh num
    getline(case_info,line);
    sph->arg->mesh_num = stoi(line.c_str());
    //case_info << sph->arg->mesh_volume << endl;    //single mesh volume
    getline(case_info,line);
    sph->arg->mesh_volume = stoi(line.c_str());

    //case_info << sph->arg->init_step << endl;  //inital time step
    getline(case_info,line);
    sph->arg->init_setp = stoi(line.c_str());
    //case_info << sph->arg->total_step << endl; //total time step
    getline(case_info,line);
    sph->arg->total_step = stoi(line.c_str());

    //current process flags
    //case_info << sph->arg->new_case_flag << endl;  // if 1 then creat a new case,or continue to run the old case
    getline(case_info,line);
    sph->arg->new_case_flag = stoi(line.c_str());
    //case_info << sph->arg->init_impac_flag << endl; //if 1 then run the init step,or run the impac step
    getline(case_info,line);
    sph->arg->init_impac_flag = stoi(line.c_str());
    //case_info << sph->arg->save_last_flag << endl; //if 1 then save the last step,or donnot save it
    getline(case_info,line);
    sph->arg->save_last_flag = stoi(line.c_str());

    //rigid info
    //case_info << sph->rigid->vx << endl;  //rigid body x-direction velocity
    getline(case_info,line);
    sph->rigid->vx = stod(line.c_str());
    //case_info << sph->rigid->vy << endl;  //rigid body y-direction velocity
    getline(case_info,line);
    sph->rigid->vy = stod(line.c_str());
    //case_info << sph->rigid->omega << endl;   //rigid body angular velocity
    getline(case_info,line);
    sph->rigid->omega = stod(line.c_str());
    //case_info << sph->rigid->accx << endl;    //rigid body x-direciton acceleration
    getline(case_info,line);
    sph->rigid->accx = stod(line.c_str());
    //case_info << sph->rigid->accy << endl;    //rigid body y-direction acceleration
    getline(case_info,line);
    sph->rigid->accy = stod(line.c_str());
    //case_info << sph->rigid->alpha << endl;   //rigid body angular acceleration
    getline(case_info,line);
    sph->rigid->alpha = stod(line.c_str());
    //case_info << sph->rigid->cogx << endl;    //x-direction center of gravity coordinate
    getline(case_info,line);
    sph->rigid->cogx = stod(line.c_str());
    //case_info << sph->rigid->cogy << endl;    //y-direction center of gravity coordinate 
    getline(case_info,line);
    sph->rigid->cogy = stod(line.c_str());
    //case_info << sph->rigid->mass << endl;    //rigid body mass 
    getline(case_info,line);
    sph->rigid->mass = stod(line.c_str());
    //case_info << sph->rigid->moi << endl;     //rigid body moment of inertia
    getline(case_info,line);
    sph->rigid->moi = stod(line.c_str());
    //case_info << sph->rigid->total << endl;   //rigid body ptc num
    getline(case_info,line);
    sph->rigid->total = stoi(line.c_str());
    
    case_info.close();

}

void sph_write_info(SPH *sph)
{
    string filename = sph->arg->case_dir;
    filename += "/case_info"
    ofstream case_info;
    case_info.open(filename.c_str());

    //sph paraments
    case_info << sph->arg->c << endl;
    case_info << sph->arg->g <<endl;   //gravity acceleration
    case_info << sph->arg->ref_rho << endl;     //reference density
    case_info << sph->arg->ptc_dx << endl;      //ptc delta spacing
    case_info << sph-arg->r << endl;       //ptc radius
    case_info << sph->arg->h << endl;   //smoothed length
    case_info << sph->arg->m << endl;   //ptc mass
    case_info << sph->arg->alpha << endl;   //kernel function's para
    case_info << sph->arg->sst << endl;    //single step time
    case_info << sph->arg->dt << endl;  //delta t

    case_info << sph->arg->fluid_x << endl;     //fluid length
    case_info << sph->arg->fluid_y << endl;     //fluid depth
    case_info << sph->arg->fluid_xnum << endl;     //fluid length direction ptc num
    case_info << sph->arg->fluid_ynum << endl;     //fluid depth direction ptc num
    
    case_info << sph->arg->domain_x << endl;    //total domain length
    case_info << sph->arg->domain_y << endl;    //total domain depth

    case_info << sph->arg->mesh_dx << endl;     //mesh delta spacing
    case_info << sph->arg->mesh_xnum << endl;      //mesh length direction num
    case_info << sph->arg->mesh_ynum << endl;      //mesh depth direction num
    case_info << sph->arg->mesh_num << endl;       //total mesh num
    case_info << sph->arg->mesh_volume << endl;    //single mesh volume

    case_info << sph->arg->init_step << endl;  //inital time step
    case_info << sph->arg->total_step << endl; //total time step

    //current process flags
    case_info << sph->arg->new_case_flag << endl;  // if 1 then creat a new case,or continue to run the old case
    case_info << sph->arg->init_impac_flag << endl; //if 1 then run the init step,or run the impac step
    case_info << sph->arg->save_last_flag << endl; //if 1 then save the last step,or donnot save it

    //rigid info
    case_info << sph->rigid->vx << endl;  //rigid body x-direction velocity
    case_info << sph->rigid->vy << endl;  //rigid body y-direction velocity
    case_info << sph->rigid->omega << endl;   //rigid body angular velocity
    case_info << sph->rigid->accx << endl;    //rigid body x-direciton acceleration
    case_info << sph->rigid->accy << endl;    //rigid body y-direction acceleration
    case_info << sph->rigid->alpha << endl;   //rigid body angular acceleration
    case_info << sph->rigid->cogx << endl;    //x-direction center of gravity coordinate
    case_info << sph->rigid->cogy << endl;    //y-direction center of gravity coordinate 
    case_info << sph->rigid->mass << endl;    //rigid body mass 
    case_info << sph->rigid->moi << endl;     //rigid body moment of inertia
    case_info << sph->rigid->total << endl;   //rigid body ptc num

    case_info.close();
}

void sph_generate_case(SPH *sph)
{
    cout << "ptc spacing is:";
    cin >> sph->arg->ptc_dx;

    cout << "fluid domain length is:";
    cin >> sph->arg->fluid_x;

    cout << "fluid domain depth is:";
    cin >> sph->arg->fluid_y;

    sph->arg->fluid_xnum = int(sph->arg->fluid_x/sph->arg->ptc_dx)+1;
    sph->arg->fluid_ynum = int(sph->arg->fluid_y/sph->arg->ptc_dx)+1;

    sph->arg->domain_x = sph->arg->fluid_x;
    sph->arg->domain_y = 1.5*sph->arg->fluid.y;
    
    sph->arg->h = 1.005*sph->arg->dx;
    sph->arg->r = 2.0*sph->arg->h;
    sph->arg->mesh_dx = sph->arg->r;
    sph->arg->mesh_xnum = int(sph->arg->domain_x/sph->arg->mesh_dx)+1;
    sph->arg->mesh_ynum = int(sph->arg->domain_y/sph->arg->mesh_dx)+1;
    sph->arg->mesh_num = sph->arg->mesh_xnum*sph->arg->mesh_ynum;
    sph->arg->mesh_volume = 32;

    sph->arg->g = 9.8;
    sph->arg->c = 10*sqrt(sph->arg->g*sph->arg->fluid_y);
    sph->arg->ref_rho = 1000.0;
    sph->arg->m = sph->arg->ref_rho*sph->arg->ptc_dx*sph->arg->ptc_dx;
    sph->arg->alpha = 15.0/(7*3.14159265358*sph->arg->h*sph->arg->h);

    cout << "the init step:";
    cin >> sph->arg->init_step;
    cout << "the total step:";
    cin >> sph->arg->total_step;

    cout << "new case flag,if 1 then creat a new case:";
    cin >> sph->arg->new_case_flag;
    cout << "init or impac simulation,if 1 then run the init step:";
    cin >> sph->arg->init_impac_flag;
    cout << "save the last time step info,if 1 to save:";
    cin >> sph->arg->save_last_flag;
}

void sph_change_case(SPH *sph)
{

}

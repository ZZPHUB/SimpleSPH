#include "SPH.cuh"
#include <iostream>
using namespace std;

void get_input(SPH *);
void fluid_ptc_generate(SPH *);
void rigid_ptc_generate(SPH *);

int main(int argc,char *argv[])
{
    SPH_ARG arg;
    SPH_RIGID rigid;
    SPH_PARTICLE particle;
    SPH sph;
    sph.host_arg = &arg;
    sph.host_rigid = &rigid;
    sph.particle = &particle;
    if(argc != 2) printf("\033[0;32;31m Error in %s:%d\n",__FILE__,__LINE__);
    arg.case_dir = argv[1];
    get_input(&sph);
    fluid_ptc_generate(&sph);
    rigid_ptc_generate(&sph);
    sph_write_info(&sph);

    return 0;
}

void get_input(SPH *sph)
{
    SPH_ARG *arg;
    arg = sph->host_arg;
    SPH_RIGID *rigid;
    rigid = sph->host_rigid;
    SPH_PARTICLE *particle;
    particle = sph->particle;

    cout << "fluid length (fluid_x) is:" << endl;
    cin >> arg->fluid_x;
    cout << "fluid depth (fluid_y) is:" << endl;
    cin >> arg->fluid_y;
    cout << "particle spacing (ptc_dx) is:" << endl;
    cin >> arg->ptc_dx;

    arg->h = 1.0005*arg->ptc_dx;
    arg->r = 2.0*arg->h;
    arg->wall_layer = 2;
    arg->fluid_xnum = (int)(arg->fluid_x/arg->ptc_dx)+1-2*arg->wall_layer;
    arg->fluid_ynum = (int)(arg->fluid_y/arg->ptc_dx)+1-arg->wall_layer;
    particle->rigid_ptc_num = 0;
    particle->fluid_ptc_num = arg->fluid_xnum*arg->fluid_ynum;
    particle->wall_ptc_num = ((int)(arg->fluid_x/arg->ptc_dx)+1)*((int)(1.1*arg->fluid_y/arg->ptc_dx)+1)- \
                             ((int)(arg->fluid_x/arg->ptc_dx)+1-2*arg->wall_layer)*((int)(1.1*arg->fluid_y/arg->ptc_dx)+1-arg->wall_layer);
    particle->total = particle->fluid_ptc_num + particle->wall_ptc_num;
    arg->ptc_num = particle->total;

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
    arg->alpha = 15.0/(7.0*3.14159265358*pow(arg->h,2));

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
    rigid->offset_x = 0.0;
    rigid->offset_y = 0.0;
    rigid->offset_angl = 0.0;
    rigid->cog_ptc_id = 0.0;
    rigid->cogx = 0.0;
    rigid->cogy = 0.0; 
    rigid->total = 0;

    //return 0;
}

void fluid_ptc_generate(SPH *sph)
{
    //return 0;
}

void rigid_ptc_generate(SPH *sph)
{
    //return 0;
}


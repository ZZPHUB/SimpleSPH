#include "SPH.cuh"
#include <assert.h>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>

using namespace std;

int main(int argc,char *argv[])
{
    SPH_ARG arg;
    SPH_RIGID rigid;
    SPH_PARTICLE particle;
    SPH sph;
    sph.host_arg = &arg;
    sph.host_rigid = &rigid;
    sph.particle = &particle;

    assert(argc == 3);
    arg.case_dir = argv[1];
    int loop_num = stoi(argv[2]);

    sph_read_info(&sph);
    particle.x = (double *)calloc(arg.ptc_num,sizeof(double));
    particle.y = (double *)calloc(arg.ptc_num,sizeof(double));
    particle.vx = (double *)calloc(arg.ptc_num,sizeof(double));
    particle.vy = (double *)calloc(arg.ptc_num,sizeof(double));
    particle.accx = (double *)calloc(arg.ptc_num,sizeof(double));
    particle.accy = (double *)calloc(arg.ptc_num,sizeof(double));
    particle.density = (double *)calloc(arg.ptc_num,sizeof(double));
    particle.pressure = (double *)calloc(arg.ptc_num,sizeof(double));
    particle.type = (int *)calloc(arg.ptc_num,sizeof(int));
    sph_read_vtk(&sph);

    vector<double> tmp_x;
    vector<double> tmp_y;
    vector<double> tmp_vx;
    vector<double> tmp_vy;
    vector<double> tmp_accx;
    vector<double> tmp_accy;
    vector<double> tmp_rho;
    vector<double> tmp_p;
    vector<int> tmp_type;

    int loop_fluid_num = 0;
    int loop_wall_num = 0;
    int offset_xnum = 0;
    for(int i=0;i<loop_num;i++)
    {
        if(i == 0)
        {
            for(int j=0;j<arg.ptc_num;j++)
            {
                if(particle.x[j] < (arg.fluid_xnum+arg.wall_layer)*arg.ptc_dx)
                {
                    if(particle.type[j] != 1)
                    {
                        tmp_x.push_back(particle.x[j]);
                        tmp_y.push_back(particle.y[j]);
                        tmp_vx.push_back(particle.vx[j]);
                        tmp_vy.push_back(particle.vy[j]);
                        tmp_accx.push_back(particle.accx[j]);
                        tmp_accy.push_back(particle.accy[j]);
                        tmp_rho.push_back(particle.density[j]);
                        tmp_p.push_back(particle.pressure[j]);
                        tmp_type.push_back(particle.type[j]);
                        if(particle.type[j] == 0) loop_fluid_num ++;
                        else loop_wall_num ++;
                    }
                }
            }
            offset_xnum = arg.wall_layer + arg.fluid_xnum;
        }
        else if (i == loop_num-1)
        {
           for(int j=0;j<arg.ptc_num;j++)
            {
                if(particle.x[j] > (arg.wall_layer-1)*arg.ptc_dx)
                {
                    if(particle.type[j] != 1)
                    {
                        tmp_x.push_back(particle.x[j]+(offset_xnum-arg.wall_layer)*arg.ptc_dx);
                        tmp_y.push_back(particle.y[j]);
                        tmp_vx.push_back(particle.vx[j]);
                        tmp_vy.push_back(particle.vy[j]);
                        tmp_accx.push_back(particle.accx[j]);
                        tmp_accy.push_back(particle.accy[j]);
                        tmp_rho.push_back(particle.density[j]);
                        tmp_p.push_back(particle.pressure[j]);
                        tmp_type.push_back(particle.type[j]);
                        if(particle.type[j] == 0) loop_fluid_num ++;
                        else loop_wall_num ++;
                    }
                }
            } 
        }
        else
        {
            for(int j=0;j<arg.ptc_num;j++)
            {
                if(particle.x[j] < (arg.fluid_xnum+arg.wall_layer)*arg.ptc_dx && particle.x[j] > (arg.wall_layer-1)*arg.ptc_dx)
                {
                    if(particle.type[j] != 1)
                    {
                        tmp_x.push_back(particle.x[j]+(offset_xnum-arg.wall_layer)*arg.ptc_dx);
                        tmp_y.push_back(particle.y[j]);
                        tmp_vx.push_back(particle.vx[j]);
                        tmp_vy.push_back(particle.vy[j]);
                        tmp_accx.push_back(particle.accx[j]);
                        tmp_accy.push_back(particle.accy[j]);
                        tmp_rho.push_back(particle.density[j]);
                        tmp_p.push_back(particle.pressure[j]);
                        tmp_type.push_back(particle.type[j]);
                        if(particle.type[j] == 0) loop_fluid_num ++;
                        else loop_wall_num ++;
                    }
                }
            }
            offset_xnum += arg.fluid_xnum;
        }
    }
    double tmp_fluid_x = 0.0;
    double tmp_fluid_y = 0.0;
    int tmp_fluid_xnum = 0;
    int tmp_fluid_ynum = 0;
    int tmp_fluid_num = 0;
    int tmp_wall_num = 0;
    tmp_fluid_xnum = loop_num*arg.fluid_xnum;
    tmp_fluid_ynum = arg.fluid_ynum;
    tmp_fluid_x = (tmp_fluid_xnum + 2*arg.wall_layer -1)*arg.ptc_dx;
    tmp_fluid_y = (arg.fluid_ynum + arg.wall_layer -1)*arg.ptc_dx;
    assert(tmp_fluid_y == arg.fluid_y);
    tmp_fluid_num = tmp_fluid_xnum*tmp_fluid_ynum;
    tmp_wall_num = ((int)(tmp_fluid_x/arg.ptc_dx) + 1)*((int)(1.1*tmp_fluid_y/arg.ptc_dx)+1) -\
                   ((int)(tmp_fluid_x/arg.ptc_dx + 1 - 2*arg.wall_layer))*((int)(1.1*tmp_fluid_y/arg.ptc_dx)+1 - arg.wall_layer);
    double rigid_offset = tmp_fluid_x/2.0 - rigid.cogx;

    int loop_rigid_num = 0;
    for(int i=0;i<arg.ptc_num;i++)
    {
        if(particle.type[i] == 1)
        {
            tmp_x.push_back(particle.x[i]+rigid_offset);
            tmp_y.push_back(particle.y[i]);
            tmp_vx.push_back(particle.vx[i]);
            tmp_vy.push_back(particle.vy[i]);
            tmp_accx.push_back(particle.accx[i]);
            tmp_accy.push_back(particle.accy[i]);
            tmp_rho.push_back(particle.density[i]);
            tmp_p.push_back(particle.pressure[i]);
            tmp_type.push_back(particle.type[i]); 
            if(rigid.cog_ptc_id == i)
            {
                rigid.cog_ptc_id = loop_fluid_num + loop_wall_num + loop_rigid_num;
            }
            loop_rigid_num++;
        }
    }

    
    assert(loop_fluid_num == tmp_fluid_num);
    assert(loop_wall_num == tmp_wall_num);
    assert(loop_rigid_num == arg.rigid_ptc_num);
    //rigid.cogx = particle.x[rigid.cog_ptc_id];
    //rigid.cogy = particle.y[rigid.cog_ptc_id];    

    arg.fluid_x = tmp_fluid_x;
    arg.fluid_y = tmp_fluid_y;
    arg.fluid_xnum = tmp_fluid_xnum;
    arg.fluid_ynum = tmp_fluid_ynum;
    arg.fluid_ptc_num = tmp_fluid_num;
    arg.wall_ptc_num = tmp_wall_num;
    arg.ptc_num = arg.fluid_ptc_num + arg.wall_ptc_num + arg.rigid_ptc_num;
   
    arg.mesh_x = arg.fluid_x;
    arg.mesh_y = arg.fluid_y * 1.5;
    arg.mesh_xnum = (int)(arg.mesh_x/arg.mesh_dx)+1;
    arg.mesh_ynum = (int)(arg.mesh_y/arg.mesh_dx)+1;
    arg.mesh_num = arg.mesh_xnum*arg.mesh_ynum;

    arg.pair_volume = (int)(32*arg.ptc_num/arg.mesh_num);

    

    free(particle.x);
    free(particle.y);
    free(particle.vx);
    free(particle.vy);
    free(particle.accx);
    free(particle.accy);
    free(particle.density);
    free(particle.pressure);
    free(particle.type);

    particle.x = (double *)calloc(arg.ptc_num,sizeof(double));
    particle.y = (double *)calloc(arg.ptc_num,sizeof(double));
    particle.vx = (double *)calloc(arg.ptc_num,sizeof(double));
    particle.vy = (double *)calloc(arg.ptc_num,sizeof(double));
    particle.accx = (double *)calloc(arg.ptc_num,sizeof(double));
    particle.accy = (double *)calloc(arg.ptc_num,sizeof(double));
    particle.density = (double *)calloc(arg.ptc_num,sizeof(double));
    particle.pressure = (double *)calloc(arg.ptc_num,sizeof(double));
    particle.type = (int *)calloc(arg.ptc_num,sizeof(int));

    
    assert( tmp_x.size() == arg.ptc_num);
    for(int i=0;i<arg.ptc_num;i++)
    {
        particle.x[i] = tmp_x.at(i);
        particle.y[i] = tmp_y.at(i);
        particle.vx[i] = tmp_vx.at(i);
        particle.vy[i] = tmp_vy.at(i);
        particle.accx[i] = tmp_accx.at(i);
        particle.accy[i] = tmp_accy.at(i);
        particle.density[i] = tmp_rho.at(i);
        particle.pressure[i] = tmp_p.at(i);
        particle.type[i] = tmp_type.at(i);
    }

    rigid.cogx = particle.x[rigid.cog_ptc_id];
    rigid.cogy = particle.y[rigid.cog_ptc_id]; 
    
    sph_save_last(&sph);    
    return 0;
}


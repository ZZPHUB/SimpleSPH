#include "IO.cuh"
#include <nlohmann/json.hpp>
using json = nlohmann::json;

void sph_read_info(SPH *sph)
{
    SPH_ARG *arg;
    SPH_RIGID *rigid;
    SPH_PARTICLE *particle;
    arg = sph->host_arg;
    rigid = sph->host_rigid;
    particle = sph->particle;

    std::string info_file = arg->case_dir;
    info_file += "/init_case_info.json";
    std::ifstream f(info_file.c_str());
    json sph_info = json::parse(f);

    arg->fluid_x = sph_info["fluid"]["x"];
    arg->fluid_y = sph_info["fluid"]["y"];
    arg->ptc_dx = sph_info["fluid"]["dx"];
    arg->fluid_xnum = sph_info["fluid"]["xnum"];
    arg->fluid_ynum = sph_info["fluid"]["ynum"];
    particle->fluid_ptc_num = sph_info["fluid"]["fluid_num"];
    particle->wall_ptc_num = sph_info["fluid"]["wall_num"];

    arg->mesh_dx = sph_info["mesh"]["dx"];
    arg->mesh_x = sph_info["mesh"]["x"];
    arg->mesh_y = sph_info["mesh"]["y"];
    arg->mesh_xnum = sph_info["mesh"]["xnum"];
    arg->mesh_ynum = sph_info["mesh"]["ynum"];
    arg->mesh_num = sph_info["mesh"]["num"];
    arg->mesh_volume = sph_info["mesh"]["volume"];

    arg->pair_volume = sph_info["pair"]["volume"];

    arg->c = sph_info["arg"]["c"];
    arg->h = sph_info["arg"]["h"];
    arg->m = sph_info["arg"]["m"];
    arg->g = sph_info["arg"]["g"];
    arg->ref_rho = sph_info["arg"]["ref_rho"];
    arg->alpha = sph_info["arg"]["alpha"];
    arg->r = sph_info["arg"]["r"];
    arg->ptc_num = sph_info["arg"]["ptc_num"];
    arg->wall_layer = sph_info["arg"]["wall_layer"];
    
    arg->dt = sph_info["time"]["dt"];
    arg->sst = sph_info["time"]["sst"];
    arg->init_step = sph_info["time"]["init_step"];
    arg->total_step = sph_info["time"]["total_step"];
    arg->print_step = sph_info["time"]["print_step"];

    arg->new_case_flag = sph_info["flag"]["new_case_flag"];
    arg->init_impac_flag = sph_info["flag"]["init_impac_flag"];
    arg->save_last_flag = sph_info["flag"]["save_last_flag"];

    rigid->vx = sph_info["rigid"]["vx"];
    rigid->vy = sph_info["rigid"]["vy"];
    rigid->omega = sph_info["rigid"]["omega"];
    rigid->accx = sph_info["rigid"]["accx"];
    rigid->accy = sph_info["rigid"]["accy"];
    rigid->alpha = sph_info["rigid"]["alpha"];
    rigid->cogx = sph_info["rigid"]["cogx"];
    rigid->cogy = sph_info["rigid"]["cogy"];
    rigid->cog_ptc_id = sph_info["rigid"]["cog_ptc_id"];
    rigid->offset_x = sph_info["rigid"]["offset_x"];
    rigid->offset_y = sph_info["rigid"]["offset_y"];
    rigid->offset_angl = sph_info["rigid"]["offset_angl"];
    rigid->mass = sph_info["rigid"]["mass"];
    rigid->moi = sph_info["rigid"]["moi"];
    rigid->total = sph_info["rigid"]["rigid_num"];
    particle->rigid_ptc_num = rigid->total;
}

void sph_write_info(SPH *sph)
{
    SPH_ARG *arg;
    SPH_RIGID *rigid;
    SPH_PARTICLE *particle;
    arg = sph->host_arg;
    rigid = sph->host_rigid;
    particle = sph->particle;
    
    std::string info_file = arg->case_dir;
    info_file += "/last_case_info.json";
    std::ofstream f(info_file.c_str());
    json sph_info;

    //arg->fluid_x = sph_info["fluid"]["x"];
    sph_info["fluid"]["x"] = arg->fluid_x;
    //arg->fluid_y = sph_info["fluid"]["y"];
    sph_info["fluid"]["y"] = arg->fluid_y;
    //arg->ptc_dx = sph_info["fluid"]["dx"];
    sph_info["fluid"]["dx"] = arg->ptc_dx;
    //arg->fluid_xnum = sph_info["fluid"]["xnum"];
    sph_info["fluid"]["xnum"] = arg->fluid_xnum;
    //arg->fluid_ynum = sph_info["fluid"]["ynum"];
    sph_info["fluid"]["ynum"] = arg->fluid_ynum;
    //particle->fulid_ptc_num = sph_info["fluid"]["fluid_num"];
    sph_info["fluid"]["fluid_num"] = particle->fluid_ptc_num;
    //particle->wall_ptc_num = sph_info["fluid"]["wall_num"];
    sph_info["fluid"]["wall_num"] = particle->wall_ptc_num;

    //arg->mesh_dx = sph_info["mesh"]["dx"];
    sph_info["mesh"]["dx"] = arg->mesh_dx;
    //arg->mesh_x = sph_info["mesh"]["x"];
    sph_info["mesh"]["x"] = arg->mesh_x;
    //arg->mesh_y = sph_info["mesh"]["y"];
    sph_info["mesh"]["y"] = arg->mesh_y;
    //arg->mesh_xnum = sph_info["mesh"]["xnum"];
    sph_info["mesh"]["xnum"] = arg->mesh_xnum;
    //arg->mesh_ynum = sph_info["mesh"]["ynum"];
    sph_info["mesh"]["ynum"] = arg->mesh_ynum;
    //arg->mesh_num = sph_info["mesh"]["num"];
    sph_info["mesh"]["num"] = arg->mesh_num;
    //arg->mesh_volume = sph_info["mesh"]["volume"];
    sph_info["mesh"]["volume"] = arg->mesh_volume;

    //arg->pair_volume = sph_info["pair"]["volume"];
    sph_info["pair"]["volume"] = arg->pair_volume;

    //arg->c = sph_info["arg"]["c"];
    sph_info["arg"]["c"] = arg->c;
    //arg->h = sph_info["arg"]["h"];
    sph_info["arg"]["h"] = arg->h;
    //arg->m = sph_info["arg"]["m"];
    sph_info["arg"]["m"] = arg->m;
    //arg->g = sph_info["arg"]["g"];
    sph_info["arg"]["g"] = arg->g;
    //arg->ref_rho = sph_info["arg"]["ref_rho"];
    sph_info["arg"]["ref_rho"] = arg->ref_rho;
    //arg->alpha = sph_info["arg"]["alpha"];
    sph_info["arg"]["alpha"] = arg->alpha;
    //arg->r = sph_info["arg"]["r"];
    sph_info["arg"]["r"] = arg->r;
    //arg->ptc_num = sph_info["arg"]["ptc_num"];
    if(arg->ptc_num == particle->total) sph_info["arg"]["ptc_num"] = arg->ptc_num;
    else printf("\033[0;32;31m Error in %s:%d\033[m\n",__FILE__,__LINE__);
    //arg->wall_layer = sph_info["arg"]["wall_layer"];
    sph_info["arg"]["wall_layer"] = arg->wall_layer;
    
    //arg->dt = sph_info["time"]["dt"];
    sph_info["time"]["dt"] = arg->dt;
    //arg->sst = sph_info["time"]["sst"];
    sph_info["time"]["sst"] = arg->sst;
    //arg->init_step = sph_info["time"]["init_step"];
    sph_info["time"]["init_step"] = arg->init_step;
    //arg->total_step = sph_info["time"]["total_step"];
    sph_info["time"]["total_step"] = arg->total_step;
    //arg->print_step = sph_info["time"]["print_step"];
    sph_info["time"]["print_step"] = arg->print_step;

    //arg->new_case_flag = sph_info["flag"]["new_case_flag"];
    sph_info["flag"]["new_case_flag"] = arg->new_case_flag;
    //arg->init_impac_flag = sph_info["flag"]["init_impac_flag"];
    sph_info["flag"]["init_impac_flag"] = arg->init_impac_flag;
    //arg->save_last_flag = sph_info["flag"]["save_last_flag"];
    sph_info["flag"]["save_last_flag"] = arg->save_last_flag;

    //rigid->vx = sph_info["rigid"]["vx"];
    sph_info["rigid"]["vx"] = rigid->vx;
    //rigid->vy = sph_info["rigid"]["vy"];
    sph_info["rigid"]["vy"] = rigid->vy;
    //rigid->omega = sph_info["rigid"]["omega"];
    sph_info["rigid"]["omega"] = rigid->omega;
    //rigid->accx = sph_info["rigid"]["accx"];
    sph_info["rigid"]["accx"] = rigid->accx;
    //rigid->accy = sph_info["rigid"]["accy"];
    sph_info["rigid"]["accy"] = rigid->accy;
    //rigid->alpha = sph_info["rigid"]["alpha"];
    sph_info["rigid"]["alpha"] = rigid->alpha;
    //rigid->cogx = sph_info["rigid"]["cogx"];
    sph_info["rigid"]["cogx"] = rigid->cogx;
    //rigid->cogy = sph_info["rigid"]["cogy"];
    sph_info["rigid"]["cogy"] = rigid->cogy;
    //rigid->cog_ptc_id = sph_info["rigid"]["cog_ptc_id"];
    sph_info["rigid"]["cog_ptc_id"] = rigid->cog_ptc_id;
    //rigid->offset_x = sph_info["rigid"]["offset_x"];
    sph_info["rigid"]["offset_x"] = rigid->offset_x;
    //rigid->offset_y = sph_info["rigid"]["offset_y"];
    sph_info["rigid"]["offset_y"] = rigid->offset_y; 
    //rigid->offset_angl = sph_info["rigid"]["offset_angl"];
    sph_info["rigid"]["offset_angl"] = rigid->offset_angl;
    //rigid->mass = sph_info["rigid"]["mass"];
    sph_info["rigid"]["mass"] = rigid->mass;
    //rigid->moi = sph_info["rigid"]["moi"];
    sph_info["rigid"]["moi"] = rigid->moi;
    //rigid->total = sph_info["rigid"]["rigid_num"];
    if(particle->rigid_ptc_num == rigid->total) sph_info["rigid"]["rigid_num"] = rigid->total;
    else printf("\033[0;32;31m Error in %s:%d\033[m\n",__FILE__,__LINE__); 

    f << sph_info << std::endl;
}

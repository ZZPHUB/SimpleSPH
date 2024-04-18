#ifndef __DATASTRUCTURE__
#define __DATASTRUCTURE__

/* Data Structure Declare */
typedef struct 
{
    /* data */
    int *ptc;
    int *count;
}SPH_MESH;

typedef  struct 
{
    /* declare the position,velosity,pressure,density,type of the particle */
    double *x;  //x coordinations of position,iterative
    double *y;  //y coordinations of position,iterative
    double *vx; //x-direction velosity,iterative
    double *vy; //y-direction velosity,iterative
    double *accx;//
    double *accy;//
    double *dif_density;//
    double *pressure;   //pressure of paritcle,non-iterative
    double *density;    //density of particle,iterative
    double *temp_x;     //for predict-correct scheme,non-iterative
    double *temp_y;     //for predict-correct scheme,non-iterative
    double *temp_vx;    //for predict-correct scheme,non-iterative
    double *temp_vy;    //for predict-correct scheme,non-iterative
    double *temp_density;   //for predict-correct scheme,non-iterative
    double *mass;   //mass of particle
    double *w; //sum of kernel value
    int *type; //particle type:0 denote fulid;1 denote rigid;-1 denote dummy particles

    //unsigned int fluid_ptc_num;  //total fluid particle number
    //unsigned int wall_ptc_num;   //total wall particle number
    //unsigned int rigid_ptc_num;  //total rigid particle number
    //unsigned int total; //total particles number
}SPH_PARTICLE;

typedef struct 
{
    /* declare the kernel and differential kernel value of each pair */
    double *w;  //kernel value
    double *dwdx;   //differential kernel value in x-direction
    double *dwdy;   //differential kernel value in y-direction
}SPH_KERNEL;

typedef struct 
{
    /* particle pare generated by NNPS algorithm */
    unsigned int total;
    unsigned int *i;
    unsigned int *j;
}SPH_PAIR;

typedef struct 
{
    /* rigid body kinematics information */
    double vx;  //rigid body x-direction velocity
    double vy;  //rigid body y-direction velocity
    double omega;   //rigid body angular velocity
    double accx;    //rigid body x-direciton acceleration
    double accy;    //rigid body y-direction acceleration
    double alpha;   //rigid body angular acceleration
    double cogx;    //x-direction center of gravity coordinate
    double cogy;    //y-direction center of gravity coordinate 
    double offset_x;    //offset in x direction
    double offset_y;    //offset in y direction
    double offset_angl;  //offset in angular
    double mass;    //rigid body mass 
    double moi;     //rigid body moment of inertia
    int cog_ptc_id;
    //int total;   //rigid body ptc num
}SPH_RIGID;

typedef struct 
{
    double *x;
    double *y;
    double *vx;
    double *vy;
    double *accx;
    double *accy;
    double *rho;
    double *drho;
    double *p;
    double *ptc_w;
    int *type;

    double *temp_x;
    double *temp_y;
    double *temp_vx;
    double *temp_vy;
    double *temp_rho;

    int *mesh;
    int *mesh_count;

    int *pair_i;
    int *pair_j;
    int *pair_count;

    double *pair_w;
    double *dwdx;
    double *dwdy;
}SPH_CUDA;

typedef struct 
{
    double c;   //sound speed
    double g;   //gravity acceleration
    double ref_rho;     //reference density
    double r;       //ptc radius
    double h;   //smoothed length
    double m;   //ptc mass
    double alpha;   //kernel function's para
    double sst;    //single step time
    double dt;  //delta t
    int wall_layer;

    double fluid_x;     //fluid length
    double fluid_y;     //fluid depth
    int fluid_xnum;     //fluid length direction ptc num
    int fluid_ynum;     //fluid depth direction ptc num
    double ptc_dx;      //ptc delta spacing
    
    double mesh_x;    //total domain length
    double mesh_y;    //total domain depth
    double mesh_dx;     //mesh delta spacing
    int mesh_xnum;      //mesh length direction num
    int mesh_ynum;      //mesh depth direction num
    int mesh_num;       //total mesh num
    int mesh_volume;    //single mesh volume

    int init_step;  //inital time step
    int total_step; //total time step
    int print_step; //print time step

     //current process flags
    int new_case_flag;  // if 1 then creat a new case,or continue to run the old case
    int init_impac_flag; //if 1 then run the init step,or run the impac step
    int save_last_flag; //if 1 then save the last step,or donnot save it

    char *case_dir;

    int ptc_num;    //ptc's total num
    int fluid_ptc_num;  //fluid's total num
    int wall_ptc_num;   //wall's total num
    int rigid_ptc_num;   //rigid_ptc_num

    int pair_num;   //pair total num
    int pair_volume; //pear mesh's pair num
    int lock;
    int tmp;   //to count some debug tmp num
}SPH_ARG;

typedef struct 
{
    /* SPH Program Struct */
    SPH_PARTICLE *particle;
    SPH_PAIR *pair;
    SPH_KERNEL *kernel;
    SPH_CUDA *cuda;
    SPH_CUDA *tmp_cuda;
    SPH_ARG *dev_arg;
    SPH_ARG *host_arg;
    SPH_RIGID *host_rigid;
    SPH_RIGID *dev_rigid;
    SPH_MESH *mesh;
}SPH;




#endif
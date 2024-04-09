#include "SPH.cuh"
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <time.h>

using namespace std;



void sph_init(SPH *sph)
{
    SPH_PARTICLE *particle;
    SPH_PAIR *pair;
    SPH_KERNEL *kernel;
    SPH_CUDA temp_cuda;
    SPH_MESH mesh;
    particle = sph->particle;
    pair = sph->pair;
    kernel = sph->kernel;

    particle->fulid_ptc_num = FLUID_PTC_NUM;  //fluid ptc num
    particle->wall_ptc_num = WALL_PTC_NUM;    //wall ptc num
    particle->rigid_ptc_num = ptc_rigid_num(); //rigid ptc num
    //get all of the particle number
    particle->total = particle->fulid_ptc_num+particle->wall_ptc_num+particle->rigid_ptc_num; 

    /************stack is too small,so init data in heap***************/
    //particle data init
    particle->x = (double *)(calloc(particle->total,sizeof(double)));
    particle->y = (double *)(calloc(particle->total,sizeof(double)));
    particle->vx = (double *)(calloc(particle->total,sizeof(double)));
    particle->vy = (double *)(calloc(particle->total,sizeof(double)));
    particle->accx = (double *)(calloc(particle->total,sizeof(double)));
    particle->accy = (double *)(calloc(particle->total,sizeof(double)));
    particle->dif_density = (double *)(calloc(particle->total,sizeof(double)));
    particle->density = (double *)(calloc(particle->total,sizeof(double)));
    particle->temp_x = (double *)(calloc(particle->total,sizeof(double)));
    particle->temp_y = (double *)(calloc(particle->total,sizeof(double)));
    particle->temp_vx = (double *)(calloc(particle->total,sizeof(double)));
    particle->temp_vy = (double *)(calloc(particle->total,sizeof(double)));
    particle->temp_density = (double *)(calloc(particle->total,sizeof(double))); 
    particle->mass = (double *)(calloc(particle->total,sizeof(double))); 
    particle->w = (double *)(calloc(particle->total,sizeof(double)));
    particle->pressure = (double *)(calloc(particle->total,sizeof(double)));
    particle->type = (int *)(calloc(particle->total,sizeof(int)));  

    //kernel data init
    kernel->w = (double *)(calloc(30*particle->total,sizeof(double)));  //this code donnot use kernel value
    kernel->dwdx = (double *)(calloc(30*particle->total,sizeof(double)));
    kernel->dwdy = (double *)(calloc(30*particle->total,sizeof(double)));
   
    //pair data init
    pair->total = 0; 
    pair->i = (unsigned int *)(calloc(30*particle->total,sizeof(unsigned int)));
    pair->j = (unsigned int *)(calloc(30*particle->total,sizeof(unsigned int)));

    //mesh data init
    mesh = (SPH_MESH)(calloc(MESH_DEEPTH_NUM*MESH_LENGTH_NUM*MESH_PTC_NUM,sizeof(int)));
    /*
    for(int i=0;i<MESH_DEEPTH_NUM;i++)
    {
        mesh[i] = (unsigned int **)(calloc(MESH_LENGTH_NUM,sizeof(unsigned int *)));
        for(int j=0;j<MESH_LENGTH_NUM;j++)
        {
            mesh[i][j] = (unsigned int *)(calloc(MESH_PTC_NUM,sizeof(unsigned int)));
        }
    }
    */
    
    sph->mesh = mesh;
    sph->d_time = DELTA_TIME;
    sph->c = ART_SOUND_VEL;
    sph->g = 0.0;
    sph->avg_time = 0.0;


    /*cuda mem alloc*/
    cudaMalloc(&(sph->cuda),sizeof(SPH_CUDA));
    cudaMalloc(&(temp_cuda.x),particle->total*sizeof(double));
    cudaMalloc(&(temp_cuda.y),particle->total*sizeof(double));
    cudaMalloc(&(temp_cuda.temp_x),particle->total*sizeof(double));
    cudaMalloc(&(temp_cuda.temp_y),particle->total*sizeof(double));
    cudaMalloc(&(temp_cuda.vx),particle->total*sizeof(double));
    cudaMalloc(&(temp_cuda.vy),particle->total*sizeof(double));
    cudaMalloc(&(temp_cuda.temp_vx),particle->total*sizeof(double));
    cudaMalloc(&(temp_cuda.temp_vy),particle->total*sizeof(double));
    cudaMalloc(&(temp_cuda.accx),particle->total*sizeof(double));
    cudaMalloc(&(temp_cuda.accy),particle->total*sizeof(double));
    cudaMalloc(&(temp_cuda.rho),particle->total*sizeof(double));
    cudaMalloc(&(temp_cuda.drho),particle->total*sizeof(double));
    cudaMalloc(&(temp_cuda.temp_rho),particle->total*sizeof(double));
    cudaMalloc(&(temp_cuda.p),particle->total*sizeof(double));
    cudaMalloc(&(temp_cuda.type),particle->total*sizeof(int));
    cudaMalloc(&(temp_cuda.ptc_w),particle->total*sizeof(double));

    cudaMalloc(&(temp_cuda.pair_w),32*particle->total*sizeof(double));
    cudaMalloc(&(temp_cuda.dwdx),32*particle->total*sizeof(double));
    cudaMalloc(&(temp_cuda.dwdy),32*particle->total*sizeof(double));
    cudaMalloc(&(temp_cuda.pair_i),32*particle->total*sizeof(int));
    cudaMalloc(&(temp_cuda.pair_j),32*particle->total*sizeof(int));
    cudaMalloc(&(temp_cuda.mesh),MESH_DEEPTH_NUM*MESH_LENGTH_NUM*MESH_PTC_NUM*sizeof(int));
    cudaMalloc(&(temp_cuda.mesh_count),MESH_DEEPTH_NUM*MESH_LENGTH_NUM*sizeof(int));

    cudaMemcpy(sph->cuda,&temp_cuda,sizeof(SPH_CUDA),cudaMemcpyHostToDevice);

    cudaMemcpy(temp_cuda.x, particle->x, particle->total*sizeof(double), cudaMemcpyHostToDevice); 
    cudaMemcpy(temp_cuda.y, particle->y, particle->total*sizeof(double), cudaMemcpyHostToDevice); 
    cudaMemcpy(temp_cuda.vx, particle->vx, particle->total*sizeof(double), cudaMemcpyHostToDevice); 
    cudaMemcpy(temp_cuda.vy, particle->vy, particle->total*sizeof(double), cudaMemcpyHostToDevice); 
    cudaMemcpy(temp_cuda.type, particle->type, particle->total*sizeof(int), cudaMemcpyHostToDevice); 
    cudaMemcpy(temp_cuda.rho, particle->density, particle->total*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(temp_cuda.accx, particle->accx, particle->total*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(temp_cuda.accy, particle->accx, particle->total*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(temp_cuda.drho, particle->dif_density, particle->total*sizeof(double), cudaMemcpyHostToDevice);
    
    cudaMemset(temp_cuda.p,0,particle->total*sizeof(double));
    cudaMemset(temp_cuda.temp_x,0,particle->total*sizeof(double));
    cudaMemset(temp_cuda.temp_y,0,particle->total*sizeof(double));
    cudaMemset(temp_cuda.temp_vx,0,particle->total*sizeof(double));
    cudaMemset(temp_cuda.temp_vy,0,particle->total*sizeof(double));
    cudaMemset(temp_cuda.temp_rho,0,particle->total*sizeof(double));
    cudaMemset(temp_cuda.ptc_w,0,particle->total*sizeof(double));

    cudaMemset(temp_cuda.pair_w,0,32*particle->total*sizeof(double));
    cudaMemset(temp_cuda.dwdx,0,32*particle->total*sizeof(double));
    cudaMemset(temp_cuda.dwdy,0,32*particle->total*sizeof(double));
    cudaMemset(temp_cuda.pair_i,0,32*particle->total*sizeof(int));
    cudaMemset(temp_cuda.pair_j,0,32*particle->total*sizeof(int));
    cudaMemset(temp_cuda.mesh,0,MESH_DEEPTH_NUM*MESH_LENGTH_NUM*MESH_PTC_NUM*sizeof(int));
    cudaMemset(temp_cuda.mesh_count,0,MESH_DEEPTH_NUM*MESH_LENGTH_NUM*sizeof(int));


    cout << "run a new case or an old case(press 1 for new,0 for old)" << endl;
    cin >> sph->new_case_flag;

    if(sph->new_case_flag == 1)
    {
        sph->current_step = 0;
        sph->total_step = INIT_TIME_STEP;
    }
    else if(sph->new_case_flag == 0)
    {
        cout << "the sph current time step is: " << endl;
        cin >> sph->current_step;
        cout << "the total sph time step is: " << endl;
        cin >> sph->total_step;
    }

    cout << "run a init case or a dynamic case(press 1 for init,0 for dynamic)" << endl;
    cin >> sph->init_impac_flag;
    
    cout << "save the last step or not(press 1 ta save,0 for not)" << endl;
    cin >> sph->save_last_flag;

    ptc_generate(sph);
    ptc_init(sph);
}
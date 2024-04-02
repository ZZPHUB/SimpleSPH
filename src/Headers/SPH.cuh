#ifndef __SPH_H__
#define __SPH_H__

#include "DataStructure.cuh"

#define VX 0
#define VY 1
#define OMEGA 2
#define ACCX 3
#define ACCY 4
#define R_ALPHA 5
#define COGX 6
#define COGY 7
#define MASS 8
#define MOI 9
/* for fluid particles,they x-direction arrange betweent(0,12-2*PTC_SPACING) */
/* for fluid particles,they y-direction arrange betweent(0,5-2*PTC_SPACING) */
/* for dummy particles,they x-direction arrange betweent (-2*PTC_SPACING,0)and(12-2*PTC_SPACOING,12) */
/* for dummy particles,they y-direction arrange betweent (-2*PTC_SACING,5) */

/* Define Symbols Here*/
#define FLUID_DOMAIN_LENGTH 1.0   //computational domain length
#define FLUID_DOMAIN_DEEPTH 1.0     //computational domain deepth
#define TOL_DOMAIN_LENGTH (FLUID_DOMAIN_LENGTH)
#define TOL_DOMAIN_DEEPTH (FLUID_DOMAIN_DEEPTH*1.50)

#define PTC_SPACING 0.005    //spacing between particles
#define PTC_SML (1.0005*PTC_SPACING)    //smoothed leng of particles
#define PTC_REGION_RADIUS (2.0*PTC_SML) //particle support region radius,which is PTC_SML*2

#define FLUID_LENGTH_NUM (int)(FLUID_DOMAIN_LENGTH/PTC_SPACING+1-8)
#define FLUID_DEEPTH_NUM (int)(FLUID_DOMAIN_DEEPTH/PTC_SPACING+1-4-12)

#define FLUID_PTC_NUM (FLUID_LENGTH_NUM)*(FLUID_DEEPTH_NUM)  //fluid particles number
//#define SOLID_PTC_NUM  0 //rigid body particles number
#define WALL_PTC_NUM (int)((FLUID_DOMAIN_DEEPTH/PTC_SPACING+1)*(FLUID_DOMAIN_LENGTH/PTC_SPACING+1)-FLUID_LENGTH_NUM*(FLUID_DEEPTH_NUM+12))  //virtual paritcles number
//#define PTC_TOL_NUM (FLUID_PTC_NUM+SOLID_PTC_NUM+VIRTUAL_PTC_NUM)   //total particles

#define MESH_SPACING PTC_REGION_RADIUS //mesh spacing
#define MESH_LENGTH_NUM (int)(TOL_DOMAIN_LENGTH/MESH_SPACING+1) //length-direction mesh number
#define MESH_DEEPTH_NUM (int)(TOL_DOMAIN_DEEPTH/MESH_SPACING+1) //deepth-direction mesh number
#define MESH_LENGTH_NUM_CUDA __double2int_rz(TOL_DOMAIN_LENGTH/MESH_SPACING+1) //length-direction mesh number
#define MESH_DEEPTH_NUM_CUDA __double2int_rz(TOL_DOMAIN_DEEPTH/MESH_SPACING+1) //deepth-direction mesh number
#define MESH_TOL_NUM (MESH_LENGTH_NUM*MESH_DEEPTH_NUM) //total mesh number
#define MESH_PTC_NUM 33 //per mesh grid contain max paticle num

#define REF_DENSITY 1000.0    //reference density for eos
#define GRAVITY_ACC 9.80     //gravity acceleration defin here
#define ART_SOUND_VEL (10*sqrt(GRAVITY_ACC*FLUID_DOMAIN_DEEPTH))  //art_sound_velosity,it's $p=10 \times (\rho-\rho_{ref})$

#define MU 0.001 //define the constant viscous

#define TH_NUM 8 //paraller threads
#define INIT_TIME_STEP 80000
#define DELTA_TIME 0.00002
#define PRINT_TIME_STEP 400  //every 50 time step to print

#define PARA (0x01|0x02|0x04|0x08)
//0x01----------->density
//0x02----------->pressure
//0x04----------->velosity
//0x08----------->acceleration

#define PTC_MASS (REF_DENSITY*pow(PTC_SPACING,2))   //every particle's mass is constant

#define T_START current_time=time(NULL);{ //get start time
#define T_END(a) }cout << a << " use time is " << time(NULL)-current_time << " s" << endl;  //echo used time

// two paritcles distance
#define PTC_DISTANCE(a,b) (sqrt(pow(particle->x[a]-particle->x[b],2)+pow(particle->y[a]-particle->y[b],2)))

#define PI 3.14159265358
#define ALPHA 15.0/(7*PI*pow(PTC_SML,2))

#define LINEAR_EOS

/*
extern __constant__ __device__ int c;
extern __constant__ __device__ int rho_0;
extern __constant__ __device__ int mesh_lnum;
extern __constant__ __device__ int mesh_dnum;
extern __constant__ __device__ int mesh_pnum;
extern __constant__ __device__ int mesh_spacing;
*/
__device__ int dev_mesh_tnum=MESH_DEEPTH_NUM*MESH_LENGTH_NUM;
__device__ int dev_mesh_lnum=MESH_LENGTH_NUM;
__device__ int dev_mesh_dnum=MESH_DEEPTH_NUM;
__device__ double dev_mesh_spacing=MESH_SPACING;
__device__ double dev_a=15.0/(7*PI*PTC_SML*PTC_SML);
__device__ double dev_h=PTC_SML;
__device__ double dev_c=10.0*(GRAVITY_ACC*FLUID_DOMAIN_DEEPTH)*(GRAVITY_ACC*FLUID_DOMAIN_DEEPTH);
__device__ double dev_m=REF_DENSITY*PTC_SPACING*PTC_SPACING;
__device__ double dev_dt=DELTA_TIME;


/* Headers Include Here*/

#include <omp.h>
#include "Lib.cuh"
#include "PreProcess.cuh"
#include "PostProcess.cuh"
#include "Equations.cuh"

/* Extern Functions Here*/
void sph_init(SPH *);
void sph_free(SPH *);
void sph_time_integral(SPH *);
void sph_rigid_integral(SPH *sph);
void ptc_dummy(SPH *);
__global__ void sph_predict_cuda(double *x,double *y,double *temp_x,double *temp_y,double *vx,double *vy,double *temp_vx,double *temp_vy,double *accx,double *accy,double *rho,double *temp_rho,double *drho,double *p,int ptc_num);
__global__ void sph_correct_cuda(double *x,double *y,double *temp_x,double *temp_y,double *vx,double *vy,double *temp_vx,double *temp_vy,double *accx,double *accy,double *rho,double *temp_rho,double *drho,double *p,int ptc_num);

//void ptc_density_correct(SPH *);


#endif

#include "SPH.cuh"

__global__ void check_ptc(SPH_CUDA *cuda,SPH_ARG *arg)
{
    const int id = threadIdx.x + threadIdx.y*blockDim.x;
    if(id >= arg->ptc_num)return;
    printf("%lf %lf\n",cuda->x[id],cuda->y[id]);
}

__global__ void check_pair(SPH_ARG *arg)
{
    printf("the pair num is:%d\n",arg->pair_num);
}

__global__ void check_mesh(SPH_CUDA *cuda,SPH_ARG *arg)
{
    const int mesh_id = blockIdx.x + blockIdx.y* gridDim.x;
    if(cuda->mesh_count[mesh_id]!=0)
    {
        printf("%d %d\n",mesh_id,cuda->mesh_count[mesh_id]);
    }
    /*
    if(cuda->mesh_count[mesh_id] != 0)
    {
        printf("mesh id is:%d ptc in mesh is:%d they are:",mesh_id,cuda->mesh_count[mesh_id]);
        for(int i=0;i<cuda->mesh_count[mesh_id];i++)
        {
            printf("%d",cuda->mesh[mesh_id+i*arg->mesh_num]);
        }
        printf("\n");
    }*/
}

int main(void)
{
    SPH_PARTICLE particle;
    SPH_KERNEL kernel;
    SPH_PAIR pair;
    SPH_RIGID wedge;
    SPH_MESH mesh = NULL;
    SPH_ARG arg;
    SPH sph;
    sph.particle = &particle;
    sph.kernel = &kernel;
    sph.pair = &pair;
    sph.host_rigid = &wedge;
    sph.host_arg = &arg;
    sph.mesh = mesh;

    cudaSetDevice(0);
    sph_init(&sph); 

    //define the seed for ptc data structure
    dim3 ptc_block(256);
    dim3 ptc_grid((int)(sph.particle->total/256)+1);
    //define the seed for mesh data structure
    dim3 mesh_block(32,32);
    dim3 mesh_grid(MESH_LENGTH_NUM,MESH_DEEPTH_NUM);
    //define the seed for pair data structre
    dim3 pair_block(512);
    dim3 pair_grid((int)(sph.particle->total/16)+1);

    //check_ptc<<<ptc_grid,ptc_block>>>(sph.cuda,sph.dev_arg);
    //cudaDeviceSynchronize();
    sph_mesh_cuda<<<ptc_grid,ptc_block>>>(sph.cuda,sph.dev_arg);
    cudaDeviceSynchronize();
    //check_mesh<<<mesh_grid,1>>>(sph.cuda,sph.dev_arg);
    //cudaDeviceSynchronize();
    //sph_nnps_cuda<<<mesh_grid,mesh_block>>>(sph.cuda,sph.dev_arg,sph.dev_rigid);
    //cudaDeviceSynchronize();
    //check_pair<<<1,1>>>(sph.dev_arg);
    //cudaDeviceSynchronize();

    int *host_mesh;
    int *host_mesh_count;
    SPH_CUDA cuda;
    cudaMemcpy(&cuda,sph.cuda,sizeof(SPH_CUDA),cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();

    host_mesh = (int *)malloc(sizeof(int)*sph.host_arg->mesh_num*sph.host_arg->mesh_volume);
    host_mesh_count = (int *)mall0c(sizeof(int)*sph.host_arg->mesh_num);

    cudaMemcpy(host_mesh,cuda.mesh,sizeof(int)*sph.host_arg->mesh_num*sph.host_arg->mesh_volume);
    cudaDeviceSynchronize();
    cudaMemcpy(host_mesh_count,cuda.mesh_count,sizeof(int)*sph.host->mesh_num);
    cudaDeviceSynchronize();

    for(int i=0;i<sph.host_arg->mesh_num;i++)
    {
        printf("mesh id is:%d mesh num is:%d\n",i,host_mesh_count[i]);
    }

    sph_free(&sph);
    cudaDeviceReset();
    return 0;
}


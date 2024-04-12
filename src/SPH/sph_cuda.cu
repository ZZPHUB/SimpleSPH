#include "SPH.cuh"

__global__ void check_ptc(SPH_CUDA *cuda,SPH_ARG *arg)
{
    const int id = threadIdx.x + threadIdx.y*blockDim.x;
    if(id >= arg->ptc_num)return;
    printf("%lf %lf\n",cuda->x[id],cuda->y[id]);
}

__global__ void check_pair(SPH_CUDA *cuda,SPH_ARG *arg)
{
    const int id = threadIdx.x + blockIdx.x * blockDim.x;
    if(id >= arg->pair_num) return;
    if(id == 0)printf("the pair num is:%d\n",arg->pair_num);
    for(int i=0;i<arg->pair_num;i++)
    {
        if(cuda->pair_i[id] == cuda->pair_i[i] && cuda->pair_j[id]==cuda->pair_j[i] && id!=i)
        {
            if(cuda->pair_i[id]!=0 && cuda->pair_j[id]!=0)
            {
                printf("type1 index_1:%d index_2:%d pair_i:%d pair_j:%d\n",id,i,cuda->pair_i[id],cuda->pair_j[id]);
            }
            atomicAdd(&(arg->tmp),1);
        }
        else if(cuda->pair_i[id] == cuda->pair_j[i] && cuda->pair_j[id]==cuda->pair_i[i])
        {
            if(id == i)
            {
                //printf("type2 nnps error !!\n");
                if(cuda->pair_i[id]!=0 && cuda->pair_j[id]!=0)
                {
                    printf("type2 index_1:%d index_2:%d pair_i:%d pair_j:%d\n",id,i,cuda->pair_i[id],cuda->pair_j[id]);
                }
            }
            else 
            {
                //printf("here is same pair\n");
                if(cuda->pair_i[id]!=0 && cuda->pair_j[id]!=0)
                {
                    printf("type3 index_1:%d index_2:%d pair_i:%d pair_j:%d\n",id,i,cuda->pair_i[id],cuda->pair_j[id]);
                }
            }
            atomicAdd(&(arg->tmp),1);
        }
    }
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

    int *host_mesh;
    int *host_mesh_count;
    SPH_CUDA cuda;
    SPH_ARG tmp_arg;
    cudaMemcpy(&cuda,sph.cuda,sizeof(SPH_CUDA),cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();

    for(int i=0;i<1;i++)
    {
        //printf("current step is:%d\n",i);
        //check_ptc<<<ptc_grid,ptc_block>>>(sph.cuda,sph.dev_arg);
        //cudaDeviceSynchronize();
        sph_mesh_cuda<<<ptc_grid,ptc_block>>>(sph.cuda,sph.dev_arg);
        cudaDeviceSynchronize();
        //check_mesh<<<mesh_grid,1>>>(sph.cuda,sph.dev_arg);
        //cudaDeviceSynchronize();
        //sph_nnps_cuda<<<mesh_grid,mesh_block>>>(sph.cuda,sph.dev_arg,sph.dev_rigid);
        //cudaDeviceSynchronize();
        //check_pair<<<(int)(250000/1024)+1,1024>>>(sph.cuda,sph.dev_arg);
        //cudaDeviceSynchronize();

        //cudaMemcpy(&tmp_arg,sph.dev_arg,sizeof(SPH_ARG),cudaMemcpyDeviceToHost);
        //printf("the total same pair num is:%d \n",tmp_arg.tmp);
        
        host_mesh = (int *)malloc(sizeof(int)*sph.host_arg->mesh_num*sph.host_arg->mesh_volume);
        host_mesh_count = (int *)malloc(sizeof(int)*sph.host_arg->mesh_num);

        cudaMemcpy(host_mesh,cuda.mesh,sizeof(int)*sph.host_arg->mesh_num*sph.host_arg->mesh_volume,cudaMemcpyDeviceToHost);
        cudaDeviceSynchronize();
        cudaMemcpy(host_mesh_count,cuda.mesh_count,sizeof(int)*sph.host_arg->mesh_num,cudaMemcpyDeviceToHost);
        cudaDeviceSynchronize();

        for(int j=0;j<sph.host_arg->mesh_num;j++)
        {
            //if(host_mesh_count[j]!=0) printf("error!!!!!!\n");
            printf("mesh id is:%d mesh num is:%d they are:",j,host_mesh_count[j]);
            for(int k=0;k<host_mesh_count[j];k++)
            {
                printf("%d,",host_mesh[j+k*sph.host_arg->mesh_num]);
            }
            printf("\n");
        }
    }

    /*
    for(int i=0;i<sph.host_arg->mesh_num;i++)
    {
        printf("mesh id is:%d mesh num is:%d\n",i,host_mesh_count[i]);
    }*/

    sph_free(&sph);
    cudaDeviceReset();
    return 0;
}


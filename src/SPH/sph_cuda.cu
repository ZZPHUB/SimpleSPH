#include "SPH.cuh"

__global__ void check_ptc(SPH_CUDA *cuda,SPH_ARG *arg)
{
    const int id = threadIdx.x + threadIdx.y*blockDim.x;
    if(id >= arg->ptc_num)return;
    printf("%lf %lf\n",cuda->x[id],cuda->y[id]);
}
__global__ void check_nnps(SPH_CUDA *cuda,SPH_ARG *arg,int *pair_i,int *pair_j,int p_num)
{
    const int mesh_id = blockIdx.x + blockIdx.y * gridDim.x;
    int id = 0;
    if( threadIdx.x >= cuda->pair_count[mesh_id]) return;
    else id = mesh_id*arg->pair_volume + threadIdx.x;
    int j=0;
    for(int i=0;i<p_num;i++)
    {   
        if(cuda->pair_i[id] == pair_i[i] && cuda->pair_j[id] == pair_j[i])
        {
            j++;
        }
        else if(cuda->pair_i[id] == pair_j[i] && cuda->pair_j[id] == pair_i[i])
        {
            j++;
        }
        else if(cuda->pair_i[id] == cuda->pair_j[id])
        {
            j++;
        }
    }
    if(j!=1)
        {
            printf("error the same pair num is:%d\n",j);
        }
}

__global__ void check_pair(SPH_CUDA *cuda,SPH_ARG *arg)
{
    //double dx=0.0;
    //double dy=0.0;
    //double q=0.0;
    //const int id = threadIdx.x + blockIdx.x * blockDim.x;
    //if(id >= arg->pair_num) return;
    const int mesh_id = blockIdx.x + blockIdx.y * gridDim.x;
    int id = 0;
    int tmp_id = 0;
    if(threadIdx.x >= cuda->pair_count[mesh_id]) return;
    else id = mesh_id*arg->pair_volume + threadIdx.x;
    for(int i=0;i<arg->mesh_num;i++)
    {
        for(int j=0;j<cuda->pair_count[i];j++)
        {
            tmp_id = i*arg->pair_volume + j;
            if(id != tmp_id)
            {
                if(cuda->pair_i[id] == cuda->pair_j[tmp_id] && cuda->pair_j[id] == cuda->pair_j[tmp_id])
                {
                    printf("type 1 errors\n");
                }
                else if(cuda->pair_i[id] == cuda->pair_j[tmp_id] && cuda->pair_j[id] == cuda->pair_i[tmp_id])
                {
                    printf("type 2 errors\n");
                }
            }
        }
    }

    /*
    dx = cuda->x[cuda->pair_i[id]] - cuda->x[cuda->pair_j[id]];
    dy = cuda->y[cuda->pair_i[id]] - cuda->y[cuda->pair_j[id]];
    q = sqrt(dx*dx+dy*dy)/arg->h;
    if(q > 2.0) printf("error !!!\n");*/

    
    //if(id == 0)printf("the pair num is:%d\n",arg->pair_num);
    /*for(int i=0;i<arg->pair_num;i++)
    {
        if(cuda->pair_i[id] == cuda->pair_i[i] && cuda->pair_j[id]==cuda->pair_j[i] && id!=i)
        {
            if(cuda->pair_i[id]!=0 && cuda->pair_j[id]!=0)
            {
                //printf("type1 index_1:%d index_2:%d pair_i:%d pair_j:%d\n",id,i,cuda->pair_i[id],cuda->pair_j[id]);
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
                    //printf("type2 index_1:%d index_2:%d pair_i:%d pair_j:%d\n",id,i,cuda->pair_i[id],cuda->pair_j[id]);
                }
            }
            else 
            {
                //printf("here is same pair\n");
                if(cuda->pair_i[id]!=0 && cuda->pair_j[id]!=0)
                {
                    //printf("type3 index_1:%d index_2:%d pair_i:%d pair_j:%d\n",id,i,cuda->pair_i[id],cuda->pair_j[id]);
                }
            }
            atomicAdd(&(arg->tmp),1);
        }
    }*/
}

__global__ void check_mesh(SPH_CUDA *cuda,SPH_ARG *arg)
{
    int mid = 0;
    int id = 0;
    const int mesh_id = blockIdx.x + blockIdx.y* gridDim.x;
    atomicAdd(&arg->tmp,cuda->pair_count[mesh_id]);
    cuda->pair_count[mesh_id]=0;
    /*
    for(int i=0;i<cuda->mesh_count[mesh_id];i++)
    {
        id = cuda->mesh[i*arg->mesh_num + mesh_id];
        if(cuda->y[id] < arg->domain_y && cuda->y[id] >= 0.0)
        {
            mid = __double2int_rz(cuda->y[id]/arg->mesh_dx)*arg->mesh_xnum;
        }
        else
        {
            mid = (arg->mesh_ynum - 1)*arg->mesh_xnum;
        }
        if(cuda->x[id] < arg->domain_x && cuda->x[id] >= 0.0)
        {
            mid += __double2int_rz(cuda->x[id]/arg->mesh_dx);
            if(mid == mesh_id) atomicAdd(&(arg->tmp),1);
            //if(mid != mesh_id) printf("mid:%d mesh_id:%d id:%d x:%lf y:%lf\n",mid,mesh_id,id,cuda->x[id],cuda->y[id]);

        }
        else
        {
            mid += arg->mesh_xnum - 1;
            if(mid == mesh_id) atomicAdd(&(arg->tmp),1);
            //if(mid != mesh_id) printf("mid:%d mesh_id:%d id:%d x:%lf y:%lf\n",mid,mesh_id,id,cuda->x[id],cuda->y[id]);
        }
    }
    cuda->mesh_count[mesh_id] = 0;*/
    /*
    if(cuda->mesh_count[mesh_id]!=0)
    {
        printf("%d %d\n",mesh_id,cuda->mesh_count[mesh_id]);
    }*/
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
    SPH_MESH mesh;
    SPH_ARG arg;
    SPH sph;
    sph.particle = &particle;
    sph.kernel = &kernel;
    sph.pair = &pair;
    sph.host_rigid = &wedge;
    sph.host_arg = &arg;
    sph.mesh = &mesh;

    cudaSetDevice(0);
    sph_init(&sph); 

    //define the seed for ptc data structure
    dim3 ptc_block(256);
    dim3 ptc_grid((int)(sph.particle->total/256)+1);
    //define the seed for mesh data structure
    dim3 mesh_block(32,32);
    dim3 mesh_grid(MESH_LENGTH_NUM,MESH_DEEPTH_NUM);
    //define the seed for pair data structre
    dim3 pair_block(128);
    dim3 pair_grid(MESH_LENGTH_NUM,MESH_DEEPTH_NUM);


    SPH_CUDA cuda;
    SPH_ARG tmp_arg;
    cudaMemcpy(&cuda,sph.cuda,sizeof(SPH_CUDA),cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();
    //int *host_pair_count = (int *)calloc(sph.host_arg->mesh_num,sizeof(int));
    //int id = 0;
    int *cpu_pair_i;
    int *cpu_pair_j;
    cudaMalloc(&cpu_pair_i,sizeof(int)*32*sph.particle->total);
    cudaMalloc(&cpu_pair_j,sizeof(int)*32*sph.particle->total);

    for(int i=0;i<100;i++)
    {
        printf("current step is:%d\n",i);
        //check_ptc<<<ptc_grid,ptc_block>>>(sph.cuda,sph.dev_arg);
        //cudaDeviceSynchronize();
        sph_mesh_cuda<<<ptc_grid,ptc_block>>>(sph.cuda,sph.dev_arg);
        cudaDeviceSynchronize();
        //check_mesh<<<mesh_grid,1>>>(sph.cuda,sph.dev_arg);
        //cudaDeviceSynchronize();
        
        cudaMemcpy(sph.mesh->ptc,cuda.mesh,sizeof(int)*sph.host_arg->mesh_num*sph.host_arg->mesh_volume,cudaMemcpyDeviceToHost);
        cudaDeviceSynchronize();
        cudaMemcpy(sph.mesh->count,cuda.mesh_count,sizeof(int)*sph.host_arg->mesh_num,cudaMemcpyDeviceToHost);
        cudaDeviceSynchronize();
        sph_nnps_cpu(&sph);

        sph_nnps_cuda<<<mesh_grid,mesh_block>>>(sph.cuda,sph.dev_arg,sph.dev_rigid);
        cudaDeviceSynchronize();

        
        /*
        for(int i=0;i<sph.host_arg->mesh_num;i++)
        {
            printf("%d\n",sph.mesh->count[i]);
        }*/
        
        cudaMemcpy(cpu_pair_i,sph.pair->i,sizeof(int)*32*sph.particle->total,cudaMemcpyHostToDevice);
        cudaDeviceSynchronize();
        cudaMemcpy(cpu_pair_j,sph.pair->j,sizeof(int)*32*sph.particle->total,cudaMemcpyHostToDevice);
        cudaDeviceSynchronize();

        check_nnps<<<pair_grid,pair_block>>>(sph.cuda,sph.dev_arg,cpu_pair_i,cpu_pair_j,sph.host_arg->pair_num);
        cudaDeviceSynchronize();
        check_pair<<<pair_grid,pair_block>>>(sph.cuda,sph.dev_arg);
        cudaDeviceSynchronize();
        check_mesh<<<mesh_grid,1>>>(sph.cuda,sph.dev_arg);
        cudaDeviceSynchronize();


        cudaError_t sph_error = cudaGetLastError();
        printf("%s\n",cudaGetErrorName(sph_error));
    /*
        cudaMemcpy(sph.pair->i,cuda.pair_i,sizeof(int)*32*sph.particle->total,cudaMemcpyDeviceToHost);
        cudaDeviceSynchronize();
        cudaMemcpy(sph.pair->j,cuda.pair_j,sizeof(int)*32*sph.particle->total,cudaMemcpyDeviceToHost);
        cudaDeviceSynchronize();
        cudaMemcpy(&tmp_arg,sph.dev_arg,sizeof(SPH_ARG),cudaMemcpyDeviceToHost);
        cudaDeviceSynchronize();
        cudaMemcpy(host_pair_count,cuda.pair_count,sph.host_arg->mesh_num*sizeof(int),cudaMemcpyDeviceToHost);
        cudaDeviceSynchronize();
        for(int i=0;i<sph.host_arg->mesh_num;i++)
        {
            for(int j=0;j<host_pair_count[i];j++)
            {
                id = i*sph.host_arg->pair_volume+j;
                printf("id:%d i:%d j:%d\n",id,sph.pair->i[id],sph.pair->j[id]); 
            }
        }
        */
        
        //printf("the total same pair num is:%d \n",tmp_arg.tmp);
        
        /*
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
        }*/
        
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


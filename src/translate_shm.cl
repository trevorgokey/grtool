
#pragma OPENCL EXTENSION cl_khr_fp64 : enable

#define SHR_SIZE (TRANS_SHR_SIZE*3)

#ifdef USE_DOUBLE
typedef double data_size;
#else
typedef float data_size;
#endif

void reduce3(__private double* sum, __local double* shm, __private size_t *a)
{
	size_t x = get_local_id(0);
	size_t x_len = get_local_size(0);
	size_t blk = get_global_id(1);
	const uint r = (x_len % (SHR_SIZE/3)) > 0;
	size_t i,j;
	const size_t shr_idx = x % (SHR_SIZE/3);
	
	char write = 0;
	
	//zero out shr mem
	barrier(CLK_LOCAL_MEM_FENCE);
	for(i = x; i < SHR_SIZE; i+=x_len)
		shm[i] = 0.0;

	
	//condense into shr mem
	for(j = 1; j <= (x_len/(SHR_SIZE/3)) + r ;  j++ )
	{
		barrier(CLK_LOCAL_MEM_FENCE);
		if(!write && x < j*(SHR_SIZE/3) && (write = 1))
		{
			shm[3*(shr_idx)    ] += sum[0];
			shm[3*(shr_idx) + 1] += sum[1];
			shm[3*(shr_idx) + 2] += sum[2];
		}
	}
	//reduce the shr mem
	for(i = SHR_SIZE/3 >> 1 ; i > 0; i >>= 1)
	{
		barrier(CLK_LOCAL_MEM_FENCE);
		if(x < i)
		{
			shm[3*x    ] += shm[3*(x + i)    ];
			shm[3*x + 1] += shm[3*(x + i) + 1];
			shm[3*x + 2] += shm[3*(x + i) + 2];
		}
	}
	barrier(CLK_LOCAL_MEM_FENCE);
	sum[0] = shm[0];
	sum[1] = shm[1];
	sum[2] = shm[2];
}

__kernel void translate(	__global size_t *sizes,
						__global data_size *master,
						__global data_size *mask,
						__global data_size *translate)
{
	size_t tid = get_global_id(0);
	size_t blk = get_global_id(1);
	__local double shm[SHR_SIZE];

    size_t j = 0,i = 0;
    double offset[] = {0.0,0.0,0.0};
    size_t st = sizes[0];
    size_t anum = sizes[0];
    double fanum = 1.0/(data_size)anum;
    size_t mnum = sizes[1];
//	uint ms = mnum*st*3;
//	uint sz = anum*st*3;
    size_t idx;


	// loop for shr mem
    for(j = tid*3 ; j < anum*3; j+=get_local_size(0)*3)
    {	
		
        idx = ( blk * 3 * anum ) + ( j );
        offset[0] += ( mask[idx   ] / anum);
        offset[1] += ( mask[idx+ 1] / anum);
        offset[2] += ( mask[idx+ 2] / anum);
		//printf("TID %zu BLK %zu PIECE IS %16.16f add from %16.16f\n",tid,blk, offset[0],( mask[idx   ] * fanum ));
    }
	reduce3(offset,shm,&anum);
    //apply to mask for future calculations
    for(j = tid*3 ; j < anum*3 ; j+=get_local_size(0)*3)
    {	
        idx = ( blk * 3 * anum ) + j;
        master[idx]   -=  offset[0];
        master[idx+1] -=  offset[1];
        master[idx+2] -=  offset[2];
		//master[idx] = 888.0;
    }
    barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);
	if(tid == 0)
	{
		idx = (blk * 3);
    	translate[idx  ] = offset[0] ;
    	translate[idx+1] = offset[1] ;
    	translate[idx+2] = offset[2] ;
	}
    //barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);
//	if((!tid))
//	{
//		printf("DIMS ARE %zu %zu and offset is %f %f %f\n",st,anum,offset[0],offset[1],offset[2]);
//	}
    //apply to data
    //for(j = 0; j < anum; j++)
    //{
    //    idx = ( blk * 3 * anum ) + ( j * 3 ) + tid;
    //    master[idx] = master[idx] - offset;
    //   // master[idx] = 99.0; //master[idx] - offset;
    //}
}

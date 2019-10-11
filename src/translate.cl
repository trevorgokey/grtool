#pragma OPENCL EXTENSION cl_khr_fp64 : enable

#if USE_DOUBLE
typedef double data_size;
#else
typedef float data_size;
#endif


//__kernel void center(	__global uint *sizes,
//						__global double *master,
//						__global double *mask,
//						__global double *translate) 
__kernel void translate(	__global size_t *sizes,
						__global data_size *master,
						__global data_size *mask,
						__global data_size *translate)
{
	size_t tid = get_global_id(0);
	size_t blk = get_global_id(1);
		
    size_t j = 0;
    data_size offset = 0.0;
	size_t st = sizes[0];
    size_t anum = sizes[0];
	data_size fanum = (data_size)anum;
    size_t mnum = sizes[1];
//	uint ms = mnum*st*3;
//	uint sz = anum*st*3;
    size_t idx;
    for(j = 0 ; j < anum; j++)
    {	
        idx = ( blk * 3 * anum ) + ( j * 3 ) + tid;
        offset += ( mask[idx] / fanum);
    }
    //apply to mask for future calculations
    for(j = 0 ; j < anum; j++)
    {	
        idx = ( blk * 3 * anum ) + ( j * 3 ) + tid;
        mask[idx] = mask[idx] - offset;
		//master[idx] = 888.0;
    }
    idx = (blk * 3 + tid);
    translate[idx] = offset ;
	return;
	//if((!tid) && (!blk))
	//{
	//	printf("DIMS ARE %zu %d and offset is %f\n",st,anum,offset);
	//}
    //apply to data
    for(j = 0; j < anum; j++)
    {
        idx = ( blk * 3 * anum ) + ( j * 3 ) + tid;
        //master[idx] = master[idx] - offset;
        master[idx] = 99.0; //master[idx] - offset;
    }
}

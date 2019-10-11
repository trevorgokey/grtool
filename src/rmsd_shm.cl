
#pragma OPENCL EXTENSION cl_khr_fp64 : enable

#define SHR_SIZE 128
#ifdef USE_DOUBLE
typedef double data_size;
#else
typedef float data_size;
#endif
//#define NT get_global_size(0)
//#define TID get_global_id(0)

void reduce(__private double* sum, __local double* shm)
{
	size_t x = get_local_id(0);
	size_t x_len = get_local_size(0);
	size_t blk = get_global_id(1);
	const size_t ss = SHR_SIZE /3 ;
	const uint r = (x_len % (ss)) > 0;
	size_t i,j;
	const size_t shr_idx = x % (ss);
	
	char write = 0;
	
	//printf("MY ID IS %zu MY LEN IS %zu\n",x,x_len);
	barrier(CLK_LOCAL_MEM_FENCE);
	//zero out shr mem
	for(i = x; i < ss; i+=x_len)
		shm[i] = 0.0;

	
	//condense into shr mem
	for(j = 1; j <= (x_len/(ss)) + r ;  j++ )
	{
		barrier(CLK_LOCAL_MEM_FENCE);
		if(!write && x < j*(ss) && (write = 1))
		{
			shm[(shr_idx)    ] += sum[0];
		}
	}
	//reduce the shr mem
	for(i = ss >> 1 ; i > 0; i >>= 1)
	{
		barrier(CLK_LOCAL_MEM_FENCE);
		if(x < i)
		{
			shm[x    ] += shm[(x + i)    ];
		}
	}
	barrier(CLK_LOCAL_MEM_FENCE);
	sum[0] = shm[0];
}

//double reduce(__private size_t len, __private double val, __local double *ssum)
//{
//	size_t j;
//	int x = 1;
//	for(j = TID; j < SSUM_SIZE ; j+=NT)
//			ssum[j] = 0.0;
//	barrier(CLK_LOCAL_MEM_FENCE);
//	for(j = 0; j < NT; j++)
//	{
//		if(TID == j)
//			ssum[0] += 1.0;
//		barrier(CLK_GLOBAL_MEM_FENCE|CLK_LOCAL_MEM_FENCE);
//	}
//	if(TID == 0)
//		printf("my TID is FOR FRAME IS %f\n",ssum[0]);
//	return 1.0;
//	//for(int p = 0; p < ((len/SSUM_SIZE) + (len%SSUM_SIZE > 0))  ; p++)
//	for(j = 0; j < len  ; j+=SSUM_SIZE)
//	{
//		if(TID == 0)
//			printf("my TID is FOR FRAME IS %f\n",ssum[0]);
//		if(x && TID < (j + SSUM_SIZE))
//		{
//			x = 0;
//			ssum[TID % SSUM_SIZE] += val;
//			//break;
//		}
//		barrier(CLK_GLOBAL_MEM_FENCE|CLK_LOCAL_MEM_FENCE);
//	}
//	return ssum[0];
//    //printf("my TID is FOR FRAME IS %f\n",ssum[TID % SSUM_SIZE]);
//	for(j = 1; j < SSUM_SIZE; j *= 2)
//	{
//		if((TID < SSUM_SIZE - j) && (TID % (j*2)) == 0)
//			ssum[TID] += ssum[TID + j];
//		barrier(CLK_LOCAL_MEM_FENCE);
//	}
//	return ssum[0];
//}

__kernel void rmsd(		__global size_t *sz,
					    __global data_size *mask,
						__global data_size *ref,
						__global data_size *rmsf,
						__global data_size *rmsd) 
{ 
	__local double ssum[SHR_SIZE];
	size_t tid = get_local_id(0);
	size_t k = get_global_id(1);
	size_t len = *sz;
	size_t stride = get_global_size(1);
	size_t p,o,j;
	// double n = (double)atoms[1];
	data_size d = 0.0;
	data_size e = 0.0;
	for(j = tid*3; j < len*3; j+=3*get_local_size(0))
	{
		const size_t base = (3*(len*k ) + j);
		data_size xr = ref[j + 0];
		data_size yr = ref[j + 1];
		data_size zr = ref[j + 2];

		data_size xm = mask[base + 0];
		data_size ym = mask[base + 1];
		data_size zm = mask[base + 2];
		 
		 e = ((xm - xr) * (xm - xr)) + 
					((ym - yr) * (ym - yr)) + 
					((zm - zr) * (zm - zr));
		 d += e;
		 //rmsf[k*len + j/3] = sqrt(e);
		// barrier(CLK_GLOBAL_MEM_FENCE|CLK_LOCAL_MEM_FENCE);
	//	 rmsf[k*len + (j/3)] = d;
	}
		 //barrier(CLK_GLOBAL_MEM_FENCE|CLK_LOCAL_MEM_FENCE);
		 double rmsdv = 0.0;
		 //barrier(CLK_GLOBAL_MEM_FENCE|CLK_LOCAL_MEM_FENCE);
		 
	//		printf("BLK %zu TID %zu D IS %f\n",k,tid,d);
		//	printf("I AM BLK %zu TID %zu GX is %zu GY is %zu IS %f\n",k,tid,len,stride,d);
//			printf("I AM BLK %zu TID %zu LX is %zu LY is %zu IS %f\n",k,tid,get_local_size(0),get_local_size(1),d);
		//	return;
		 reduce(&d,ssum);
		 if(tid == 0)
		 {
		 	d = sqrt(d / (data_size)len);
			rmsd[k] = d;
			//printf("BLK %zu RMSD IS %f\n",k,d);
		 }
		 //if(tid == 0)
		 //{
		 //	rmsdv = sqrt(rmsdv / ((double)len));
		 //	rmsd[k] = rmsdv;
		 //}
		 //printf("block %zu tid %zu RMSD FOR FRAME IS %6.4f %6.4f %6.4f vs %6.4f %6.4f %6.4f d = %f\n", k,tid, xr,yr,zr,xm,ym,zm,d);
		 
		 //for(p = tid; p < SSUM_SIZE; p+=len)
		 //{
		 //   ssum[p] = 0.0;
		 //}
		 //barrier(CLK_LOCAL_MEM_FENCE);
		 //for(p = tid ; tid < SSUM_SIZE && p < len; p+=SSUM_SIZE) 
		 //{
		 //   rmsdv = rmsdv + rmsf[k*len + p]	;
		 //	if(tid == 0) printf("SUM IS %f ( %f )\n",rmsdv, rmsf[k*len + p]);
		 //}
		 //if( tid < SSUM_SIZE)
		 //   ssum[tid] = rmsdv;
		 //barrier(CLK_LOCAL_MEM_FENCE);
		 //if(tid == 0) printf("SUM IS %f ( %f )\n",rmsdv, rmsf[k*len + p]);
		 //for(o = 2; o <= SSUM_SIZE; o *= 2)
		 //{
		 //    for(p = tid, rmsdv = 0.0; p < SSUM_SIZE/o; p+=o)
		 //    {
		 //   	rmsdv = rmsdv + ssum[p];
		 //    }
		 //    if(tid < SSUM_SIZE/o)
		 //   	 ssum[tid] = rmsdv;
		 //	barrier(CLK_LOCAL_MEM_FENCE);
		 //}
		 //if(tid == 0)
		 //{
		 //   rmsdv = sqrt(ssum[0] / ((double)len));
		 //   rmsd[k] = rmsdv;
		 //}
		 
//		 if(tid == 0)
//		 {
//		    for(p = 0; p < len; p++)
//			{
//		    	rmsdv = rmsdv + rmsf[k*len + p];
//	//	 	printf("BLK %zu SUM IS %f ( %f )\n",k,rmsdv, rmsf[k*len + p]);
//			}
//		 	//printf("BLK %zu FINAL IS %f ( %f )\n",k,rmsdv, rmsf[k*len + p]);
//		    rmsdv = sqrt(rmsdv / ((double)len));
//		    rmsd[k] = rmsdv;
//		 }
}

// vim:ft=c

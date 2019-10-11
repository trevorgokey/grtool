#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#define EndianSwap(n) (rotate(n & 0x00FF00FF, 24U)|(rotate(n, 8U) & 0x00FF00FF))

#ifdef USE_DOUBLE
typedef double data_size;
#else
typedef float data_size;
#endif

inline data_size convert(float*d)
{
	return *d;
}
inline uchar4 swap(uchar4 *b)
{
	return (*b).wzyx; 
}

__kernel void netcdf(__global uchar4* in,__global data_size *out, __global size_t* size)
{
	uchar4 h;
	size_t idx;
	const size_t total = size[0]; 
    for(idx = get_global_id(0); idx < total; idx += get_global_size(0))
	{
		h = in[idx];
		h = swap(&h);
		out[idx] = convert(&h);
	}
}

// vi:ft=c

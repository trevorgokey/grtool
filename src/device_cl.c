
#include <stdio.h>
#include "device_cl.h"
#include "buffer_cl.h"

void device_cl_dty(struct device*);
static unsigned long long get_memory_size(struct device*);
int device_cl_buffer_bld(struct device*, struct device*, struct buffer*, int,size_t, size_t*, size_t);
void device_cl_buffer_dty(struct device*, struct device*, struct buffer*);
void device_cl_deallocator(uintptr_t);
void* device_cl_map(struct device*, struct buffer*);
void device_cl_unmap(struct device*, struct buffer*);

void device_cl_mem_dealloc(struct device* self, void* mem, size_t memsize)
{
#if (USE_RECYCLE == 1)
	memspace__push__bytes_(&self->memspace,mem,memsize);
#else
	clReleaseMemObject(mem);
#endif
}
uintptr_t device_cl_mem_alloc(struct device* self,size_t memsize)
{
#if (USE_RECYCLE == 1)
	return memspace__pop__bytes_(&self->memspace,memsize);
#else
	return 0;
#endif
}

struct device_interface device_cl_interface =
{
	.buffer_bld = device_cl_buffer_bld,
	.buffer_dty = device_cl_buffer_dty,
	.mem_alloc = device_cl_mem_alloc,
	.mem_dealloc = device_cl_mem_dealloc,
	.map = device_cl_map,
	.unmap = device_cl_unmap,
	.dty = device_cl_dty,
	.deallocator = device_cl_deallocator,
	.get_memory_size = get_memory_size
};

void* device_cl_map(struct device* self, struct buffer* buffer)
{
	struct device_cl_properties* props = self->properties;
	struct buffer_cl_properties* buf_props = buffer->properties;
//	double *d = malloc(buffer->size);	
//	clEnqueueReadBuffer(props->out,
//				buffer->data,
//				CL_TRUE,
//				0,
//				buffer->size,
//				d,
//				0, NULL, NULL);
//	buf_props->map_ptr = d;
//	return d;
	cl_int ret ;
	 buf_props->map_ptr = clEnqueueMapBuffer(props->in,
			buffer->data,
			buf_props->async,
			buf_props->map_flags,
			0,
			buffer->size,
			0,
			NULL,
			(cl_event*)&buffer->lock,
			&ret);
	if(ret != CL_SUCCESS)
		printf("MAP ERROR: %d\n",ret);
	 return buf_props->map_ptr;
}
void device_cl_unmap(struct device* self, struct buffer* buffer)
{
	struct device_cl_properties* props = self->properties;
	struct buffer_cl_properties* buf_props = buffer->properties;
	cl_int ret ;
//		ret = clEnqueueWriteBuffer(props->in,
//				buffer->data,
//				CL_TRUE,
//				0,
//				buffer->size,
//				buf_props->map_ptr,
//				0, NULL, NULL);
//	free(buf_props->map_ptr);
//	buf_props->map_ptr = NULL;
//	return;
	ret = clEnqueueUnmapMemObject(props->in,
			buffer->data,
			buf_props->map_ptr,
			0,
			NULL,
			NULL);
        if(buffer->lock) clReleaseEvent(buffer->lock) ; buffer->lock = NULL;
	if(ret != CL_SUCCESS)
		printf("UNMAP ERROR: %d\n",ret);
}

static unsigned long long get_memory_size(struct device* self)
{
	cl_ulong mem;
	clGetDeviceInfo( ((struct device_cl_properties*)self->properties)->id ,
			CL_DEVICE_GLOBAL_MEM_SIZE,sizeof(cl_ulong),&mem,0);
	return mem;
}

int device_cl_bld(struct device* self, cl_device_id id, cl_context ctx)
{
	cl_int ret = 0;
	struct device_cl_properties* props = malloc(sizeof(struct device_cl_properties));
	
	if(id)  props->id = id; else goto fail;
	if(ctx) 
	{
		props->ctx = ctx;
		clRetainContext(ctx);
	}
	else goto fail;
	//self->mem_size = get_memory_size(id)
	
	cl_command_queue_properties prop = 0;
#if STATS
	prop |= CL_QUEUE_PROFILING_ENABLE;
	prop |= CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE;
#endif
#if USE_ASYNC
	prop |= CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE;
#endif
	props->in = props->out = props->exe = NULL;
	props->exe  = clCreateCommandQueue(ctx,id,prop,&ret); 
	if(ret != CL_SUCCESS) goto fail;
	props->in  = clCreateCommandQueue(ctx,id,prop,&ret); 
	if(ret != CL_SUCCESS) goto fail;
	props->out = clCreateCommandQueue(ctx,id,prop,&ret); 
	if(ret != CL_SUCCESS) goto fail;
	props->exe = clCreateCommandQueue(ctx,id,prop,&ret); 
	if(ret != CL_SUCCESS) goto fail;
	
	self->properties = props;
	ret = device_bld(self, &device_cl_interface, NULL ); if(ret) goto fail;
	goto done;
	
fail:
	printf("FAIL: device\n");
	device_dty(self);
done:
	return ret;
}

void device_cl_print_kernel_info(cl_kernel k,cl_device_id i)
{
	size_t val,val2;
	cl_ulong val3,val4;
	clGetKernelWorkGroupInfo (k,
			i,
			CL_KERNEL_WORK_GROUP_SIZE,
			sizeof(size_t),
			&val,
			0);
	clGetKernelWorkGroupInfo (k,
			i,
		 CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE	,
			sizeof(size_t),
			&val2,
			0);
	clGetKernelWorkGroupInfo (k,
			i,
			CL_KERNEL_LOCAL_MEM_SIZE,
			sizeof(cl_ulong),
			&val3,
			0);
	clGetKernelWorkGroupInfo (k,
			i,
		 CL_KERNEL_PRIVATE_MEM_SIZE,
			sizeof(cl_ulong),
			&val4,
			0);
			printf("INFO: CL_KERNEL_WORK_GROUP_SIZE is                    %zu\n"
				   "      CL_KERNEL_LOCAL_MEM_SIZE  is                    %lu\n"
				   "      CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE is %zu\n"
				   "      CL_KERNEL_PRIVATE_MEM_SIZE is                   %lu\n"
				   ,val,val3,val2,val4);

}
void device_cl_properties_dty(struct device* device)
{
	struct device_cl_properties* props = device->properties;
}

double device_cl_get_exe_time_seconds(cl_event e)
{
	cl_ulong s,t,*ptr;
	ptr = &s;
	clGetEventProfilingInfo(e, CL_PROFILING_COMMAND_START, sizeof(*ptr), ptr, NULL);
	ptr = &t;
	clGetEventProfilingInfo(e, CL_PROFILING_COMMAND_END, sizeof(*ptr), ptr, NULL);
	return (t - s)/1e9;
}
double device_cl_get_delay_time_seconds(cl_event e)
{
	cl_ulong s,t,*ptr;
	ptr = &s;
	clGetEventProfilingInfo(e, CL_PROFILING_COMMAND_SUBMIT, sizeof(*ptr), ptr, NULL);
	ptr = &t;
	clGetEventProfilingInfo(e, CL_PROFILING_COMMAND_START, sizeof(*ptr), ptr, NULL);
	return (t - s)/1e9;
}

unsigned long long device_cl_get_time_ns(cl_event e, cl_profiling_info s)
{
	cl_ulong t;
	clGetEventProfilingInfo(e, s, sizeof(t), &t, NULL);
	return t;
}
double device_cl_get_queued_time_seconds(cl_event e)
{
	cl_ulong s,t,*ptr;
	ptr = &s;
	clGetEventProfilingInfo(e, CL_PROFILING_COMMAND_QUEUED, sizeof(*ptr), ptr, NULL);
	ptr = &t;
	clGetEventProfilingInfo(e, CL_PROFILING_COMMAND_SUBMIT, sizeof(*ptr), ptr, NULL);
	return (t - s)/1e9 ;
}

void device_cl_dty(struct device* self)
{
	struct device_cl_properties* props = self->properties;
	if(props->in)  clReleaseCommandQueue(props->in );
	if(props->out) clReleaseCommandQueue(props->out);
	if(props->exe) clReleaseCommandQueue(props->exe);
	if(props->ctx) clReleaseContext     (props->ctx);
	props->in = props->out = props->exe = NULL;
	props->id = NULL;
}

int device_cl_buffer_bld(struct device* self, struct device* host, struct buffer* buf, 
		int mode,size_t bytes_per_elem, size_t* dim, size_t dim_size)
{
	buffer_bld(self,host,buf,bytes_per_elem,dim,dim_size);
	struct device_cl_properties* device_cl = self->properties;
	struct buffer_cl_properties* buf_props =  
		(struct buffer_cl_properties*)host->fn->mem_alloc(host,
				sizeof(struct buffer_cl_properties));
	buf->data = (void*)device_cl_mem_alloc(self,buf->size);
	cl_int ret;
	if(buf->data == NULL && buf->size)
	{
		cl_mem_flags flags;
		flags = CL_MEM_READ_WRITE;
		if(mode == 'r') 
			flags = CL_MEM_READ_ONLY;
		if(mode == 'w') 
			flags = CL_MEM_WRITE_ONLY;
		if(mode == ('r' + 'w')) 
			flags = CL_MEM_READ_WRITE;
		if(mode == ('r' + 'p')) 
			flags = CL_MEM_READ_ONLY|CL_MEM_ALLOC_HOST_PTR;
		if(mode == ('w' + 'p')) 
			flags = CL_MEM_WRITE_ONLY|CL_MEM_ALLOC_HOST_PTR;
		if(mode == ('w' + 'p' + 'r')) 
			flags = CL_MEM_READ_ONLY|CL_MEM_ALLOC_HOST_PTR;
		//if(mode == ('r' + 'm')) flags = CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR;
		buf->data = clCreateBuffer
			(device_cl->ctx
			 ,flags
			 ,buf->size
			 ,0
			 ,&ret);
		if(ret != CL_SUCCESS)
			printf("OPENCL MALLOC ERROR: %d\n",ret);
		buf_props->flags = flags;
	}
	buf_props->map_flags = CL_MAP_WRITE | CL_MAP_READ;
	buf_props->host_buffer = NULL;
	buf_props->async = CL_TRUE;
	buf->properties = buf_props;

	return ret;
}

void device_cl_deallocator(uintptr_t mem)
{
	clReleaseMemObject((cl_mem)mem);
}


void device_cl_buffer_dty(struct device* self, struct device* host, struct buffer* b)
{
	//device_cl_mem_dealloc(b->device, (uintptr_t)b->data, b->size);
	struct buffer_cl_properties* buf_props =  b->properties;
	if(buf_props->host_buffer)
		buffer_dty(self,buf_props->host_buffer);
	host->fn->mem_dealloc(host,b->properties,sizeof(struct buffer_cl_properties));
	if(b->lock)
		clReleaseEvent((cl_event)b->lock);
	b->lock = NULL;
	//clReleaseMemObject((cl_mem)b->data);
}

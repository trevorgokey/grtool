#include <string.h> //memcpy
#include <stddef.h>
#include <stdio.h>
#include "sys/prctl.h"
#include "defines.h"
#include "system_cl.h"
#include "device_cl.h"
#include "buffer.h"

#define MAX_DEVICES 1

struct system_interface system_cl_interface =
{
	.dty = system_cl_dty
};

int system_cl_properties_bld(struct system* self)
{
	
}

void getIDs()
{
	//cl_uint num_platforms = 0;
	//cl_platform_id platforms[16];
	//ret = clGetPlatformIDs(16, platforms, &num_platforms);
	//printf("Notified %d platforms\n",num_platforms);
	//if(ret != CL_SUCCESS)
	//{
	//	printf("Error: OpenCl returned %d\n",ret);
	//}
}

int system_cl_bld(struct system* self, cl_platform_id id)
{
	cl_int ret;
	cl_device_id d[MAX_DEVICES];
	size_t i;
	cl_uint device_size;
	unsigned int num_bytes;

	struct system_cl_properties* props = malloc(sizeof(struct system_cl_properties));
	props->id = id;
	self->properties = props;
	
	prctl(PR_SET_NAME,"opencl_worker");
	
	ret = clGetDeviceIDs(id, CL_DEVICE_TYPE_ALL, MAX_DEVICES, d, &device_size);
	
	if( unlikely (device_size > MAX_DEVICES) )
	{
		printf("WARNING: only %d/%d devices will be used\n", MAX_DEVICES,device_size);
		device_size = MAX_DEVICES;
	}
	cl_context_properties ctx_props[] = {
		CL_CONTEXT_PLATFORM, (cl_context_properties)id, 0 };
	props->ctx = clCreateContext(ctx_props,device_size,d,0,0,&ret);
	if( unlikely (ret != CL_SUCCESS) )
	{
		printf("Error in clCreateContext: %d\n",ret);
		return ret;
	}
	ret |= system_bld(self,&system_cl_interface,device_size);
	
	for(i = 0, self->device_size = 0; i < device_size; ++i)
	{
		if(! device_cl_bld(&self->device[self->device_size],d[i],props->ctx) )
		{
			printf("CL: Device %zu: %p online\n",i,d[i]);
			self->device_size++;
		}
		else
		{
			printf("CL: Device %zu: %p FAILED\n",i,d[i]);
			device_dty(&self->device[self->device_size]);
		}
	}
	
	prctl(PR_SET_NAME,"a.out");
}

void system_cl_dty(struct system* self)
{
	struct system_cl_properties* props = self->properties;
	if(props->ctx != NULL)
		clReleaseContext(props->ctx);
	props->ctx = NULL;
}


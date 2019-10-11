
#include <stdio.h>
#include <string.h>
#include "defines.h"
#include "driver_cl.h"

#define LOG_SIZE (1024*1024)

struct driver_interface driver_cl_interface =
{
	.dty = driver_cl_dty
};

static char* get_build_log(cl_device_id id, cl_program program, char* log, size_t* lsize)
{
	cl_int ret = clGetProgramBuildInfo(program,
			id,
			CL_PROGRAM_BUILD_LOG,
			LOG_SIZE,
			log,
			lsize);
	return log;
}

int  driver_cl_bld(struct driver* self,cl_program program, char* kname, struct device* device)
{
	driver_bld(self,device);
	struct driver_cl_properties* props = malloc(sizeof(struct driver_cl_properties));
	struct device_cl_properties* device_props = device->properties;
	cl_int ret = CL_SUCCESS;
	cl_kernel kernel;
	int r;
	
	
	self->properties = props;
	props->program = NULL;
	props->kernel = NULL;
	self->fn = &driver_cl_interface;
    char opts[128];
    memset(opts,0,128);
    size_t off = 0;
	//const char opts[] = "-cl-opt-disable";
	//const char opts[] = "-cl-fast-relaxed-math  -cl-mad-enable";
#if USE_DOUBLE
	const char dbl[] = " -DUSE_DOUBLE ";
    memcpy(opts+off,dbl,strlen(dbl));
    off += strlen(dbl);
#endif
    sprintf(opts+off," -DALIGN_SHR_SIZE=%d -DTRANS_SHR_SIZE=%d ",ALIGN_SHR_SIZE,TRANS_SHR_SIZE);
//	cl_command_queue_properties prop = 0;
//#if STATS
//	prop |= CL_QUEUE_PROFILING_ENABLE;
//	//prop |= CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE;
//#endif
//	props->exe  = clCreateCommandQueue(device_props->ctx,device_props->id,prop,&ret); 
	ret = clBuildProgram(program,1,&device_props->id,opts,0,0);
	if( unlikely (ret != CL_SUCCESS) )
	{	
		char log[LOG_SIZE];
		printf("ERROR in clBuildProgram ret %d\n%s\n****\n", ret,
				get_build_log(device_props->id, program, log, 0) );
		r = ret;
	}
	else
	{
		//clRetainProgram(program);
		props->program = program;
		kernel = clCreateKernel(program,kname,&ret);
		if( likely (ret == CL_SUCCESS) )
			props->kernel = kernel;	
		r = ret;
	}
	return r;
}
void driver_cl_dty(struct driver* self)
{
	struct driver_cl_properties* props = self->properties;
	if(props->kernel ) clReleaseKernel (props->kernel );
	if(props->program) clReleaseProgram(props->program);
	//if(props->exe) clReleaseCommandQueue(props->exe);
	props->kernel = NULL;
	props->program = NULL;
}

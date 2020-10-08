
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "module_cl.h"
#include "driver_cl.h"


static struct module_interface interface = 
{
	.start = NULL,
	.stop = NULL,
	.dty = module_cl_dty
};

int module_cl_bld(struct module* self, struct system* system, char* filename, char* kname)
{
	module_bld(self,&interface,system);
	size_t i,j;
	FILE* fp = fopen(filename,"r");
	const int len = 1024;
	char buffer[len];
	char**  ptr;
	size_t* slen;
	char**  ksource = 0;
	size_t  ksource_lines = 0;
	size_t* ksource_lines_size = 0;
	size_t  lnum = 0;
	struct module_cl_properties* module_props = malloc(sizeof(struct module_cl_properties));
	struct system_cl_properties* system_props = system->properties;
	if(module_props == NULL) return -1;
	module_props->program = NULL;
	module_props->kernel = NULL;
	self->properties = module_props;
	while(fgets(buffer,len,fp) != NULL)
	{
		int strln = 1024;
		if (strlen(buffer) < strln) strln = strlen(buffer);
		ptr = realloc(ksource,sizeof(char*) * (lnum+1));
		slen = realloc(ksource_lines_size,sizeof(size_t)*(lnum+1));
		if(ptr && slen)
		{
			ptr[lnum] = malloc(strln+1);
			slen[lnum] = strln;
			strncpy(ptr[lnum],buffer,strln+1);
			lnum++;
			ksource = ptr;
			ksource_lines = lnum;
			ksource_lines_size = slen;
		}
		else
		{
			printf("TRANSlATE: COULD NOT ReAD CL FILE\n");
			break;
		}
	}
	cl_int ret = 0;
	fclose(fp);



	for(i = 0; i < self->driver_size; i++)
	{
		cl_program program = clCreateProgramWithSource(system_props->ctx
				,ksource_lines
				,(const char**)ksource
				,ksource_lines_size
				,&ret);
		if(ret == CL_SUCCESS)
		{
			if( driver_cl_bld(&self->driver[i], program, kname, 
					&system->device[i]) != CL_SUCCESS)
			{
				for(j = 0; j < ksource_lines; j++)
					printf("%6zu: %s",j,ksource[j]);
				printf("\n");
				
			}
			clReleaseProgram(module_props->program);
		}
	}
	
	for(i = 0; i < ksource_lines; i++)
	{
		free(ksource[i]);
	}
	free(ksource);
	free(ksource_lines_size);
	assert(ret == CL_SUCCESS);
	//s->sizes[0] = 1;
	//s->sizes[1] = 0;
	//s->sizes[2] = 0;
	//s->clsizes = clCreateBuffer(stx->cunit[1].system.cl.ctx,
	//		CL_MEM_READ_ONLY|CL_MEM_USE_HOST_PTR,
	//		sizeof(cl_uint)*3,
	//		&s->sizes,
	//		&ret);
	//assert(ret == CL_SUCCESS);
	//assert(ret == CL_SUCCESS);
	//ret = clSetKernelArg(s->kernel,0,sizeof(cl_mem),&s->clsizes);
	//assert(ret == CL_SUCCESS);
	
	
	return 0;
}

void module_cl_dty(struct module* self)
{
	struct module_cl_properties* props = self->properties;
	if(props->program) clReleaseProgram(props->program);
	props->program = NULL;
}

#undef base

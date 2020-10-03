
#ifndef MODULE_CL_H
#define MODULE_CL_H

#include "module.h"
#include "system_cl.h"


struct module_cl_properties
{
	size_t* dim_global;
	cl_uint  dim_global_size;
	
	size_t* dim_local;
	cl_uint  dim_local_size;
	
	cl_program program;
	cl_kernel  kernel;
};


int  module_cl_bld(struct module*, struct system*, char*, char*);
void module_cl_dty(struct module*);

#endif /* MODULE_CL_H */

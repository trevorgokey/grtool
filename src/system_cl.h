
#ifndef SYSTEM_CL_H
#define SYSTEM_CL_H

#include <CL/cl.h>
#include "buffer.h"
#include "system.h"
#include "buffer_cl.h"


struct system_cl_properties
{
	cl_platform_id id;
	cl_context ctx;
	cl_int cl_err;
};


int  system_cl_bld(struct system*, cl_platform_id);
void system_cl_dty(struct system*);

int system_cl_allocator(struct buffer*);
void system_cl_deallocator(struct buffer*);
#endif

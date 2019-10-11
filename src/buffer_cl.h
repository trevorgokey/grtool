
#ifndef BUFFER_CL_H
#define BUFFER_CL_H

#include <CL/cl.h>
#include "buffer.h"

struct buffer_cl_properties
{
	cl_mem_flags flags;
	cl_map_flags map_flags;
	cl_bool async;
	void* map_ptr;
	struct buffer* host_buffer;
};

void buffer_cl_dty(struct buffer*);

#endif

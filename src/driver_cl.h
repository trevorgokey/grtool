
#ifndef DRIVER_CL_H
#define DRIVER_CL_H

#include "driver.h"
#include "device_cl.h"

struct driver_cl_properties
{
	cl_command_queue exe;
	cl_program program;
	cl_kernel kernel;
};

int  driver_cl_bld(struct driver*,cl_program, char*, struct device*);
void driver_cl_dty(struct driver*);

#endif /* DRIVER_CL_H */

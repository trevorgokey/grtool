
#ifndef DEVICE_CL_H
#define DEVICE_CL_H

#include <CL/cl.h>
#include "device.h"


struct device_cl_properties
{
	cl_command_queue in,out,exe;
	cl_device_id id;
	cl_context ctx;
};


int  device_cl_bld(struct device*, cl_device_id,cl_context);
void device_cl_dty(struct device*);
double device_cl_get_exe_time_seconds(cl_event);
double device_cl_get_queued_time_seconds(cl_event);
double device_cl_get_delay_time_seconds(cl_event);

void device_cl_print_kernel_info(cl_kernel,cl_device_id);
unsigned long long device_cl_get_time_ns(cl_event, cl_profiling_info);
#endif

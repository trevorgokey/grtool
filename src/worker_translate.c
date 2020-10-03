#include <stdio.h>
#include "defines.h"
#include "worker_translate.h"
#include "system_cl.h"
#include "module_cl.h"
#include "module_linux_simple.h"
#include "device_cl.h"
#include "driver_cl.h"


static int worker_translate_start(struct worker*);
static void worker_translate_process(struct worker*, struct module*, struct frame*); 
static int worker_translate_stop(struct worker*);
static int worker_translate_module_cl_bld(struct worker*, struct module*, struct system*);
static int worker_translate_module_cl_dty(struct worker*, struct module*); 
static int worker_translate_module_linux_simple_bld(struct worker*, struct module*, 
		struct system*);
static int worker_translate_module_linux_simple_dty(struct worker*, struct module*); 
static void worker_translate_prepare_buffers(struct worker*, struct module*, struct frame**);
static void worker_translate_process_buffers(struct worker*);
struct worker_interface worker_translate_interface = 
{
	.start = worker_translate_start,
	.stop = worker_translate_stop,
	.prepare_buffers = worker_translate_prepare_buffers,
	.process_buffers = worker_translate_process_buffers,
	.dty = worker_translate_dty,
	.module_linux_simple_bld = worker_translate_module_linux_simple_bld,
	.module_linux_simple_dty = worker_translate_module_linux_simple_dty,
	.module_cl_bld = worker_translate_module_cl_bld,
	.module_cl_dty = worker_translate_module_cl_dty
};

int worker_translate_bld(struct worker* self, struct device* host)
{
	return worker_bld(self, &worker_translate_interface, host);
}

void worker_translate_dty(struct worker* self)
{
	int in_used[] = {2,-1};
	worker_clean_buffers(self,self->input,in_used);
	int out_used[] = {0,-1};
	worker_clean_buffers(self,self->output,out_used);

}
static int worker_translate_start(struct worker* self) { return WORKER_UNIMPLEMENTED; }
static void worker_translate_process(struct worker* self, struct module* module, struct frame* frame) 
{
	return; 
}
static int worker_translate_stop(struct worker* self) { return WORKER_UNIMPLEMENTED; }


static void worker_translate_module_linux_simple_process(struct worker* self, struct module* module, 
		struct buffer** in, struct buffer** out)
{
	size_t i,j,atoms;
	data_size* main = in[0]->data;
	data_size* mask = in[0]->data;
	data_size* mean = out[0]->data;
	atoms = in[0]->dim[1];
#if STATS
	struct timespec start,finish;
		clock_gettime(CLOCK_REALTIME, &start);
#endif
	for(j = 0; j < in[0]->dim[0] ; j++)
	{
		mean[3*j + 0] = 0.0;
		mean[3*j + 1] = 0.0;
		mean[3*j + 2] = 0.0;
		for(i = 0; i < atoms ; i++)
		{
			mean[3*j + 0] += mask[j*atoms*3 + i*3 + 0] / atoms;
			mean[3*j + 1] += mask[j*atoms*3 + i*3 + 1] / atoms;
			mean[3*j + 2] += mask[j*atoms*3 + i*3 + 2] / atoms;
	//		printf("ID %zu SUM IS %16.16f add from %16.16f\n",i,mean[3*j + 0],  mask[j*atoms*3 + i*3 + 0] / atoms);
		}
		for(i = 0; i < atoms ; i++)
		{
			main[j*atoms*3 + i*3 + 0] -= mean[3*j + 0];
			main[j*atoms*3 + i*3 + 1] -= mean[3*j + 1];
			main[j*atoms*3 + i*3 + 2] -= mean[3*j + 2];
		}
	}
#if STATS
	clock_gettime(CLOCK_REALTIME, &finish);
	self->exe_tot += (finish.tv_sec - start.tv_sec) + 
		((finish.tv_nsec - start.tv_nsec) / 1e9);
#endif
}
static void worker_translate_module_cl_process(struct worker* self, struct module* module, 
		struct buffer** in, struct buffer** out)
{
	struct worker_translate_properties* props = 
		(struct worker_translate_properties*)self->properties;
	struct device_cl_properties* device_props = 
		(struct device_cl_properties*)module->driver[0].device->properties;
	struct driver_cl_properties* driver_props = 
		(struct driver_cl_properties*)module->driver[0].properties;
	size_t l = TRANS_BLK;
	size_t dims[] = { l,in[0]->dim[0] };
	size_t ldims[] = { l,1 };
	size_t dim[] = { in[0]->dim[1] };
	cl_uint num_dims = 2;
	static int first = 0;
#if STATS
	struct timespec start,finish;
		clock_gettime(CLOCK_REALTIME, &start);
	unsigned long long time;
#endif
	cl_bool blk = !CL_ASYNC; 
	cl_int ret  = 0;
	if(first == 0)
	{
		ret = clEnqueueWriteBuffer(device_props->in,
			in[2]->data,
			blk,
			0,
			sizeof(size_t),
			dim,
			0, NULL,
#if CL_ASYNC
			(cl_event*)&in[2]->lock
#else
			NULL
#endif
			);
	//	clFlush(device_props->in);
	if(unlikely(ret != CL_SUCCESS))
	{
		printf("TRANS: CPY FAIL %d\n",ret);
	}
	ret = clSetKernelArg(driver_props->kernel,
			0,
			sizeof(cl_mem),
			&in[2]->data);
	}
	ret |= clSetKernelArg(driver_props->kernel,
			1,
			sizeof(cl_mem),
			&in[0]->data);
	ret |= clSetKernelArg(driver_props->kernel,
			2,
			sizeof(cl_mem),
			&in[0]->data);
	ret |= clSetKernelArg(driver_props->kernel,
			3,
			sizeof(cl_mem),
			&out[0]->data);
	if(unlikely(ret != CL_SUCCESS))
	{
		printf("ARG FAIL %d\n",ret);
	}
#if CL_ASYNC 
	if(first == 0)
	{
	    clWaitForEvents(1, (cl_event*)&in[2]->lock);
#if CLEAN_UP
	    clReleaseEvent(in[2]->lock); in[2]->lock = 0;
#endif
	}
#endif
    first = 1;
#if STATS
	clock_gettime(CLOCK_REALTIME, &finish);
	self->memin_tot += (finish.tv_sec - start.tv_sec) + 
	    ((finish.tv_nsec - start.tv_nsec) / 1e9);
	if(!first)
	{
	    self->memin_exe += device_cl_get_exe_time_seconds(in[2]->lock);
	    self->memin_wait += device_cl_get_queued_time_seconds(in[2]->lock);
	    self->memin_delay += device_cl_get_delay_time_seconds(in[2]->lock);
	}
	clock_gettime(CLOCK_REALTIME, &start);
	register_profile_event(self,"proc_ocl_kern_wait_start",self->frame_cnt,0);
#endif
	ret = clEnqueueNDRangeKernel(device_props->exe,
			driver_props->kernel,
			num_dims,
			NULL,
			dims,
			ldims,
			0,
			NULL,
			(cl_event*)&in[0]->lock
			);
	if(unlikely(ret != CL_SUCCESS))
	{
		printf("EXEC ERROR: %d\n",ret);
	}
	clWaitForEvents(1,(cl_event*)&in[0]->lock);
#if STATS
	//clock_gettime(CLOCK_REALTIME, &finish);
	//	self->exe_tot += (finish.tv_sec - start.tv_sec) + 
	//		((finish.tv_nsec - start.tv_nsec) / 1e9);
	//self->exe_exe += device_cl_get_exe_time_seconds(in[0]->lock);
	//self->exe_wait += device_cl_get_queued_time_seconds(in[0]->lock);
	//self->exe_delay += device_cl_get_delay_time_seconds(in[0]->lock);
	time = (start.tv_sec * 1e9) + (start.tv_nsec);
	time +=  (unsigned long long)
		(device_cl_get_queued_time_seconds(in[0]->lock)*1e9);
	register_profile_event(self,"proc_ocl_kern_delay_start",self->frame_cnt,&time);
	time +=  (unsigned long long)
		(device_cl_get_delay_time_seconds(in[0]->lock)*1e9);
	register_profile_event(self,"proc_ocl_kern_exe_start",self->frame_cnt,&time);
	time +=  (unsigned long long)
		(device_cl_get_exe_time_seconds(in[0]->lock)*1e9);
	//time = (finish.tv_sec * 1e9) + (finish.tv_nsec);
	register_profile_event(self,"proc_fin",self->frame_cnt,&time);
	//register_profile_event(self,"proc_fin",self->frame_cnt,0);
#endif
#if CLEAN_UP
	if(in[0]->lock) clReleaseEvent(in[0]->lock); in[0]->lock = 0;
#endif
}

static int worker_translate_module_linux_simple_bld(struct worker* self, 
		struct module* module, struct system* system) 
{ 
	if(! module_linux_simple_bld(module,system) )
	{
		worker_add_module(self,module,&worker_translate_module_linux_simple_process);
		return  0;
	}
	else
		return -1;
}

static int worker_translate_module_linux_simple_dty(struct worker* self,struct module* module) 
{ 
	return 0; 
}
static int worker_translate_module_cl_bld(struct worker* self, struct module* module,
		struct system* system) 
{ 
	char* filename = "translate_shm.cl";
	char* kname    = "translate";
	if(! module_cl_bld(module,system,filename,kname) )
	{
		worker_add_module(self,module,&worker_translate_module_cl_process);
		return 0;
	}
	else
		return -1;
}

static int worker_translate_module_cl_dty(struct worker* self,struct module* module) 
{ 
	return 0; 
}

static void worker_translate_prepare_buffers(struct worker* self, struct module* module, 
		struct frame** f)
{
	struct frame* frame = f[0];
	struct device* device = module->driver[0].device;
	if(unlikely (self->input == NULL) ) 
	{
		self->input = malloc(sizeof(struct buffer*) * 4);
		if(likely   (self->input != NULL) )
		{
			self->input[0] = NULL;
			self->input[1] = NULL;
			self->input[2] = NULL;
			self->input[3] = NULL;
		}
	}

	if(likely   (self->input != NULL) )
	{
		self->input[0] = &frame->buffer[0];
		self->input[1] = &frame->buffer[1];
		if(unlikely(self->input[2] == NULL))
		{
	                size_t dims[] = { 1 };
	                size_t num_dims = 1;
			self->input[2] = (struct buffer*)self->host->fn->mem_alloc(self->host,
					sizeof(struct buffer) );
			device->fn->buffer_bld(device,self->host,self->input[2], 
					'r',sizeof(size_t), dims,num_dims);
		}
		self->input[3] = NULL;
	}
	if(unlikely (self->output == NULL) ) 
	{
		self->output = malloc(sizeof(struct buffer*) * 2);
		self->output[0] = NULL;
		self->output[1] = NULL;
	}

	if(likely   (self->output != NULL) )
	{
	//	const int num_dims = 2;
	//	size_t d[] = { 3, frame->buffer[1].dim[0] };
	//	self->output[0] = (struct buffer*)self->host->fn->mem_alloc(self->host,sizeof(struct buffer));
	//	struct device* device = module->driver[module->driver_next_idx].device;
	//	device->fn->buffer_bld(device,self->host,self->output[0], 'w',sizeof(double),
	//			d,num_dims);
		if(unlikely(self->output[0] == NULL))
		{
			size_t dims[] = { self->input[0]->dim[0],3 };
			size_t num_dims = 2;
			self->output[0] = (struct buffer*)self->host->fn->mem_alloc(self->host,
					sizeof(struct buffer) );
			device->fn->buffer_bld(device,self->host,self->output[0], 
					'w' ,sizeof(data_size), dims,num_dims);
		}
		self->output[1] = NULL;
	}
}
static void worker_translate_process_buffers(struct worker* self)
{
	FILE* fout = NULL;
	//struct device_cl_properties* device_props = self->input[0]->device->properties;
	data_size* d;
	size_t i;
#if 0
	fout = fopen("trans.dat","a");
	if(fout)
	{
		d = malloc(self->output[0]->size / self->output[0]->dim[0]);
		//clEnqueueReadBuffer(device_props->out,
		//		self->output[0]->data,
		//		CL_TRUE,
		//		0,
		//		self->output[0]->size / self->output[0]->dim[0],
		//		d,
		//		0, NULL, NULL);
		//free(d);
	//((struct buffer_cl_properties*)self->input[0]->properties)->map_flags = CL_MAP_READ;
	d = self->input[0]->device->fn->map(self->input[0]->device,
			self->input[0]);
	for(i = 0; i < self->input[0]->dim[0] * self->input[0]->dim[1] * 3; i+=3)
		//fprintf(fout,"%p %p %p\n",&d[i],&d[i+1],&d[i+2]);
		//fprintf(fout,"%x %x %x\n",(uint)d[i],(uint)d[i+1],(uint)d[i+2]);
		fprintf(fout,"%zu: %f %f %f\n",i,d[i],d[i+1],d[i+2]);
	fprintf(fout,"**\n");
	self->input[0]->device->fn->unmap(self->input[0]->device,self->input[0]);
	fclose(fout);
	}
#endif
#if WRITE_TRANS
	fout = fopen("off.dat","a");
	if(fout)
	{
//		d = malloc(self->output[0]->size);
//		clEnqueueReadBuffer(device_props->out,
//				self->output[0]->data,
//				CL_TRUE,
//				0,
//				self->output[0]->size,
//				d,
//				0, NULL, NULL);
	//((struct buffer_cl_properties*)self->output[0]->properties)->map_flags = CL_MAP_READ;
	d = self->output[0]->device->fn->map(self->output[0]->device,
			self->output[0]);
	for(i = 0; i < self->output[0]->dim[0] * self->output[0]->dim[1] ; i+=3)
		//fprintf(fout,"%p %p %p\n",&d[i],&d[i+1],&d[i+2]);
		//fprintf(fout,"%x %x %x\n",(uint)d[i],(uint)d[i+1],(uint)d[i+2]);
		fprintf(fout,"%zu: %16.12f %16.12f %16.12f\n",i,d[i],d[i+1],d[i+2]);
	fprintf(fout,"**\n");
	self->output[0]->device->fn->unmap(self->output[0]->device,self->output[0]);
	//	free(d);
	fclose(fout);
	}
#endif

}

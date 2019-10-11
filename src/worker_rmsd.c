#include <stdio.h>
#include <math.h>
#include <string.h>
#include "defines.h"
#include "worker_rmsd.h"
#include "system_cl.h"
#include "module_cl.h"
#include "device_cl.h"
#include "driver_cl.h"


static int worker_rmsd_start(struct worker*);
static void worker_rmsd_process(struct worker*, struct module*, struct frame*); 
static int worker_rmsd_stop(struct worker*);
static int worker_rmsd_module_cl_bld(struct worker*, struct module*, struct system*);
static int worker_rmsd_module_cl_dty(struct worker*, struct module*); 
static int worker_rmsd_module_linux_simple_bld(struct worker*, struct module*, 
		struct system*);
static int worker_rmsd_module_linux_simple_dty(struct worker*, struct module*); 
static void worker_rmsd_prepare_buffers(struct worker*, struct module*, struct frame**);
static void worker_rmsd_process_buffers(struct worker*);
struct worker_interface worker_rmsd_interface = 
{
	.start = worker_rmsd_start,
	.stop = worker_rmsd_stop,
	.prepare_buffers = worker_rmsd_prepare_buffers,
	.process_buffers = worker_rmsd_process_buffers,
	.dty = worker_rmsd_dty,
	.module_linux_simple_bld = worker_rmsd_module_linux_simple_bld,
	.module_linux_simple_dty = worker_rmsd_module_linux_simple_dty,
	.module_cl_bld = worker_rmsd_module_cl_bld,
	.module_cl_dty = worker_rmsd_module_cl_dty
};

int worker_rmsd_bld(struct worker* self, struct device* host)
{
	int ret = worker_bld(self, &worker_rmsd_interface, host);
	self->properties = malloc(sizeof(struct worker_rmsd_properties));
	((struct worker_rmsd_properties*)self->properties)->frmsd = fopen("rmsd.out","w");
	return ret;
}

void worker_rmsd_dty(struct worker* self)
{
	fclose( ((struct worker_rmsd_properties*)self->properties)->frmsd );
	((struct worker_rmsd_properties*)self->properties)->frmsd = NULL;
	int in_used[] = {1,-1};
	worker_clean_buffers(self,self->input,in_used);
	int out_used[] = {0,1,-1};
	worker_clean_buffers(self,self->output,out_used);

}
static int worker_rmsd_start(struct worker* self) { return WORKER_UNIMPLEMENTED; }
static void worker_rmsd_process(struct worker* self, struct module* module, struct frame* frame) 
{
	return; 
}
static int worker_rmsd_stop(struct worker* self) { return WORKER_UNIMPLEMENTED; }


static void worker_rmsd_module_linux_simple_process(struct worker* self, struct module* module, 
		struct buffer** in, struct buffer** out)
{
	int i = 0;
	size_t k,p;
	unsigned int index = 0;
	unsigned int stride3 = 0;
	size_t stride = in[0]->dim[0];
	size_t atoms = in[0]->dim[1];
	double* data = in[0]->data;
	double* mask = in[0]->data;
	size_t maskSize = atoms;
	double* maskRef = in[1]->data;
	size_t counter;
	const double fms = maskSize;

	double *rmsf = out[0]->data;
	double *rmsd = out[1]->data;
	memset(rmsf,0,out[0]->size);
	memset(rmsd,0,out[1]->size);
	
#if STATS
	struct timespec start,finish;
	clock_gettime(CLOCK_REALTIME,&start);
#endif
	for(counter = 0; counter < stride; counter++)
	{
		for(k = 0; k < maskSize*3; k+=3)
		{
			double d = 0.0;
			double x; 
			x =  mask[k + 0] - maskRef[k + 0];
			d += x * x;
			x =  mask[k + 1] - maskRef[k + 1];
			d += x * x;
			x =  mask[k + 2] - maskRef[k + 2];
			d += x * x;
			rmsf[k/3] = d;
		}
		double rmsdval = 0.0;
		for(p = 0; p < maskSize; p++){

			rmsdval += ((rmsf[p]) / fms);
	//	printf("rmsf for %zu is %f\n",p,rmsf[p]);
		}
	//	printf("rmsd for %zu is %f\n",p,rmsdval * fms);
		rmsdval = sqrt(rmsdval);
		rmsd[counter] = rmsdval;
		mask += maskSize * 3;
		rmsf += maskSize;
	}
#if STATS
	clock_gettime(CLOCK_REALTIME, &finish);
	self->exe_tot += (finish.tv_sec - start.tv_sec) + 
		((finish.tv_nsec - start.tv_nsec) / 1e9);
#endif
}

static void worker_rmsd_module_cl_process(struct worker* self, struct module* module, 
		struct buffer** in, struct buffer** out)
{
	struct worker_rmsd_properties* props = 
		(struct worker_rmsd_properties*)self->properties;
	struct device_cl_properties* device_props = 
		(struct device_cl_properties*)module->driver[0].device->properties;
	struct driver_cl_properties* driver_props = 
		(struct driver_cl_properties*)module->driver[0].properties;
	
	unsigned long long time;
	size_t l = 128;
	size_t ldims[] = {l,1};
	size_t dims[] = { l,in[0]->dim[0]};
	//size_t dim[] = { in[0]->dim[0], in[0]->dim[1], in[0]->dim[1] };
	cl_uint num_dims = 2;
	static int first = 0;
	cl_int ret = 0;

#if STATS
	struct timespec start,finish;
	clock_gettime(CLOCK_REALTIME, &start);
#endif
	ret = clEnqueueWriteBuffer(device_props->in,
			in[0]->data,
			CL_TRUE,
			0,
			sizeof(size_t) * 1,
			&in[0]->dim[1],
			//0, NULL, NULL);
			0, NULL, (cl_event*)&in[2]->lock);
	if(unlikely(ret != CL_SUCCESS))
	{
		printf("CPY FAIL %d\n",ret);
	}
	ret = clSetKernelArg(driver_props->kernel,
			1,
			sizeof(cl_mem),
			&in[0]->data);
	if(!first)
	{
	ret |= clSetKernelArg(driver_props->kernel,
			0,
			sizeof(cl_mem),
			&in[2]->data);
	ret |= clSetKernelArg(driver_props->kernel,
			2,
			sizeof(cl_mem),
			&in[1]->data);
	}

	first = 1;
	ret |= clSetKernelArg(driver_props->kernel,
			3,
			sizeof(cl_mem),
			&out[0]->data);
	ret |= clSetKernelArg(driver_props->kernel,
			4,
			sizeof(cl_mem),
			&out[1]->data);
	if(unlikely(ret != CL_SUCCESS))
	{
		printf("ARG FAIL %d\n",ret);
	}
	clWaitForEvents(1, (cl_event*)&in[2]->lock);
#if STATS
		clock_gettime(CLOCK_REALTIME, &finish);
			self->memin_tot += (finish.tv_sec - start.tv_sec) + 
				((finish.tv_nsec - start.tv_nsec) / 1e9);
	clock_gettime(CLOCK_REALTIME,&start);
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
			//NULL);
			(cl_event*)&in[0]->lock);
		//clFlush(device_props->exe);
	if(unlikely(ret != CL_SUCCESS))
	{
		printf("EXEC ERROR: %d\n",ret);
	}
	clWaitForEvents(1,(cl_event*)&in[0]->lock);
#if STATS
	clock_gettime(CLOCK_REALTIME, &finish);
			self->exe_tot += (finish.tv_sec - start.tv_sec) + 
				((finish.tv_nsec - start.tv_nsec) / 1e9);
		self->exe_exe += device_cl_get_exe_time_seconds(in[0]->lock);
		self->exe_wait += device_cl_get_queued_time_seconds(in[0]->lock);
		self->exe_delay += device_cl_get_delay_time_seconds(in[0]->lock);
		time = (start.tv_sec * 1e9) + (start.tv_nsec);
		time +=  (unsigned long long)
			(device_cl_get_queued_time_seconds(in[0]->lock)*1e9);
		register_profile_event(self,"proc_ocl_kern_delay_start",self->frame_cnt,&time);
		time +=  (unsigned long long)
			(device_cl_get_delay_time_seconds(in[0]->lock)*1e9);
		register_profile_event(self,"proc_ocl_kern_exe_start",self->frame_cnt,&time);
		time +=  (unsigned long long)
			(device_cl_get_exe_time_seconds(in[0]->lock)*1e9);
		time = (finish.tv_sec * 1e9) + (finish.tv_nsec);
	clock_gettime(CLOCK_REALTIME, &start);
	register_profile_event(self,"proc_fin",self->frame_cnt,&time);
#endif
#if CLEAN_UP
	if(in[0]->lock) clReleaseEvent(in[0]->lock); in[0]->lock = 0;
	if(in[2]->lock) clReleaseEvent(in[2]->lock); in[2]->lock = 0;
#endif
	//clock_gettime(CLOCK_REALTIME, &finish);
	//		self->proc_fin += (finish.tv_sec - start.tv_sec) + 
	//			((finish.tv_nsec - start.tv_nsec) / 1e9);
}

static int worker_rmsd_module_linux_simple_bld(struct worker* self, 
		struct module* module, struct system* system) 
{ 
	if(! module_linux_simple_bld(module,system) )
	{
		worker_add_module(self,module,&worker_rmsd_module_linux_simple_process);
		return  0;
	}
	else
		return -1;
}

static int worker_rmsd_module_linux_simple_dty(struct worker* self,struct module* module) 
{ 
	return 0; 
}
static int worker_rmsd_module_cl_bld(struct worker* self, struct module* module,
		struct system* system) 
{ 
	char* filename = "rmsd_shm.cl";
	char* kname    = "rmsd";
	if(! module_cl_bld(module,system,filename,kname) )
	{
		worker_add_module(self,module,&worker_rmsd_module_cl_process);
		return 0;
	}
	else
		return -1;
}

static int worker_rmsd_module_cl_dty(struct worker* self,struct module* module) 
{ 
	return 0; 
}

static void worker_rmsd_prepare_buffers(struct worker* self, struct module* module, 
		struct frame** f)
{
	struct frame* frame = f[0];
	size_t dims[] = { 3 };
	size_t num_dims = 1;
	struct device* device = module->driver[0].device;
	if(unlikely (self->input == NULL) ) 
	{
		self->input = malloc(sizeof(struct buffer*) * 4);
		if(likely   (self->input != NULL) )
		{
			self->input[0] = NULL;
			self->input[1] = NULL;
			self->input[2] = NULL;
		}
	}

	if(likely   (self->input != NULL) )
	{
		self->input[0] = &frame->buffer[0];
		if(unlikely(self->input[1] == NULL))
		{
			size_t dims[] = { self->input[0]->dim[1], self->input[0]->dim[2] };
			size_t num_dims = 2;
			self->input[1] = (struct buffer*)self->host->fn->mem_alloc(self->host,
					sizeof(struct buffer) );
			device->fn->buffer_bld(device,self->host,self->input[1], 
					'r',sizeof(data_size), dims,num_dims);
			double* tdata = device->fn->map(device,self->input[0]);
			double* sdata = device->fn->map(device,self->input[1]);
			memcpy(sdata,tdata,self->input[1]->size);
			device->fn->unmap(device,self->input[1]);
			device->fn->unmap(device,self->input[0]);
		}
		if(unlikely(self->input[2] == NULL))
		{
			size_t dims[] = { self->input[0]->dim[1] };
			size_t num_dims = 1;
			self->input[2] = (struct buffer*)self->host->fn->mem_alloc(self->host,
					sizeof(struct buffer) );
			device->fn->buffer_bld(device,self->host,self->input[2], 
					'r' + 'p',sizeof(data_size), dims,num_dims);
			data_size* tdata = device->fn->map(device,self->input[0]);
			data_size* sdata = device->fn->map(device,self->input[2]);
			memcpy(sdata,tdata,self->input[2]->size);
			device->fn->unmap(device,self->input[0]);
			device->fn->unmap(device,self->input[2]);
		}
		self->input[3] = NULL;
	}
	if(unlikely (self->output == NULL) ) 
	{
		self->output = malloc(sizeof(struct buffer*) * 3);
		self->output[0] = NULL;
		self->output[1] = NULL;
		self->output[2] = NULL;
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
			size_t dims[] = { self->input[0]->dim[0], self->input[0]->dim[1] };
			size_t num_dims = 2;
			self->output[0] = (struct buffer*)self->host->fn->mem_alloc(self->host,
					sizeof(struct buffer) );
			device->fn->buffer_bld(device,self->host,self->output[0], 
					'w' ,sizeof(data_size), dims,num_dims);
		}
		if(unlikely(self->output[1] == NULL))
		{
			size_t dims[] = { self->input[0]->dim[0] };
			size_t num_dims = 1;
			self->output[1] = (struct buffer*)self->host->fn->mem_alloc(self->host,
					sizeof(struct buffer) );
			device->fn->buffer_bld(device,self->host,self->output[1], 
					'w',sizeof(data_size), dims,num_dims);
		}
		self->output[2] = NULL;
	}
}
static void worker_rmsd_process_buffers(struct worker* self)
{
	FILE* fout = ((struct worker_rmsd_properties*)self->properties)->frmsd;
	size_t i;
#if 0
	double* d;
	fout = fopen("ormsd.dat","a");
	if(fout)
	{
		//double*d = malloc(self->output[0]->size / self->output[0]->dim[0]);
		//clEnqueueReadBuffer(device_props->out,
		//		self->output[0]->data,
		//		CL_TRUE,
		//		0,
		//		self->output[0]->size / self->output[0]->dim[0],
		//		d,
		//		0, NULL, NULL);
		//free(d);
	((struct buffer_cl_properties*)self->input[0]->properties)->map_flags = CL_MAP_READ;
	double* d = self->input[0]->device->fn->map(self->input[0]->device,
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
#if 1
	//fout = fopen("rmsd.dat","a");
	if(fout)
	{
	struct device_cl_properties* device_props = self->output[1]->device->properties;
		//double*d = malloc(self->output[1]->size);
		//clEnqueueReadBuffer(device_props->out,
		//		self->output[1]->data,
		//		CL_TRUE,
		//		0,
		//		self->output[1]->size,
		//		d,
		//		0, NULL, NULL);
	//((struct buffer_cl_properties*)self->output[0]->properties)->map_flags = CL_MAP_READ;
	data_size* d = self->output[1]->device->fn->map(self->output[1]->device,
			self->output[1]); 
	for(i = 0; i < self->output[1]->dim[0] ; i++)
		//fprintf(fout,"%p %p %p\n",&d[i],&d[i+1],&d[i+2]);
		//fprintf(fout,"%x %x %x\n",(uint)d[i],(uint)d[i+1],(uint)d[i+2]);
		fprintf(fout,"%12.4f\n", d[i]);
	//fprintf(fout,"**\n");
	self->output[1]->device->fn->unmap(self->output[1]->device,self->output[1]);
		//free(d);
	//fclose(fout);
	}
#endif
#if 0
	FILE* pout = fopen("rmsf.dat","a");
	if(pout)
	{
		size_t j;
	//struct device_cl_properties* device_props = self->output[0]->device->properties;
		//double*d = malloc(self->output[0]->size);
		//clEnqueueReadBuffer(device_props->out,
		//		self->output[0]->data,
		//		CL_TRUE,
		//		0,
		//		self->output[0]->size,
		//		d,
		//		0, NULL, NULL);
	//((struct buffer_cl_properties*)self->output[0]->properties)->map_flags = CL_MAP_READ;
	double* d = self->output[0]->device->fn->map(self->output[0]->device,
			self->output[0]); 
	for(i = 0; i < self->output[0]->dim[0] ; i++)
	{
		fprintf(pout,"%zu: ",i);
		for(j = 0; j < self->output[0]->dim[1]; j++)
		//fprintf(pout,"%p %p %p\n",&d[i],&d[i+1],&d[i+2]);
		//fprintf(pout,"%x %x %x\n",(uint)d[i],(uint)d[i+1],(uint)d[i+2]);
			fprintf(pout,"%12.8f ", d[i*self->output[0]->dim[1] + j]);
		fprintf(pout,"\n");
	}
	self->output[0]->device->fn->unmap(self->output[0]->device,self->output[0]);
	//	free(d);
	fclose(pout);
	}
#endif

}

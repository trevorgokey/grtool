
#include <string.h>
#include <netcdf.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include <stdio.h>
#include "defines.h"
#include "worker_write.h"
#include "device_cl.h"
#include "driver_cl.h"
#include "buffer_cl.h"


static int worker_write_start(struct worker*);
static int worker_write_stop(struct worker*);
static void worker_write_process(struct worker*, struct module*, struct frame*);
static int worker_write_module_cl_bld(struct worker*, struct module*, struct system*);
static int worker_write_module_cl_dty(struct worker*, struct module*);
static int worker_write_module_linux_simple_bld(struct worker*, struct module*, 
		struct system*);
static int worker_write_module_linux_simple_dty(struct worker*, struct module*);
static void worker_write_prepare_buffers(struct worker*, struct module*, struct frame**);
static void worker_write_process_buffers(struct worker*);


struct worker_interface worker_write_interface =
{
	.start = worker_write_start,
	.process = worker_write_process,
	.stop = worker_write_stop,
	.prepare_buffers = worker_write_prepare_buffers,
	.process_buffers = worker_write_process_buffers,
	.dty = worker_write_dty,
	.module_linux_simple_bld = worker_write_module_linux_simple_bld,
	.module_linux_simple_dty = worker_write_module_linux_simple_dty,
	.module_cl_bld = worker_write_module_cl_bld,
	.module_cl_dty = worker_write_module_cl_dty
};

static void swap_bytes(void *bytes,void* swapped,size_t size)
{
	char swap;
	char* b = bytes;
	char* s = swapped;
	size_t idx;
	size_t siz = size -1;
	for(idx = 0; idx < size; idx++,siz--)
	{
		//s[idx] = b[siz];
		s[siz] = b[idx];
	}
	
}

int worker_write_bld(struct worker* self, struct device* host, char *filename)
{
	worker_bld(self, &worker_write_interface, host);
	struct worker_write_properties* props = malloc(sizeof(struct worker_write_properties));
	props->fd = 0;
	props->fp = 0;
	props->crd = 0;
	props->beg[0] = 0;
	props->beg[1] = 0;
	props->beg[2] = 0;
	props->end[0] = 0;
	props->end[1] = 0;
	props->end[2] = 0;
	props->filename = malloc(strlen(filename) + 1);
	strcpy(props->filename,filename);
	self->properties = props;
}
void worker_write_dty(struct worker* self)
{
	struct worker_write_properties* props = self->properties;
	nc_close(props->fd);
	if( likely (props->filename != NULL) ) free(props->filename);
	props->filename = NULL;
	int in_used[] = {-1};
	worker_clean_buffers(self,self->input,in_used);
	int out_used[] = {-1};
	worker_clean_buffers(self,self->output,out_used);
}

static void worker_write_module_linux_simple_process(struct worker* self, 
		struct module* module, struct buffer** in, struct buffer** out)
{
	struct worker_write_properties* props = self->properties;
	float* fdata;
	unsigned char* map;
	double* data;
	double* mask;
	size_t i;
	int fp = props->fp;
	long pagesize = sysconf(_SC_PAGE_SIZE);
	off_t off = (props->offset + (props->beg[0] * props->end[1] * props->end[2] * 4)) / pagesize ;
	size_t off_pad = (props->offset + (props->beg[0] * props->end[1] * props->end[2] * 4)) % pagesize ;
#if 0
		map = mmap(NULL,(props->end[0] * props->end[1] * props->end[2] * 4), 
				PROT_READ, MAP_SHARED, fp, off*pagesize);
		data = out[0]->device->fn->map(out[0]->device, out[0]);
		mask = out[1]->device->fn->map(out[1]->device, out[1]);
		memset(data,0,out[0]->size);
		memset(mask,0,out[0]->size);
		map += off_pad ;
		
		for(i = 0; i < out[0]->size / out[0]->bytes_per_elem ; i+=3)
		{
			float d = 0;
			swap_bytes(&map[4*i + 0],&d,4);
			data[i + 0] = d;
			swap_bytes(&map[4*i + 4],&d,4);
			data[i + 1] = d;
			swap_bytes(&map[4*i + 8],&d,4);
			data[i + 2] = d;
		}
		for(i = 0; i < out[1]->size / out[1]->bytes_per_elem ; i+=3)
		{
			float d = 0;
			swap_bytes(&map[4*i + 0],&d,4);
			mask[i + 0] = d;
			swap_bytes(&map[4*i + 4],&d,4);
			mask[i + 1] = d;
			swap_bytes(&map[4*i + 8],&d,4);
			mask[i + 2] = d;
		}
		munmap(map - off_pad,props->end[0] * props->end[1] * props->end[2] * 4);
#else
#if STATS
		struct timespec start,finish;
	clock_gettime(CLOCK_REALTIME, &start);	
#endif
		fdata = (float*)self->host->fn->mem_alloc(self->host,
					sizeof(float) * in[0]->size / in[0]->bytes_per_elem);
		data = in[0]->device->fn->map(in[0]->device, in[0]);
		for(i = 0; i < in[0]->dim[0] * in[0]->dim[1] * in[0]->dim[2] ; i++)
		{
			fdata[i] = data[i];
		}
#if STATS
		clock_gettime(CLOCK_REALTIME, &finish);
		self->memin_tot += (finish.tv_sec - start.tv_sec) + 
			((finish.tv_nsec - start.tv_nsec) / 1e9);
		self->memin_exe += (finish.tv_sec - start.tv_sec) + 
			((finish.tv_nsec - start.tv_nsec) / 1e9);
		clock_gettime(CLOCK_REALTIME, &start);	
#endif
		int retval = nc_put_vara_float(props->fd, props->crd, props->beg,in[0]->dim,fdata);
		nc_sync(props->fd);
#if STATS
		clock_gettime(CLOCK_REALTIME, &finish);
		self->exe_tot += (finish.tv_sec - start.tv_sec) + 
			((finish.tv_nsec - start.tv_nsec) / 1e9);
		self->exe_exe += (finish.tv_sec - start.tv_sec) + 
			((finish.tv_nsec - start.tv_nsec) / 1e9);
		//if(NC_NOERR != nc_get_vara_float(props->fd, 
		//			props->crd,props->beg,props->end,fdata))
		clock_gettime(CLOCK_REALTIME, &start);	
#endif
		in[0]->device->fn->unmap(in[0]->device,in[0]);
		//self->host->fn->mem_dealloc(self->host,fdata,
		//		sizeof(float) * out[0]->size / out[0]->bytes_per_elem);

		self->host->fn->mem_dealloc(self->host,fdata,
				sizeof(float) * in[0]->size / in[0]->bytes_per_elem);
#endif
	props->beg[0] += in[0]->dim[0];
#if STATS
		clock_gettime(CLOCK_REALTIME, &finish);
		self->proc_fin += (finish.tv_sec - start.tv_sec) + 
			((finish.tv_nsec - start.tv_nsec) / 1e9);
#endif
}

void worker_write_module_cl_process(struct worker* self, struct module* module, struct buffer** in, struct buffer** out)
{
	worker_write_module_linux_simple_process(self,module,in,out);
	return;
	/*
	struct worker_write_properties* props = 
		(struct worker_write_properties*)self->properties;
	struct device_cl_properties* device_props = 
		(struct device_cl_properties*)module->driver[0].device->properties;
	struct driver_cl_properties* driver_props = 
		(struct driver_cl_properties*)module->driver[0].properties;
	float* fdata;
	unsigned char* map;
	unsigned char* cdata;
	double* data;
	double* mask;
	size_t i;
	cl_int ret;
	size_t* sz;
	size_t dims[] = { 1024 };
	cl_uint num_dims = 1;
	long pagesize = sysconf(_SC_PAGE_SIZE);
	size_t record_pad = 13*4;
	size_t record_len = record_pad + 4 * props->dim_len[1] * props->dim_len[2];
	off_t off = (props->offset + props->beg[0] * record_len) / pagesize ;
	off_t off_pad = (props->offset + props->beg[0] * record_len) % pagesize ;
	size_t total_pad = off_pad;
	double mib, sec;
	
	struct timespec start,finish;
	cl_bool blk = CL_FALSE;
	if(likely (props->beg[0] < props->dim_len[0]))
	{
		size_t siz = props->end[0] * props->end[1] * props->end[2];
		ret = clEnqueueWriteBuffer(device_props->in,
				in[0]->data,
				blk,
				0,
				sizeof(size_t),
				&siz,
				//0, NULL, NULL);
				0, NULL, (cl_event*)&in[0]->lock);
		if(unlikely(ret != CL_SUCCESS))
		{
			printf("CPY FAIL %d\n",ret);
		}
		map = mmap(NULL,off_pad + (props->end[0] * record_len), 
				PROT_READ, MAP_SHARED|MAP_POPULATE, props->fp, off*pagesize);
		if(unlikely(map == (unsigned char*)-1))
		{
			perror("map");
		}
		map += off_pad;
		unsigned char* cpy = (unsigned char*)self->host->fn->mem_alloc(self->host,
				props->end[0]*(record_len - record_pad));
		for(i = 0; i < props->end[0]; i++)
			memcpy(cpy + (i*(record_len - record_pad)),
					map + (i*(record_len)),
					(record_len - record_pad));
		clock_gettime(CLOCK_REALTIME,&start);
#if 0

		clEnqueueWriteBuffer(device_props->in,
				in[1]->data,
				CL_FALSE,
				0,
				 props->end[0]*(record_len - record_pad),
				cpy,
				//0, NULL, NULL);
				0, NULL, (cl_event*)&in[1]->lock);
//		for(i = 0; i < props->end[0]; i++)
//		clEnqueueWriteBuffer(device_props->in,
//				in[1]->data,
//				CL_TRUE,
//				(i*(record_len - record_pad)),
//				record_len - record_pad,
//				map + (i*record_len),
//				0, NULL, NULL);
	//	clEnqueueWriteBuffer(device_props->in,
	//			in[1]->data,
	//			CL_TRUE,
	//			0,
	//			props->end[0]*(record_len - record_pad),
	//			map,
	//			0, NULL, NULL);
#else
		((struct buffer_cl_properties*)in[1]->properties)->map_flags = CL_MAP_WRITE;
		cdata = in[1]->device->fn->map(in[1]->device, in[1]);
	
		//memcpy(cdata , map, props->end[0]*(record_len - record_pad));
		
		
		//for(i = 0; i < props->end[0]; i++)
		//	memcpy(cdata + (i*(record_len - record_pad)),
		//			map + (i*(record_len)),
		//			(record_len - record_pad));


			memcpy(cdata , cpy, props->end[0]*(record_len - record_pad));
		in[1]->device->fn->unmap(in[1]->device, in[1]);
#endif
		//clock_gettime(CLOCK_REALTIME, &finish);
		
	//	FILE* fout = fopen("off.dat","a");
	//	fprintf(fout,"%zu to %zu off_pad is %zu\n",off*pagesize,off_pad + props->end[0] * props->end[1] * props->end[2] * 4,off_pad);
	//	fclose(fout);
		
		//double mib = (props->end[0] * props->end[1] * props->end[2] * 4)	/ 1e6;
		//double sec = (finish.tv_nsec - start.tv_nsec) / 1e9;
		//printf("Memory Bandwidth is %f MiB at %f = %f MiB/s\n", mib, sec , mib/sec);
		//((struct buffer_cl_properties*)in[0]->properties)->map_flags = CL_MAP_WRITE;
		
		
//		sz = in[0]->device->fn->map(in[0]->device, in[0]);
//		sz[0] = props->end[0] * props->end[1] * props->end[2];
//		//sz[1] = props->end[1];
//		//sz[2] = props->end[2]; 
//		in[0]->device->fn->unmap(in[0]->device, in[0]);
		//sz = out[0]->device->fn->map(out[0]->device, out[0]);
		//memset(sz,0,out[0]->size);
		//out[0]->device->fn->unmap(out[0]->device, out[0]);

		if(props->end[0] * props->end[1] * props->end[2] < dims[0])
			dims[0] = props->end[0] * props->end[1] * props->end[2];
		ret = clSetKernelArg(driver_props->kernel,
				0,
				sizeof(cl_mem),
				&in[1]->data);
		ret |= clSetKernelArg(driver_props->kernel,
				1,
				sizeof(cl_mem),
				&out[0]->data);
		ret |= clSetKernelArg(driver_props->kernel,
				2,
				sizeof(cl_mem),
				&in[0]->data);
		if(unlikely(ret != CL_SUCCESS))
		{
			printf("ARG FAIL %d\n",ret);
		}
		
		clWaitForEvents(1, (cl_event*)&in[0]->lock);
		//clWaitForEvents(1, (cl_event*)&in[1]->lock);
		clock_gettime(CLOCK_REALTIME,&start);
		ret = clEnqueueNDRangeKernel(device_props->exe,
				driver_props->kernel,
				num_dims,
				NULL,
				dims,
				NULL,
				0,
				NULL,
		//		NULL);
				(cl_event*)&out[0]->lock);
		if(unlikely(ret != CL_SUCCESS))
		{
			printf("EXEC ERROR: %d\n",ret);
		}
		self->host->fn->mem_dealloc(self->host,cpy,props->end[0]*(record_len - record_pad));
		munmap(map - off_pad, off_pad + props->end[0] * record_len);
		//clFinish(device_props->exe);
		clWaitForEvents(1,(cl_event*)&out[0]->lock);
		clock_gettime(CLOCK_REALTIME, &finish);
		//clReleaseEvent(in[0]->lock);
		//clReleaseEvent(out[0]->lock);
		mib = (props->end[0] * props->end[1] * props->end[2] * 4)	/ 1e6;
		sec = (finish.tv_sec - start.tv_sec) + ((finish.tv_nsec - start.tv_nsec) / 1e9);
#ifdef WRITE_BANDWIDTH
			printf("Kernel Bandwidth is %f MiB at %f = %f MiB/s\n", mib, sec , mib/sec);
#endif
		clReleaseEvent(out[0]->lock) ;out[0]->lock = 0;
		clReleaseEvent(in[0]->lock) ;in[0]->lock = 0;
	}
	props->beg[0] += props->end[0];
	out[0]->dim[0] = props->end[0];
	out[1]->dim[0] = props->end[0];
	*/
}

static int worker_write_module_linux_simple_bld(struct worker* self, struct module* module, 
		struct system* system)
{
	if(! module_linux_simple_bld(module,system))
	{
		worker_add_module(self,module,&worker_write_module_linux_simple_process);
		return 0;
	}
	else 
		return -1;
}
static int worker_write_module_linux_simple_dty(struct worker* self, struct module* module)
{

}

static int worker_write_module_cl_bld(struct worker* self, struct module* module,
		struct system* system)
{
	//decide the format to pick the right module
	char *kname = "netcdf";
	char* filename = "read.cl";
	if(! module_cl_bld(module,system,filename,kname) )
	{
		worker_add_module(self,module,&worker_write_module_cl_process);
		return 0;
	}
	else return -1;
}

static int worker_write_module_cl_dty(struct worker* self, struct module* module)
{

	return 0;
}

static int worker_write_start(struct worker* self) 
{ 
	struct worker_write_properties* props = self->properties;
	int ret = 0;
	int c = 0,i;
	int64_t o = 0;
	char offset[9];
	int second = 0;
	char swap;
	int idx = 0;
	uint64_t size = 0;
	int elem_size;
	
	//nc_close(props->fd);

    //props->dim_len[0] = 10;
	//props->beg[0] = props->dim_len[0] -10;
	return 0;
}
static void worker_write_process(struct worker* self, struct module* module, 
		struct frame* frame) { return; }
static int worker_write_stop(struct worker* self) { return WORKER_UNIMPLEMENTED; }
static void worker_write_prepare_buffers(struct worker* self, struct module* module, 
		struct frame** f)
{
	struct worker_write_properties* props = self->properties;
	struct frame* frame = f[0]; 
	struct device* device = module->driver[0].device;
	size_t d[] = { props->end[0], props->dim_len[1], props->dim_len[2] };
	const int num_dims = 3;

	if(unlikely (self->input == NULL) ) 
	{
		self->input = malloc(sizeof(struct buffer*) * 2);
		self->input[0] = NULL;
		self->input[1] = NULL;
	}
	self->input[0] = &frame->buffer[0];
	int ret;
	if(props->fd == 0)
	{
		ret = nc_create(props->filename,NC_64BIT_OFFSET|NC_WRITE|NC_FORMAT_CLASSIC,&props->fd);

		//ret |= nc_redef(props->fd);
		ret |= nc_def_dim (props->fd, "frame", NC_UNLIMITED, &props->dim[0]);
		ret |= nc_def_dim (props->fd, "atom", self->input[0]->dim[1], &props->dim[1]);
		ret |= nc_def_dim (props->fd, "spatial", 3, &props->dim[2]);
		ret |= nc_def_var (props->fd, "coordinates", NC_FLOAT, 3, props->dim, &props->crd);
		ret |= nc_put_att_text (props->fd, NC_GLOBAL, "application",5, "AMBER");
		ret |= nc_put_att_text (props->fd, NC_GLOBAL, "program",5, "AMBER");
		ret |= nc_put_att_text (props->fd, NC_GLOBAL, "programVersion",3, "0.1");
		ret |= nc_put_att_text (props->fd, NC_GLOBAL, "Conventions",5, "AMBER");
		ret |= nc_put_att_text (props->fd, NC_GLOBAL, "ConventionVersion",3, "1.0");
		ret |= nc_enddef(props->fd);
		if(ret != NC_NOERR)
		{
			printf("TRAJAMBER: Could not create file properly\n");
		}
	}

	if(unlikely (self->output == NULL) ) self->output = malloc(sizeof(struct buffer*) * 1);
	if(likely   (self->output != NULL) )
	{
		self->output[0] = NULL;
	}

	//frame[0]->beg[0] = props->beg[0];
	//frame[0]->beg[1] = props->beg[1];
	//frame[0]->beg[2] = props->beg[2];
	//frame[0]->end[0] = props->end[0];
	//frame[0]->end[1] = props->end[1];
	//frame[0]->end[2] = props->end[2];
	
}

static void worker_write_process_buffers(struct worker* self)
{

#if 0
	size_t i;
	//struct device_cl_properties* device_props = self->output[0]->device->properties;
	FILE* fout = fopen("write.dat","a");
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
	//((struct buffer_cl_properties*)self->output[0]->properties)->map_flags = CL_MAP_READ;
	double* d = self->input[0]->device->fn->map(self->input[0]->device,
			self->input[0]);
	for(i = 0; i < self->input[0]->dim[0] * self->input[0]->dim[1] * 3; i+=3)
		//fprintf(fin,"%p %p %p\n",&d[i],&d[i+1],&d[i+2]);
		//fprintf(fin,"%x %x %x\n",(uint)d[i],(uint)d[i+1],(uint)d[i+2]);
		fprintf(fout,"%zu: %f %f %f\n",i/3,d[i],d[i+1],d[i+2]);
	fprintf(fout,"**\n");
	self->input[0]->device->fn->unmap(self->input[0]->device,self->input[0]);
	fclose(fout);
	}
#endif
}

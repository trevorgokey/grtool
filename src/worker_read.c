
#include <string.h>
#include <netcdf.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include <stdio.h>
#include "defines.h"
#include "worker_read.h"
#include "driver_linux_simple.h"
#include "driver_cl.h"


#if (MAX_FRAMES > 0)
extern const int max_frames;
extern volatile int frames;
extern pthread_mutex_t lock;
#endif

static int worker_read_start(struct worker*);
static int worker_read_stop(struct worker*);
static void worker_read_process(struct worker*, struct module*, struct frame*);
static int worker_read_module_cl_bld(struct worker*, struct module*, struct system*);
static int worker_read_module_cl_dty(struct worker*, struct module*);
static int worker_read_module_linux_simple_bld(struct worker*, struct module*, 
		struct system*);
static int worker_read_module_linux_simple_dty(struct worker*, struct module*);
static void worker_read_prepare_buffers(struct worker*, struct module*, struct frame**);
static void worker_read_process_buffers(struct worker*);


struct worker_interface worker_read_interface =
{
	.start = worker_read_start,
	.process = worker_read_process,
	.stop = worker_read_stop,
	.prepare_buffers = worker_read_prepare_buffers,
	.process_buffers = worker_read_process_buffers,
	.dty = worker_read_dty,
	.module_linux_simple_bld = worker_read_module_linux_simple_bld,
	.module_linux_simple_dty = worker_read_module_linux_simple_dty,
	.module_cl_bld = worker_read_module_cl_bld,
	.module_cl_dty = worker_read_module_cl_dty
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

int worker_read_bld(struct worker* self, struct device* host, char *filename)
{
	worker_bld(self, &worker_read_interface, host);
	struct worker_read_properties* props = malloc(sizeof(struct worker_read_properties));
	props->filename = malloc(strlen(filename) + 1);
	strcpy(props->filename,filename);
	self->properties = props;
}

void worker_read_dty(struct worker* self)
{
	struct worker_read_properties* props = self->properties;
	close(props->fp);
	nc_close(props->fd);
	if( likely (props->filename != NULL) ) free(props->filename);
	props->filename = NULL;
	int in_used[] = {0,1,-1};
	worker_clean_buffers(self,self->input,in_used);
}

static void worker_read_module_linux_simple_process(struct worker* self, 
		struct module* module, struct buffer** in, struct buffer** out)
{
	struct worker_read_properties* props = self->properties;
	float* fdata;
	unsigned char* map;
	data_size* data;
	data_size* mask;
	size_t i;
	int fp = props->fp;
#if STATS
	struct timespec start,finish;
#endif
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
		clock_gettime(CLOCK_REALTIME, &start);
#endif
		if(sizeof(data_size) == 8)
		{
			fdata = (float*)self->host->fn->mem_alloc(self->host,
						sizeof(float) * out[0]->size / out[0]->bytes_per_elem);
			data = out[0]->device->fn->map(out[0]->device, out[0]);
		}
		else
		{
			fdata = out[0]->device->fn->map(out[0]->device, out[0]);
		}
		//memset(fdata,0,	sizeof(float) * out[0]->size / out[0]->bytes_per_elem);
		//memset(data,0,out[0]->size );
#if STATS
		clock_gettime(CLOCK_REALTIME, &finish);
		self->memin_tot += (finish.tv_sec - start.tv_sec) + 
			((finish.tv_nsec - start.tv_nsec) / 1e9);
		clock_gettime(CLOCK_REALTIME, &start);
#endif
		if(NC_NOERR != nc_get_vara_float(props->fd, 
					props->crd,props->beg,out[0]->dim,fdata))
			printf("FAIL TO READ\n");
		if (sizeof(data_size) == 8)
		{
			for(i = 0; i < out[0]->dim[0] * out[0]->dim[1] * out[0]->dim[2] ; i++)
			{
				data[i] = fdata[i];
			}
		}
#if STATS
		clock_gettime(CLOCK_REALTIME, &finish);
		self->exe_tot += (finish.tv_sec - start.tv_sec) + 
			((finish.tv_nsec - start.tv_nsec) / 1e9);
		clock_gettime(CLOCK_REALTIME, &start);
#endif
		out[0]->device->fn->unmap(out[0]->device,out[0]);
		//self->host->fn->mem_dealloc(self->host,fdata,
		//		sizeof(float) * out[0]->size / out[0]->bytes_per_elem);

		//data = out[1]->device->fn->map(out[1]->device, out[1]);
		////fdata = (float*)self->host->fn->mem_alloc(self->host,
		////			sizeof(float) * out[1]->size / out[1]->bytes_per_elem);
		//if(NC_NOERR != nc_get_vara_double(props->fd, 
		//			props->crd,props->beg,props->end,data))
		//	printf("FAIL TO READ\n");
		//else
		//	for(i = 0; i < out[1]->size / out[1]->bytes_per_elem ; i+=3)
		//	{
		//		data[i + 0] = fdata[i + 0];
		//		data[i + 1] = fdata[i + 1];
		//		data[i + 2] = fdata[i + 2];
		//	}
		//out[1]->device->fn->unmap(out[1]->device,out[1]);
		if(sizeof(data_size) == 8)
		self->host->fn->mem_dealloc(self->host,fdata,
				sizeof(float) * out[0]->size / out[0]->bytes_per_elem);
#endif
	props->beg[0] += props->end[0];
	out[0]->dim[0] = props->end[0];
	out[1]->dim[0] = props->end[0];
#if STATS	
	clock_gettime(CLOCK_REALTIME, &finish);
		self->proc_fin += (finish.tv_sec - start.tv_sec) + 
			((finish.tv_nsec - start.tv_nsec) / 1e9);
#endif
}

void worker_read_module_cl_process(struct worker* self, struct module* module, struct buffer** in, struct buffer** out)
{
	//worker_read_module_linux_simple_process(self,module,in,out);
	//return;
	struct worker_read_properties* props = 
		(struct worker_read_properties*)self->properties;
	struct device_cl_properties* device_props = 
		(struct device_cl_properties*)module->driver[0].device->properties;
	struct driver_cl_properties* driver_props = 
		(struct driver_cl_properties*)module->driver[0].properties;
	float* fdata;
	static int first = 0;
	unsigned char* map;
	unsigned char* cdata;
	data_size* data;
	data_size* mask;
	size_t i;
	cl_int ret;
	size_t* sz;
	size_t siz = props->end[0] * props->end[1] * props->end[2];
	size_t tlen = READ_BLK;
	size_t dims[] = { tlen*(siz/tlen) } ;
	size_t ldims[] = { tlen };
	cl_uint num_dims = 1;
#if STATS
	unsigned long long time = 0;
	struct timespec s,f;
	struct timespec start,finish;
	double mib, sec;
#endif
	long pagesize = sysconf(_SC_PAGE_SIZE);
	size_t record_pad = 13*4;
	size_t record_len = record_pad + 4 * props->dim_len[1] * props->dim_len[2];
	off_t off = (props->offset + props->beg[0] * record_len) / pagesize ;
	off_t off_pad = (props->offset + props->beg[0] * record_len) % pagesize ;
	size_t total_pad = off_pad;
	cl_ulong res;
	//clGetDeviceInfo(device_props->id, CL_DEVICE_PROFILING_TIMER_RESOLUTION, sizeof(cl_ulong), &res, NULL);
	cl_bool blk = !CL_ASYNC;
	register_profile_event(self,"memin_ocl_start",self->frame_cnt,0);
		if(first == 0)
		{
#if USE_MAP
			sz = in[0]->device->fn->map(in[0]->device, in[0]);
			sz[0] = props->end[0] * props->end[1] * props->end[2];
			in[0]->device->fn->unmap(in[0]->device, in[0]);
#else
			ret = clEnqueueWriteBuffer(device_props->in,
					in[0]->data,
					blk,
					0,
					sizeof(size_t),
					&siz,
					0, NULL, 
#if CL_ASYNC
					(cl_event*)&in[0]->lock
#else
					NULL
#endif
					);
			if(unlikely(ret != CL_SUCCESS))
			{
				printf("READ: CPY FAIL %d\n",ret);
			}
#endif
		ret |= clSetKernelArg(driver_props->kernel,
				2,
				sizeof(cl_mem),
				&in[0]->data);
		}
#if STATS
		//register_profile_event(self,"memin_ocl_wait_start",self->frame_cnt,0);
		clock_gettime(CLOCK_REALTIME, &s);
#endif
		((struct buffer_cl_properties*)in[1]->properties)->map_flags = CL_MAP_WRITE;
#if USE_MAP
		//((struct buffer_cl_properties*)in[1]->properties)->async = CL_TRUE;
		cdata = in[1]->device->fn->map(in[1]->device, in[1]);
		//((struct buffer_cl_properties*)in[1]->properties)->async = CL_FALSE;
#endif
		
		map = mmap(NULL,off_pad + (props->end[0] * record_len), 
				PROT_READ, MAP_SHARED|MAP_POPULATE, props->fp, off*pagesize);
		if(unlikely(map == (unsigned char*)-1))
		{
			perror("map");
		}
		map += off_pad;
#if (USE_MAP==0)
		unsigned char* cpy = (unsigned char*)self->host->fn->mem_alloc(self->host,
				props->end[0]*(record_len - record_pad));
#if SINGLE_TX
		    register_profile_event(self,"memin_local_map_start",self->frame_cnt,0);
#endif
		for(i = 0; i < props->end[0]; i++)
        {
#if !SINGLE_TX
#if STATS
		clock_gettime(CLOCK_REALTIME, &s);	
#endif
            char str[64];
            //snprintf(str,64,"%05zu_memin_ocl_wait_start",i);
            snprintf(str,64,"memin_ocl_wait_start",i);
		    register_profile_event(self,str,self->frame_cnt,0);
		    clEnqueueWriteBuffer(device_props->in,
		    		in[1]->data,
		    		blk,
				    (i*(record_len - record_pad)),
				    record_len - record_pad,
				    map + (i*record_len),
		    		0,NULL,
#if CL_ASYNC || STATS
				(cl_event*)&in[1]->lock
#else
				NULL
#endif
				);
            
		    time = (s.tv_sec * 1e9) + (s.tv_nsec);
#if CL_ASYNC || STATS
	    	clWaitForEvents(1, (cl_event*)&in[1]->lock);
#endif
		    self->memin_wait += device_cl_get_queued_time_seconds(in[1]->lock);
		    time +=  (unsigned long long)
		    	(device_cl_get_queued_time_seconds(in[1]->lock)*1e9);
            //snprintf(str,64,"%05zu_memin_ocl_delay_start",i);
            snprintf(str,64,"memin_ocl_delay_start",i);
		    register_profile_event(self,str,self->frame_cnt,&time);
		    time +=  (unsigned long long)
		    	(device_cl_get_delay_time_seconds(in[1]->lock)*1e9);
            //snprintf(str,64,"%05zu_memin_ocl_exe_start",i);
            snprintf(str,64,"memin_ocl_exe_start",i);
		    register_profile_event(self,str,self->frame_cnt,&time);
#if (CL_ASYNC || STATS) && CLEAN_UP
        	if(in[1]->lock) clReleaseEvent(in[1]->lock) ;in[1]->lock = 0;
#endif
#else
			memcpy(cpy + (i*(record_len - record_pad)),
					map + (i*(record_len)),
					(record_len - record_pad));
#endif
        }
#if (SINGLE_TX == 1)
#if STATS
		clock_gettime(CLOCK_REALTIME, &s);	
#endif
		register_profile_event(self,"memin_ocl_wait_start",self->frame_cnt,0);
		clEnqueueWriteBuffer(device_props->in,
				in[1]->data,
				blk,
				0,
				 props->end[0]*(record_len - record_pad),
				cpy,
				0,NULL,
#if CL_ASYNC || STATS
				(cl_event*)&in[1]->lock
#else
				NULL
#endif
				);
#else // SINGLE_TX
		    register_profile_event(self,"proc_ocl_kern_pre_start",self->frame_cnt,0);
#endif
		//clFlush(device_props->in);
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
#if STATS
		clock_gettime(CLOCK_REALTIME, &f);
		self->memin_exe += (f.tv_sec - s.tv_sec) + 
			((f.tv_nsec - s.tv_nsec) / 1e9);
#endif
#else // MAP == 0
	
		//memcpy(cdata , map, props->end[0]*(record_len - record_pad));
		

#if CL_ASYNC || STATS
		clWaitForEvents(1, (cl_event*)&in[1]->lock);
#endif
		//self->memin_exe += device_cl_get_exe_time_seconds(in[1]->lock);
		register_profile_event(self,"memin_ocl_exe_start",self->frame_cnt,0);
		for(i = 0; i < props->end[0]; i++)
			memcpy(cdata + (i*(record_len - record_pad)),
					map + (i*(record_len)),
					(record_len - record_pad));
			//memcpy(cdata , cpy, props->end[0]*(record_len - record_pad));
		in[1]->device->fn->unmap(in[1]->device, in[1]);
		register_profile_event(self,"proc_ocl_kern_pre_start",self->frame_cnt,0);
#endif //USE_MAP==0
#if STATS && (USE_MAP == 0) && SINGLE_TX
		clock_gettime(CLOCK_REALTIME, &f);
		self->memin_exe += (f.tv_sec - s.tv_sec) + 
			((f.tv_nsec - s.tv_nsec) / 1e9);
		time = (s.tv_sec * 1e9) + (s.tv_nsec);
		self->memin_wait += device_cl_get_queued_time_seconds(in[1]->lock);
		time +=  (unsigned long long)
			(device_cl_get_queued_time_seconds(in[1]->lock)*1e9);
		register_profile_event(self,"memin_ocl_delay_start",self->frame_cnt,&time);
		time +=  (unsigned long long)
			(device_cl_get_delay_time_seconds(in[1]->lock)*1e9);
		register_profile_event(self,"memin_ocl_exe_start",self->frame_cnt,&time);
		time +=  (unsigned long long)
			(device_cl_get_exe_time_seconds(in[1]->lock)*1e9);
		register_profile_event(self,"proc_ocl_kern_pre_start",self->frame_cnt,&time);

#endif //STATS
		//clock_gettime(CLOCK_REALTIME, &finish);
		
	//	FILE* fout = fopen("off.dat","a");
	//	fprintf(fout,"%zu to %zu off_pad is %zu\n",off*pagesize,off_pad + props->end[0] * props->end[1] * props->end[2] * 4,off_pad);
	//	fclose(fout);
		
		//double mib = (props->end[0] * props->end[1] * props->end[2] * 4)	/ 1e6;
		//double sec = (finish.tv_nsec - start.tv_nsec) / 1e9;
		//printf("Memory Bandwidth is %f MiB at %f = %f MiB/s\n", mib, sec , mib/sec);
		//((struct buffer_cl_properties*)in[0]->properties)->map_flags = CL_MAP_WRITE;
		
		
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
		if(unlikely(ret != CL_SUCCESS))
		{
			printf("ARG FAIL %d\n",ret);
		}
#if CL_ASYNC
		if(!first) clWaitForEvents(1, (cl_event*)&in[0]->lock);
#endif
#if STATS
		self->memin_tot += (finish.tv_sec - start.tv_sec) + 
			((finish.tv_nsec - start.tv_nsec) / 1e9);
		if(!first)
		{
			//self->memin_exe += device_cl_get_exe_time_seconds(in[0]->lock);
			self->memin_wait += device_cl_get_queued_time_seconds(in[0]->lock);
			self->memin_delay += device_cl_get_delay_time_seconds(in[0]->lock);
		}
		//device_cl_print_kernel_info(driver_props->kernel,device_props->id);

		clock_gettime(CLOCK_REALTIME, &start);
#endif
		register_profile_event(self,"proc_ocl_kern_wait_start",self->frame_cnt,0);
		ret = clEnqueueNDRangeKernel(device_props->exe,
				driver_props->kernel,
				num_dims,
				NULL,
				dims,
				ldims,
				0, NULL,
				(cl_event*)&out[0]->lock
                                );
		//clFinish(device_props->exe);
		if(unlikely(ret != CL_SUCCESS))
		{
			printf("EXEC ERROR: %d\n",ret);
		}
		first = 1;
		//res = 1;
	clWaitForEvents(1,(cl_event*)&out[0]->lock);
#if STATS
	clock_gettime(CLOCK_REALTIME, &finish);
	time = (start.tv_sec * 1e9) + (start.tv_nsec);
	time +=  (unsigned long long)
		(device_cl_get_queued_time_seconds(out[0]->lock)*1e9);
	register_profile_event(self,"proc_ocl_kern_delay_start",self->frame_cnt,&time);
	time +=  (unsigned long long)
		(device_cl_get_delay_time_seconds(out[0]->lock)*1e9);
	register_profile_event(self,"proc_ocl_kern_exe_start",self->frame_cnt,&time);
	time +=  (unsigned long long)
		(device_cl_get_exe_time_seconds(out[0]->lock)*1e9);
	time = (finish.tv_sec * 1e9) + (finish.tv_nsec);
	register_profile_event(self,"proc_fin",self->frame_cnt,&time);
	

	self->exe_exe += device_cl_get_exe_time_seconds(out[0]->lock);
	self->exe_wait += device_cl_get_queued_time_seconds(out[0]->lock);
	self->exe_delay += device_cl_get_delay_time_seconds(out[0]->lock);
		self->exe_tot += (finish.tv_sec - start.tv_sec) + 
			((finish.tv_nsec - start.tv_nsec) / 1e9);
	clock_gettime(CLOCK_REALTIME,&start);
#endif
	munmap(map - off_pad, off_pad + props->end[0] * record_len);
	self->host->fn->mem_dealloc(self->host,cpy,props->end[0]*(record_len - record_pad));
	props->beg[0] += props->end[0];
	out[0]->dim[0] = props->end[0];
	out[1]->dim[0] = props->end[0];
		
#if CLEAN_UP
	if(out[0]->lock)clReleaseEvent(out[0]->lock); out[0]->lock = 0;
#endif
#if CL_ASYNC && CLEAN_UP
	if(in[0]->lock) clReleaseEvent(in[0]->lock) ;in[0]->lock = 0;
#endif
#if (CL_ASYNC || STATS) && CLEAN_UP
	if(in[1]->lock) clReleaseEvent(in[1]->lock) ;in[1]->lock = 0;
#endif
#if STATS
	clock_gettime(CLOCK_REALTIME, &finish);
	self->proc_fin += (finish.tv_sec - start.tv_sec) + 
		((finish.tv_nsec - start.tv_nsec) / 1e9);
#endif
}

static int worker_read_module_linux_simple_bld(struct worker* self, struct module* module, 
		struct system* system)
{
	if(! module_linux_simple_bld(module,system))
	{
		worker_add_module(self,module,&worker_read_module_linux_simple_process);
		return 0;
	}
	else 
		return -1;
}
static int worker_read_module_linux_simple_dty(struct worker* self, struct module* module)
{

}

static int worker_read_module_cl_bld(struct worker* self, struct module* module,
		struct system* system)
{
	//decide the format to pick the right module
	char *kname = "netcdf";
	char* filename = "read.cl";
	if(! module_cl_bld(module,system,filename,kname) )
	{
		worker_add_module(self,module,&worker_read_module_cl_process);
		return 0;
	}
	else return -1;
}

static int worker_read_module_cl_dty(struct worker* self, struct module* module)
{

	return 0;
}

static int worker_read_start(struct worker* self) 
{ 
	struct worker_read_properties* props = self->properties;
	int ret = 0;
	FILE* fp = fopen(props->filename, "rb");
	int c = 0,i;
	int64_t o = 0;
	char offset[9];
	int second = 0;
	char swap;
	int idx = 0;
	uint64_t size = 0;
	int elem_size;
	if(fp != NULL)	
	{
		do
		{
			
			if((fgetc(fp) == 'r') && (fgetc(fp) == 'o') && (fgetc(fp) == 'm'))
			{
				second++;
				if( second == 1 )
				{
					fseek(fp,4,SEEK_CUR);
					char disp[8];
					elem_size = 4;
					ret = fread(&offset,elem_size,1,fp);
					swap_bytes(offset,&size,elem_size);
					elem_size = 8;
					ret = fread(&offset,elem_size,1,fp);
					swap_bytes(offset,&o,elem_size);
					break;
				}
			}
			fseek(fp,-1,SEEK_CUR);
		}
		while( (c = fgetc(fp)) != EOF);
		fclose(fp);
		props->offset = o;
	}
	
	props->fp = open(props->filename,O_RDONLY);
	
	ret = nc_open(props->filename,NC_NOWRITE,&props->fd);

	ret |= nc_inq_dimid  ( props->fd, "frame"  , props->dim + 0);	
	ret |= nc_inq_dimid  ( props->fd, "atom"   , props->dim + 1);	
	ret |= nc_inq_dimid  ( props->fd, "spatial", props->dim + 2);	
	ret |= nc_inq_dimlen ( props->fd, props->dim[0], props->dim_len + 0);
	ret |= nc_inq_dimlen ( props->fd, props->dim[1], props->dim_len + 1);
	ret |= nc_inq_dimlen ( props->fd, props->dim[2], props->dim_len + 2);

	ret |= nc_inq_varid  ( props->fd, "coordinates", &props->crd);
	//nc_close(props->fd);

#ifdef LAST
	if(LAST <= props->dim_len[0])
		props->dim_len[0] = LAST;
#endif
	props->beg[0] = props->beg[1] = props->beg[2] = 0;
	props->end[0] = (STRIDE > props->dim_len[0]) ? props->dim_len[0] : STRIDE;
	props->end[1] = props->dim_len[1];
	props->end[2] = props->dim_len[2];
	//props->beg[0] = props->dim_len[0] -10;
	return 0;
}
static void worker_read_process(struct worker* self, struct module* module, 
		struct frame* frame) { return; }
static int worker_read_stop(struct worker* self) { return WORKER_UNIMPLEMENTED; }
static void worker_read_prepare_buffers(struct worker* self, struct module* module, 
		struct frame** frame)
{
	static size_t id = 0;
	struct worker_read_properties* props = self->properties;
	*frame = (struct frame*)self->host->fn->mem_alloc(self->host,
				sizeof(struct frame));
#if (MAX_FRAMES > 0)
#if STATS
	register_profile_event(self,"frame_wait",id+1,NULL);
#endif
	while(frames == max_frames);
#if STATS
	register_profile_event(self,"prepare_cont",id+1,NULL);
#endif
	pthread_mutex_lock(&lock);
	frames++;
	pthread_mutex_unlock(&lock);
#endif
	frame[0]->id = id;
	id++;

	frame[0]->buffer = (struct buffer*)self->host->fn->mem_alloc(self->host,
				sizeof(struct buffer) * 2);
	frame[0]->buffer_size = 1;
	struct device* device = module->driver[0].device;
	if(unlikely (self->input == NULL) ) 
	{
		self->input = malloc(sizeof(struct buffer*) * 3);
		self->input[0] = NULL;
		self->input[1] = NULL;
		self->input[2] = NULL;
	}

	size_t d[] = { props->end[0], props->dim_len[1], props->dim_len[2] };
	const int num_dims = 3;
	if(likely   (self->input != NULL) )
	{
		size_t dims[] = { 3 };
		if(unlikely(self->input[0] == NULL))
		{
			self->input[0] = (struct buffer*)self->host->fn->mem_alloc(self->host,
					sizeof(struct buffer) );
			device->fn->buffer_bld(device,self->host,self->input[0], 
					'r' + 'p',sizeof(size_t), dims,1);
		}
		if(unlikely(self->input[1] == NULL))
		{
			int mode = 'r';
#if USE_PIN
			mode += 'p';
#endif
			self->input[1] = (struct buffer*)self->host->fn->mem_alloc(self->host,
					sizeof(struct buffer) );
			device->fn->buffer_bld(device,self->host,self->input[1], 
					mode  ,sizeof(data_size), d,num_dims);
		}
		self->input[2] = NULL;
	}
	if(unlikely (self->output == NULL) ) self->output = malloc(sizeof(struct buffer*) * 3);
	if(likely   (self->output != NULL) )
	{
		self->output[0] = &frame[0]->buffer[0];
		self->output[1] = &frame[0]->buffer[0];
		self->output[2] = NULL;
	}
	if (props->beg[0] + props->end[0] >= props->dim_len[0])
		props->end[0] = props->dim_len[0] - props->beg[0];
	frame[0]->len[0] = props->dim_len[0];
	frame[0]->len[1] = props->dim_len[1];
	frame[0]->len[2] = props->dim_len[2];
	frame[0]->beg[0] = props->beg[0];
	frame[0]->beg[1] = props->beg[1];
	frame[0]->beg[2] = props->beg[2];
	frame[0]->end[0] = props->end[0];
	frame[0]->end[1] = props->end[1];
	frame[0]->end[2] = props->end[2];
	
	int mode = 'r' + 'w';
	device->fn->buffer_bld(device,self->host,self->output[0], 
			mode,sizeof(data_size), d,num_dims);
	//device->fn->buffer_bld(device,self->host,self->output[1], 
	//		mode,sizeof(double), d,num_dims);
	if(props->end[0] != 0)
	{
		frame[0]->status = 0;
	}
	else
	{
		frame[0]->status = 1;
	}
}

static void worker_read_process_buffers(struct worker* self)
{
#if 0
	size_t i;
	//struct device_cl_properties* device_props = self->output[0]->device->properties;
	FILE* fout = fopen("out.dat","a");
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
	data_size *d = self->output[0]->device->fn->map(self->output[0]->device,
			self->output[0]);
	for(i = 0; i < self->output[0]->dim[0] * self->output[0]->dim[1] * 3; i+=3)
		//fprintf(fout,"%p %p %p\n",&d[i],&d[i+1],&d[i+2]);
		//fprintf(fout,"%x %x %x\n",(uint)d[i],(uint)d[i+1],(uint)d[i+2]);
		fprintf(fout,"%zu: %f %f %f\n",i/3,d[i],d[i+1],d[i+2]);
	fprintf(fout,"**\n");
	self->output[0]->device->fn->unmap(self->output[0]->device,self->output[0]);
	fflush(fout);
	fclose(fout);
	}
#endif
}

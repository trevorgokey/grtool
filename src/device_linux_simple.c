
//#include <sys/resource.h>
#include <assert.h>
#include<stdio.h>
#include <unistd.h>
#include "device_linux_simple.h"

void device_linux_simple_dty(struct device*);
size_t get_memory_size(struct device*);

uintptr_t device_linux_simple_mem_alloc(struct device*, size_t);
void device_linux_simple_mem_dealloc(struct device*, void*,size_t);
int   device_linux_simple_buffer_bld(struct device*, struct device*, struct buffer*, int,
		size_t, size_t*, size_t);
void*  device_linux_simple_map(struct device*, struct buffer*);
void   device_linux_simple_unmap(struct device*, struct buffer*);
void device_linux_simple_buffer_dty(struct device*, struct device*,struct buffer*);
void device_linux_simple_deallocator(uintptr_t);

struct device_interface device_linux_simple_interface =
{
	.deallocator = device_linux_simple_deallocator,
	.buffer_bld = device_linux_simple_buffer_bld,
	.buffer_dty = device_linux_simple_buffer_dty,
	.map = device_linux_simple_map,
	.unmap = device_linux_simple_unmap,
	.mem_alloc = device_linux_simple_mem_alloc,
	.mem_dealloc = device_linux_simple_mem_dealloc,
	.dty = device_linux_simple_dty
};

void* device_linux_simple_allocator(size_t memsize)
{
	return malloc(memsize);
}
void device_linux_simple_deallocator(uintptr_t mem)
{
	free((void*)mem);
}

void*  device_linux_simple_map(struct device* self, struct buffer* buffer)
{
	return buffer->data;
}
void  device_linux_simple_unmap(struct device* self, struct buffer* buffer)
{

}

int device_linux_simple_bld(struct device* self, unsigned long long *memsize)
{
	unsigned long long mem;
	if(memsize == NULL)
	{
			mem = sysconf(_SC_PAGESIZE) * sysconf(_SC_PHYS_PAGES);
	}
	else mem = *memsize;
	self->properties = NULL;
	device_bld(self,&device_linux_simple_interface,&mem);
}

void device_linux_simple_dty(struct device* self)
{

}
size_t get_memory_size(struct device* self)
{
	unsigned long long ps = sysconf(_SC_PAGESIZE);
	//unsigned long long pn = sysconf(_SC_AVPHYS_PAGES);
	unsigned long long max = ps * sysconf(_SC_PHYS_PAGES);
	return max;
}

uintptr_t device_linux_simple_mem_alloc(struct device* self,size_t memsize)
{
	void* mem = NULL;
#if (USE_RECYCLE==1)
        mem = (void*)memspace__pop__bytes_(&self->memspace,memsize);
#endif
	if(mem == NULL)
		return (uintptr_t)malloc(memsize);
	//void* m = malloc(memsize);
	//if(m == NULL)
	//{
	//	perror("memory alloc error");
	//	assert(0);
	//}
	//return (uintptr_t)m;
}
void device_linux_simple_mem_dealloc(struct device* self,void* mem,size_t memsize)
{
#if (USE_RECYCLE==1)
	memspace__push__bytes_(&self->memspace,mem,memsize);
#else
	free(mem);
#endif
}
int  device_linux_simple_buffer_bld(
	struct device* self,
	struct device* host,
	struct buffer* buf,
	int mode,
	size_t bytes_per_elem,
	size_t* dim,
	size_t dim_size)
{
	buffer_bld(
		self,
		host,
		buf,
		bytes_per_elem,
		dim,
		dim_size);
	struct device_linux_simple_properties* device_linux_simple = self->properties;
	//struct buffer_properties* buf_props =
	//	host->fn->mem_alloc(host,sizeof(struct buffer_cl_properties));
	buf->data = (void*)self->fn->mem_alloc(self,buf->size);
	if(buf->data == NULL)
		buf->data = malloc(buf->size);
	buf->device = self;
	buf->properties = NULL;

	if(buf->data != NULL) return 0;
	else return -1;
}

void device_linux_simple_buffer_dty(struct device* self, struct device* host, struct buffer* buffer)
{

	//memspace__push__bytes_(&self->memspace,buffer->data,buffer->size);
}

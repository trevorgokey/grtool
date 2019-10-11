
#include <stdio.h>
#include <string.h>
#include "defines.h"
#include "device.h"

int buffer_bld(struct device* self,struct device* host, struct buffer* buffer,
		size_t bytes_per_element,size_t* dim,size_t dim_size)
{
	size_t i, size;
	buffer->dim = (size_t*)host->fn->mem_alloc(host,sizeof(size_t) * dim_size);
	buffer->device = self;
	for(i = 0; i < dim_size; i++)
	{
		buffer->dim[i] = dim[i];
	}
	buffer->dim_size = dim_size;
	buffer->valid = 0;
	buffer->lock = NULL;
	buffer->properties = NULL;
	buffer->data = NULL;
	for(i = 0,size = 1; i < dim_size; i++)
		size *= dim[i];
	buffer->size = size * bytes_per_element; 
	buffer->bytes_per_elem = bytes_per_element;
	return 0;
}

void buffer_dty(struct device* device, struct buffer* buffer)
{
	device->fn->mem_dealloc(device,buffer->dim,sizeof(size_t) * buffer->dim_size);
	buffer->device->fn->buffer_dty(buffer->device,device,buffer);
	buffer->device->fn->mem_dealloc(buffer->device,buffer->data,buffer->size);
	memset(buffer,0,sizeof(struct buffer));
}

int device_bld(struct device* self, struct device_interface* interface, unsigned long long *mem)
{
	//printf("HELLO: device\n");
	self->fn = interface;
	//self->device_properties_bld(self);
	unsigned long long memsize = (mem == NULL ? self->fn->get_memory_size(self) : *mem);
	self->memspace = memspace__bld__with_size_(memsize);
	self->memspace.deallocator = interface->deallocator; 
	return 0;
}

void device_dty(struct device* self)
{
	memspace__dty_(&self->memspace);
	if( likely(self->fn != NULL) ) self->fn->dty(self);
	self->fn = NULL;
	if(self->properties) free(self->properties);
	self->properties = NULL;
	//printf("BYE: device\n");
}

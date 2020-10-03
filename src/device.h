
#ifndef DEVICE_H
#define DEVICE_H

#include <stddef.h>
#include <stdlib.h>
#include "memspace_.h"
#include "buffer.h"
#include "defines.h"

struct device
{
	memspace_ memspace;
	size_t mem_size;
	struct device_interface* fn;
	void* properties;
};
struct device_interface
{
	void   (*deallocator)(uintptr_t);
	void   (*dty)(struct device*);
	int    (*buffer_bld)(struct device*, struct device*, struct buffer*, int, 
			size_t, size_t*, size_t);
	void*  (*map)(struct device*, struct buffer*);
	void   (*unmap)(struct device*, struct buffer*);

	void   (*buffer_dty)(struct device*, struct device*, struct buffer*);
	uintptr_t  (*mem_alloc)  (struct device*, size_t);
	void   (*mem_dealloc)(struct device*, void*,size_t);
	unsigned long long (*get_memory_size)(struct device*);
};

struct buffer_interface
{
	void (*dty)(struct device*, struct device*, struct buffer*);
};

int buffer_bld(struct device*, struct device*, struct buffer*, 
			size_t, size_t*, size_t);

void buffer_dty(struct device*, struct buffer*);

int  device_bld(struct device*, struct device_interface*, unsigned long long*);
void device_dty(struct device*);

#endif /* DEVICE_H */

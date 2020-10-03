#ifndef BUFFER_H
#define BUFFER_H

#include <stddef.h>
#include <stdlib.h>

struct device;

struct buffer
{
	void* data;
	size_t  size;
	unsigned int bytes_per_elem;
	size_t* dim;
	size_t  dim_size;
	char valid;
	struct device* device;
	void* lock;
	struct buffer_interface* fn;
	void* properties;
};



#endif //BUFFER_H

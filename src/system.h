
#ifndef SYSTEM_H
#define SYSTEM_H

#include <stddef.h>
#include "device.h"
#include "buffer.h"


typedef enum { NOERR, INVALID, ERR_NO_SYSTEM_AVAILABLE } status;


struct system
{
	struct system_interface* fn;
	
	struct device* device;
	size_t device_size;
	void*  properties;
};
struct system_interface
{
	void   (*dty)(struct system*);
};


int  system_bld(struct system*, struct system_interface*, size_t);
void system_dty(struct system*);

#endif /* SYSTEM_H */


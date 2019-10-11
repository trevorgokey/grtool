
#ifndef DRIVER_H
#define DRIVER_H

#include "device.h"

struct driver
{
	struct device* device;
	struct driver_interface* fn;
	void*  properties;

};
struct driver_interface
{
	void (*dty)(struct driver*);
};


int  driver_bld(struct driver*, struct device*);
void driver_dty(struct driver*);

#endif /* DRIVER_H */

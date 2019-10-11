
#ifndef DRIVER_LINUX_SIMPLE_H
#define DRIVER_LINUX_SIMPLE_H

#include "driver.h"
#include "device_linux_simple.h"

struct driver_linux_simple_properties
{

};

int  driver_linux_simple_bld(struct driver*, struct device*);
void driver_linux_simple_dty(struct driver*);

#endif /* DRIVER_LINUX_SIMPLE_H */

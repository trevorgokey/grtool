
#include "driver_linux_simple.h"

static struct driver_interface interface = 
{
	.dty = driver_linux_simple_dty
};

int  driver_linux_simple_bld(struct driver* self, struct device* host)
{
	driver_bld(self,host);
	self->fn = &interface;
	return 0;
}
void driver_linux_simple_dty(struct driver* self)
{

}


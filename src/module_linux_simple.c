
#include "module_linux_simple.h"
#include "driver_linux_simple.h"


static struct module_interface interface = 
{
	.dty = module_linux_simple_dty
};

int  module_linux_simple_bld(struct module* self, struct system* system)
{
	size_t i;
	module_bld(self,&interface,system);
	for(i = 0; i < self->driver_size; i++)
	{
		driver_linux_simple_bld(&self->driver[i],&system->device[i]);
	}
	return 0;
}
void module_linux_simple_dty(struct module* self)
{
	
}


#include <unistd.h>
#include <stdio.h>
#include "system_linux_simple.h"

void system_linux_simple_dty(struct system*);

struct system_interface system_linux_simple_interface = 
{
	.dty = system_linux_simple_dty
};


int system_linux_simple_bld(struct system* self, size_t* memsize)
{
	size_t device_size = 1;
	system_bld(self,&system_linux_simple_interface,device_size);
	self->properties = NULL;
	size_t i;
	for(i = 0, self->device_size = 0; i < device_size; ++i)
	{
		if(! device_linux_simple_bld(&self->device[self->device_size], memsize) )
		{
			printf("LINUX_SIMPLE: Device %zu: %p online\n",i,&self->device[self->device_size]);
			self->device_size++;
		}
		else
		{
			printf("LINUX_SIMPLE: Device %zu: %p FAILED\n",i,&self->device[self->device_size]);
			device_dty(&self->device[self->device_size]);
		}
	}
}

void system_linux_simple_dty(struct system* self)
{

}




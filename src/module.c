
#include <string.h>
#include <stdio.h>
#include "defines.h"
#include "module.h"

int  module_bld(struct module* self, struct module_interface* interface, struct system* system)
{
	memset(self,0,sizeof(struct module));
	//printf("HELLO: module\n");
	self->fn = interface;
	self->driver = malloc(system->device_size * sizeof(struct driver));
	if( unlikely(self->driver == NULL) ) return -1;
	self->driver_size = system->device_size;
	self->properties = NULL;
	self->driver_next_idx = 0;
	return 0;
}
void module_dty(struct module* self)
{
	size_t i;
	//printf("BYE: module\n");
	if(self->driver)
	{
		for(i = 0; i < self->driver_size; i++)
		{
			driver_dty(&self->driver[i]);
		}
		free(self->driver);
		self->driver = NULL;
		self->driver_size = 0;
	}
	self->fn->dty(self);
	if(self->properties) free(self->properties);
	self->properties = NULL;
}

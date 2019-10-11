
#include <stdio.h>
#include "system.h"
#include "defines.h"


int  system_bld(struct system* self, struct system_interface* interface, size_t num_devices)
{
	printf("HELLO: system\n");
	self->fn = interface;
	
	size_t num_bytes = sizeof(struct device) * num_devices;
	self->device = malloc(num_bytes);
	
	if(self->device == NULL)
		perror("ERROR: could not allocate memory for devices\n");
	return 0;
}
void system_dty(struct system* self)
{
	size_t i;
	if(self->fn != NULL)
		self->fn->dty(self);
	self->fn = NULL;
	if(self->properties) free(self->properties);
	self->properties = NULL;
	if(self->device)
	{
		for(i = 0; i < self->device_size; i++)
			device_dty(&self->device[i]);
		free(self->device);
		self->device = NULL;
		self->device_size = 0;
	}
	printf("BYE: system\n");
}

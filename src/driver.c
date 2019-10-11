
#include "driver.h"
#include <stdio.h>

int  driver_bld(struct driver* self, struct device* device)
{	
	//printf("HELLO: driver\n");
	self->device = device;
	self->properties = NULL;
	return 0;
}

void driver_dty(struct driver* self)
{
	self->device = NULL;
	if(self->fn != NULL) 
		self->fn->dty(self);
	self->fn = NULL;
	if(self->properties != NULL) free(self->properties);
	self->properties = NULL;
	//printf("BYE: driver\n");
}

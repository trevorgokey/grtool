
#include <stdio.h>
#include <string.h>
#include "defines.h"
#include "frame.h"

void frame_dty(struct frame* self, struct device* device)
{
	size_t i;
	for(i = 0; i < self->buffer_size; i++)
	{
		buffer_dty(device,&self->buffer[i]);
	}
	device->fn->mem_dealloc(device, self->buffer, sizeof(struct buffer) * self->buffer_size);
	memset(self,0,sizeof(struct frame));
}

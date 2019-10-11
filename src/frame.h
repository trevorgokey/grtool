
#ifndef FRAME_H
#define FRAME_H

#include <stddef.h>
#include "buffer.h"
#include "device.h"

struct frame
{
	struct buffer* buffer;
	size_t buffer_size;
	size_t id;
	unsigned int beg[3];
	unsigned int end[3];
	unsigned int len[3];
	int status;
};

void frame_dty(struct frame*, struct device*);

#endif

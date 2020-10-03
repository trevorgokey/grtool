
#ifndef MODULE_H
#define MODULE_H

#include "system.h"
#include "driver.h"
#include "buffer.h"

struct worker;

struct module
{
	struct module_interface* fn;	
	void (*process) (struct worker*, struct module*, struct buffer**, struct buffer**);
	struct system* system;
	struct driver* driver;
	size_t driver_size;
	unsigned int   driver_next_idx;
	void*  properties;
};


struct module_interface
{
	void (*start)  (struct module*);
	void (*stop)   (struct module*);
	void (*dty)    (struct module*);
};


int  module_bld(struct module*, struct module_interface*, struct system*);
void module_dty(struct module*);

#endif /* MODULE_H */

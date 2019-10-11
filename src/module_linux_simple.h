
#ifndef MODULE_LINUX_SIMPLE_H
#define MODULE_LINUX_SIMPLE_H

#include "module.h"
#include "system_linux_simple.h"


struct module_linux_simple_properties
{

};


int  module_linux_simple_bld(struct module*, struct system*);
void module_linux_simple_dty(struct module*);

#endif /* MODULE_LINUX_SIMPLE_H */

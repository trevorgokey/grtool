
#ifndef SYSTEM_LINUX_SIMPLE_H
#define SYSTEM_LINUX_SIMPLE_H

#include "system.h"

struct system_linux_simple_properties
{

};

int system_linux_simple_bld(struct system*, size_t*);
size_t system_host_mem_alloc(size_t);
void   system_host_mem_dealloc(void*);

#endif

#ifndef MEMSPACE_
#define MEMSPACE_

#include <stdint.h>
#include <stddef.h>
#include "buffer.h"
#include "queue_.h"


typedef struct
{
	queue_ queue;
	size_t size;
#ifdef STATS_INSTRUMENTATION
	size_t min, max, mode, pops, pushes;
#endif
} 
size_bin_;

typedef struct
{
	size_bin_* size_bin;
	uintptr_t (*allocator)(size_t);
	void  (*deallocator)(uintptr_t);
	pthread_mutex_t lock;
	size_t num_bins;
	size_t max;
	size_t used;
	size_t avail;
	size_t total;
#ifdef STATS_INSTRUMENTATION
	size_t min, max, mode, pops, pushes;
#endif
}
memspace_;

memspace_  memspace__bld__with_size_(size_t m);
void       memspace__dty_(memspace_ *m);
size_bin_* memspace__add_bin_(memspace_ *m,size_bin_ *bin);
size_bin_* memspace__find_bin_of_size(memspace_* m,size_t sz);
void       memspace__push__bytes_(memspace_* m, void*, size_t);
uintptr_t      memspace__pop__bytes_(memspace_* m, size_t);

#endif

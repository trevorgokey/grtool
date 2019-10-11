#ifndef QUEUE_
#define QUEUE_

#include <pthread.h>
#include <stdint.h>

typedef struct
{
	uintptr_t* item;
	pthread_mutex_t lock;
	pthread_cond_t not_empty;
	pthread_cond_t not_full;
	size_t max,cty,size,front,end,err;
#ifdef STATS_INSTRUMENTATION
	size_t min, max, mode, pops, pushes;
#endif
}
queue_ ;

queue_ queue__bld_();
queue_ queue__bld__with_size_(unsigned int m);

void   queue__dty_(queue_*);
void   queue__dty__n_items_(queue_ *, size_t);
void   queue__dty__all_items_(queue_ *, void (*)() );
void   queue__fill_n_items_(queue_ *, size_t, size_t);

void   queue__extend_(queue_*, size_t);
uintptr_t queue__pop_(queue_*,char);
uintptr_t queue__pop__wait_(queue_*);
uintptr_t queue__pop__nowait_(queue_*);
void   queue__push_(queue_*, uintptr_t );

#endif

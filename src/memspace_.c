#include <string.h>
#include <stdio.h>
#include <stddef.h>
#include<assert.h>
#include "memspace_.h"

memspace_ memspace__bld__with_size_(size_t m)
{
	printf("HELLO: Memspace with %f MiB\n",m/1000000.0);
	return (memspace_) 
	{ 
		  .size_bin = 0
		, .lock = PTHREAD_MUTEX_INITIALIZER
		, .num_bins = 0
		, .avail = m
		, .total = 0
		, .used = 0
		, .max = m };
}

void memspace__dty_(memspace_ *m)
{
	size_t i;
	size_t bytes = 0;
	size_t arrays = 0;
	pthread_mutex_lock(&m->lock);
	if(m->size_bin)
	{
		for(i = 0; i < m->num_bins ; i++)
		{
			arrays += m->size_bin[i].queue.size;
			bytes += m->size_bin[i].queue.size * m->size_bin[i].size;
			queue__dty__all_items_(&m->size_bin[i].queue,m->deallocator);
			queue__dty_(&m->size_bin[i].queue);
		}
		free(m->size_bin);
	}
	printf("BYE: Memspace %p (%f) MiB cleared %zu arrays totalling %f of %f MiB bytes, reuse ratio %.2e\n",
			m,m->max/1000000.0,arrays,bytes/ 1000000.0, m->total / 1000000.0,(double)m->total/(double)bytes * 100.);
	pthread_mutex_unlock(&m->lock);
	pthread_mutex_destroy(&m->lock);
	memset(m,0,sizeof(memspace_));
}
size_bin_* memspace__add_bin_(memspace_ *m,size_bin_ *bin)
{
	pthread_mutex_lock(&m->lock);
	size_t new_sz = sizeof(size_bin_)*(m->num_bins + 1);
	size_bin_* nptr = realloc(m->size_bin,new_sz);
	if(nptr != 0)
	{
		m->size_bin = nptr;
		m->size_bin[m->num_bins] = *bin;
		free(bin);
		bin = &m->size_bin[m->num_bins];
		m->num_bins++;
	}
	pthread_mutex_unlock(&m->lock);
	return bin;
}

size_bin_* memspace__find_bin_of_size(memspace_* m,size_t sz)
{
	size_t i;
	size_bin_* b = 0;
	pthread_mutex_lock(&m->lock);
	for(i = 0; i < m->num_bins; i++)
		if(m->size_bin[i].size == sz)
		{
			b = &m->size_bin[i];
			break;
		}
	pthread_mutex_unlock(&m->lock);
	return b;
}

void memspace__push__bytes_(memspace_* m, void* mem, size_t size)
{
	//m->deallocator((uintptr_t)mem); return;
	size_bin_ *bin = memspace__find_bin_of_size(m,size);
	if(bin == 0)
	{
		bin = malloc(sizeof(size_bin_));
		//bin->queue = queue__bld__with_size_(3);
		bin->queue = queue__bld_();
		//bin->queue.max = 0;
		bin->size = size;
		bin = memspace__add_bin_(m,bin);
	}
	queue__push_(&bin->queue,(uintptr_t)mem);

}

uintptr_t memspace__pop__bytes_(memspace_* m, size_t size)
{
	//return 0;
	//return (size_t)m->allocator(flgs);
	uintptr_t mem = 0;
	size_bin_ *bin = memspace__find_bin_of_size(m,size);
	if(bin == 0)
	{
		bin = malloc(sizeof(size_bin_));
		bin->queue = queue__bld_();
		//bin->queue = queue__bld__with_size_(3);
		bin->size = size;
		bin = memspace__add_bin_(m,bin);
	}
	mem = queue__pop__nowait_(&bin->queue);
	if(mem == 0)
	{
		if (size <= m->avail)
		{
			//mem = m->allocator(b);
			m->used += size;
			m->avail -= size;
		}
		else if( bin->queue.cty == 0 )
		{
			//printf("DEBUG: ADDED %zu BYTES\n",size - m->avail);
			//m->avail = 0;
			//m->used += size;
			//m->max += size - m->avail;
			
			printf("FATAL: CANNOT ALLOCATE %zu BYTES\n",size);
            assert(0);
			//mem = queue__pop__wait_(&bin->queue);
		}
		else
			mem = queue__pop__wait_(&bin->queue);
	}
	m->total += size;
	return mem;
}

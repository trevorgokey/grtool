#include "queue_.h"
#include "stdlib.h"
#include "stdio.h"
#include "assert.h"

queue_ queue__bld_()
{
	//printf("HELLO: Queue with unlimited capacity\n");
	return (queue_) { 
		  .item = NULL
		, .lock = PTHREAD_MUTEX_INITIALIZER
		, .not_full = PTHREAD_COND_INITIALIZER
		, .not_empty = PTHREAD_COND_INITIALIZER
		, .max = 0
		, .cty = 0
		, .size = 0
		, .front = 0
		, .end = 0
		, .err = 0
	};
}
queue_ queue__bld__with_size_(unsigned int m)
{
	//printf("HELLO: Queue with %d capacity\n",m);
	void *ptr = malloc(sizeof(void*) * m);
	char err = 0;
	if(ptr == 0)
		err = 1;
	return (queue_) { 
		  .item = ptr
		, .lock = PTHREAD_MUTEX_INITIALIZER
		, .not_full = PTHREAD_COND_INITIALIZER
		, .not_empty = PTHREAD_COND_INITIALIZER
		, .max = m
		, .cty = m
		, .size = 0
		, .front = 0
		, .end = 0
		, .err = err
	};
}
void queue__dty_(queue_ *q)
{
	if(q->item)
		free(q->item);
	q->item = 0;
	q->max = 0;
	q->cty = 0;
	q->size = 0;
	q->err = 1;
	pthread_mutex_destroy(&q->lock);
	pthread_cond_destroy(&q->not_full);
	pthread_cond_destroy(&q->not_empty);
}
void queue__extend_(queue_* q, size_t cty)
{
	assert(q->err == 0);
	//pthread_mutex_lock(&q->lock);
	//if(q->max > 0)
	//	cty = cty > q->max ? q->max : cty;
	if(cty > q->cty)
	{
		size_t* nptr = malloc(sizeof(size_t)*cty);
		if(nptr != 0)
		{
			//printf("Notice: extending queue to size %zu items (%p to %p)\n",cty, q->item, nptr);
			size_t i = 0;
			for(; i < q->size; i++)
			{
				nptr[i] = q->item[(q->front + i) % q->cty];
			}
			q->cty = cty;
			if(q->item) free(q->item);
			q->item = nptr;
			q->front = 0;
			q->end = q->size;
		}
		else
		{
			printf("FAIL\n");
			assert(0);
		}
	}
	//pthread_mutex_unlock(&q->lock);
}
uintptr_t queue__pop_(queue_* q,char wait)
{
	uintptr_t val = 0;
	pthread_mutex_lock(&q->lock);
	if(q->size == 0)
	{
		if(wait)
		{
			while(q->size == 0)
				pthread_cond_wait(&q->not_empty,&q->lock);
		}
		else
			goto finish;
	}
	//printf("POP: took spot %d\n", q->front);
	val = q->item[q->front];
	//printf("POP %f on %p (sz %zu) HAS PTRS: ", *(double*)val,q->item,q->size);
	//fflush(stdout);
	//size_t i = 0;
	//for(i = 0; i < q->size; i++)
	//{
	//	printf(" %zu: %p",i,(void*)q->item[(q->front + i) % q->cty]);
	//	fflush(stdout);
	//}
	//printf(" REMOVING %p\n", (void*)val);
	//fflush(stdout);
	q->size--;
	q->front = (q->front + 1) % q->cty;
	//printf("POP: front is %d end is %d size is %d cty is %d\n",
	//		q->front,
	//		q->end,
	//		q->size,
	//		q->cty);
finish:
	pthread_cond_signal(&q->not_full);
	pthread_mutex_unlock(&q->lock);
	return val;
}
uintptr_t queue__pop__wait_(queue_* q)
{
	return queue__pop_(q,1);
}
uintptr_t queue__pop__nowait_(queue_* q)
{
	return queue__pop_(q,0);
}
void queue__push_(queue_* q, uintptr_t val)
{
	pthread_mutex_lock(&q->lock);
//restart:
	if(q->size == q->cty)
	{
		if(q->max > 0)
		{
			while(q->size == q->cty)
				pthread_cond_wait(&q->not_full,&q->lock);
		}
		else
		{
			queue__extend_(q,(q->cty+1)*2);
			//goto restart;
		}
	}
	//printf("PUSH: put on spot %d\n,", q->end);
	//size_t i = 0;
	//printf("PUS %f on %p (sz %zu) HAS PTRS: ", *(double*)val,q->item,q->size);
	//fflush(stdout);
	//for(; i < q->size; i++)
	//{
	//	printf(" %zu: %p",i,(void*)q->item[(q->front + i) % q->cty]);
	//	fflush(stdout);
	//}
	//printf(" PUSHED %p\n", val);
	//fflush(stdout);
	q->item[q->end] = (size_t)val;
	q->end = (q->end+1)%q->cty;
	q->size++;
	//printf("PUSH: front is %zu end is %zu size is %zu cty is %zu\n",
	//		q->front,
	//		q->end,
	//		q->size,
	//		q->cty);
	//if( q->cty < q->max )
	//{
	//	size_t new_size = (q->cty * 2) > q->max ? q->max : (q->cty * 2);
	//	queue__extend_(q, new_size);
	//}
	pthread_cond_signal(&q->not_empty);
	pthread_mutex_unlock(&q->lock);
}
void queue__dty__n_items_(queue_ *q, size_t num)
{
	size_t i;
	for(i = 0; i < num; i++)
	{
		size_t v = queue__pop__nowait_(q);
		if(v != 0)
			free((void*)v);
	}
}
void queue__dty__all_items_(queue_ *q, void (*deallocator)() )
{
	size_t v = 0;
	while( (v = queue__pop__nowait_(q)) != 0)
	{
		deallocator((void*)v);
	}
}
void queue__fill_n_items_(queue_ *q, size_t num, size_t bytes)
{
	size_t i = 0;
	for(i = 0; i < num; i++)
	{
		uintptr_t v = (uintptr_t)malloc(bytes);
		if(v != 0)
			queue__push_(q,v);
	}
}

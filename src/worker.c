
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "defines.h"
#include "worker.h"

static int worker_counter = 1;
#if (MAX_FRAMES > 0)
const int max_frames = MAX_FRAMES;
volatile int frames = 0;
pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER;
#endif

int worker_bld(struct worker* self, struct worker_interface* interface, struct device* host)
{
	//printf("HELLO: worker\n");
	self->id = worker_counter++;
	register_profile_event(self,"bld_start",0,NULL);
	self->frame_cnt = 0;
	self->host = host;
	self->fn = interface;
	self->name = NULL;
	self->name_size = 0;
	self->in = self->out = NULL;
	self->err = 0;
	self->module = NULL;
	self->module_size = 0;
	self->properties = NULL;
	self->thread_launched = false;
#if STATS
	self->events = NULL;
	   self->tot=
	   self->wait_in= self->wait_out= 
	   self->prep= self->proc= self->proc_buf= self->proc_fin=
	   self->memin_tot= self->memin_wait= self->memin_exe=
	   self->exe_tot= self->exe_wait= self->exe_exe= 
	   self->memout_tot= self->memout_wait= self->memout_exe = 
	   self->memin_delay = self->exe_delay = 0.0;
#endif
	//register_profile_event(self,"bld_end",0,NULL);
	return 0;
}

struct module* worker_pick_module(struct worker* self)
{
	return &self->module[0];
}

static void print_iter(size_t beg, size_t len, size_t end)
{
	unsigned int tenpercent = (unsigned int)((float)end*0.1);
	unsigned int twopercent = (unsigned int)((float)end*0.01);
	unsigned int fivepercent = (unsigned int)((float)end*0.05);
	size_t i = 0;
	while(i < len)
	{
		if((beg + len) == end)
		{
			printf("100%%\n");
			break;
		}
		else if(tenpercent && ((i + beg) % tenpercent) == 0 )
			printf("%.f%%",(float)(beg + i)/end * 100);
		else if(fivepercent && ((i + beg) % fivepercent) == 0 )
			printf("%.f%%",(float)(beg + i)/end * 100);
		else if(twopercent && ((i + beg) % twopercent) == 0 )
			printf(".");
		fflush(stdout);
		i++;
	}
}

void* worker_start(struct worker* self)
{
#if STATS
	struct timespec s,t,start,finish;
#endif
	unsigned long long i = 0;
	int done = 0;
	//register_profile_event(self,"start",0,NULL);
	//printf("WORKER %d:LAUNCHED\n", self->id);

#if STATS
	clock_gettime(CLOCK_REALTIME, &s);	
#endif
	struct frame* frame;
	if( likely(self->fn->start != NULL) ) self->fn->start(self);
	while(! done )
	{
		i = ++self->frame_cnt;
		
#if STATS
		register_profile_event(self,"pop_start",i,NULL);
#endif
		if( likely(self->in != NULL) )
		{
#if STATS
			clock_gettime(CLOCK_REALTIME, &start);	
#endif
			frame = (struct frame*)queue__pop__wait_(self->in);
#if STATS
			clock_gettime(CLOCK_REALTIME, &finish);
			self->wait_in += (finish.tv_sec - start.tv_sec) + 
				((finish.tv_nsec - start.tv_nsec) / 1e9);
#endif
		}
		struct module* module = worker_pick_module(self);
		
#if STATS
		register_profile_event(self,"prepare_start",i,NULL);
		clock_gettime(CLOCK_REALTIME, &start);	
#endif
		self->fn->prepare_buffers(self,module,&frame);
#if STATS
		clock_gettime(CLOCK_REALTIME, &finish);
		self->prep += (finish.tv_sec - start.tv_sec) + 
			((finish.tv_nsec - start.tv_nsec) / 1e9);
#endif
		done = frame->status;
#if STATS
		register_profile_event(self,"proc_start",i,NULL);
#endif
		if(likely(done == 0))
		{
#if STATS
			clock_gettime(CLOCK_REALTIME, &start);	
#endif
			module->process(self, module, self->input, self->output);
#if STATS
			clock_gettime(CLOCK_REALTIME, &finish);
			self->proc += (finish.tv_sec - start.tv_sec) + 
				((finish.tv_nsec - start.tv_nsec) / 1e9);
			clock_gettime(CLOCK_REALTIME, &start);	
#endif
			self->fn->process_buffers(self);
#if STATS
			register_profile_event(self,"proc_buf_start",i,NULL);
			clock_gettime(CLOCK_REALTIME, &finish);
			self->proc_buf += (finish.tv_sec - start.tv_sec) + 
				((finish.tv_nsec - start.tv_nsec) / 1e9);
#endif
		}
#if STATS
		register_profile_event(self,"push_start",i,NULL);
#endif
		if( likely(self->out != NULL))
		{
#if STATS
			clock_gettime(CLOCK_REALTIME, &start);	
#endif
			queue__push_(self->out,(uintptr_t)frame);
#if STATS
			clock_gettime(CLOCK_REALTIME, &finish);
			self->wait_out += (finish.tv_sec - start.tv_sec) + 
				((finish.tv_nsec - start.tv_nsec) / 1e9);
#endif
		}
		else
		{
			if(likely(done == 0) )
				print_iter(frame->beg[0], frame->end[0], frame->len[0]);
#if (MAX_FRAMES > 0)
			while(frames == 0);
			pthread_mutex_lock(&lock);
			frames--;
			pthread_mutex_unlock(&lock);
#endif
			frame_dty(frame,self->host);
			self->host->fn->mem_dealloc(self->host,frame,sizeof(struct frame));
		}
#if STATS
		register_profile_event(self,"push_end",i,NULL);
#endif
	}
	if( likely(self->fn->stop != NULL) ) self->fn->stop(self);
	
#if STATS
	clock_gettime(CLOCK_REALTIME, &t);
	self->tot = (t.tv_sec - s.tv_sec) + 
		((t.tv_nsec - s.tv_nsec) / 1e9);
#endif

	return NULL;
}

void worker_launch(struct worker* self)
{
	self->thread_launched = true;
#if STATS
	register_profile_event(self,"worker_start",0,NULL);
#endif
	pthread_create(&self->thread,NULL,(void* (*)(void*))worker_start,self);
}

void worker_add_module(struct worker* self, struct module* module, 
		void (*process)(struct worker*, struct module*, struct buffer**, struct buffer**))
{
	struct module* modules = realloc(self->module, 
			sizeof(struct module) * (self->module_size + 1));
	if( likely( modules != NULL ))
	{
		self->module = modules;
		self->module[self->module_size] = *module;
		self->module[self->module_size].process = process;
		self->module_size++;
	}
}

void worker_join(struct worker* self)
{
	if(self->thread_launched) pthread_join(self->thread,NULL);
	self->thread_launched = false;
}

static void* clean_array(struct buffer** b)
{
	if(b)
	{
		free(b);
	}
	return NULL;
}

int worker_connect(struct worker* out,struct worker*in , int* size)
{
	int sz;
	if(size == NULL)
		sz = 2;
	else
		sz = *size;

	queue_* q = malloc(sizeof(queue_));
	*q = queue__bld__with_size_(sz);
	out->out = in->in = q;
}

void worker_clean_buffers(struct worker* self, struct buffer** b, int* idx)
{
	int i = 0;
	if(idx == NULL) return;
	struct buffer* buffer;
	while(idx[i] != -1)
	{
		buffer = b[idx[i]];
		if(buffer)
			buffer_dty(self->host,buffer);
		buffer = NULL;
		self->host->fn->mem_dealloc(self->host, buffer, sizeof(struct buffer));
		i++;
	}
}

void print_profile_event(struct profile_event* e, FILE *fp)
{
	fprintf(fp,"%-16s %08u %-32s %16llu\n", e->worker_name,e->frame_id,e->event_name,
			e->time);
}

void worker_dty(struct worker* self)
{
	//register_profile_event(self,"dty_start",0,NULL);
	static int header = 0;
	if( likely(self->fn != NULL) ) self->fn->dty(self);
	self->fn = NULL;
	
	if(self->in) queue__dty_(self->in); free(self->in);
	self->in = NULL;
	self->out = NULL;

	if( likely(self->module != NULL) ) free(self->module); 
	self->module = NULL;
	if( likely(self->name != NULL) ) free(self->name); 
	self->name = 0;
	self->name_size = 0;
	if(self->properties) free(self->properties);
	self->properties = NULL;
	
	self->input  = clean_array(self->input);
	self->output = clean_array(self->output);
#if STATS
	if(header == 0 && (header = 1))
	printf("DONE -- %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s\n",
	       "tot",
		   "wait_in", "wait_out", 
		   "prep", "proc", "proc_buf", "proc_fin",
		   "memin_tot", "memin_wait", "memin_delay","memin_exe",
		   "exe_tot", "exe_wait", "exe_delay","exe_exe", 
		   "memout_tot", "memout_wait", "memout_exe");
	printf("DONE -- %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f\n",
		   self->tot,
		   self->wait_in, self->wait_out, 
		   self->prep, self->proc, self->proc_buf, self->proc_fin,
		   self->memin_tot, self->memin_wait, self->memin_delay, self->memin_exe,
		   self->exe_tot, self->exe_wait, self->exe_delay, self->exe_exe, 
		   self->memout_tot, self->memout_wait, self->memout_exe);
	int i = 0;
	//register_profile_event(self,"dty_end",0,NULL);
	if(self->events != NULL)
	{
		FILE *fp = fopen("events.gant","a");
		while(self->events[i] != NULL)
		{
			print_profile_event(self->events[i],fp);
			free(self->events[i]);
			self->events[i] = NULL;
			i++;
		}
		free(self->events);
		fclose(fp);
	}
#endif
}

void register_profile_event(struct worker* self, char* event_name,unsigned long id, unsigned long long* time)
{
#if STATS
	struct timespec start;
	unsigned long t = 0;
	clock_gettime(CLOCK_REALTIME,&start);
	int i = 0;
	if(self->events == NULL)
	{
		self->events = malloc(sizeof(void*));
		self->events[0] = NULL;
	}

	while(self->events[i] != NULL) i++;
	struct profile_event** ptr = realloc(self->events,sizeof(void*)*(i+2));
	if(ptr != NULL) self->events = ptr;
	self->events[i+1] = NULL;
	self->events[i] = malloc(sizeof(struct profile_event));
	int len;
	len = strlen(event_name);
	if(len > sizeof(self->events[i]->event_name))
		len = sizeof(self->events[i]->event_name);
	strncpy(self->events[i]->event_name, event_name,31);
	self->events[i]->event_name[31] = '\0';
	char name[16];
	snprintf(name,15,"%u",self->id);
	strncpy(self->events[i]->worker_name, name,15);
	self->events[i]->worker_name[15] = '\0';
	self->events[i]->frame_id = id;
	if(time == NULL)
		t = (start.tv_sec *1e9) + (start.tv_nsec );
	else
		t = *time;
	self->events[i]->is_end = 0;

	self->events[i]->time = t;
#endif
}




#ifndef WORKER_H
#define WORKER_H

//#include "config_list.h"
#include "frame.h"
#include "buffer.h"
#include "queue_.h"
#include "device.h"
#include "module.h"

#define WORKER_UNIMPLEMENTED -1
#define WORKER_NORMAL         0

struct profile_event
{
	int is_end;
	unsigned int frame_id;
	char event_name[32];
	char worker_name[16];
	unsigned long long time;
};



struct worker
{
	char*  name;
	size_t name_size;
	struct worker_interface* fn;
	int    id;
	int    err;
	struct device* host;
	struct buffer** input;
	struct buffer** output;
	unsigned int frame_cnt;
	queue_* in;
	queue_* out;
	pthread_t thread;
	int thread_launched;

#if STATS
	double tot,
		   wait_in, wait_out, 
		   prep, proc, proc_buf, proc_fin,
		   memin_tot, memin_wait, memin_delay, memin_exe,
		   exe_tot, exe_wait, exe_delay,exe_exe, 
		   memout_tot, memout_wait, memout_exe;
	struct profile_event** events;
#endif
	struct module* module;
	size_t module_size;
	void*  properties;
};


struct worker_interface
{
	int  (*start)(struct worker*);
	void (*process)(struct worker*, struct module*, struct frame*);
	int  (*stop)(struct worker*); 
	void (*dty)(struct worker*);
	void (*prepare_buffers)(struct worker*, struct module*, struct frame**);
	void (*process_buffers)(struct worker*);
	int  (*module_linux_simple_bld)  (struct worker*, struct module*, struct system*);
	int  (*module_linux_simple_dty)  (struct worker*, struct module*);
	int  (*module_cl_bld)  (struct worker*, struct module*, struct system*);
	int  (*module_cl_dty)  (struct worker*, struct module*);
};

int  worker_bld(struct worker*, struct worker_interface*, struct device*);
void worker_dty(struct worker*);
void worker_launch(struct worker*);
void worker_join(struct worker*);
void worker_add_module(struct worker*, struct module*, void (*)(struct worker*, struct module*, struct buffer**, struct buffer**));

void worker_clean_buffers(struct worker*, struct buffer**, int*);

void register_profile_event(struct worker*, char*,unsigned long, unsigned long long*);

int worker_connect(struct worker*,struct worker*, int*);
//int worker_launch(struct worker*);
//int worker_status(struct worker*);
//
//void worker_print_messages(struct worker*);
//void worker_print_finish(struct worker*);

#endif

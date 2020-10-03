#include <stdio.h>
#include "defines.h"
#include "system_cl.h"
#include "system_linux_simple.h"
#include "worker.h"
#include "worker_translate.h"
#include "worker_read.h"
#include "worker_write.h"
#include "worker_align.h"


int main(int argc, char** argv)
{
	int use_cl;
	size_t i, k;
	size_t num_workers = 5;
	struct system host;
	unsigned long long mem = 700000;
	unsigned long long* memptr = &mem;
	memptr = NULL;
	system_linux_simple_bld(&host,memptr);
	struct system system[16];
	struct worker* worker = calloc(num_workers, sizeof(struct worker));
	queue_* queue = calloc(num_workers, sizeof(queue_));
	struct module* module;
	cl_uint num_cl_platforms = 0;
	int total_platforms;
	cl_platform_id platform[2];
	struct timespec start;
	clock_gettime(CLOCK_REALTIME, &start);
#if STATS
	FILE *fp = fopen("events.gant","w");
	fprintf(fp,"%-16s %08u %-32s %16llu\n", "0",0,"prog_start",
			(unsigned long long)(start.tv_sec *1e9 + start.tv_nsec));
#endif
	if(argc > 1)
		use_cl = atoi(argv[1]);
	else
		use_cl = USE_CL;
	if(use_cl)
	{
		cl_int ret = 0;
		//ret = clGetPlatformIDs(1, 0, &num_cl_platforms);
		num_cl_platforms = 2;
		ret |= clGetPlatformIDs(num_cl_platforms, platform, &num_cl_platforms);
		printf("Notified %d CL platforms\n",num_cl_platforms);
		if(num_cl_platforms > 1) num_cl_platforms = 1;
		if( unlikely (ret != CL_SUCCESS) )
		{
			printf("Error: clGetPlatformIDs returned %d\n",ret);
		}
		for(i = 0; i < num_cl_platforms; i++)
		{
			system_cl_bld(&system[i],platform[i]);
			printf("SYSTEM %zu: %zu devices\n", i, system[i].device_size);
		}
	}
	total_platforms = num_cl_platforms + 1;
	module = calloc((1 + num_workers) * (total_platforms), sizeof(struct module));
	i = 0;
	printf("Starting reader\n");
	worker_read_bld(&worker[i++],&host.device[0], TRAJECTORY_LOCATION);
#if TRANS 
	printf("Starting translate worker\n");
	worker_translate_bld(&worker[i++],&host.device[0]);
#endif
#if ALIGN
	printf("Starting align worker\n");
	worker_align_bld(&worker[i++],&host.device[0]);
#endif
#if RMSD
	printf("Starting RMSD worker\n");
	worker_rmsd_bld(&worker[i++],&host.device[0]);
#endif
#if WRITE
	printf("Starting writer\n");
	worker_write_bld(&worker[i++],&host.device[0],"test.nc");
#endif
	num_workers = i;
	int size[] = {WORKER_QUEUE_SIZE};
	for(i = 0; i < num_workers - 1; i++)
		worker_connect(&worker[i],&worker[i+1],size);
	for(i = 0; i < num_workers; i++)
	{
		for(k = 0; k < num_cl_platforms; k++)
		{
			worker[i].fn->module_cl_bld(&worker[i], &module[i*total_platforms + k],
					&system[k]);
		}
		worker[i].fn->module_linux_simple_bld(&worker[i], &module[i*total_platforms + k],
				&host);
	}

	for(i = 0; i < num_workers; i++)
		worker_launch(&worker[i]);
	for(i = 0; i < num_workers; i++)
		worker_join(&worker[i]);
	
	for(i = 0; i < num_workers; i++)
	{
		worker_dty(&worker[i]);
		//queue__dty_(&queue[i]);
		for(k = 0; k < total_platforms; k++)
		{
			module_dty(&module[i*total_platforms + k]);
		}
	}
	
#if STATS
	fclose(fp);
#endif

	free(module);
	free(worker);
	free(queue);
	for(i = 0; i < num_cl_platforms; i++)
	{
		system_dty(&system[i]);
	}
	system_dty(&host);
//	printf("SYSTEMCL is status %d with cl_status %d\n" scl->base.err,scl->cl_err);
	//if( ! scl.err )
	//	systems.add(scl);
	
//	system_native snative = build_system_native();
//	if( ! scl.err )
//		systems.add(scl);
//
//	if( systems.empty() ) 
//	{
//		ret = ERR_NO_SYSTEM_AVAILABLE;
//		goto clean_systems;
//	}
//	
//	ret = 0;
//	if( ret = parse_config(filename,config) ) 
//		goto clean_systems;
//	
//	if( ret = build_workers(config,systems,workers)
//		goto clean_workers;
//	
//	if ( ret = connect_workers(workers,config) )
//		goto clean_systems;
//	
//	if ( ret = build_modules(workers,systems) )
//		goto clean_modules;
//	
//	if ( ret = launch_workers(workers) )
//		goto clean_modules
//
//	while( ret = collect_worker_status(workers) )
//		print_worker_messages(workers);
//
//	if( ret == SUCCESS )
//		foreach worker in workers
//			print_worker_finish(worker);
//
//clean_modules:
//	release_modules(workers);
//clean_workers:
//	release_workers(workers);
//clean_systems:
//	release_systems(systems);
//	release_system_cl(scl);	
//	print_status(ret);
	return 0;
}

//build(name,in,out,opts)
//{
//
//}
//
//rmsd_module_build_all(rmsd,systems,mask)
//{
//	foreach system in systems:
//		switch(system.type)
//		{
//			case SYSTEM_CL:
//				
//			break;
//		}
//}
//
//build_workers(config)
//{
//	foreach unit in config
//	{
//		switch(unit.type)
//		{
//		case WORKER_RMSD:
//			worker = rmsd_worker_build(unit);
//			ret |=   rmsd_module_build(rmsd,systems);
//			break;
//		case WORKER_READ:
//			worker = read_worker_build(unit);
//			ret |=   read_module_build(rmsd,systems);
//			break;
//		}
//		if(ret)
//			break;
//	}
//	return ret;
//}
//
//build_modules(workers,systems)
//{
//	foreach worker in workers
//	{
//		switch(worker.type)
//		{
//			case WORKER_RMSD:
//			break;
//		}
//	}
//	return ret;
//}


//struct system system;
//system_cl_bld(&system,id);
//struct module module;
//
////associates devices to drivers
//module_cl_bld(&module,&system);
//
//struct worker worker;
//worker_translate_bld(&worker);
//worker->fn->module_cl_bld(&worker,&module);
//
//worker->fn->launch();
//
//worker->fn->dty(&worker);
//module->fn->dty(&module);
//system->fn->dty(&system);


	//const int num_workers = 5;
	//struct system_cl scl[2];
	//struct module_cl** mcl;
	//cl_uint num_platforms = 0;
	//cl_platform_id platform[16];
	//cl_int ret = clGetPlatformIDs(16, platform, &num_platforms);
	//printf("Notified %d platforms\n",num_platforms);
	//if(ret != CL_SUCCESS)
	//{
	//	printf("Error: OpenCl returned %d\n",ret);
	//}
	//size_t i, k;
	//for(i = 0; i < num_platforms; i++)
	//{
	//  system_cl_bld(&scl[i],platform[i]);
	//}

	//struct worker** worker = malloc(num_workers * sizeof(struct worker*));
	//mcl = malloc(num_workers * sizeof(struct module_cl*));
	//for(i = 0; i < num_workers; i++)
	//{
	//	worker[i] = malloc(sizeof(struct worker_translate));
	//	worker_translate_bld((struct worker_translate*)worker[i]);
	//	
	//	mcl[i] = malloc(num_platforms * sizeof(struct module_cl));
	//	worker[i]->module_cl = mcl[i];
	//	worker[i]->module_cl_size = num_platforms;
	//	for(k = 0; k < num_platforms; k++)
	//	{
	//		module_bld(&mcl[i][k],&scl[k]);
	//		worker[i]->fn->module_cl_bld(worker[i],&mcl[i][k],&scl[k]);
	//	}
	//}
	//for(i = 0; i < num_workers; i++)
	//{
	//	for(k = 0; k < num_platforms; k++)
	//	{
	//		worker[i]->fn->module_cl_dty(worker[i],&mcl[i][k]);
	//	}
	//	free(mcl[i]);
	//	worker_translate_dty((struct worker_translate*)worker[i]);
	//	free(worker[i]);
	//}
	//free(mcl);
	//free(worker);
	//for(i = 0; i < num_platforms; i++)
	//{
	//  system_cl_dty(&scl[i]);
	//}
// vim: ft=c


#ifndef WORKER_WRITE_H
#define WORKER_WRITE_H

#include "worker.h"
#include "system.h"

struct worker_write_properties
{
	char* filename;
	int fd;
	int fp;
	int crd;
	int dim[3];
	uint64_t offset;
	size_t dim_len[3];
	size_t beg[3];
	size_t end[3];
};

int worker_write_bld(struct worker*, struct device*, char*);
void worker_write_dty(struct worker*);

#endif /* WORKER_WRITE_H */

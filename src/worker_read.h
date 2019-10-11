
#ifndef WORKER_READ_H
#define WORKER_READ_H

#include "worker.h"
#include "system.h"

struct worker_read_properties
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

int worker_read_bld(struct worker*, struct device*, char*);
void worker_read_dty(struct worker*);

#endif /* WORKER_READ_H */

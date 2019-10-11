
#ifndef WORKER_ALIGN_H
#define WORKER_ALIGN_H

#include "worker.h"

struct worker_align_properties
{
	FILE* frmsd;
};

int worker_align_bld(struct worker*, struct device*);
void worker_align_dty(struct worker*);

#endif /* WORKER_ALIGN_H */

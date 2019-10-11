
#ifndef WORKER_RMSD_H
#define WORKER_RMSD_H

#include "worker.h"

struct worker_rmsd_properties
{
	FILE* frmsd;
};

int worker_rmsd_bld(struct worker*, struct device*);
void worker_rmsd_dty(struct worker*);

#endif /* WORKER_RMSD_H */


#ifndef WORKER_TRANSLATE_H
#define WORKER_TRANSLATE_H

#include "worker.h"

struct worker_translate_properties
{

};

int worker_translate_bld(struct worker*, struct device*);
void worker_translate_dty(struct worker*);

#endif /* WORKER_TRANSLATE_H */

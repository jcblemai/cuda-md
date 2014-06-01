#include "alloc.h"
#include <stdio.h>

void _AllocFailure(size_t size)
{
	fprintf(stderr,"Alloc failure: size=%d\n",size);
	abort();
}

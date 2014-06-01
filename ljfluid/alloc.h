/*
 * alloc.h -- Basic fail-safe allocation templates. 
 * 
 * Copyright (c) 2004 by Wolfgang Wieser, email: > wwieser -a- gmx -*- de <
 * 
 * This file may be distributed and/or modified under the terms of the 
 * GNU General Public License version 2 as published by the Free Software 
 * Foundation. 
 * 
 * This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
 * WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 * 
 */

#ifndef _ALLOC_H_
#define _ALLOC_H_ 1

#include <stdlib.h>

extern void _AllocFailure(size_t size);

// Use that only for simple (POD) types like int or double. 
// Allocate an array of nelem elements: 
template<class T>inline T *ALLOC(size_t nelem)
{
	nelem*=sizeof(T);
	T *ptr=(T*)malloc(nelem);
	if(!ptr && nelem)  _AllocFailure(nelem);
	return(ptr);
}
// Free an array as allocated by ALLOC(): 
template<class T>inline T *FREE(T *ptr)
	{  if(ptr) free(ptr);  return(NULL);  }
// Reallocation, just as uaual...
template<class T>inline T *REALLOC(T *ptr,size_t nelem)
{
	nelem*=sizeof(T);
	T *nptr=(T*)realloc(ptr,nelem);
	if(!nptr && nelem)  _AllocFailure(nelem);
	return(nptr);
}

#endif  /* _ALLOC_H_ */

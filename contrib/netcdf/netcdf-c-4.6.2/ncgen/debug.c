/*
Copyright (c) 1998-2017 University Corporation for Atmospheric Research/Unidata
See LICENSE.txt for license information.
*/

#include "includes.h"

#define TRACE

extern char* ncclassname(nc_class);

#ifdef TRACE
#define T(fcn,mem) {if(trace) {fprintf(stderr,"X: %s: %p\n", fcn,mem);}}
#else
#define T(fcn,mem)
#endif

int debug = 0;
static int trace = 0;

int
settrace(int tf)
{
    trace = tf;
    return 1;
}


void fdebug(const char *fmt, ...)
{
    va_list argv;
    if(debug == 0) return;
    va_start(argv,fmt);
    (void)vfprintf(stderr,fmt,argv) ;
}

/**************************************************/

/* Support debugging of memory*/

void
chkfree(void* memory)
{
    if(memory == NULL) {
	panic("free: null memory");
    }
T("free",memory);
    free(memory);
}

void*
chkmalloc(size_t size)
{
    void* memory = malloc(size);
    if(memory == NULL) {
	panic("malloc:out of memory");
    }
T("malloc",memory);
    return memory;
}

void*
chkcalloc(size_t size)
{
    void* memory = calloc(size,1); /* use calloc to zero memory*/
    if(memory == NULL) {
	panic("calloc:out of memory");
    }
T("calloc",memory);
    return memory;
}

void*
chkrealloc(void* ptr, size_t size)
{
    void* memory = realloc(ptr,size);
    if(memory == NULL) {
	panic("realloc:out of memory");
    }
if(ptr != memory) {T("free",memory); T("realloc",memory);}
    return memory;
}

char*
chkstrdup(const char* s)
{
    char* dup;
    if(s == NULL) {
	panic("strdup: null argument");
    }
    dup = strdup(s);
    if(dup == NULL) {
	panic("strdup: out of memory");
    }
T("strdup",dup);
    return dup;
}

int
panic(const char* fmt, ...)
{
    va_list args;
    if(fmt != NULL) {
      va_start(args, fmt);
      vfprintf(stderr, fmt, args);
      fprintf(stderr, "\n" );
      va_end( args );
    } else {
      fprintf(stderr, "panic" );
    }
    fprintf(stderr, "\n" );
    fflush(stderr);
    abort();
    return 0;
}

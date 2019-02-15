#ifndef NCGEN_DEBUG_H
#define NCGEN_DEBUG_H

/*********************************************************************
 *   Copyright 1993, UCAR/Unidata
 *   See netcdf/COPYRIGHT file for copying and redistribution conditions.
 *   $Header: /upc/share/CVS/netcdf-3/ncgen/debug.h,v 1.2 2010/03/31 18:18:34 dmh Exp $
 *********************************************************************/

#include <stdarg.h>
#include <assert.h>
#include "generr.h"
#include "bytebuffer.h"

#if 0
#define GENDEBUG 2
#endif

#ifdef GENDEBUG
#  if GENDEBUG > 0
#    define GENDEBUG1
#  endif
#  if GENDEBUG > 1
#    define GENDEBUG2
#  endif
#  if GENDEBUG > 2
#    define GENDEBUG3
#  endif
#endif

extern int settrace(int);

extern int debug;
extern int ncgdebug;

extern void fdebug(const char *fmt, ...);

#define PANIC(msg) assert(panic(msg))
#define PANIC1(msg,arg) assert(panic(msg,arg))
#define PANIC2(msg,arg1,arg2) assert(panic(msg,arg1,arg2))
#define ASSERT(expr) {if(!(expr)) {panic("assertion failure: %s",#expr);}}
extern int panic(const char* fmt, ...);

/*
Provide wrapped versions of XXXalloc for debugging/
The wrapped version:
1. fails if size is zero or memory is NULL
2. fails if memory is exhausted.
3. zeros all allocated memory.
*/

#define emalloc(x) chkmalloc(x) /*note only single arg */
#define ecalloc(x) chkcalloc(x) /*note only single arg */
#define erealloc(p,x)   chkrealloc(p,x)
#define efree(x) chkfree(x)
#define estrdup(x) chkstrdup(x)
extern void* chkmalloc(size_t);
extern void* chkcalloc(size_t);
extern void* chkrealloc(void*,size_t);
extern void  chkfree(void*);
extern char* chkstrdup(const char* s);

#define MEMCHECK(var,throw) {if((var)==NULL) return (throw);}

#endif /*NCGEN_DEBUG_H*/

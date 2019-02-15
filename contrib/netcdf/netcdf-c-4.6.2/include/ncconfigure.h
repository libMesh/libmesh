/*
 * Copyright 2010 University Corporation for Atmospheric
 * Research/Unidata. See COPYRIGHT file for more info.
 *
 * This header file is for the parallel I/O functions of netCDF.
 *
 */
/* "$Id: netcdf_par.h,v 1.1 2010/06/01 15:46:49 ed Exp $" */

#ifndef NCCONFIGURE_H
#define NCCONFIGURE_H 1

#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif

/*
This is included in bottom
of config.h. It is where,
typically, alternatives to
missing functions should be
defined and missing types defined.
*/

#ifndef HAVE_STRDUP
extern char* strdup(const char*);
#endif

/*
#ifndef HAVE_SSIZE_T
typedef long ssize_t;
#define HAVE_SSIZE_T
#endif
*/
/* handle null arguments */
#ifndef nulldup
#ifdef HAVE_STRDUP
#define nulldup(s) ((s)==NULL?NULL:strdup(s))
#else
char *nulldup(const char* s);
#endif
#endif

#ifdef _MSC_VER
#ifndef HAVE_SSIZE_T
#include <basetsd.h>
typedef SSIZE_T ssize_t;
#define HAVE_SSIZE_T 1
#endif
#endif

#ifndef HAVE_STRLCAT
#ifdef _MSC_VER
/* Windows strlcat_s is equivalent to strlcat, but different arg order */
#define strlcat(d,s,n) strcat_s((d),(n),(s))
#else
extern size_t strlcat(char* dst, const char* src, size_t dsize);
#endif
#endif

#ifndef nulldup
#define nulldup(s) ((s)==NULL?NULL:strdup(s))
#endif
#ifndef nulllen
#define nulllen(s) ((s)==NULL?0:strlen(s))
#endif
#ifndef nullfree
#define nullfree(s) {if((s)!=NULL) {free(s);} else {}}
#endif

#ifndef HAVE_UCHAR
typedef unsigned char uchar;
#endif

#ifndef HAVE_LONGLONG
typedef long long longlong;
typedef unsigned long long ulonglong;
#endif

#ifndef HAVE_USHORT
typedef unsigned short ushort;
#endif

#ifndef HAVE_UINT
typedef unsigned int uint;
#endif

#endif /* NCCONFIGURE_H */

/*********************************************************************
 *   Copyright 2009, UCAR/Unidata
 *   See netcdf/COPYRIGHT file for copying and redistribution conditions.
 *********************************************************************/
/* $Id: generr.h,v 1.2 2010/05/24 19:59:57 dmh Exp $ */
/* $Header: /upc/share/CVS/netcdf-3/ncgen/generr.h,v 1.2 2010/05/24 19:59:57 dmh Exp $ */

#ifndef GENERR_H
#define GENERR_H

#include <stdarg.h>

extern int error_count;

#if 0
#define vastart(argv,fmt) va_start(argv,fmt)
#define vaend(argv,fmt) va_end(argv)
#endif

extern void vderror(const char *fmt, va_list argv);
extern void vdwarn(const char *fmt, va_list argv);
extern void derror(const char *fmt, ...);
extern int panic(const char* fmt, ...);
extern void nprintf(char* buffer, size_t size, const char *fmt, ...);
extern  void semerror(const int, const char *fmt, ...);
extern  void semwarn(const int, const char *fmt, ...);

#endif /*GENERR_H*/

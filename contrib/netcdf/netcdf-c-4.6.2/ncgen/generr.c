/*********************************************************************
 *   Copyright 1993, UCAR/Unidata
 *   See netcdf/COPYRIGHT file for copying and redistribution conditions.
 *   $Header: /upc/share/CVS/netcdf-3/ncgen/generr.c,v 1.1 2009/09/25 18:22:22 dmh Exp $
 *********************************************************************/

#include "includes.h"
#include <ctype.h>	/* for isprint() */

int error_count;

#if 0
#define vastart(argv,fmt) va_start(argv,fmt)
#define vaend(argv,fmt) va_end(argv)
#endif

/*
 * For logging error conditions.
 * Designed to be called by other vararg procedures
 */
void
vderror(const char *fmt, va_list argv)
{
    (void) vdwarn(fmt,argv);
    error_count++;
}

/*
 * For logging error conditions.
 * Designed to be called by other vararg procedures
 */
void
vdwarn(const char *fmt, va_list argv)
{
    (void) vfprintf(stderr,fmt,argv) ;
    (void) fputc('\n',stderr) ;
    (void) fflush(stderr);	/* to ensure log files are current */
}

void
derror(const char *fmt, ...)
{
    va_list argv;
    va_start(argv,fmt);
    vderror(fmt,argv);
}

/* Report version errors */
void
verror(const char *fmt, ...)
{
    char newfmt[2048];
    va_list argv;
    va_start(argv,fmt);
    strcpy(newfmt,"netCDF classic: not supported: ");
    strncat(newfmt,fmt,2000);
    vderror(newfmt,argv);
    va_end(argv);
}

void
semwarn(const int lno, const char *fmt, ...)
{
    va_list argv;
    va_start(argv,fmt);
    (void)fprintf(stderr,"%s: %s line %d: ", progname, cdlname, lno);
    vdwarn(fmt,argv);
}

void
semerror(const int lno, const char *fmt, ...)
{
    va_list argv;
    va_start(argv,fmt);
    (void)fprintf(stderr,"%s: %s line %d: ", progname, cdlname, lno);
    vderror(fmt,argv);
    finalize_netcdf(1); /* immediately fatal */
}

/* Capture potential version errors */
static char* markcdf4_msg = NULL;
void
markcdf4(const char* msg)
{
    enhanced_flag = 1;
    if(markcdf4_msg == NULL)
        markcdf4_msg = (char*)msg;
}

char*
getmarkcdf4(void)
{
    return markcdf4_msg;
}

/* Capture potential version errors */
static char* markcdf5_msg = NULL;
void
markcdf5(const char* msg)
{
    cdf5_flag = 1;
    if(markcdf5_msg == NULL)
        markcdf5_msg = (char*)msg;
}

char*
getmarkcdf5(void)
{
    return markcdf5_msg;
}

/**************************************************/
/* Provide a version of snprintf that panics if*/
/* the buffer is overrun*/

void
nprintf(char* buffer, size_t size, const char *fmt, ...)
{
    long written;
    va_list args;
    va_start(args,fmt);
    written = vsnprintf(buffer,size,fmt,args);
    if(written < 0 || written >= size) {
	PANIC("snprintf failure");
    }
}

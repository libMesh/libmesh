/*********************************************************************
 *   Copyright 1993, UCAR/Unidata
 *   See netcdf/COPYRIGHT file for copying and redistribution conditions.
 *********************************************************************/
#include "config.h"
#include <stdarg.h>
#include <stdio.h>

#include "nclog.h"
#include "ncdap.h"

int ncdap3debug = 0;

#ifdef CATCHERROR
/* Place breakpoint here to catch errors close to where they occur*/
int
dapbreakpoint(int err) {return err;}

int
dapthrow(int err)
{
    if(err == 0) return err;
    return dapbreakpoint(err);
}
#endif

int
dappanic(const char* fmt, ...)
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
    return 0;
}

/*
Provide a way to print the oc full name for
an ocnode
*/

char*
ocfqn(OCddsnode node)
{
    OClink conn;
    oc_get_connection(node,&conn);
    return makeocpathstring(conn,node,".");
}

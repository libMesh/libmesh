/*********************************************************************
 *   Copyright 2010, UCAR/Unidata
 *   See netcdf/COPYRIGHT file for copying and redistribution conditions.
 *********************************************************************/

#include "config.h"
#include "ncdispatch.h"

extern int NC3_initialize(void);

extern int NCSUBSTRATE_initialize(void);

#ifdef USE_NETCDF4
extern int NC4_initialize(void);
#endif

#ifdef USE_DAP
extern int NCD2_initialize(void);
#endif

#ifdef USE_CDMREMOTE
extern int NCCR_initialize(void);
#endif

#ifdef USE_PNETCDF
extern int NC5_initialize(void);
#endif

/**
This procedure invokes all defined
initializers, and there is an initializer
for every known dispatch table.
So if you modify the format of NC_Dispatch,
then you need to fix it everywhere.
*/

int
NC_initialize(void)
{
    int stat = NC_NOERR;

    /* Allow libdispatch to do initialization */
    if((stat = NCDISPATCH_initialize())) return stat;

    /* Initialize each active protocol */

    if((stat = NC3_initialize())) return stat;

#ifdef USE_DAP
    if((stat = NCD2_initialize())) return stat;
#endif

#ifdef USE_PNETCDF
    if((stat = NC5_initialize())) return stat;
#endif

#ifdef USE_NETCDF4
    if((stat = NC4_initialize())) return stat;

    /* if((stat = NCD_initialize())) return stat; */

#ifdef USE_DAP
#ifdef NOTUSED
    if((stat = NCD4_initialize())) return stat;
#endif
#endif

#ifdef USE_CDMREMOTE
    if((stat = NCCR_initialize())) return stat;
#endif

#endif /* USE_NETCDF4 */

    /* Finally, initialize the SUBSTRATE table (dsubstrate.c) */
    if((stat = NCSUBSTRATE_initialize())) return stat;

    return NC_NOERR;
}



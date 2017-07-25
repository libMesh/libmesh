/*********************************************************************
  *   Copyright 1993, UCAR/Unidata
  *   See netcdf/COPYRIGHT file for copying and redistribution conditions.
  *********************************************************************/
#ifndef NCDAP_H
#define NCDAP_H 1

#include "config.h"
#include <stdlib.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#include <fcntl.h>
#include <sys/stat.h>
#include <string.h>
#include <stdio.h>

#include "ncbytes.h"
#include "nclist.h"
#include "nchashmap.h"
#include "nclog.h"
#include "ncuri.h"

#include "fbits.h"
#include "dceconstraints.h"

#include "netcdf.h"
#include "ncdispatch.h"
#include "nc.h"
#include "nc3internal.h"
 /* netcdf overrides*/
#include "dapnc.h"

#include "oc.h"

#include "dapdebug.h"
#include "daputil.h"

/**************************************************/
/* sigh, do the forwards */
struct NCprojection;
struct NCselection;
struct Getvara;
struct NCcachenode;
struct NCcache;
struct NCslice;
struct NCsegment;

/**************************************************/

#include "nccommon.h"
#include "getvara.h"
#include "constraints.h"

/**************************************************/

#ifndef USE_NETCDF4
#define	NC_UBYTE 	7	/* unsigned 1 byte int */
#define	NC_USHORT 	8	/* unsigned 2-byte int */
#define	NC_UINT 	9	/* unsigned 4-byte int */
#define	NC_INT64 	10	/* signed 8-byte int */
#define	NC_UINT64 	11	/* unsigned 8-byte int */
#define	NC_STRING 	12	/* string */

#define NC_MAX_BYTE 127
#define NC_MIN_BYTE (-NC_MAX_BYTE-1)
#define NC_MAX_CHAR 255
#define NC_MAX_SHORT 32767
#define NC_MIN_SHORT (-NC_MAX_SHORT - 1)
#define NC_MAX_INT 2147483647
#define NC_MIN_INT (-NC_MAX_INT - 1)
#define NC_MAX_FLOAT 3.402823466e+38f
#define NC_MIN_FLOAT (-NC_MAX_FLOAT)
#define NC_MAX_DOUBLE 1.7976931348623157e+308 
#define NC_MIN_DOUBLE (-NC_MAX_DOUBLE)
#define NC_MAX_UBYTE NC_MAX_CHAR
#define NC_MAX_USHORT 65535U
#define NC_MAX_UINT 4294967295U
#define NC_MAX_INT64 (9223372036854775807LL)
#define NC_MIN_INT64 (-9223372036854775807LL-1)
#define NC_MAX_UINT64 (18446744073709551615ULL)
#define X_INT64_MAX     (9223372036854775807LL)
#define X_INT64_MIN     (-X_INT64_MAX - 1)
#define X_UINT64_MAX    (18446744073709551615ULL)
#endif /*USE_NETCDF4*/


/**************************************************/

extern struct NCTMODEL nctmodels[];

/**************************************************/
/* Import some internal procedures from libsrc*/

/* Internal, but non-static procedures */
extern NCerror computecdfvarnames(NCDAPCOMMON*,CDFnode*,NClist*);
extern NCerror computecdfnodesets(NCDAPCOMMON* nccomm, CDFtree* tree);
extern NCerror computevarnodes(NCDAPCOMMON*, NClist*, NClist*);
extern NCerror collectvardefdims(NCDAPCOMMON* drno, CDFnode* var, NClist* dimset);
extern NCerror fixgrids(NCDAPCOMMON* drno);
extern NCerror dapmerge(NCDAPCOMMON* drno, CDFnode* node, OCobject dasroot);
extern NCerror sequencecheck(NCDAPCOMMON* drno);
extern NCerror computecdfdimnames(NCDAPCOMMON*);
extern NCerror attachdatadds(struct NCDAPCOMMON*);
extern NCerror detachdatadds(struct NCDAPCOMMON*);
extern void dapdispatch3init(void);

/*
extern void dereference(NCconstraint* constraint);
extern NCerror rereference(NCconstraint*, NClist*);
*/

extern NCerror dapbuildvaraprojection(CDFnode*,
		     const size_t* startp, const size_t* countp, const ptrdiff_t* stridep,
		     struct DCEprojection** projectionlist);

extern NCerror nc3d_getvarx(int ncid, int varid,
	    const size_t *startp,
	    const size_t *countp,
	    const ptrdiff_t *stridep,
	    void *data,
	    nc_type dsttype0);

/**************************************************/

/* From: ncd2dispatch.c*/
extern size_t dap_one[NC_MAX_VAR_DIMS];
extern size_t dap_zero[NC_MAX_VAR_DIMS];

extern NCerror nc3d_open(const char* path, int mode, int* ncidp);
extern int nc3d_close(int ncid);
extern NCerror restruct(NCDAPCOMMON*, CDFnode* ddsroot, CDFnode* pattern, NClist*);
extern void setvisible(CDFnode* root, int visible);
extern NCerror mapnodes(CDFnode* dstroot, CDFnode* srcroot);
extern void unmap(CDFnode* root);

#if 0
extern NCerror fetchpatternmetadata(NCDAPCOMMON* nccomm);
extern NCerror fetchconstrainedmetadata(NCDAPCOMMON* nccomm);
extern void applyclientparamcontrols(NCDAPCOMMON*);
extern NCerror suppressunusablevars(NCDAPCOMMON*);
extern void estimatevarsizes(NCDAPCOMMON*);
extern NCerror showprojection(NCDAPCOMMON*, CDFnode* var);
#endif

/* From: dapcvt.c*/
extern NCerror dapconvert(nc_type, nc_type, char*, char*, size_t);
extern int dapcvtattrval(nc_type, void*, NClist*);

#endif /*NCDAP_H*/

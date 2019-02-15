/*********************************************************************
 *   Copyright 2016, UCAR/Unidata
 *   See netcdf/COPYRIGHT file for copying and redistribution conditions.
 *********************************************************************/
#include "config.h"
#include <stdarg.h>
#include <stdio.h>

#include "d4includes.h"
#include "ncdispatch.h"
#include "netcdf_aux.h"

int ncdap4debug = 0;

#ifdef D4CATCH
/* Place breakpoint here to catch errors close to where they occur*/
int
d4breakpoint(int err) {return err;}

int
d4throw(int err)
{
    if(err == 0) return err;
    return d4breakpoint(err);
}
#endif

int
d4panic(const char* fmt, ...)
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

const char*
NCD4_sortname(NCD4sort sort)
{
    switch (sort) {
    case NCD4_NULL: return "NCD4_NULL";
    case NCD4_ATTR: return "NCD4_ATTR";
    case NCD4_ATTRSET: return "NCD4_ATTRSET";
    case NCD4_XML: return "NCD4_XML";
    case NCD4_DIM: return "NCD4_DIM";
    case NCD4_GROUP: return "NCD4_GROUP";
    case NCD4_TYPE: return "NCD4_TYPE";
    case NCD4_VAR: return "NCD4_VAR";
    case NCD4_ECONST: return "NCD4_ECONST";
    default: break;
    }
    return "unknown";
}

const char*
NCD4_subsortname(nc_type subsort)
{
    switch (subsort) {
    case NC_NAT: return "NC_NAT";
    case NC_BYTE: return "NC_BYTE";
    case NC_CHAR: return "NC_CHAR";
    case NC_SHORT: return "NC_SHORT";
    case NC_INT: return "NC_INT";
    case NC_FLOAT: return "NC_FLOAT";
    case NC_DOUBLE: return "NC_DOUBLE";
    case NC_UBYTE: return "NC_UBYTE";
    case NC_USHORT: return "NC_USHORT";
    case NC_UINT: return "NC_UINT";
    case NC_INT64: return "NC_INT64";
    case NC_UINT64: return "NC_UINT64";
    case NC_STRING: return "NC_STRING";
    case NC_VLEN: return "NC_VLEN";
    case NC_OPAQUE: return "NC_OPAQUE";
    case NC_ENUM: return "NC_ENUM";
    case NC_COMPOUND: return "NC_COMPOUND";
    default: break;
    }
    return "unknown";
}

/*
For debugging purposes, it is desirable to fake an nccopy
bv inserting the data into the substrate and then writing it out.
*/

int
NCD4_debugcopy(NCD4INFO* info)
{
    int i,ret=NC_NOERR;
    NCD4meta* meta = info->substrate.metadata;
    NClist* topvars = nclistnew();
    NC* ncp = info->controller;

    /* Walk each top level variable, read all of it and write it to the substrate */
    if((ret=NCD4_getToplevelVars(meta, NULL, topvars)))
	goto done;
    /* Read from the dap data by going thru the dap4 interface */
    for(i=0;i<nclistlength(topvars);i++) {
	NCD4node* var = nclistget(topvars,i);
	NCD4node* type = var->basetype;
	NCD4node* grp = NCD4_groupFor(var);
	int grpid = grp->meta.id;
	int varid = var->meta.id;
	d4size_t varsize;
	void* memory = NULL;
	size_t dimprod = NCD4_dimproduct(var);
	int ncid = info->substrate.nc4id;

	varsize = type->meta.memsize * dimprod;
	memory = d4alloc(varsize);
        if(memory == NULL)
	    {ret = NC_ENOMEM; goto done;}		
	{
	    /* We need to read via NCD4 */
	    int d4gid = makedap4id(ncp,grpid);
            if((ret=nc_get_var(d4gid,varid,memory)))
	        goto done;
 	}
	/* Now, turn around and write it to the substrate.
	  WARNING: we have to specify the shape ourselves
          because, if unlimited is involved then there is
          potentially a difference between the substrate unlimited
          size and the dap4 data specified size. In fact,
          the substrate will always be zero unless debugcopy is used.
	*/
	{	
	    size_t edges[NC_MAX_VAR_DIMS];
	    int d;
	    for(d=0;d<nclistlength(var->dims);d++) {
		NCD4node* dim = (NCD4node*)nclistget(var->dims,d);
		edges[d] = (size_t)dim->dim.size;
	    }
            if((ret=nc_put_vara(grpid,varid,nc_sizevector0,edges,memory)))
	        goto done;
	}
	if((ret=ncaux_reclaim_data(ncid,type->meta.id,memory,dimprod)))
	    goto done;
	free(memory);
	memory = NULL;
    }	    
done:
    if(topvars)
        nclistfree(topvars);
    if(ret != NC_NOERR) {
        fprintf(stderr,"debugcopy: %d %s\n",ret,nc_strerror(ret));
    }
    return THROW(ret);
}

/* Provide a string printer that can be called from gdb */
void
NCD4_printstring(const char* s)
{
    fprintf(stderr,"%s\n",s);
    fflush(stderr);
}

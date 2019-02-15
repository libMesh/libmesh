/*********************************************************************
 *   Copyright 2016, UCAR/Unidata
 *   See netcdf/COPYRIGHT file for copying and redistribution conditions.
 *********************************************************************/

#include "config.h"
#include "ncdispatch.h"
#include "ncd4dispatch.h"
#include "nc4internal.h"
#include "d4includes.h"
#include "d4odom.h"

/* Forward */
static int getvarx(int ncid, int varid, NCD4INFO**, NCD4node** varp, nc_type* xtypep, size_t*, nc_type* nc4typep, size_t*);

int
NCD4_get_vara(int ncid, int varid,
	    const size_t *start, const size_t *edges,
            void *value,
	    nc_type memtype)
{
    int ret;
    /* TODO: optimize since we know stride is 1 */
    ret = NCD4_get_vars(ncid,varid,start,edges,nc_ptrdiffvector1,value,memtype);
    return ret;
}

int
NCD4_get_vars(int ncid, int varid,
	    const size_t *start, const size_t *edges, const ptrdiff_t* stride,
            void *memoryin, nc_type xtype)
{
    int i,ret;
    NCD4INFO* info;
    NCD4meta* meta;
    NCD4node* ncvar;
    NCD4node* nctype;
    D4odometer* odom = NULL;
    nc_type nc4type;
    size_t nc4size, xsize;
    void* instance = NULL; /* Staging area in case we have to convert */
    NClist* blobs = NULL;
    int rank;
    size_t dimsizes[NC_MAX_VAR_DIMS];
    d4size_t dimproduct;
    size_t dstcount;
    
    if((ret=getvarx(ncid, varid, &info, &ncvar, &xtype, &xsize, &nc4type, &nc4size)))
	{THROW(ret); goto done;}

    meta = info->substrate.metadata;
    nctype = ncvar->basetype;
    rank = nclistlength(ncvar->dims);
    blobs = nclistnew();

    instance = malloc(nc4size);
    if(instance == NULL)
	{ret = THROW(NC_ENOMEM); goto done;}	

    dimproduct = NCD4_dimproduct(ncvar);
    /* build size vector */
    for(i=0;i<rank;i++)  {
	NCD4node* dim = nclistget(ncvar->dims,i);
	dimsizes[i] = (size_t)dim->dim.size;
    }
	
    /* Extract and desired subset of data */
    if(rank > 0)
        odom = d4odom_new(rank,start,edges,stride,dimsizes);
    else
        odom = d4scalarodom_new();
    dstcount = 0; /* We always write into dst starting at position 0*/
    for(;d4odom_more(odom);dstcount++) {
	void* xpos;
	void* offset;
	void* dst;
	d4size_t count;
        count = d4odom_next(odom);
	if(count >= dimproduct) {
	    ret = THROW(NC_EINVALCOORDS);
	    goto done;
	}
        xpos = INCR(memoryin,(xsize * dstcount)); /* ultimate destination */
	/* We need to compute the offset in the dap4 data of this instance;
	   for fixed size types, this is easy, otherwise we have to walk
	   the variable size type
	*/
	if(nctype->meta.isfixedsize) {
	    offset = INCR(ncvar->data.dap4data.memory,(nc4size * count));
	} else {
            offset = ncvar->data.dap4data.memory;
	    /* We have to walk to the count'th location in the data */
	    if((ret=NCD4_moveto(meta,ncvar,count,&offset)))
	        {THROW(ret); goto done;}		    
	}
	dst = instance;
	if((ret=NCD4_fillinstance(meta,nctype,&offset,&dst,blobs)))
	    {THROW(ret); goto done;}
	if(xtype == nc4type) {
	    /* We can just copy out the data */
	    memcpy(xpos,instance,nc4size);
	} else { /* Need to convert */
	    if((ret=NCD4_convert(nc4type,xtype,xpos,instance,1)))
	        {THROW(ret); goto done;}
	}
    }

done:
    /* cleanup */
    if(odom != NULL)
	d4odom_free(odom);
    if(instance != NULL)
	free(instance);
    if(ret != NC_NOERR) { /* reclaim all malloc'd data if there is an error*/
	for(i=0;i<nclistlength(blobs);i++) {
	    nullfree(nclistget(blobs,i));
	}
    }
    if(blobs) nclistfree(blobs);
    return (ret);
}

static int
getvarx(int ncid, int varid, NCD4INFO** infop, NCD4node** varp,
	nc_type* xtypep, size_t* xsizep, nc_type* nc4typep, size_t* nc4sizep)
{
    int ret = NC_NOERR;
    NC* ncp;
    NCD4INFO* info;
    NCD4meta* meta;
    NCD4node* group;
    NCD4node* var;
    NCD4node* type;
    nc_type xtype, actualtype;
    size_t instancesize, xsize;
    int grp_id;

    if((ret = NC_check_id(ncid, (NC**)&ncp)) != NC_NOERR)
	{THROW(ret); goto done;}

    info = getdap(ncp);
    if(info == NULL)
	{ret = THROW(NC_EBADID); goto done;}
    meta = info->substrate.metadata;
    if(meta == NULL)
	{ret = THROW(NC_EBADID); goto done;}

    /* Locate var node via (grpid,varid) */
    grp_id = GROUPIDPART(ncid);
        
    group = nclistget(meta->groupbyid,grp_id);
    if(group == NULL)
	return THROW(NC_EBADID);
    var = nclistget(group->group.varbyid,varid);
    if(var == NULL)
	return THROW(NC_EBADID);
    type = var->basetype;
    actualtype = type->meta.id;
    instancesize = type->meta.memsize;

    /* Figure out the type conversion, if any */
    xtype = *xtypep;
    if(xtype == NC_NAT)
	xtype = actualtype;
    if(xtype != actualtype && xtype > NC_MAX_ATOMIC_TYPE)
	return THROW(NC_EBADTYPE);
    if((xtype == NC_CHAR || xtype == NC_STRING)
	&& (actualtype != NC_CHAR && actualtype != NC_STRING))
	return THROW(NC_ECHAR);
    if(xtype <= NC_MAX_ATOMIC_TYPE)
        xsize = NCD4_typesize(xtype);
    else
	xsize = instancesize;

    /* Return relevant info */
    if(infop) *infop = info;
    if(xtypep) *xtypep = xtype;
    if(xsizep) *xsizep = xsize;
    if(nc4typep) *nc4typep = actualtype;
    if(nc4sizep) *nc4sizep = instancesize;
    if(varp) *varp = var;
done:
    return THROW(ret);    
}


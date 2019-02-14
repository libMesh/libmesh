/*********************************************************************
 *   Copyright 2016, UCAR/Unidata
 *   See netcdf/COPYRIGHT file for copying and redistribution conditions.
 *********************************************************************/

#include "d4includes.h"
#include <stdarg.h>
#include <assert.h>
#include "ezxml.h"
#include "d4includes.h"
#include "d4odom.h"

/**
This code serves two purposes
1. Preprocess the dap4 serialization wrt endianness, etc.
   (NCD4_processdata)
2. Walk a specified variable instance to convert to netcdf4
   memory representation.
   (NCD4_fillinstance)

*/

/***************************************************/
/* Forwards */
static int fillstring(NCD4meta*, void** offsetp, void** dstp, NClist* blobs);
static int fillopfixed(NCD4meta*, d4size_t opaquesize, void** offsetp, void** dstp);
static int fillopvar(NCD4meta*, NCD4node* type, void** offsetp, void** dstp, NClist* blobs);
static int fillstruct(NCD4meta*, NCD4node* type, void** offsetp, void** dstp, NClist* blobs);
static int fillseq(NCD4meta*, NCD4node* type, void** offsetp, void** dstp, NClist* blobs);

/***************************************************/
/* Macro define procedures */

#ifdef D4DUMPCSUM
static unsigned int debugcrc32(unsigned int crc, const void *buf, size_t size)
{
    int i;
    fprintf(stderr,"crc32: ");
    for(i=0;i<size;i++) {fprintf(stderr,"%02x",((unsigned char*)buf)[i]);}
    fprintf(stderr,"\n");
    return NCD4_crc32(crc,buf,size);
}
#define CRC32 debugcrc32
#else
#define CRC32 NCD4_crc32
#endif

#define ISTOPLEVEL(var) ((var)->container == NULL || (var)->container->sort == NCD4_GROUP)

/***************************************************/
/* API */

int
NCD4_processdata(NCD4meta* meta)
{
    int ret = NC_NOERR;
    int i;
    NClist* toplevel = NULL;
    NCD4node* root = meta->root;
    void* offset;

    /* Recursively walk the tree in prefix order 
       to get the top-level variables; also mark as unvisited */
    toplevel = nclistnew();
    NCD4_getToplevelVars(meta,root,toplevel);

    /* If necessary, byte swap the serialized data */
    /* Do we need to swap the dap4 data? */
    meta->swap = (meta->serial.hostlittleendian != meta->serial.remotelittleendian);

    /* Compute the  offset and size of the toplevel vars in the raw dap data. */
    offset = meta->serial.dap;
    for(i=0;i<nclistlength(toplevel);i++) {
	NCD4node* var = (NCD4node*)nclistget(toplevel,i);
        if((ret=NCD4_delimit(meta,var,&offset)))
	    FAIL(ret,"delimit failure");
    }

    /* Swap the data for each top level variable,
	including the checksum (if any)
    */
    if(meta->swap) {
        if((ret=NCD4_swapdata(meta,toplevel)))
	    FAIL(ret,"byte swapping failed");
    }

    /* Compute the checksums of the top variables */
    if(meta->localchecksumming) {
	for(i=0;i<nclistlength(toplevel);i++) {
	    unsigned int csum = 0;
	    NCD4node* var = (NCD4node*)nclistget(toplevel,i);
            csum = CRC32(csum,var->data.dap4data.memory,var->data.dap4data.size);
            var->data.localchecksum = csum;
	}
    }

    /* verify checksums */
    if(!meta->ignorechecksums && meta->serial.remotechecksumming) {
        for(i=0;i<nclistlength(toplevel);i++) {
	    NCD4node* var = (NCD4node*)nclistget(toplevel,i);
	    if(var->data.localchecksum != var->data.remotechecksum) {
		nclog(NCLOGERR,"Checksum mismatch: %s\n",var->name);
		ret = NC_EDAP;
		goto done;
	    }
        }
    }
done:
    if(toplevel) nclistfree(toplevel);
    return THROW(ret);
}

/*
Build a single instance of a type. The blobs
argument accumulates any malloc'd data so we can
reclaim it in case of an error.

Activity is to walk the variable's data to
produce a copy that is compatible with the
netcdf4 memory format.

Assumes that NCD4_processdata has been called.
*/

int
NCD4_fillinstance(NCD4meta* meta, NCD4node* type, void** offsetp, void** dstp, NClist* blobs)
{
    int ret = NC_NOERR;
    void* offset = *offsetp;
    void* dst = *dstp;
    d4size_t memsize = type->meta.memsize;
    d4size_t dapsize = type->meta.dapsize;

    /* If the type is fixed size, then just copy it  */
    if(type->subsort <= NC_UINT64 || type->subsort == NC_ENUM) {
	/* memsize and dapsize are the same */
	assert(memsize == dapsize);
	memcpy(dst,offset,dapsize);
	offset = INCR(offset,dapsize);
    } else switch(type->subsort) {
        case NC_STRING: /* oob strings */
	    if((ret=fillstring(meta,&offset,&dst,blobs)))
	        FAIL(ret,"fillinstance");
	    break;
	case NC_OPAQUE:
	    if(type->opaque.size > 0) {
	        /* We know the size and its the same for all instances */
	        if((ret=fillopfixed(meta,type->opaque.size,&offset,&dst)))
	            FAIL(ret,"fillinstance");
   	    } else {
	        /* Size differs per instance, so we need to convert each opaque to a vlen */
	        if((ret=fillopvar(meta,type,&offset,&dst,blobs)))
	            FAIL(ret,"fillinstance");
	    }
	    break;
	case NC_STRUCT:
	    if((ret=fillstruct(meta,type,&offset,&dst,blobs)))
                FAIL(ret,"fillinstance");
	    break;
	case NC_SEQ:
	    if((ret=fillseq(meta,type,&offset,&dst,blobs)))
                FAIL(ret,"fillinstance");
	    break;
	default:
	    ret = NC_EINVAL;
            FAIL(ret,"fillinstance");
    }
    *dstp = dst;
    *offsetp = offset; /* return just past this object in dap data */
done:
    return THROW(ret);
}

static int
fillstruct(NCD4meta* meta, NCD4node* type, void** offsetp, void** dstp, NClist* blobs)
{
    int i,ret = NC_NOERR;
    void* offset = *offsetp;
    void* dst = *dstp;
 
#ifdef CLEARSTRUCT
    /* Avoid random data within aligned structs */
    memset(dst,0,type->meta.memsize);
#endif

    /* Walk and read each field taking alignments into account */
    for(i=0;i<nclistlength(type->vars);i++) {
	NCD4node* field = nclistget(type->vars,i);
	NCD4node* ftype = field->basetype;
	void* fdst = INCR(dst,field->meta.offset);
	if((ret=NCD4_fillinstance(meta,ftype,&offset,&fdst,blobs)))
            FAIL(ret,"fillstruct");
    }
    dst = INCR(dst,type->meta.memsize);
    *dstp = dst;
    *offsetp = offset;    
done:
    return THROW(ret);
}

static int
fillseq(NCD4meta* meta, NCD4node* type, void** offsetp, void** dstp, NClist* blobs)
{
    int ret = NC_NOERR;
    d4size_t i,recordcount;    
    void* offset;
    nc_vlen_t* dst;
    NCD4node* vlentype;
    d4size_t recordsize;

    offset = *offsetp;

    dst = (nc_vlen_t*)*dstp;
    vlentype = type->basetype;    
    recordsize = vlentype->meta.memsize;

    /* Get record count (remember, it is already properly swapped) */
    recordcount = GETCOUNTER(offset);
    SKIPCOUNTER(offset);
    dst->len = (size_t)recordcount;

    /* compute the required memory */
    dst->p = d4alloc(recordsize*recordcount);
    if(dst->p == NULL) 
	FAIL(NC_ENOMEM,"fillseq");

    for(i=0;i<recordcount;i++) {
	/* Read each record instance */
	void* recdst = INCR((dst->p),(recordsize * i));
	if((ret=NCD4_fillinstance(meta,vlentype,&offset,&recdst,blobs)))
	    FAIL(ret,"fillseq");
    }
    dst++;
    *dstp = dst;
    *offsetp = offset;
done:
    return THROW(ret);
}

/*
Extract and oob a single string instance
*/
static int
fillstring(NCD4meta* meta, void** offsetp, void** dstp, NClist* blobs)
{
    int ret = NC_NOERR;
    d4size_t count;    
    void* offset = *offsetp;
    char** dst = *dstp;
    char* q;

    /* Get string count (remember, it is already properly swapped) */
    count = GETCOUNTER(offset);
    SKIPCOUNTER(offset);
    /* Transfer out of band */
    q = (char*)d4alloc(count+1);
    if(q == NULL)
	{FAIL(NC_ENOMEM,"out of space");}
    memcpy(q,offset,count);
    q[count] = '\0';
    /* Write the pointer to the string */
    *dst = q;
    dst++;
    *dstp = dst;    
    offset = INCR(offset,count);
    *offsetp = offset;
#if 0
    nclistpush(blobs,q);
#else
    q = NULL;
#endif
done:
    return THROW(ret);
}

static int
fillopfixed(NCD4meta* meta, d4size_t opaquesize, void** offsetp, void** dstp)
{
    int ret = NC_NOERR;
    d4size_t count, actual;
    int delta;
    void* offset = *offsetp;
    void* dst = *dstp;

    /* Get opaque count */
    count = GETCOUNTER(offset);
    SKIPCOUNTER(offset);
    /* verify that it is the correct size */
    actual = count;
    delta = actual - opaquesize;
    if(delta != 0) {
#ifdef FIXEDOPAQUE
	nclog(NCLOGWARN,"opaque changed from %lu to %lu",actual,opaquesize);
	memset(dst,0,opaquesize); /* clear in case we have short case */
	count = (delta < 0 ? actual : opaquesize);
#else
        FAIL(NC_EVARSIZE,"Expected opaque size to be %lld; found %lld",opaquesize,count);
#endif
    }
    /* move */
    memcpy(dst,offset,count);
    dst = INCR(dst,count);
    *dstp = dst;
    offset = INCR(offset,count);
    *offsetp = offset;
#ifndef FIXEDOPAQUE
done:
#endif
    return THROW(ret);
}

/*
Move a dap4 variable length opaque out of band.
We treat as if it was (in cdl) ubyte(*).
*/

static int
fillopvar(NCD4meta* meta, NCD4node* type, void** offsetp, void** dstp, NClist* blobs)
{
    int ret = NC_NOERR;
    d4size_t count;
    nc_vlen_t* vlen;
    void* offset = *offsetp;
    void* dst = *dstp;
    char* q;

    /* alias dst format */
    vlen = (nc_vlen_t*)dst;

    /* Get opaque count */
    count = GETCOUNTER(offset);
    SKIPCOUNTER(offset);
    /* Transfer out of band */
    q = (char*)d4alloc(count);
    if(q == NULL) FAIL(NC_ENOMEM,"out of space");
    memcpy(q,offset,count);
    vlen->p = q;
    vlen->len = (size_t)count;
    q = NULL; /*nclistpush(blobs,q);*/
    dst = INCR(dst,sizeof(nc_vlen_t));
    *dstp = dst;
    offset = INCR(offset,count);
    *offsetp = offset;
done:
    return THROW(ret);
}


/**************************************************/
/* Utilities */
int
NCD4_getToplevelVars(NCD4meta* meta, NCD4node* group, NClist* toplevel)
{
    int ret = NC_NOERR;
    int i;

    if(group == NULL)
	group = meta->root;

    /* Collect vars in this group */
    for(i=0;i<nclistlength(group->vars);i++) {
        NCD4node* node = (NCD4node*)nclistget(group->vars,i);
        nclistpush(toplevel,node);
        node->visited = 0; /* We will set later to indicate written vars */
#ifdef D4DEBUGDATA
fprintf(stderr,"toplevel: var=%s\n",node->name);
#endif
    }
    /* Now, recurse into subgroups; will produce prefix order */
    for(i=0;i<nclistlength(group->groups);i++) {
        NCD4node* g = (NCD4node*)nclistget(group->groups,i);
	if((ret=NCD4_getToplevelVars(meta,g,toplevel))) goto done;
    }
done:
    return THROW(ret);
}

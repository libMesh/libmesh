/*********************************************************************
 *   Copyright 2016, UCAR/Unidata
 *   See netcdf/COPYRIGHT file for copying and redistribution conditions.
 *********************************************************************/

#include <stdarg.h>
#include <assert.h>

#include "d4includes.h"
#include "ezxml.h"

/*
The primary purpose of this code is to provide node and data walkers
to do a variety of things.

* (topo)
  a. topologically sort the set of allnodes in postfix order.

* (delimit)
  a. Delimit the top-level vars in the dap4data and also compute checksum

* (move)
  a. Walk the (toplevel) var's data to get to the count'th instance.

*/

/**************************************************/

#define COUNTSIZE 8

/**************************************************/
/* Forward */

static int delimitAtomicVar(NCD4meta*, NCD4node* var, void** offsetp);
static int delimitOpaqueVar(NCD4meta*,  NCD4node* var, void** offsetp);
static int delimitSeq(NCD4meta*, NCD4node* var, void** offsetp);
static int delimitSeqArray(NCD4meta*, NCD4node* var,  void** offsetp);
static int delimitStruct(NCD4meta*, NCD4node* var, void** offsetp);
static int delimitStructArray(NCD4meta*, NCD4node* var,  void** offsetp);
static int skipInstance(NCD4meta*, NCD4node* type, void** offsetp);
static int skipAtomicInstance(NCD4meta*, NCD4node* type, void** offsetp);
static int skipStructInstance(NCD4meta*, NCD4node* type,  void** offsetp);
static int skipSeqInstance(NCD4meta*, NCD4node* vlentype, void** offsetp);

static void walk(NCD4node* node, NClist* sorted);

#ifdef D4DEBUGDATA
static void
ADDNODE(NClist* list, NCD4node* node)
{
    fprintf(stderr,"addnode: %s: %s %s (%llu)\n",
	node->name,
	NCD4_sortname(node->sort),
	NCD4_subsortname(node->subsort),
	node);
    fflush(stderr);
    nclistpush(list,node);
}
#else
#define ADDNODE(list,node) nclistpush((list),(node))
#endif


	
/**************************************************/
/* Topo sort in postfix order */

int
NCD4_toposort(NCD4meta* compiler)
{
    int i, ret=NC_NOERR;
    size_t len = nclistlength(compiler->allnodes);
    NCD4node** list = (NCD4node**)nclistcontents(compiler->allnodes);
    NCD4node** p;
    NClist* sorted = nclistnew();
    nclistsetalloc(sorted,len);
    for(i=0,p=list;i<len;p++,i++) {
	NCD4node* node = *p;
	switch (node->sort) { /* Collect things known to have no dependencies */
	case NCD4_DIM:
	    node->visited = 1;
	    ADDNODE(sorted,node);
	    break;
	case NCD4_TYPE:
	    if(node->subsort <= NC_MAX_ATOMIC_TYPE || node->subsort == NC_OPAQUE) {
		node->visited = 1;
		ADDNODE(sorted,node);
	    } else
		node->visited = 0;
	    break;
	default:
	    node->visited = 0;
	}
    }
    walk(compiler->root,sorted);
    /* Last step is to add in any remaining unvisited nodes, but report them */
    for(i=0,p=list;i<len;p++,i++) {
	NCD4node* node = *p;
	if(node->visited) continue;
#ifdef D4DEBUGDATA
fprintf(stderr,"unvisited node: %s\n",node->name); fflush(stderr);
#endif
	node->visited = 1;
	ADDNODE(sorted,node);
    }
    nclistfree(compiler->allnodes);
    compiler->allnodes = sorted;
#ifdef D4DEBUGDATA
{int i;
for(i=0;i<nclistlength(sorted);i++)
fprintf(stderr,"sorted: %s\n",((NCD4node*)nclistget(sorted,i))->name);
fflush(stderr);
}
#endif
    return THROW(ret);
}

/*
Do depth first search
*/
static void
walk(NCD4node* node, NClist* sorted)
{
    int i;

    if(node->visited) return;
    node->visited = 1;

    switch (node->sort) {
    case NCD4_GROUP: /* depends on its elements and attributes and subgroups */
	for(i=0;i<nclistlength(node->group.elements);i++) {
	    NCD4node* elem = (NCD4node*)nclistget(node->group.elements,i);
	    walk(elem,sorted);
	}		
	break;	
    case NCD4_TYPE: /* Need to discriminate on the subsort */
	switch (node->subsort) {
	case NC_SEQ:
	    /* Depends on its basetype */
	    walk(node->basetype,sorted);
	    break;
	case NC_STRUCT: /* Depends on its fields */
            for(i=0;i<nclistlength(node->vars);i++) {
	        NCD4node* f = (NCD4node*)nclistget(node->vars,i);
		walk(f,sorted);
	    }
	    break;	
	case NC_ENUM: /* Depends on its basetype, but since that is atomic, we can ignore */
	    /* fall thru */
	default: /* Atomic or opaque, so already marked */
	   break;
	}
	break;

    case NCD4_VAR: /* Depends on: dimensions and basetype and maps */
        for(i=0;i<nclistlength(node->dims);i++) {
	    NCD4node* d = (NCD4node*)nclistget(node->dims,i);
  	    walk(d,sorted);
	}
	walk(node->basetype,sorted);
        for(i=0;i<nclistlength(node->maps);i++) {
	    NCD4node* m = (NCD4node*)nclistget(node->maps,i);
  	    walk(m,sorted);
	}
	break;
	
    case NCD4_ATTR: /* Depends on its base type */
	walk(node->basetype,sorted);
	break;

    case NCD4_ATTRSET: /* Depends on its contained attributes, but handled after switch */
	/* fall thru */
    default: /* depends on nothing else */
	break;
    }

    /* Do Attributes last */
    for(i=0;i<nclistlength(node->attributes);i++) {
	NCD4node* a = (NCD4node*)nclistget(node->attributes,i);
	walk(a,sorted);
    }		
    ADDNODE(sorted,node);    
}


/**************************************************/
/* Mark the offset and length of each var/field
   inside the raw dapdata.
   Assumes it is called before byte swapping, so we
   need to do swapping of counts and the final remote checksum.
*/

int
NCD4_delimit(NCD4meta* compiler, NCD4node* topvar, void** offsetp)
{
    int ret = NC_NOERR;
    void* offset;

    offset = *offsetp;
    ASSERT((ISTOPLEVEL(topvar)));
    topvar->data.dap4data.memory = offset;
    if(topvar->sort == NCD4_VAR) {
	switch (topvar->subsort) {
        case NC_STRUCT:
            if((ret=delimitStructArray(compiler,topvar,&offset))) goto done;
            break;
        case NC_SEQ:
            if((ret=delimitSeqArray(compiler,topvar,&offset))) goto done;
            break;
	default:
            if((ret=delimitAtomicVar(compiler,topvar,&offset))) goto done;
            break;
	}
    }
    /* Track the variable size (in the dap4 data) but do not include
       any checksum */
    topvar->data.dap4data.size = (d4size_t)(((char*)offset) - ((char*)*offsetp));
    /* extract the dap4 data checksum, if present */
    if(compiler->serial.remotechecksumming) {
	union ATOMICS csum;
        memcpy(csum.u8,offset,CHECKSUMSIZE);
        topvar->data.remotechecksum = csum.u32[0];
	if(compiler->swap) swapinline32(&topvar->data.remotechecksum);
	offset = INCR(offset,CHECKSUMSIZE);
    }
    *offsetp = offset;
done:
    return THROW(ret);
}

/* Includes opaque and enum */
static int
delimitAtomicVar(NCD4meta* compiler, NCD4node* var, void** offsetp)
{
    int ret = NC_NOERR;
    void* offset;
    int typesize;
    d4size_t i;
    d4size_t dimproduct;
    nc_type tid;
    NCD4node* basetype;
    NCD4node* truetype;

    assert(var->sort == NCD4_VAR);
    dimproduct = NCD4_dimproduct(var);
    basetype = var->basetype;
    if(basetype->subsort == NC_OPAQUE)
        return delimitOpaqueVar(compiler,var,offsetp);

    truetype = basetype;
    if(truetype->subsort == NC_ENUM)
        truetype = basetype->basetype;

    offset = *offsetp;
    tid = truetype->subsort;
    typesize = NCD4_typesize(tid);
    if(tid != NC_STRING) {
	offset = INCR(offset,(typesize*dimproduct));
    } else if(tid == NC_STRING) { /* walk the counts */
        unsigned long long count;
        for(i=0;i<dimproduct;i++) {
            /* Get string count */
            count = GETCOUNTER(offset);
            SKIPCOUNTER(offset);
  	    if(compiler->swap) swapinline64(&count);
            /* skip count bytes */
  	    offset = INCR(offset,count);
        }
    }
    *offsetp = offset;
    return THROW(ret);
}

static int
delimitOpaqueVar(NCD4meta* compiler,  NCD4node* var, void** offsetp)
{
    int ret = NC_NOERR;
    void* offset;
    d4size_t i;
    unsigned long long count;
    d4size_t dimproduct = NCD4_dimproduct(var);

    offset = *offsetp;
    for(i=0;i<dimproduct;i++) {
        /* Walk the instances */
        count = GETCOUNTER(offset);
        SKIPCOUNTER(offset);
        if(compiler->swap) swapinline64(&count);
        offset = INCR(offset,count);
    }
    *offsetp = offset;
    return THROW(ret);
}

static int
delimitStructArray(NCD4meta* compiler, NCD4node* varortype,  void** offsetp)
{
    int ret = NC_NOERR;
    void* offset;
    d4size_t i;
    d4size_t dimproduct;
    NCD4node* type;

   if(varortype->sort == NCD4_VAR) {
	dimproduct = NCD4_dimproduct(varortype);	
	type = varortype->basetype;
   } else {
	dimproduct = 1;
	type = varortype;	
   }

    offset = *offsetp;
    for(i=0;i<dimproduct;i++) {
        if((ret=delimitStruct(compiler,type,&offset))) goto done;
    }
    *offsetp = offset;
done:
    return THROW(ret);
}

static int
delimitStruct(NCD4meta* compiler, NCD4node* basetype, void** offsetp)
{
    int ret = NC_NOERR;
    int i;
    void* offset;

    offset = *offsetp;
    /* The fields are associated with the basetype struct */
    for(i=0;i<nclistlength(basetype->vars);i++) {
        NCD4node* field = (NCD4node*)nclistget(basetype->vars,i);
        switch (field->subsort) {
	default:
            if((ret=delimitAtomicVar(compiler,field,&offset))) goto done;
            break;
        case NC_STRUCT: /* recurse */
            if((ret=delimitStructArray(compiler,field,&offset))) goto done;
            break;
        case NC_SEQ:
            if((ret=delimitSeqArray(compiler,field,&offset))) goto done;
            break;
        }
    }
    *offsetp = offset;
done:
    return THROW(ret);
}

static int
delimitSeqArray(NCD4meta* compiler, NCD4node* varortype,  void** offsetp)
{
    int ret = NC_NOERR;
    void* offset;
    d4size_t i;
    d4size_t dimproduct;
    NCD4node* type;

    if(varortype->sort == NCD4_VAR) {
	dimproduct = NCD4_dimproduct(varortype);	
	type = varortype->basetype;
    } else {
	dimproduct = 1;
	type = varortype;	
    }

    offset = *offsetp;
    for(i=0;i<dimproduct;i++) {
        if((ret=delimitSeq(compiler,type,&offset))) goto done;
    }
    *offsetp = offset;
done:
    return THROW(ret);
}

static int
delimitSeq(NCD4meta* compiler, NCD4node* vlentype, void** offsetp)
{
    int ret = NC_NOERR;
    int i;
    void* offset;
    d4size_t recordcount;
    NCD4node* recordtype;

    /* The true type of the record is the basetype of the vlen,
       where the vlen type is the basetype of the var
     */ 
    assert(vlentype->subsort == NC_VLEN);
    recordtype = vlentype->basetype;

    offset = *offsetp;

    /* Get he record count */
    recordcount = GETCOUNTER(offset);
    SKIPCOUNTER(offset);
    if(compiler->swap) swapinline64(&recordcount);

    for(i=0;i<recordcount;i++) {
	switch (recordtype->subsort) {
        case NC_STRUCT:
            if((ret=delimitStructArray(compiler,recordtype,&offset))) goto done;
            break;
        case NC_SEQ:
            if((ret=delimitSeqArray(compiler,recordtype,&offset))) goto done;
            break;
	default:
            if((ret=delimitAtomicVar(compiler,recordtype,&offset))) goto done;
            break;
	}
    }
    *offsetp = offset;
done:
    return THROW(ret);
}

/**************************************************/
/*
Walk the (toplevel) var's data to get to the count'th instance.
For efficiency, it can be supplied with a previous case.
Assumes it is called after byte swapping and offsetting.
Assumes that var is not fixed size.
*/

int 
NCD4_moveto(NCD4meta* compiler, NCD4node* var, d4size_t count, void** offsetp)
{
    int ret = NC_NOERR;
    void* offset = NULL;
    d4size_t startcount = 0;
    NCD4node* basetype = NULL;

    ASSERT((ISTOPLEVEL(var)));

    offset = *offsetp;
    startcount = 0;
    basetype = var->basetype;
    for(;startcount < count;startcount++) {
	if((ret=skipInstance(compiler,basetype,&offset)))
	    goto done;
    }
    *offsetp = offset;
done:
    return THROW(ret);
}

static int
skipInstance(NCD4meta* compiler, NCD4node* type, void** offsetp)
{
    int ret = NC_NOERR;
    void* offset = NULL;
 
    offset = *offsetp;
    switch (type->subsort) {
    case NC_STRUCT:
        if((ret=skipStructInstance(compiler,type,&offset))) goto done;
        break;
    case NC_SEQ:
        if((ret=skipSeqInstance(compiler,type,&offset))) goto done;
        break;
    default:
        if((ret=skipAtomicInstance(compiler,type,&offset))) goto done;
        break;
    }
    *offsetp = offset;
done:
    return THROW(ret);
}

/* Includes opaque and enum */
static int
skipAtomicInstance(NCD4meta* compiler, NCD4node* type, void** offsetp)
{
    int ret = NC_NOERR;
    void* offset;
    d4size_t count;
    int typesize;

    offset = *offsetp;

    switch (type->subsort) {
    default: /* fixed size atomic type */
        typesize = NCD4_typesize(type->meta.id);
	offset = INCR(offset,typesize);
	break;
    case NC_STRING:
        /* Get string count */
        count = GETCOUNTER(offset);
        SKIPCOUNTER(offset);
        /* skip count bytes */
	offset = INCR(offset,count);
	break;
    case NC_OPAQUE:
        /* get count */
        count = GETCOUNTER(offset);
        SKIPCOUNTER(offset);
	offset = INCR(offset,count);
	break;
    case NC_ENUM:
        return skipAtomicInstance(compiler,type->basetype,offsetp);
    }
    *offsetp = offset;
    return THROW(ret);
}

static int
skipStructInstance(NCD4meta* compiler, NCD4node* type,  void** offsetp)
{
    int ret = NC_NOERR;
    d4size_t i,j;
    void* offset;

    offset = *offsetp;
    /* Skip each field */
    for(i=0;i<nclistlength(type->vars);i++) {
        NCD4node* field = (NCD4node*)nclistget(type->vars,i);
	NCD4node* ftype = field->basetype;
        d4size_t dimproduct = NCD4_dimproduct(field);	
	for(j=0;j<dimproduct;j++) {
	    if((ret=skipInstance(compiler,ftype,&offset)))
		goto done;
	}
    }
    *offsetp = offset;
done:
    return THROW(ret);
}

static int
skipSeqInstance(NCD4meta* compiler, NCD4node* vlentype, void** offsetp)
{
    int ret = NC_NOERR;
    d4size_t i;
    void* offset;
    NCD4node* structtype;
    d4size_t recordcount;

    offset = *offsetp;

    structtype = vlentype->basetype;
    ASSERT((structtype->subsort == NC_STRUCT));

    /* Get record count */
    recordcount = GETCOUNTER(offset);
    SKIPCOUNTER(offset);

    for(i=0;i<recordcount;i++) {
	/* Skip a record instance */
	if((ret=skipStructInstance(compiler,structtype,&offset)))
	    goto done;
    }
    *offsetp = offset;
done:
    return THROW(ret);
}


/*********************************************************************
 *   Copyright 2009, UCAR/Unidata
 *   See netcdf/COPYRIGHT file for copying and redistribution conditions.
 *********************************************************************/

#include "includes.h"
#include "nc_iter.h"
#include "nclog.h"

#ifdef ENABLE_BINARY

/* Forward */
static void alignto(int alignment, Bytebuffer* buf, ptrdiff_t base);

static int bin_uid = 0;

static int
bin_charconstant(Generator* generator, Symbol* sym, Bytebuffer* buf, ...)
{
    /* Just transfer charbuf to codebuf */
    Bytebuffer* charbuf;
    va_list ap;
    va_start(ap,buf);
    charbuf = va_arg(ap, Bytebuffer*);
    va_end(ap);
    bbNull(charbuf);
    bbCatbuf(buf,charbuf);
    return 1;
}

static int
bin_constant(Generator* generator, Symbol* sym, NCConstant* con, Bytebuffer* buf,...)
{
    if(con->nctype != NC_ECONST) {
        alignbuffer(con,buf);
    }
    switch (con->nctype) {
    case NC_OPAQUE: {
        unsigned char* bytes = NULL;
        size_t len;
	/* Assume the opaque string has been normalized */
        bytes=makebytestring(con->value.opaquev.stringv,&len);
        bbAppendn(buf,(void*)bytes,len);
	efree(bytes);
    } break;
    case NC_CHAR:
        bbAppendn(buf,&con->value.charv,sizeof(con->value.charv));
        break;
    case NC_BYTE:
        bbAppendn(buf,(void*)&con->value.int8v,sizeof(con->value.int8v));
        break;
    case NC_SHORT:
        bbAppendn(buf,(void*)&con->value.int16v,sizeof(con->value.int16v));
        break;
    case NC_INT:
        bbAppendn(buf,(void*)&con->value.int32v,sizeof(con->value.int32v));
        break;
    case NC_FLOAT:
        bbAppendn(buf,(void*)&con->value.floatv,sizeof(con->value.floatv));
        break;
    case NC_DOUBLE:
        bbAppendn(buf,(void*)&con->value.doublev,sizeof(con->value.doublev));
        break;
    case NC_UBYTE:
        bbAppendn(buf,(void*)&con->value.uint8v,sizeof(con->value.uint8v));
        break;
    case NC_USHORT:
        bbAppendn(buf,(void*)&con->value.uint16v,sizeof(con->value.uint16v));
        break;
    case NC_UINT:
        bbAppendn(buf,(void*)&con->value.uint32v,sizeof(con->value.uint32v));
        break;
    case NC_INT64: {
        union SI64 { char ch[8]; long long i64;} si64;
        si64.i64 = con->value.int64v;
        bbAppendn(buf,(void*)si64.ch,sizeof(si64.ch));
        } break;
    case NC_UINT64: {
        union SU64 { char ch[8]; unsigned long long i64;} su64;
        su64.i64 = con->value.uint64v;
        bbAppendn(buf,(void*)su64.ch,sizeof(su64.ch));
        } break;
    case NC_NIL:
    case NC_STRING: {
        int len = (size_t)con->value.stringv.len;
	if(len == 0 && con->value.stringv.stringv == NULL) {
	    char* nil = NULL;
            bbAppendn(buf,(void*)&nil,sizeof(nil));
	} else {
            char* ptr = (char*)ecalloc(len+1);
	    memcpy(ptr,con->value.stringv.stringv,len);
	    ptr[len] = '\0';
            bbAppendn(buf,(void*)&ptr,sizeof(ptr));
	    ptr = NULL;
        }
	} break;

    default: PANIC1("bin_constant: unexpected type: %d",con->nctype);
    }
    return 1;
}

static int
bin_listbegin(Generator* generator, Symbol* tsym, void* liststate, ListClass lc, size_t size, Bytebuffer* buf, int* uidp, ...)
{
    if(uidp) *uidp = ++bin_uid;
    if(lc == LISTCOMPOUND)
        *((int*)liststate) = bbLength(buf);
    return 1;
}

static int
bin_list(Generator* generator, Symbol* tsym, void* liststate, ListClass lc, int uid, size_t count, Bytebuffer* buf, ...)
{
    if(lc == LISTCOMPOUND) {
        int offsetbase = *((int*)liststate);
        /* Pad for the alignment */
	alignto(tsym->typ.alignment,buf,offsetbase);		
    }
    return 1;
}

static int
bin_listend(Generator* generator, Symbol* tsym, void* liststate, ListClass lc, int uid, size_t count, Bytebuffer* buf, ...)
{
    if(lc == LISTCOMPOUND) {
        int offsetbase = *((int*)liststate);
        /* Pad out the whole instance */
	alignto(tsym->typ.cmpdalign,buf,offsetbase);		
    }
    return 1;
}


static int
bin_vlendecl(Generator* generator, Symbol* tsym, Bytebuffer* buf, int uid, size_t count,...)
{
    va_list ap;
    Bytebuffer* vlenmem;
    nc_vlen_t ptr;
    va_start(ap,count);
    vlenmem = va_arg(ap, Bytebuffer*);
    va_end(ap);
    ptr.len = count;
    ptr.p = bbExtract(vlenmem);
    bbAppendn(buf,(char*)&ptr,sizeof(ptr));
    return 1;
}

static int
bin_vlenstring(Generator* generator, Symbol* sym, Bytebuffer* codebuf, int* uidp, size_t* sizep,...)
{
    Bytebuffer* vlenmem;
    nc_vlen_t ptr;
    va_list ap;
    if(uidp) *uidp = ++bin_uid;
    va_start(ap,sizep);
    vlenmem = va_arg(ap, Bytebuffer*);
    va_end(ap);
    ptr.len = bbLength(vlenmem);
    ptr.p = bbDup(vlenmem);
    bbAppendn(codebuf,(char*)&ptr,sizeof(ptr));
    return 1;
}

static const char zeros[] =
    "\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0";

static void
alignto(int alignment, Bytebuffer* buf, ptrdiff_t base)
{
    int pad = 0;
    ptrdiff_t offset = bbLength(buf);
    offset -= base; /* Need to actually align wrt to the base */
    pad = getpadding(offset,alignment);
    if(pad > 0) {
	bbAppendn(buf,(void*)zeros,pad);
    }
}

/* Define the single static bin data generator  */
static Generator bin_generator_singleton = {
    NULL,
    bin_charconstant,
    bin_constant,
    bin_listbegin,
    bin_list,
    bin_listend,
    bin_vlendecl,
    bin_vlenstring
};
Generator* bin_generator = &bin_generator_singleton;

/**************************************************/

static int bin_generate_data_r(NCConstant* instance, Symbol* tsym, Datalist* fillvalue, Bytebuffer* databuf);

static void
write_alignment(int alignment, Bytebuffer* buf)
{
    int pad = 0;
    ptrdiff_t offset = bbLength(buf);
    pad = getpadding(offset,alignment);
    if(pad > 0) {
	bbAppendn(buf,(void*)zeros,pad);
    }
}

/**
Alternate binary data generator.
Inputs:
	Datalist* data - to use to generate the binary data
	Symbol* tsym - the top-level type for which instances
	               are to be generated
	Datalist* fillvalue - the fillvalue for the toplevel type
	Bytebuffer* databuf - the buffer into which instances are to be stored
*/

int
binary_generate_data(Datalist* data, Symbol* tsym, Datalist* fillvalue, Bytebuffer* databuf)
{
    int stat = NC_NOERR;
    size_t count = data->length;
    size_t i;

    bbClear(databuf);
    for(i=0;i<count;i++) {
	NCConstant* instance = datalistith(data,i);
	if((stat = bin_generate_data_r(instance, tsym, fillvalue, databuf))) goto done;
    }
done:
    return stat;
}

/* Recursive helper that does the bulk of the work */
static int
bin_generate_data_r(NCConstant* instance, Symbol* tsym, Datalist* fillvalue, Bytebuffer* databuf)
{
    int stat = NC_NOERR;

    if(instance->nctype == NC_FILLVALUE) {
        /* replace with fillvalue for the type */
	Datalist* filllist = (fillvalue == NULL ? getfiller(tsym) : fillvalue);
	ASSERT(datalistlen(filllist)==1)
	instance = datalistith(filllist,0);
    }

    switch (tsym->subclass) {
    case NC_PRIM: {
	switch (tsym->nc_id) {
        case NC_CHAR: {
            char* p = NULL;
            NCConstant* tmp = nullconst();
            tmp->nctype = NC_CHAR;
            convert1(instance,tmp);
            p = &tmp->value.charv;;
            bbAppendn(databuf,p,sizeof(char));
            reclaimconstant(tmp);
            } break;
        case NC_BYTE: {
            signed char* p = NULL;
            NCConstant* tmp = nullconst();
            tmp->nctype = NC_BYTE;
            convert1(instance,tmp);
            p = &tmp->value.int8v;
            bbAppendn(databuf,p,sizeof(signed char));
            reclaimconstant(tmp);
            } break;
        case NC_UBYTE: {
            unsigned char* p = NULL;
            NCConstant* tmp = nullconst();
            tmp->nctype = NC_UBYTE;
            convert1(instance,tmp);
            p = &tmp->value.uint8v;
            bbAppendn(databuf,p,sizeof(unsigned char));
            reclaimconstant(tmp);
            } break;
        case NC_SHORT: {
            short* p = NULL;
            NCConstant* tmp = nullconst();
            tmp->nctype = NC_SHORT;
            convert1(instance,tmp);
            p = &tmp->value.int16v;
            bbAppendn(databuf,p,sizeof(short));
            reclaimconstant(tmp);
            } break;
        case NC_USHORT: {
            unsigned short* p = NULL;
            NCConstant* tmp = nullconst();
            tmp->nctype = NC_USHORT;
            convert1(instance,tmp);
            p = &tmp->value.uint16v;
            bbAppendn(databuf,p,sizeof(unsigned short));
            reclaimconstant(tmp);
            } break;
        case NC_INT: {
            int* p = NULL;
            NCConstant* tmp = nullconst();
            tmp->nctype = NC_INT;
            convert1(instance,tmp);
            p = &tmp->value.int32v;
            bbAppendn(databuf,p,sizeof(int));
            reclaimconstant(tmp);
            } break;
        case NC_UINT: {
            unsigned int* p = NULL;
            NCConstant* tmp = nullconst();
            tmp->nctype = NC_UINT;
            convert1(instance,tmp);
            p = &tmp->value.uint32v;
            bbAppendn(databuf,p,sizeof(unsigned int));
            reclaimconstant(tmp);
            } break;
        case NC_INT64: {
            long long* p = NULL;
            NCConstant* tmp = nullconst();
            tmp->nctype = NC_INT64;
            convert1(instance,tmp);
            p = &tmp->value.int64v;
            bbAppendn(databuf,p,sizeof(long long));
            reclaimconstant(tmp);
            } break;
        case NC_UINT64: {
            unsigned long long* p = NULL;
            NCConstant* tmp = nullconst();
            tmp->nctype = NC_UINT64;
            convert1(instance,tmp);
            p = &tmp->value.uint64v;
            bbAppendn(databuf,p,sizeof(unsigned long long));
            reclaimconstant(tmp);
            } break;
        case NC_FLOAT: {
            float* p = NULL;
            NCConstant* tmp = nullconst();
            tmp->nctype = NC_FLOAT;
            convert1(instance,tmp);
            p = &tmp->value.floatv;
            bbAppendn(databuf,p,sizeof(float));
            reclaimconstant(tmp);
            } break;
        case NC_DOUBLE: {
            double* p = NULL;
            NCConstant* tmp = nullconst();
            tmp->nctype = NC_DOUBLE;
            convert1(instance,tmp);
            p = &tmp->value.doublev;
            bbAppendn(databuf,p,sizeof(double));
            reclaimconstant(tmp);
            } break;
        case NC_STRING: {
            char* p = NULL;
            NCConstant* tmp = nullconst();
            tmp->nctype = NC_STRING;
            convert1(instance,tmp);
            p = emalloc(tmp->value.stringv.len+1);
	    memcpy(p,tmp->value.stringv.stringv,tmp->value.stringv.len);
	    p[tmp->value.stringv.len] = '\0';
            bbAppendn(databuf,&p,sizeof(char*));
            reclaimconstant(tmp);
            } break;
	default: stat = NC_EINTERNAL; goto done; /* Should never happen */
	} break; /*switch*/
	} break; /*NC_PRIM*/
	
    case NC_ENUM: {
	Symbol* basetype = tsym->typ.basetype;
	/* Pretend */
	stat = bin_generate_data_r(instance,basetype,fillvalue,databuf);
        } break;
    case NC_OPAQUE: {
	unsigned char* bytes = NULL;
	size_t len = 0;
	if(instance->nctype != NC_OPAQUE)
	    {stat = NC_EBADTYPE; goto done;}
	/* Assume the opaque string has been normalized */
        bytes=makebytestring(instance->value.opaquev.stringv,&len);
	if(bytes == NULL) {stat = NC_ENOMEM; goto done;}
        bbAppendn(databuf,(void*)bytes,len);
	free(bytes);
        } break;
    case NC_VLEN: {
	Datalist* sublist = NULL;
	Bytebuffer* vlendata = NULL;
	nc_vlen_t p;
	if(instance->nctype != NC_COMPOUND) {
	    nclog(NCLOGERR,"Translating vlen: expected sublist");
	    stat = NC_EBADTYPE; goto done;
	}
	sublist = instance->value.compoundv;
	vlendata = bbNew();
	if((stat = binary_generate_data(sublist,tsym->typ.basetype,NULL,vlendata))) goto done;
	p.len = datalistlen(sublist);
	p.p = bbContents(vlendata);
        bbAppendn(databuf,(char*)&p,sizeof(nc_vlen_t));
        } break;
    case NC_COMPOUND: { /* The really hard one */
	size_t nfields, fid, i;
	Datalist* cmpd = instance->value.compoundv;
        write_alignment(tsym->typ.cmpdalign,databuf);
        /* Get info about each field in turn and build it*/
        nfields = listlength(tsym->subnodes);
        for(fid=0;fid<nfields;fid++) {
	    Symbol* field = listget(tsym->subnodes,fid);
	    NCConstant* fieldinstance = datalistith(cmpd,fid);
	    int ndims = field->typ.dimset.ndims;
	    size_t arraycount;
	    if(ndims == 0) {
	        ndims=1; /* fake the scalar case */
	    }
  	    /* compute the total number of elements in the field array */
	    arraycount = 1;
	    for(i=0;i<ndims;i++) arraycount *= field->typ.dimset.dimsyms[i]->dim.declsize;
	    write_alignment(field->typ.alignment,databuf);
	    /* Write the instances */
	    for(i=0;i<arraycount;i++) {
	        if((stat = bin_generate_data_r(fieldinstance, field->typ.basetype, NULL, databuf))) goto done;
	    }
	}		
        } break;

    default: stat = NC_EINTERNAL; goto done; /* Should never happen */
    }
done:
    return stat;
}

/**
Internal equivalent of ncaux_reclaim_data.
*/

/* It is helpful to have a structure that contains memory and an offset */
typedef struct Reclaim {char* memory; ptrdiff_t offset;} Reclaim;

static ptrdiff_t read_alignment(ptrdiff_t offset, unsigned long alignment);
static int bin_reclaim_datar(Symbol* tsym, Reclaim* reclaim);
static int bin_reclaim_usertype(Symbol* tsym, Reclaim* reclaim);
static int bin_reclaim_compound(Symbol* tsym, Reclaim* reclaim);
static int bin_reclaim_vlen(Symbol* tsym, Reclaim* reclaim);
static int bin_reclaim_enum(Symbol* tsym, Reclaim* reclaim);
static int bin_reclaim_opaque(Symbol* tsym, Reclaim* reclaim);

int
binary_reclaim_data(Symbol* tsym, void* memory, size_t count)
{
    int stat = NC_NOERR;
    size_t i;
    Reclaim reclaimer;
    
    if(tsym == NULL
       || (memory == NULL && count > 0))
        {stat = NC_EINVAL; goto done;}
    if(memory == NULL || count == 0)
        goto done; /* ok, do nothing */
    reclaimer.offset = 0;
    reclaimer.memory = memory;
    for(i=0;i<count;i++) {
	if((stat=bin_reclaim_datar(tsym,&reclaimer))) /* reclaim one instance */
	    break;
    }
done:
    return stat;
}

/* Recursive type walker: reclaim a single instance */
static int
bin_reclaim_datar(Symbol* tsym, Reclaim* reclaimer)
{
    int stat = NC_NOERR;
    
    switch  (tsym->subclass) {
    case NC_CHAR: case NC_BYTE: case NC_UBYTE:
    case NC_SHORT: case NC_USHORT:
    case NC_INT: case NC_UINT: case NC_FLOAT:
    case NC_INT64: case NC_UINT64: case NC_DOUBLE:
        reclaimer->offset += tsym->typ.size;
	break;
#ifdef USE_NETCDF4
    case NC_STRING: {
	char** sp = (char**)(reclaimer->memory+reclaimer->offset);
        /* Need to reclaim string */
	if(*sp != NULL) efree(*sp);
	reclaimer->offset += tsym->typ.size;
	} break;
    default:
    	/* reclaim a user type */
	stat = bin_reclaim_usertype(tsym,reclaimer);
#else
    default:
	stat = NC_ENOTNC4;
#endif
	break;
    }
    return stat;
}
	
static int
bin_reclaim_usertype(Symbol* tsym, Reclaim* reclaimer)
{
    int stat = NC_NOERR;

    /* Get info about the xtype */
    switch (tsym->subclass) {
    case NC_OPAQUE: stat = bin_reclaim_opaque(tsym,reclaimer); break;
    case NC_ENUM: stat = bin_reclaim_enum(tsym,reclaimer); break;
    case NC_VLEN: stat = bin_reclaim_vlen(tsym,reclaimer); break;
    case NC_COMPOUND: stat = bin_reclaim_compound(tsym,reclaimer); break;
    default:
        stat = NC_EINVAL;
	break;
    }
    return stat;
}

static ptrdiff_t
read_alignment(ptrdiff_t offset, unsigned long alignment)
{
    size_t delta = (offset % alignment);
    if(delta == 0) return offset;
    return offset + (alignment - delta);
}


static int
bin_reclaim_vlen(Symbol* tsym, Reclaim* reclaimer)
{
    int stat = NC_NOERR;
    size_t i;
    Symbol* basetype = tsym->typ.basetype;
    nc_vlen_t* vl = (nc_vlen_t*)(reclaimer->memory+reclaimer->offset);

    /* Free up each entry in the vlen list */
    if(vl->p != NULL) {
	Reclaim vreclaimer;
	vreclaimer.memory = vl->p;
	vreclaimer.offset = 0;
        for(i=0;i<vl->len;i++) {
	    vreclaimer.offset = read_alignment(vreclaimer.offset,basetype->typ.alignment);
	    if((stat = bin_reclaim_datar(basetype,&vreclaimer))) goto done;
	    vreclaimer.offset += basetype->typ.size;
	}
	reclaimer->offset += tsym->typ.size;
	efree(vl->p);
    }
done:
    return stat;
}

static int
bin_reclaim_enum(Symbol* tsym, Reclaim* reclaimer)
{
    return bin_reclaim_datar(tsym->typ.basetype,reclaimer);
}

static int
bin_reclaim_opaque(Symbol* tsym, Reclaim* reclaimer)
{
    /* basically a fixed size sequence of bytes */
    reclaimer->offset += tsym->typ.size;
    return NC_NOERR;
}

static int
bin_reclaim_compound(Symbol* tsym, Reclaim* reclaimer)
{
    int stat = NC_NOERR;
    int nfields;
    size_t fid, i, arraycount;
    ptrdiff_t saveoffset;

    reclaimer->offset = read_alignment(reclaimer->offset,tsym->typ.cmpdalign);
    saveoffset = reclaimer->offset;

    /* Get info about each field in turn and reclaim it */
    nfields = listlength(tsym->subnodes);
    for(fid=0;fid<nfields;fid++) {
	Symbol* field = listget(tsym->subnodes,fid);
	int ndims = field->typ.dimset.ndims;
	/* compute the total number of elements in the field array */
	for(i=0;i<ndims;i++) arraycount *= field->typ.dimset.dimsyms[i]->dim.declsize;
	reclaimer->offset = read_alignment(reclaimer->offset,field->typ.alignment);
	for(i=0;i<arraycount;i++) {
	    if((stat = bin_reclaim_datar(field->typ.basetype, reclaimer))) goto done;
	}		
    }
    reclaimer->offset = saveoffset;
    reclaimer->offset += tsym->typ.size;
done:
    return stat;
}

#endif /*ENABLE_BINARY*/


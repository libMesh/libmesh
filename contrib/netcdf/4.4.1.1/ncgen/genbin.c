/*********************************************************************
 *   Copyright 1993, UCAR/Unidata
 *   See netcdf/COPYRIGHT file for copying and redistribution conditions.
 *   $Header: /upc/share/CVS/netcdf-3/ncgen/genbin.c,v 1.4 2010/05/27 21:34:17 dmh Exp $
 *********************************************************************/

#include "includes.h"
#include <ctype.h>	/* for isprint() */

#ifdef ENABLE_BINARY

#undef TRACE

/* Forward*/
static void genbin_defineattr(Symbol* asym);
static void genbin_definevardata(Symbol* vsym);
static int  genbin_write(Generator*,Symbol*,Bytebuffer*,int,size_t*,size_t*);
static int genbin_writevar(Generator*,Symbol*,Bytebuffer*,int,size_t*,size_t*);
static int genbin_writeattr(Generator*,Symbol*,Bytebuffer*,int,size_t*,size_t*);
#ifdef USE_NETCDF4
static void genbin_deftype(Symbol* tsym);
static void genbin_definespecialattributes(Symbol* var);
#endif

/*
 * Generate C code for creating netCDF from in-memory structure.
 */
void
gen_netcdf(const char *filename)
{
    int stat, ncid;
    int idim, ivar, iatt;
    int ndims, nvars, natts, ngatts;

#ifdef USE_NETCDF4
    int ntyps, ngrps, igrp;
#endif

    Bytebuffer* databuf = bbNew();

    ndims = listlength(dimdefs);
    nvars = listlength(vardefs);
    natts = listlength(attdefs);
    ngatts = listlength(gattdefs);
#ifdef USE_NETCDF4
    ntyps = listlength(typdefs);
    ngrps = listlength(grpdefs);
#endif /*USE_NETCDF4*/

    /* Turn on logging */
#ifdef LOGGING
    nc_set_log_level(ncloglevel);
#endif

    /* create netCDF file, uses NC_CLOBBER mode */
    cmode_modifier |= NC_CLOBBER;
#ifdef USE_NETCDF4
    if(!usingclassic)
        cmode_modifier |= NC_NETCDF4;
#endif

    stat = nc_create(filename, cmode_modifier, &ncid);
    check_err(stat,__LINE__,__FILE__);

    /* ncid created above is also root group*/
    rootgroup->ncid = ncid;

#ifdef USE_NETCDF4
    /* Define the group structure */
    /* walking grdefs list will do a preorder walk of all defined groups*/
    for(igrp=0;igrp<ngrps;igrp++) {
	Symbol* gsym = (Symbol*)listget(grpdefs,igrp);
	if(gsym == rootgroup) continue; /* ignore root group*/
	stat = nc_def_grp(gsym->container->ncid,gsym->name,&gsym->ncid);
	check_err(stat,__LINE__,__FILE__);
    }
#endif

#ifdef USE_NETCDF4
    /* Define the types*/
    if (ntyps > 0) {
	int ityp;
	for(ityp = 0; ityp < ntyps; ityp++) {
	    Symbol* tsym = (Symbol*)listget(typdefs,ityp);
	    genbin_deftype(tsym);
	}
    }
#endif

    /* define dimensions from info in dims array */
    if (ndims > 0) {
        for(idim = 0; idim < ndims; idim++) {
            Symbol* dsym = (Symbol*)listget(dimdefs,idim);
	    stat = nc_def_dim(dsym->container->ncid,
			      dsym->name,
			      (dsym->dim.isunlimited?NC_UNLIMITED:dsym->dim.declsize),
			      &dsym->ncid);
	    check_err(stat,__LINE__,__FILE__);
       }
    }

    /* define variables from info in vars array */
    if (nvars > 0) {
	for(ivar = 0; ivar < nvars; ivar++) {
            Symbol* vsym = (Symbol*)listget(vardefs,ivar);
	    if (vsym->typ.dimset.ndims > 0) {	/* a dimensioned variable */
		/* construct a vector of dimension ids*/
		int dimids[NC_MAX_VAR_DIMS];
		for(idim=0;idim<vsym->typ.dimset.ndims;idim++)
		    dimids[idim] = vsym->typ.dimset.dimsyms[idim]->ncid;
		stat = nc_def_var(vsym->container->ncid,
				  vsym->name,
			          vsym->typ.basetype->ncid,
		        	  vsym->typ.dimset.ndims,
				  dimids,
				  &vsym->ncid);
	    } else { /* a scalar */
		stat = nc_def_var(vsym->container->ncid,
				  vsym->name,
			          vsym->typ.basetype->ncid,
		        	  vsym->typ.dimset.ndims,
				  NULL,
				  &vsym->ncid);
	    }
	    check_err(stat,__LINE__,__FILE__);
	}
    }

#ifdef USE_NETCDF4
    /* define special variable properties */
    if(nvars > 0) {
	for(ivar = 0; ivar < nvars; ivar++) {
            Symbol* var = (Symbol*)listget(vardefs,ivar);
	    genbin_definespecialattributes(var);
	}
    }
#endif /*USE_NETCDF4*/

/* define global attributes */
    if(ngatts > 0) {
	for(iatt = 0; iatt < ngatts; iatt++) {
	    Symbol* gasym = (Symbol*)listget(gattdefs,iatt);
	    genbin_defineattr(gasym);
	}
    }

    /* define per-variable attributes */
    if(natts > 0) {
	for(iatt = 0; iatt < natts; iatt++) {
	    Symbol* asym = (Symbol*)listget(attdefs,iatt);
	    genbin_defineattr(asym);
	}
    }

    if (nofill_flag) {
	stat = nc_set_fill(rootgroup->ncid, NC_NOFILL, 0);
	check_err(stat,__LINE__,__FILE__);
    }

    /* leave define mode */
    stat = nc_enddef(rootgroup->ncid);
    check_err(stat,__LINE__,__FILE__);

    if(!header_only) {
        /* Load values into those variables with defined data */
        if(nvars > 0) {
            for(ivar = 0; ivar < nvars; ivar++) {
                Symbol* vsym = (Symbol*)listget(vardefs,ivar);
                if(vsym->data != NULL) {
                    bbClear(databuf);
                    genbin_definevardata(vsym);
                }
            }
        }
    }
    bbFree(databuf);
}

#ifdef USE_NETCDF4

#if 0
Turn off for now.

static void
genbin_defineglobalspecials(void)
{
    int stat = NC_NOERR;
    const char* format = NULL;
    if(usingclassic) return;
    if(!/*Main.*/format_attribute) return;
    /* Watch out, this is a global Attribute */
    format = kind_string(/*Main.*/format_flag);
    stat = nc_put_att_text(rootgroup->ncid,NC_GLOBAL,"_Format",strlen(format),format);
    check_err(stat,__LINE__,__FILE__);
}
#endif /*0*/

static void
genbin_definespecialattributes(Symbol* var)
{
    int stat;
    Specialdata* special = &var->var.special;
    if(special->flags & _STORAGE_FLAG) {
        int storage = special->_Storage;
        size_t* chunks = special->_ChunkSizes;
        if(special->nchunks == 0 || chunks == NULL) chunks = NULL;
        stat = nc_def_var_chunking(var->container->ncid,
                                   var->ncid,
                                   (storage == NC_CONTIGUOUS?NC_CONTIGUOUS
                                                            :NC_CHUNKED),
                                   chunks);
        check_err(stat,__LINE__,__FILE__);
    }
    if(special->flags & _FLETCHER32_FLAG) {
        stat = nc_def_var_fletcher32(var->container->ncid,
                                     var->ncid,
                                     special->_Fletcher32);
        check_err(stat,__LINE__,__FILE__);
    }
    if(special->flags & (_DEFLATE_FLAG | _SHUFFLE_FLAG)) {
        stat = nc_def_var_deflate(var->container->ncid,
                                  var->ncid,
                                  (special->_Shuffle == 1?1:0),
                                  (special->_DeflateLevel >= 0?1:0),
                                  (special->_DeflateLevel >= 0?special->_DeflateLevel
                                                              :0));
        check_err(stat,__LINE__,__FILE__);
    }
    if(special->flags & _ENDIAN_FLAG) {
        stat = nc_def_var_endian(var->container->ncid,
                                 var->ncid,
                                 (special->_Endianness == NC_ENDIAN_LITTLE?
                                        NC_ENDIAN_LITTLE
                                       :NC_ENDIAN_BIG));
        check_err(stat,__LINE__,__FILE__);
    }
    if(special->flags & _NOFILL_FLAG) {
        stat = nc_def_var_fill(var->container->ncid,
                                 var->ncid,
		                 (special->_Fill?NC_FILL:NC_NOFILL),
                                 NULL);
        check_err(stat,__LINE__,__FILE__);
    }
}
#endif /*USE_NETCDF4*/


void
cl_netcdf(void)
{
    int stat;
    stat = nc_close(rootgroup->ncid);
    check_err(stat,__LINE__,__FILE__);
}

#ifdef USE_NETCDF4
/*
Generate type definitions
*/
static void
genbin_deftype(Symbol* tsym)
{
    unsigned long i;
    int stat;

    ASSERT(tsym->objectclass == NC_TYPE);
    switch (tsym->subclass) {
    case NC_PRIM: break; /* these are already taken care of*/
    case NC_OPAQUE:
	stat = nc_def_opaque(tsym->container->ncid,
			     tsym->typ.size,
			     tsym->name,
			     &tsym->ncid);
        check_err(stat,__LINE__,__FILE__);
	break;
    case NC_ENUM:
      {
        Bytebuffer* datum;
        Datalist* ecdl;
	stat = nc_def_enum(tsym->container->ncid,
			   tsym->typ.basetype->typ.typecode,
			   tsym->name,
			   &tsym->ncid);
        check_err(stat,__LINE__,__FILE__);
	datum = bbNew();
        ecdl = builddatalist(1);
        dlextend(ecdl); /* make room for one constant*/
	ecdl->length = 1;
	for(i=0;i<listlength(tsym->subnodes);i++) {
	  Symbol* econst = (Symbol*)listget(tsym->subnodes,i);
	  ASSERT(econst->subclass == NC_ECONST);
	  generator_reset(bin_generator,NULL);
	  bbClear(datum);
	  generate_basetype(econst->typ.basetype,&econst->typ.econst,datum,NULL,bin_generator);
	  stat = nc_insert_enum(tsym->container->ncid,
				tsym->ncid,
				econst->name,
				bbContents(datum));
	  check_err(stat,__LINE__,__FILE__);
	}
	bbFree(datum);
    dlfree(&ecdl);
      }
      break;
    case NC_VLEN:
	stat = nc_def_vlen(tsym->container->ncid,
			   tsym->name,
			   tsym->typ.basetype->ncid,
			   &tsym->ncid);
        check_err(stat,__LINE__,__FILE__);
	break;
    case NC_COMPOUND:
	stat = nc_def_compound(tsym->container->ncid,
			       tsym->typ.size,
			       tsym->name,
			       &tsym->ncid);
        check_err(stat,__LINE__,__FILE__);
	for(i=0;i<listlength(tsym->subnodes);i++) {
	    Symbol* efield = (Symbol*)listget(tsym->subnodes,i);
	    ASSERT(efield->subclass == NC_FIELD);
	    if(efield->typ.dimset.ndims == 0){
	        stat = nc_insert_compound(
				tsym->container->ncid,
				tsym->ncid,
				efield->name,
			        efield->typ.offset,
				efield->typ.basetype->ncid);
	    } else {
		int j;
		int dimsizes[NC_MAX_VAR_DIMS]; /* int because inside compound */
		/* Generate the field dimension constants*/
		for(j=0;j<efield->typ.dimset.ndims;j++) {
		     unsigned int size = efield->typ.dimset.dimsyms[j]->dim.declsize;
		     dimsizes[j] = size;
		}
	        stat = nc_insert_array_compound(
				tsym->container->ncid,
				tsym->ncid,
				efield->name,
			        efield->typ.offset,
				efield->typ.basetype->ncid,
				efield->typ.dimset.ndims,
				dimsizes);
	    }
            check_err(stat,__LINE__,__FILE__);
	}
	break;
    default: panic("definectype: unexpected type subclass");
    }
}
#endif /*USE_NETCDF4*/

static void
genbin_defineattr(Symbol* asym)
{
    Bytebuffer* databuf = bbNew();
    generator_reset(bin_generator,NULL);
    generate_attrdata(asym,bin_generator,(Writer)genbin_write,databuf);
}


/* Following is patterned after the walk functions in semantics.c */
static void
genbin_definevardata(Symbol* vsym)
{
    Bytebuffer* databuf;
    if(vsym->data == NULL) return;
    databuf = bbNew();
    generator_reset(bin_generator,NULL);
    generate_vardata(vsym,bin_generator,(Writer)genbin_write,databuf);
}

static int
genbin_write(Generator* generator, Symbol* sym, Bytebuffer* memory,
             int rank, size_t* start, size_t* count)
{
    if(sym->objectclass == NC_ATT)
	return genbin_writeattr(generator,sym,memory,rank,start,count);
    else if(sym->objectclass == NC_VAR)
	return genbin_writevar(generator,sym,memory,rank,start,count);
    else
	PANIC("illegal symbol for genbin_write");
    return NC_EINVAL;
}

static int
genbin_writevar(Generator* generator, Symbol* vsym, Bytebuffer* memory,
                int rank, size_t* start,
#ifdef USE_NOFILL
	        size_t* indices
#else
	        size_t* count
#endif
	       )
{
    int stat = NC_NOERR;
    char* data = bbContents(memory);
#ifdef USE_NOFILL
    size_t count[NC_MAX_VAR_DIMS];
    { int i; for(i=0;i<rank;i++) count[i] = indices[i] - start[i];}
#endif

#ifdef GENDEBUG
    {
    int i;
    fprintf(stderr,"startset = [");
    for(i=0;i<rank;i++)
	fprintf(stderr,"%s%lu",(i>0?", ":""),(unsigned long)start[i]);
    fprintf(stderr,"] ");
    fprintf(stderr,"countset = [");
    for(i=0;i<rank;i++)
	fprintf(stderr,"%s%lu",(i>0?", ":""),(unsigned long)count[i]);
    fprintf(stderr,"]\n");
    fflush(stderr);
    }
#endif

    if(rank == 0) {
	size_t count[1] = {1};
        stat = nc_put_var1(vsym->container->ncid, vsym->ncid, count, data);
    } else {
        stat = nc_put_vara(vsym->container->ncid, vsym->ncid, start, count, data);
    }
    check_err(stat,__LINE__,__FILE__);
    bbClear(memory);
    return stat;
}

static int
genbin_writeattr(Generator* generator, Symbol* asym, Bytebuffer* databuf,
           int rank, size_t* start, size_t* count)
{
    int stat;
    size_t len;
    Datalist* list;
    int varid, grpid, typid;
    Symbol* basetype = asym->typ.basetype;

    grpid = asym->container->ncid;
    varid = (asym->att.var == NULL?NC_GLOBAL : asym->att.var->ncid);
    typid = basetype->ncid;

    list = asym->data;
    len = list->length;

    /* Use the specialized put_att_XX routines if possible*/
    if(isprim(basetype->typ.typecode)) {
      switch (basetype->typ.typecode) {
      case NC_BYTE: {
        signed char* data = (signed char*)bbContents(databuf);
        stat = nc_put_att_schar(grpid,varid,asym->name,typid,len,data);
        check_err(stat,__LINE__,__FILE__);
      } break;
      case NC_CHAR: {
        char* data = (char*)bbContents(databuf);
		size_t slen = bbLength(databuf);
		/* Revise length if slen == 0 */
		if(slen == 0) {
          bbAppend(databuf,'\0');
          /* bbAppend frees the memory pointed to by char* data,
             so re-assign.  See Coverity issue: 1265731.*/
          data = (char*)bbContents(databuf);
          slen++;
        }
        stat = nc_put_att_text(grpid,varid,asym->name,slen,data);
        check_err(stat,__LINE__,__FILE__);
      } break;
      case NC_SHORT: {
        short* data = (short*)bbContents(databuf);
        stat = nc_put_att_short(grpid,varid,asym->name,typid,len,data);
        check_err(stat,__LINE__,__FILE__);
            } break;
            case NC_INT: {
                int* data = (int*)bbContents(databuf);
                stat = nc_put_att_int(grpid,varid,asym->name,typid,len,data);
                check_err(stat,__LINE__,__FILE__);
            } break;
            case NC_FLOAT: {
                float* data = (float*)bbContents(databuf);
                stat = nc_put_att_float(grpid,varid,asym->name,typid,len,data);
                check_err(stat,__LINE__,__FILE__);
            } break;
            case NC_DOUBLE: {
                double* data = (double*)bbContents(databuf);
                stat = nc_put_att_double(grpid,varid,asym->name,typid,len,data);
                check_err(stat,__LINE__,__FILE__);
            } break;
#ifdef USE_NETCDF4
	    case NC_STRING: {
		const char** data;
	        data = (const char**)bbContents(databuf);
                stat = nc_put_att_string(grpid,varid,asym->name,
					 bbLength(databuf)/sizeof(char*),
				         data);
		} break;
            case NC_UBYTE: {
                unsigned char* data = (unsigned char*)bbContents(databuf);
                stat = nc_put_att_uchar(grpid,varid,asym->name,typid,len,data);
                check_err(stat,__LINE__,__FILE__);
            } break;
            case NC_USHORT: {
                unsigned short* data = (unsigned short*)bbContents(databuf);
                stat = nc_put_att_ushort(grpid,varid,asym->name,typid,len,data);
                check_err(stat,__LINE__,__FILE__);
            } break;
            case NC_UINT: {
                unsigned int* data = (unsigned int*)bbContents(databuf);
                stat = nc_put_att_uint(grpid,varid,asym->name,typid,len,data);
                check_err(stat,__LINE__,__FILE__);
            } break;
            case NC_INT64: {
                long long* data = (long long*)bbContents(databuf);
                stat = nc_put_att_longlong(grpid,varid,asym->name,typid,len,data);
                check_err2(stat,asym->lineno,__LINE__,__FILE__);
            } break;
            case NC_UINT64: {
                unsigned long long* data = (unsigned long long*)bbContents(databuf);
                stat = nc_put_att_ulonglong(grpid,varid,asym->name,typid,len,data);
                check_err(stat,__LINE__,__FILE__);
            } break;
#endif
            default: PANIC1("genbin_defineattr: unexpected basetype: %d",basetype->typ.typecode);
	}
    } else { /* use the generic put_attribute for user defined types*/
	const char* data;
	data = (const char*)bbContents(databuf);
        stat = nc_put_att(grpid,varid,asym->name,typid,
			        len,(void*)data);
        check_err(stat,__LINE__,__FILE__);
#ifdef GENDEBUG
	{
	    char out[4096];
	    memset(out,0x77,sizeof(out));
	    stat = nc_get_att(grpid,varid,asym->name,&out);
            check_err(stat,__LINE__,__FILE__);
	}
#endif
    }
    return stat;
}

#endif /*ENABLE_BINARY*/

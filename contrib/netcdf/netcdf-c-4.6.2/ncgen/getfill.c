/*********************************************************************
 *   Copyright 1993, UCAR/Unidata
 *   See netcdf/COPYRIGHT file for copying and redistribution conditions.
 *   $Header: /upc/share/CVS/netcdf-3/ncgen/getfill.c,v 1.8 2010/04/14 22:04:56 dmh Exp $
 *********************************************************************/

#include "includes.h"
#include "dump.h"

/* mnemonic*/
#define TOPLEVEL 1

/*Forward*/
static void fill(Symbol* tsym, Datalist*);
static void fillarray(Symbol* tsym, Dimset* dimset, int index, Datalist*);
static void filllist(Symbol* tvsym, Datalist* dl);

/* Construct a Datalist representing a complete fill value.
   for a specified variable or type. Cache if needed later.
   The rules are as follows
   1. If the tvsym argument is a variable and it has a _FillValue
      attribute, then use that value.
   2. If the tvsym argunment is a variable with no specified
      _FillValue, then create one based on the variable's basetype.
   3. If the tvsym argument is a type, then if it already has a fill value,
      then use that, otherwise build and cache a fillvalue appropriate
      for that type.
*/

Datalist*
getfiller(Symbol* tvsym)
{
    Datalist* filler = NULL;
    ASSERT(tvsym->objectclass == NC_VAR || tvsym->objectclass == NC_TYPE);
    if(tvsym->objectclass == NC_VAR) {
        if(tvsym->var.special->_Fillvalue !=  NULL) {
	    /* We have a _FillValue Attribute specified */
	    filler = tvsym->var.special->_Fillvalue;
        } else { /* otherwise create a fillvalue for the base type */
	    filler = getfiller(tvsym->typ.basetype);
	}
    } else { /* (tvsym->objectclass == NC_TYPE) */
        if(tvsym->typ._Fillvalue ==  NULL) {
	    /* create and cache */
            filler = builddatalist(0);
	    fill(tvsym,filler);
	    tvsym->typ._Fillvalue = filler;
	}
	filler = tvsym->typ._Fillvalue;
    }
#ifdef GENDEBUG2
    dumpdatalist(filler,"getfiller");
#endif
    return filler;
}

static void
fill(Symbol* tsym, Datalist* filler)
{
    unsigned long i;
    NCConstant* con = NULL;

    ASSERT(tsym->objectclass == NC_TYPE);
    switch (tsym->subclass) {
    case NC_ENUM: case NC_OPAQUE: case NC_PRIM:
        con = nullconst();
        con->nctype = tsym->typ.typecode;
        nc_getfill(con,tsym);
        dlappend(filler,con);
	break;
    case NC_COMPOUND: {
	/* Given a compound type, the fill will be a sublist
	   consisting itself of N constants, where N is the number of fields.
	   Non-array fields will be direct, array fields will be sublists.
	*/
	Datalist* cmpdlist = builddatalist(listlength(tsym->subnodes)); /* list of field constants */
        for(i=0;i<listlength(tsym->subnodes);i++) {
	    NCConstant* fieldinstance = NULL;
	    Symbol* field = (Symbol*)listget(tsym->subnodes,i);
	    if(field->typ.dimset.ndims > 0) {	
  	        /*  Build a sublist for this field */
	        Datalist* arraydata = builddatalist(0);
                fillarray(field->typ.basetype,&field->typ.dimset,0,arraydata);
		/* Convert to a compound constant */
	        fieldinstance = list2const(arraydata);
		dlappend(cmpdlist,fieldinstance);
		fieldinstance = NULL;
	    } else
		/* Append directly to cmpdlist*/
		filllist(field->typ.basetype,cmpdlist);
        }	  
        /* Add compound instance to the filler list */
	con = list2const(cmpdlist);
        dlappend(filler,con);
	con = NULL;
	} break;
    case NC_VLEN: {
        Datalist* vlensublist = NULL;
	/* The vlist will have a sublist hanging off of fille, with
	   a single compound constant in filler */
	vlensublist = builddatalist(0);
	filllist(tsym->typ.basetype,vlensublist); /* generate a single instance*/
	con = builddatasublist(vlensublist);
        dlappend(filler,con);
	} break;
    default: PANIC1("fill: unexpected subclass %d",tsym->subclass);
    }
}

static void
filllist(Symbol* tsym, Datalist* dl)
{
    int i;
    Datalist* sublist;
    NCConstant* con = NULL;

    ASSERT(tsym->objectclass == NC_TYPE);
    switch (tsym->subclass) {
    case NC_ENUM: case NC_OPAQUE: case NC_PRIM:
        con = nullconst();
        con->nctype = tsym->typ.typecode;
        nc_getfill(con,tsym);
	dlappend(dl,con);
	break;
    case NC_COMPOUND:
	sublist = builddatalist(listlength(tsym->subnodes));
        for(i=0;i<listlength(tsym->subnodes);i++) {
	    Symbol* field = (Symbol*)listget(tsym->subnodes,i);
	    filllist(field->typ.basetype,sublist);
        }	  
	con = builddatasublist(sublist);
	dlappend(dl,con);
	break;
    case NC_VLEN:
	sublist = builddatalist(0);
	filllist(tsym->typ.basetype,sublist); /* generate a single instance*/
	con = builddatasublist(sublist);
	dlappend(dl,con);
	break;
    default: PANIC1("fill: unexpected subclass %d",tsym->subclass);
    }
}

/* Create an array of fill values of basetype given the dimset. This
is recursive over the dimensions as specified by the index argumen
*/
static void
fillarray(Symbol* basetype, Dimset* dimset, int index, Datalist* arraylist)
{
    int i;
    Symbol* dim = dimset->dimsyms[index];
    unsigned int size = dim->dim.declsize;
    int isunlimited = (size == 0);
    int lastdim = (index == (dimset->ndims - 1));

    if(isunlimited) {
	/* do a single entry to satisfy*/
        if(lastdim) {
	    filllist(basetype,arraylist);
	} else {
	    fillarray(basetype->typ.basetype,dimset,index+1,arraylist);
	}
    } else { /* bounded*/
        if(lastdim) {
	    for(i=0;i<size;i++) filllist(basetype,arraylist);
	} else {
	    for(i=0;i<size;i++) {
	        fillarray(basetype->typ.basetype,dimset,index+1,arraylist);
	    }
	}
    }
}

/*
 * Given primitive netCDF type, return a default fill_value appropriate for
 * that type.
 */
void
nc_getfill(NCConstant* value, Symbol* tsym)
{
    switch(value->nctype) {
      case NC_CHAR: value->value.charv = NC_FILL_CHAR; break;
      case NC_BYTE: value->value.int8v = NC_FILL_BYTE; break;
      case NC_SHORT: value->value.int16v = NC_FILL_SHORT; break;
      case NC_INT: value->value.int32v = NC_FILL_INT; break;
      case NC_FLOAT: value->value.floatv = NC_FILL_FLOAT; break;
      case NC_DOUBLE: value->value.doublev = NC_FILL_DOUBLE; break;
      case NC_UBYTE: value->value.uint8v = NC_FILL_UBYTE; break;
      case NC_USHORT: value->value.uint16v = NC_FILL_USHORT; break;
      case NC_UINT: value->value.uint32v = NC_FILL_UINT; break;
      case NC_INT64: value->value.int64v = NC_FILL_INT64; break;
      case NC_UINT64: value->value.uint64v = NC_FILL_UINT64; break;
      case NC_STRING:
        value->value.stringv.stringv = nulldup(NC_FILL_STRING);
        value->value.stringv.len = (int)strlen(NC_FILL_STRING);
	/* Exception: if string is null, then make it's length be 1 */
	if(value->value.stringv.len == 0)
	    value->value.stringv.len = 1;
	break;
      case NC_OPAQUE:
        value->value.opaquev.len = 2;
        value->value.opaquev.stringv = nulldup("00");
	break;
      case NC_ENUM: {
	Symbol* econst;
	NCConstant* eccon;
	if(tsym == NULL)
	    derror("nc_getfill: no enum type specified");
	/* Get the first value */
	if(tsym->subclass != NC_ENUM)
	    derror("nc_getfill: expected enum type");
	if(listlength(tsym->subnodes) == 0)
	    derror("nc_getfill: empty enum type");
	econst = listget(tsym->subnodes,0);
	eccon = econst->typ.econst;
	switch (eccon->nctype) {
	case NC_BYTE:
	case NC_SHORT:
	case NC_INT:
	case NC_UBYTE:
	case NC_USHORT:
	case NC_UINT:
	case NC_INT64:
	case NC_UINT64:
	    value->value = eccon->value;
	    break;
	default:
	    derror("nc_getfill: illegal enum basetype");
	}
      } break;
      default:
	derror("nc_getfill: unrecognized type: %d",value->nctype);
    }
}

char* 
nc_dfaltfillname(nc_type nctype)
{
    switch (nctype) {
    case NC_BYTE: return "NC_FILL_BYTE";
    case NC_CHAR: return "NC_FILL_CHAR";
    case NC_SHORT: return "NC_FILL_SHORT";
    case NC_INT: return "NC_FILL_INT";
    case NC_FLOAT: return "NC_FILL_FLOAT";
    case NC_DOUBLE: return "NC_FILL_DOUBLE";
    case NC_UBYTE: return "NC_FILL_UBYTE";
    case NC_USHORT: return "NC_FILL_USHORT";
    case NC_UINT: return "NC_FILL_UINT";
    case NC_INT64: return "NC_FILL_INT64";
    case NC_UINT64: return "NC_FILL_UINT64";
    case NC_STRING: return "NC_FILL_STRING";
    default: PANIC("unexpected default fill name");
    }
    return NULL;
}


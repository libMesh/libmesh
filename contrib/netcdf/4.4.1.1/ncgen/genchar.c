/*********************************************************************
 *   Copyright 2009, UCAR/Unidata
 *   See netcdf/COPYRIGHT file for copying and redistribution conditions.
 *********************************************************************/

#include "includes.h"
#include "odom.h"

/******************************************************/
/* Code for generating char variables etc; mostly
   language independent */
/******************************************************/

/*Forward*/
static size_t gen_charconstant(NCConstant*, Bytebuffer*, int fillchar);
static int getfillchar(Datalist* fillsrc);
static void gen_leafchararray(Dimset*,int,Datalist*,Bytebuffer*, int);
#if 0
static void gen_chararrayr(Dimset*,int,int,Bytebuffer*,Datalist*,int,int,int);
#endif
static NCConstant* makeconst(int lineno, int len, char* str);

/*
Matching strings to char variables, attributes, and vlen
constants is challenging because it is desirable to mimic
the original ncgen(3). The "algorithms" used there have no
simple characterization (such as "abc" == {'a','b','c'}).
So, this rather ugly code is kept in this file
and a variety of heuristics are used to mimic ncgen3.

The core algorithm is as follows.
1. Assume we have a set of dimensions D1..Dn.
   Any of the Di may be unlimited,
   but it is assumed that the sizes of the Di are all known.
2. Given a sequence of string or character constants
   C1..Cm, our goal is to construct a single string
   whose length is the cross product of D1 thru Dn.
3. For purposes of this algorithm, character constants
   are treated as strings of size 1.
4. Construct Dx = cross product of D1 thru D(n-1).
5. For each constant Ci, add fill characters, if necessary,
   so that its length is a multiple of Dn.
6. Concatenate the modified C1..Cm to produce string S.
7. Add fill characters to S to make its size be a multiple of
   Dn.
8. If S is longer than the Dx * Dn, then truncate
   and generate a warning.

Two special cases:
1. character vlen: char(*) vlen_t.
    For this case, we simply concat all the elements.
2. character attribute.
    For this case, we simply concat all the elements.
*/

void
gen_chararray(Dimset* dimset, int dimindex, Datalist* data, Bytebuffer* charbuf, Datalist* fillsrc)
{
    int fillchar = getfillchar(fillsrc);
    int rank = rankfor(dimset);
    int firstunlim = findunlimited(dimset,0);
    int nunlim = countunlimited(dimset);
    int nc3unlim = (nunlim <= 1 && (firstunlim == 0 || firstunlim == rank)); /* netcdf-3 case of at most 1 unlim in 0th dimension */

    /* Case: netcdf3 case */
    if(nc3unlim) {
	gen_leafchararray(dimset,0,data,charbuf,fillchar);
	return;
    }

    /* else generate should have done all the hard work */
    gen_leafchararray(dimset,dimindex,data,charbuf,fillchar);
}

#if 0
/* Recursive helper */
static void
gen_chararrayr(Dimset* dimset, int dimindex,
               Bytebuffer* databuf, Datalist* data, int fillchar,
	       int unitsize, int expectedsize)
{
    int i;
    size_t dimsize = declsizefor(dimset,dimindex);
    int rank = dimset->ndims;
    int firstunlim = findunlimited(dimset,0);
    int lastunlimited = findlastunlimited(dimset);
    int nextunlimited = findunlimited(dimset,dimindex+1);
    int islastgroup = (lastunlimited == rank || dimindex >= lastunlimited || dimindex == rank-1);
    Odometer* subodom = NULL;

    ASSERT(rank > 0);
    ASSERT((islastgroup));

    /* we should be at a list of simple constants */
    for(i=0;i<data->length;i++) {
        NCConstant* c = datalistith(data,i);
        ASSERT(!islistconst(c));
        if(isstringable(c->nctype)) {
            int j;
            size_t constsize;
            constsize = gen_charconstant(c,databuf,fillchar);
            if(constsize % unitsize > 0) {
                size_t padsize = unitsize - (constsize % unitsize);
                for(j=0;j<padsize;j++) bbAppend(databuf,fillchar);
            }
        } else {
            semwarn(constline(c),
                   "Encountered non-string and non-char constant in datalist; ignored");
        }
    }/* for */

    /* If |databuf| > expectedsize, complain: exception is zero length */
    if(bbLength(databuf) == 0 && expectedsize == 1) {
        /* this is okay */
    } else if(bbLength(databuf) > expectedsize) {
        semwarn(data->data[0].lineno,"character data list too long; expected %d character constant, found %d: ",expectedsize,bbLength(databuf));
    } else {
        size_t bufsize = bbLength(databuf);
        /* Pad to size dimproduct size */
        if(bufsize % expectedsize > 0) {
            size_t padsize = expectedsize - (bufsize % expectedsize);
            for(i=0;i<padsize;i++) bbAppend(databuf,fillchar);
        }
    }
}
#endif

void
gen_charattr(Datalist* data, Bytebuffer* databuf)
{
    gen_charvlen(data,databuf);
}

void
gen_charvlen(Datalist* data, Bytebuffer* databuf)
{
    int i;
    NCConstant* c;

    ASSERT(bbLength(databuf) == 0);

    for(i=0;i<data->length;i++) {
        c = datalistith(data,i);
        if(isstringable(c->nctype)) {
            (void)gen_charconstant(c,databuf,NC_FILL_CHAR);
        } else {
            semerror(constline(c),
                     "Encountered non-string and non-char constant in datalist");
            return;
        }
    }
}

static size_t
gen_charconstant(NCConstant* con, Bytebuffer* databuf, int fillchar)
{
    /* Following cases should be consistent with isstringable */
    size_t constsize = 1;
    switch (con->nctype) {
    case NC_CHAR:
        bbAppend(databuf,con->value.charv);
        break;
    case NC_BYTE:
        bbAppend(databuf,con->value.int8v);
        break;
    case NC_UBYTE:
        bbAppend(databuf,con->value.uint8v);
        break;
    case NC_STRING:
        constsize = con->value.stringv.len;
        bbAppendn(databuf,con->value.stringv.stringv,
                         con->value.stringv.len);
        bbNull(databuf);
        break;
    case NC_FILL:
        bbAppend(databuf,fillchar);
        break;
    default:
        PANIC("unexpected constant type");
    }
    return constsize;
}

static int
getfillchar(Datalist* fillsrc)
{
    /* Determine the fill char */
    int fillchar = 0;
    if(fillsrc != NULL && fillsrc->length > 0) {
        NCConstant* ccon = fillsrc->data;
        if(ccon->nctype == NC_CHAR) {
            fillchar = ccon->value.charv;
        } else if(ccon->nctype == NC_STRING) {      
            if(ccon->value.stringv.len > 0) {
                fillchar = ccon->value.stringv.stringv[0];
            }
        }
    }
    if(fillchar == 0) fillchar = NC_FILL_CHAR; /* default */
    return fillchar;
}

/* I think there is a flaw in the ncgen manual
   about handling something like this:
    dimensions:
      n = 8 ;
    variables:
      char cdata2(n) ;
    data:
      cdata2 = '\000','\001','\002','\177','\200','\201','\376','\377';
    The rules would say that each of the 8 char constants must
    be padded to length 8. I think this is only true if the dimension
    is unlimited, and even then, I am not sure.
 */

static void
gen_leafchararray(Dimset* dimset, int dimindex, Datalist* data,
                   Bytebuffer* charbuf, int fillchar)
{
    int i;
    size_t expectedsize,xproduct,unitsize;
    int rank = rankfor(dimset);

    ASSERT(bbLength(charbuf) == 0);
    ASSERT((findlastunlimited(dimset) == rank
            || findlastunlimited(dimset) == dimindex));

    /*
    There are a number of special cases that must be
    considered, mostly driven by the need to keep consistent
    with ncgen3.  These cases are driven by the number of
    dimensions, which dimensions are unlimited (if any), etc.

    The general rule is based on the size of the last
    dimension, we compute the required size (after padding)
    of each string constant. Expected size is then the size
    of concat of the string constants after padding.

    There is another special case used for back compatibility with ncgen3.
    In the datalist, all sequences of character constants (i.e. 'X')
    are concatenated into a single string; the result, however is not
    concatenated with any trailing or leading string (with double quotes).
    */

    /* Rebuild the datalist to merge 'x' constants */
    {
	int i,cccount = 0;
	/* Do initial walk */
	for(i=0;i<datalistlen(data);i++) {
	    NCConstant* con = datalistith(data,i);
	    if(consttype(con) == NC_CHAR || consttype(con) == NC_BYTE) {
		cccount++;
	    }
	}
	if(cccount > 1) {
	    char* accum = (char*)malloc(cccount+1);
	    int len = 0;
	    Datalist* newlist = builddatalist(datalistlen(data));
	    int lineno = 0;
            NCConstant* con;
	    for(i=0;i<datalistlen(data);i++) {
	        con = datalistith(data,i);
	        if(consttype(con) == NC_CHAR || consttype(con) == NC_BYTE) {
		    if(len == 0)
			lineno = constline(con);
		    accum[len] = con->value.charv;
		    len++;
		} else {
		    if(len > 0) {
			con = makeconst(lineno,len,accum);
			len = 0;
			lineno = 0;
		    }
		    dlappend(newlist,con);
		}
	    }
	    /* deal with any unclosed strings */
	    if(len > 0) {
		con = makeconst(lineno,len,accum);
		len = 0;
		lineno = 0;
	        dlappend(newlist,con);
	    }
	    free(accum);
	    data = newlist;
	}
    }

    /* Compute crossproduct up to (but not including) the last dimension */
    xproduct = crossproduct(dimset,dimindex,rank-1);

    /* Start casing it out */
    if(rank == 0) {
        unitsize = 1;
        expectedsize = (xproduct * unitsize);
    } else if(rank == 1) {
	unitsize = 1;
        expectedsize = (xproduct * declsizefor(dimset,rank-1));
    } else if(isunlimited(dimset,rank-1)) {/* last dimension is unlimited */
        unitsize = 1;
        expectedsize = (xproduct*declsizefor(dimset,rank-1));
    } else { /* rank > 0 && last dim is not unlimited */
	unitsize =  declsizefor(dimset,rank-1);
        expectedsize = (xproduct * unitsize);
    }

    for(i=0;i<data->length;i++) {
        NCConstant* c = datalistith(data,i);
        ASSERT(!islistconst(c));
        if(isstringable(c->nctype)) {
            int j;
            size_t constsize;
            constsize = gen_charconstant(c,charbuf,fillchar);
            if(constsize == 0 || constsize % unitsize > 0) {
                size_t padsize = unitsize - (constsize % unitsize);
                for(j=0;j<padsize;j++) bbAppend(charbuf,fillchar);
            }
        } else {
            semwarn(constline(c),"Encountered non-string and non-char constant in datalist; ignored");
        }
    }
    /* If |databuf| > expectedsize, complain: exception is zero length */
    if(bbLength(charbuf) == 0 && expectedsize == 1) {
        /* this is okay */
    } else if(bbLength(charbuf) > expectedsize) {
        semwarn(data->data[0].lineno,"character data list too long; expected %d character constant, found %d: ",expectedsize,bbLength(charbuf));
    } else {
        size_t bufsize = bbLength(charbuf);
        /* Pad to size dimproduct size */
        if(bufsize % expectedsize > 0) {
            size_t padsize = expectedsize - (bufsize % expectedsize);
            for(i=0;i<padsize;i++) bbAppend(charbuf,fillchar);
        }
    }
}

/* Create a new string constant */
static NCConstant*
makeconst(int lineno, int len, char* str)
{
    NCConstant* con = (NCConstant*)malloc(sizeof(NCConstant));
    con->nctype = NC_STRING;
    con->lineno = lineno;
    con->filled = 0;
    con->value.stringv.len = len;
    /* We cannot use strdup because str might have embedded nuls */
    con->value.stringv.stringv = (char*)malloc(len+1);
    memcpy((void*)con->value.stringv.stringv,(void*)str,len);
    con->value.stringv.stringv[len] = '\0';
    return con;
}

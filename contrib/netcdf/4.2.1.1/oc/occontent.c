/* Copyright 2009, UCAR/Unidata and OPeNDAP, Inc.
   See the COPYRIGHT file for more information. */

#include "config.h"
#include "ocinternal.h"
#include "occontent.h"
#include "ocdebug.h"

/* Mnemonic*/
#define ISPACKED 1

/* Define the skipstate flags */
#ifdef OCDEBUG
typedef enum Skipstate {
SKIPFIELDS = 1, /* skip instance without leading tag or arraycount */
SKIPINSTANCE = 2, /*skip leading sequence tag or array counts */
SKIPWHOLE = 4 /* skip complete object */
} Skipstate;
#else
#define SKIPFIELDS	1 /* skip instance without leading tag or arraycount */
#define SKIPINSTANCE	2 /*skip leading sequence tag or array counts */
#define SKIPWHOLE	4 /* skip complete object */
#endif

/* Forward*/
static OCmode modetransition(OCnode*,OCmode);
static int ocgetsequencetag(XXDR* xdrs);
static int ocskipcounts(XXDR* xdrs, OCnode*, off_t expected);

static int ocarrayith(OCstate*, OCcontent*, OCcontent*, size_t);
static int ocsequenceith(OCstate*, OCcontent*, OCcontent*, size_t);
static int ocfieldith(OCstate*, OCcontent*, OCcontent*, size_t);

static size_t ocarraycount(OCstate*, struct OCcontent*);
static size_t ocsequencecount(OCstate*, struct OCcontent*);
static size_t ocfieldcount(OCstate*, struct OCcontent*);
static size_t ocprimcount(OCstate*, struct OCcontent*);

static OCcontent* occlearcontent(struct OCstate* state, OCcontent* content);


#ifdef OCDEBUG
static void
report(size_t memsize, char* memory)
{
#ifdef OCIGNORE
    switch(memsize) {
    case 1:
        oc_log(LOGNOTE,"reading xdr: %lu bytes = |%2x|",(unsigned long)memsize,memory[0]);
	break;
    case 4:
        oc_log(LOGNOTE,"reading xdr: %lu bytes = |%4lu|",(unsigned long)memsize,*(unsigned int*)memory);
	break;
    case 8:
        oc_log(LOGNOTE,"reading xdr: %lu bytes = |%8lu|",(unsigned long)memsize,*(unsigned long long*)memory);
	break;
    default:
        oc_log(LOGNOTE,"reading xdr: %lu bytes",(unsigned long)memsize);
	break;
    }    
    oc_log(LOGNOTE,"reading xdr: %lu bytes",(unsigned long)memsize);
#endif
}

static char*
modestring(OCmode mode)
{
    switch(mode) {
    case OCFIELDMODE: return "FIELD";
    case OCSEQUENCEMODE: return "SEQUENCE";
    case OCARRAYMODE: return "ARRAY";
    case OCPRIMITIVEMODE: return "PRIMITIVE";
    case OCNULLMODE: return "NULL";
    case OCEMPTYMODE: return "EMPTY";
    }
    return "?";
}

void
octrace1(char* proc, OCstate* state, OCcontent* content, int start, int count)
{
    unsigned long pos = xxdr_getpos(content->tree->data.xdrs);
    fprintf(stderr,"trace: %s mode=%s node=%s (%s)",
	    proc,modestring(content->mode),content->node->fullname,
            octypetostring(content->node->octype));
    if(content->packed)
        fprintf(stderr," packed=%d",content->packed);
    if(count >= 0)
        fprintf(stderr," start=%lu count=%lu",
	    (unsigned long)start,(unsigned long)count);
    else
        fprintf(stderr," index=%lu",(unsigned long)start);
    fprintf(stderr," xdrs.pos=%lu",pos);
    fprintf(stderr,"\n");
    if(content->cache.valid) {
	fprintf(stderr,"\tcache{index=%lu maxindex=%lu offset=%lu}\n",
		(unsigned long)content->cache.index,
		(unsigned long)content->cache.maxindex,
		(unsigned long)content->cache.offset);
    }
    fflush(stderr);
}

void
octrace(char* proc, OCstate* state, OCcontent* content, int index)
{
    octrace1(proc,state,content,index,-1);
}


#else
#define octrace(proc,state,content,index)
#define octrace1(proc,state,content,start,count)
#define report(memsize,memory)
#endif /*OCDEBUG*/

OCmode
ocgetmode(OCcontent* content)
{
    return (content == NULL ? OCNULLMODE : content->mode);
}

OCcontent*
ocnewcontent(OCstate* state)
{
    OCcontent* content;
    if(state == NULL) return NULL;
    content = state->contentlist;
    /* Search for an unused content node*/
    while(content != NULL && content->mode != OCEMPTYMODE) {
	content = content->next;
    }
    if(content == NULL) {
	content = (OCcontent*)ocmalloc(sizeof(OCcontent));
	MEMCHECK(content,(OCcontent*)NULL);
	content->magic = OCMAGIC;
        content->next = state->contentlist;
        state->contentlist = content;
    }
    return occlearcontent(state,content);
}

void
ocfreecontent(OCstate* state, OCcontent* content)
{
    if(content != NULL) {content->mode = OCEMPTYMODE;}
}

static OCcontent*
occlearcontent(struct OCstate* state, OCcontent* content)
{
    /* save fields that should not be cleared */
    unsigned int magic = content->magic;
    OCcontent* next = content->next;
    memset((void*)content,sizeof(OCcontent),0);
    /* set/restore non-null fields */
    content->magic = magic;
    content->next = next;
    content->state = state;
    content->mode = OCNULLMODE;
    return content;
}

static OCcontent*
ocsetcontent(OCcontent* childcontent, OCcontent* parent, OCnode* node, int packed)
{
    childcontent->state = parent->state;
    childcontent->cache.valid = 0;
    childcontent->node = node;
    childcontent->tree = node->root->tree;
    childcontent->mode = modetransition(node,parent->mode);
    childcontent->packed = packed;
    return childcontent;
}

#ifdef OCIGNORE
static OCcontent*
occlonecontent(OCstate* state, OCcontent* content)
{
    OCcontent* clone = ocnewcontent(state);
    clone->mode = content->mode;
    clone->node = content->node;
    clone->cache = content->indexcache;
    return clone;
}
#endif

OCerror
ocdataith(OCstate* state, OCcontent* parent, size_t index, OCcontent* child)
{
    OCerror ocerr = OC_NOERR;
    switch (parent->mode) {
    case OCARRAYMODE:
        ocerr = ocarrayith(state,parent,child,index);
	break;
    case OCSEQUENCEMODE:
	ocerr = ocsequenceith(state,parent,child,index);
	break;
    case OCFIELDMODE:
	ocerr = ocfieldith(state,parent,child,index);
	break;
    default: return OC_EINVAL;
    }
    if(ocerr == OC_EDATADDS)
	ocdataddsmsg(state,parent->tree);
    return ocerr;
}

OCerror
ocdatacount(OCstate* state, OCcontent* current, size_t* sizep)
{
    OCerror ocerr = OC_NOERR;
    size_t count = 0;
    switch(current->mode) {
    case OCARRAYMODE:
	count = ocarraycount(state,current);
	break;
    case OCSEQUENCEMODE:
	count = ocsequencecount(state,current);
	break;
    case OCFIELDMODE:
	count = ocfieldcount(state,current);
	break;
    case OCPRIMITIVEMODE:
	count = ocprimcount(state,current);        
	break;
    default: 
	return OC_EINVAL;
    }
    current->cache.maxindex = ocmax(count,current->cache.maxindex);
    if(sizep) *sizep = count;
    return ocerr;
}

OCerror
ocrootdata(OCstate* state, OCnode* root, OCcontent* content)
{
    OCtree* tree;
    if(state == NULL || root == NULL || content == NULL)
	return OCTHROW(OC_EINVAL);
    if(root->tree == NULL) return OCTHROW(OC_EINVAL);
    tree = root->tree;
    if(tree->dxdclass != OCDATADDS) return OCTHROW(OC_ENODATA);
    if(tree->nodes == NULL) return OCTHROW(OC_EINVAL);
    if(tree->data.xdrs == NULL)
	return OCTHROW(OC_EXDR);

    occlearcontent(state,content);
    content->mode = OCFIELDMODE;
    content->node = root;
    content->tree = tree;

    content->cache.index = 0;
    content->cache.maxindex = oclistlength(content->node->subnodes);
    content->cache.offset = 0;
    content->cache.valid = 1;

    return OCTHROW(OC_NOERR);
}

/* Remember: we are operating wrt the datadds count, not the dds count */
static OCerror
ocarrayith(OCstate* state, OCcontent* content, OCcontent* elemcontent, size_t index)
{
    unsigned int i;
    int stat = OC_NOERR;
    XXDR* xdrs;
    int packed, scalar;
    OCtype etype,octype;
    OCnode* node;
    int startindex = 0;

    octrace("ocarrayith", state, content, index);

    if(state == NULL || content == NULL) return OCTHROW(OC_EINVAL);
    if(content->mode != OCARRAYMODE) return OCTHROW(OC_EINVAL);

    etype = content->node->etype;
    octype = content->node->octype;
    node = content->node;
    scalar = (node->array.rank == 0 ? 1 : 0);
    packed = (!scalar && octype == OC_Primitive &&
              (etype == OC_Byte || etype == OC_UByte || etype == OC_Char));

    xdrs = content->tree->data.xdrs;
    if(xdrs == NULL) return OCTHROW(OC_EXDR);

    if(!content->cache.valid) {
        content->cache.index = 0; /* because we will have to walk to index'th data */
        content->cache.maxindex = totaldimsize(node);
	content->cache.valid = 1;
        /* skip past the initial counts, if any */
	if(!scalar) {
	    if(!ocskipcounts(xdrs,node,node->skip.count)) return OCTHROW(OC_EDATADDS);
    	}
	/* checkpoint xdr position */
	content->cache.offset = xxdr_getpos(xdrs);
    }

    /* move to the checkpoint position */
    startindex = content->cache.index;
    if(!xxdr_setpos(xdrs,content->cache.offset)) return xdrerror();

    /* skip to the index'th item */
    if(packed) {
        content->cache.index = 0; /* keep at beginning */
    } else {
        for(i=startindex;i<index;i++) {
            stat = ocskipinstance(node,xdrs,SKIPINSTANCE,NULL);
            if(stat != OC_NOERR) return OCTHROW(stat);
	}
        content->cache.index = index; /* now we are at the index'th item */
    }

    /* update cache */
    content->cache.index = index;
    content->cache.offset = xxdr_getpos(xdrs);

    /* set up the content for the current item in the array */
    ocsetcontent(elemcontent,content,node,packed); /*keep same node */
    if (index == content->cache.maxindex) {
        /* mark eod */
	elemcontent->mode = OCNULLMODE;	
    }

    return OCTHROW(stat);
}

static int
ocsequenceith(OCstate* state, OCcontent* content, OCcontent* structcontent, size_t index)
{
    unsigned int i;
    int stat = OC_NOERR;
    XXDR* xdrs;
    OCtype octype,etype;
    int packed,scalar;
    OCnode* node = content->node;
    int startindex, tag;

    octrace("ocsequenceith", state, content, index);

    if(state == NULL || content == NULL) goto einval;
    if(content->mode != OCSEQUENCEMODE) goto einval;
    if(node->octype != OC_Sequence) goto einval;

    octype = node->octype;
    etype = node->etype;
    scalar = (node->array.rank == 0 ? 1 : 0);
    packed = (!scalar && octype == OC_Primitive &&
              (etype == OC_Byte || etype == OC_UByte || etype == OC_Char));

    xdrs = content->tree->data.xdrs;
    if(xdrs == NULL) goto exdr;

    if(!content->cache.valid) {
	content->cache.valid = 1;
        content->cache.index = 0;
        content->cache.maxindex = 0;
	content->cache.offset = xxdr_getpos(xdrs);
    }

    /* move to checkpoint position*/
    startindex = content->cache.index;
    if(!xxdr_setpos(xdrs,content->cache.offset)) goto exdr;

    /* Walk past the first (index-1) records */
    for(tag=StartOfSequence,i=startindex;i<index;i++) {
        /* skip instance, including tag, but leave xdr at next tag */
        stat = ocskipinstance(node,xdrs,SKIPINSTANCE,&tag);
	if(stat != OC_NOERR) goto done;
        if(tag == EndOfSequence)
            break;
        if(tag != StartOfSequence) {
            oc_log(LOGERR,"missing/invalid begin/end record marker\n");
            goto einvalcoords;
        }
    }
    if(stat != OC_NOERR) goto done;

    /* update cache; should be pointing to index'th record tag */
    content->cache.index = index;
    content->cache.maxindex = ocmax(index,content->cache.index);

    /* this is a bit (too) tricky */
    /* move to point to fields */
    tag = ocgetsequencetag(xdrs);
    if(tag == EndOfSequence) {
	/* point past end of sequence */
        content->cache.offset = xxdr_getpos(xdrs);
    } else {/*tag == StartOfSequence*/
	/* point to (next) start of sequence */
        content->cache.offset = xxdr_getpos(xdrs) - XDRUNIT;
    }
    /* at this point, xdrs should point unconditionally
       past the index'th tag */

    /* Set state of new content: keep same node */
    ocsetcontent(structcontent,content,node,packed);
    if(tag == EndOfSequence) {
        /* mark eod */
	structcontent->mode = OCNULLMODE;	
    }
done:
    return OCTHROW(stat);
einval:
    stat = OC_EINVAL;
    goto done;
exdr:
    stat = OC_EXDR;
    goto done;
einvalcoords:
    stat = OC_EINVALCOORDS;
    goto done;
}

/*
The ocfieldcontent procedure has to deal with the fact
that the dap constraints may have removed some fields
from the datadds and hence some fields may have no
representation in the xdr data (or compiled data).
Assume that xdr points to start of 0th field.
*/
static int
ocfieldith(OCstate* state, OCcontent* content, OCcontent* fieldcontent, size_t index)
{
    unsigned int i;
    int stat = OC_NOERR;
    XXDR* xdrs;
    OCtype octype,etype;
    int packed;
    int isscalar;
    OCnode* node;
    int startindex;

    octrace("ocfieldith", state, content, index);

    if(state == NULL || content == NULL) return OCTHROW(OC_EINVAL);
    if(content->mode != OCFIELDMODE) return OCTHROW(OC_EINVAL);

    node = content->node;
    octype = node->octype;
    etype = node->etype;
    isscalar = (node->array.rank == 0 ? 1 : 0);
    packed = (!isscalar && octype == OC_Primitive
              && (etype == OC_Byte || etype == OC_UByte || etype == OC_Char));

    xdrs = content->tree->data.xdrs;
    if(xdrs == NULL) return OCTHROW(OC_EXDR);

    if(!content->cache.valid) {
        content->cache.index = 0;
        content->cache.maxindex = oclistlength(node->subnodes);
	content->cache.valid = 1;
	/* checkpoint xdr position */
	content->cache.offset = xxdr_getpos(xdrs);
    }

    /* move to the checkpoint position */
    startindex = content->cache.index;
    if(!xxdr_setpos(xdrs,content->cache.offset)) return xdrerror();

    switch (octype) {
    case OC_Sequence: /* assume xdrs points past sequence tag */
    case OC_Grid: /* Note that the Grid array is field 0 and the maps are 1..nsubnodes*/
    case OC_Dataset:
    case OC_Structure:
	/* walk to (i-1)'th field */
        for(i=startindex;i<index;i++) { /* walk field by field */
  	    OCnode* ithfield = (OCnode*)oclistget(node->subnodes,i);
	    stat = ocskipinstance(ithfield,xdrs,SKIPWHOLE,NULL);
	    if(stat != OC_NOERR) return OCTHROW(stat);
	}
	break;

    default: return OCTHROW(OC_EINVAL);
    }

    /* update cache */
    content->cache.index = index;
    content->cache.offset = xxdr_getpos(xdrs);
    /* Set state of new content: node changes to field node */
    ocsetcontent(fieldcontent,content,
		     (OCnode*)oclistget(node->subnodes,index),
		     packed);
    if(index >= content->cache.maxindex) {
        /* mark eod */
        fieldcontent->mode = OCNULLMODE;	
    }

    return OCTHROW(stat);
}

/*
In order to actually extract data,
one must move to the specific primitive typed
field containing the data of interest by using
ocfieldcontent().
Then, oc_getcontent() is invoked to extract
some subsequence of items from the field.
Note that oc_getcontent() will also work for scalars,
but the start must be zero and the count must be one.
*/

int
ocgetcontent(OCstate* state, OCcontent* content, void* memory, size_t memsize,
                 size_t start, size_t count)
{
    int stat = OC_NOERR;
    XXDR* xdrs;
    OCtype etype, octype;
    int isscalar, packed;
    size_t elemsize, totalsize;
    OCnode* node = content->node;

    octrace1("ocgetcontent", state, content, start, count);

    if(state == NULL || content == NULL || memory == NULL)
	{OCTHROWCHK(stat=OC_EINVAL); goto done;}
    if(content->mode != OCPRIMITIVEMODE || node->octype != OC_Primitive)
	{OCTHROWCHK(stat=OC_EINVAL); goto done;}

    octype = node->octype;
    etype = node->etype;

    isscalar = (node->array.rank == 0);
    packed = (!isscalar && octype == OC_Primitive
              && (etype == OC_Byte || etype == OC_UByte || etype == OC_Char));

    if(isscalar && (start != 0 || count != 1))
	{OCTHROWCHK(stat=OC_EINVALCOORDS); goto done;}

    /* validate memory space*/
    elemsize = octypesize(etype);
    totalsize = elemsize*count;
    if(memsize < totalsize) return OCTHROW(OC_EINVAL);

    xdrs = content->tree->data.xdrs;
    if(xdrs == NULL) return OCTHROW(OC_EXDR);

    /* Need to setup the cache */
    if(!content->cache.valid) {
        content->cache.valid = 1;
        content->cache.index = 0;
        content->cache.maxindex = totaldimsize(content->node);
	if(!ocskipcounts(xdrs,content->node,content->cache.maxindex))
	    return OCTHROW(OC_EXDR);
        content->cache.offset = xxdr_getpos(xdrs);
    }

    if(content->cache.valid && content->cache.maxindex < (start+count))
	return OCTHROW(OC_ENODATA);

    /* utilize the cache */
    if(!xxdr_setpos(xdrs,content->cache.offset)) return OCTHROW(OC_EXDR);

    /* Extract the data */
    stat = ocxdrread(content,xdrs,(char*)memory,memsize,start,count);

#ifdef OCDEBUG
    report(memsize,memory+start);
#endif

    /* Update the cache */
    if(!packed) {
        content->cache.index = (start+count);
        content->cache.offset = xxdr_getpos(xdrs);
    }

done:
    return OCTHROW(stat);
}

static size_t
ocfieldcount(OCstate* state, OCcontent* content)
{
    OCnode* node = content->node;
    size_t count;
    OCASSERT((node != NULL));
    count = oclistlength(node->subnodes);
    return count;
}

static size_t
ocarraycount(OCstate* state, OCcontent* content)
{
    unsigned int count;
    OCnode* node = content->node;

    OCASSERT((node != NULL));
    OCASSERT((content->mode == OCARRAYMODE));

    count = totaldimsize(node);

#ifdef VERIFY
    if(node->array.rank > 0) {
	off_t checkpoint;
	XXDR* xdrs;
	unsigned int xdrcount;
        /* verify against xdr */
	xdrs = content->tree->data.xdrs;
        OCASSERT((xdrs != NULL));
        /* checkpoint current location */
        checkpoint = xxdr_getpos(xdrs);
        /* extract the count*/
        if(!xxdr_uint(xdrs,&xdrcount)) return 0;
	if(xdrcount != count) return 0;
        /* return to checkpoint position*/
        if(!xxdr_setpos(xdrs,checkpoint)) return 0;
    }
#endif /*VERIFY*/
    return (size_t)count;
}

/* Counting records actually requires walking the xdr packet
   so it is not necessarily cheap*/
static size_t
ocsequencecount(OCstate* state, OCcontent* content)
{
    size_t count;
    OCnode* node = content->node;
    XXDR* xdrs;
    off_t checkpoint;

    OCASSERT((node != NULL));
    OCASSERT((node->octype == OC_Sequence));
    OCASSERT((content->mode == OCSEQUENCEMODE));

    xdrs = content->tree->data.xdrs;
    OCASSERT((xdrs != NULL));

    /* checkpoint location */
    checkpoint = xxdr_getpos(xdrs);

    for(count=0;;count++) {
	int tag;
        OCerror stat = ocskipinstance(node,xdrs,SKIPINSTANCE,&tag);
	if(stat != OC_NOERR) {count = 0; break;}
        if(tag == EndOfSequence) {
            break; /* done with the count*/
        } else if(tag != StartOfSequence) {
            oc_log(LOGERR,"missing/invalid begin/end record marker\n");
	    return 0;
	}
    }

    /* move back to checkpoint position*/
    if(!xxdr_setpos(xdrs,checkpoint)) return 0;

    return count;
}

static size_t
ocprimcount(OCstate* state, OCcontent* content)
{
    unsigned int count;
    OCnode* node = content->node;

    OCASSERT((node != NULL));
    OCASSERT((content->mode == OCPRIMITIVEMODE));

    count = totaldimsize(node);

#ifdef VERIFY
    if(node->array.rank > 0) {
	off_t checkpoint;
	XXDR* xdrs;
	unsigned int xdrcount;
        /* verify against xdr */
	xdrs = content->tree->data.xdrs;
        OCASSERT((xdrs != NULL));
        /* checkpoint current location */
        checkpoint = xxdr_getpos(xdrs);
        /* extract the count*/
        if(!xxdr_uint(xdrs,&xdrcount)) return 0;
	if(xdrcount != count) return 0;
        /* return to checkpoint position*/
        if(!xxdr_setpos(xdrs,checkpoint)) return 0;
    }
#endif /*VERIFY*/
    return (size_t)count;
}

static OCmode
modetransition(OCnode* node, OCmode srcmode)
{
    OCmode  newmode = OCNULLMODE;
    switch (srcmode) {
    case OCARRAYMODE:
	switch (node->octype) {
	case OC_Sequence:
	    newmode = OCSEQUENCEMODE;
	    break;
	case OC_Grid:
	case OC_Structure:
	    newmode = OCFIELDMODE;
	    break;
	default:
	    break;
	}
	break;

    case OCSEQUENCEMODE:
	switch (node->octype) {
	default:
	    newmode = OCFIELDMODE;
	    break;
	}
	break;

    case OCFIELDMODE:
	switch (node->octype) {
	case OC_Sequence:
	case OC_Grid:
	case OC_Structure:
	    newmode = OCARRAYMODE;
	    break;
	case OC_Primitive:
	    newmode = OCPRIMITIVEMODE;
	    break;
	default:
	    break;
	}
	break;

    case OCPRIMITIVEMODE:
    case OCNULLMODE:
    case OCEMPTYMODE:
    default:
	newmode = OCNULLMODE;
	break;
    }
    if(newmode == OCNULLMODE)
        OCPANIC1("No defined mode transition: %d",(int)srcmode);
    return newmode;
}

/* get the presumed current sequence tag */
static int
ocgetsequencetag(XXDR* xdrs)
{
    char tag[XDRUNIT];
    if(!xxdr_getbytes(xdrs,tag,sizeof(tag))) return 0;
    return tag[0];
}

static int
ocskipcounts(XXDR* xdrs, OCnode* node, off_t expected)
{
    if(node->array.rank == 0) return 1; /* simple scalar */
#ifdef VERIFY
    unsigned int xdrcount0,xdrcount1;
    /* Collect the dimension count from the xdr data packet*/
    if(!xxdr_uint(xdrs,&xdrcount0)) OCGOTO(shortxdr);
    if(expected >= 0 && xdrcount0 != expected) return 0;
    /* pull out redundant second count*/
    /* (note that String/URL do not have redundant count)*/
    if(node->octype == OC_Primitive
       && node->etype != OC_String && node->etype != OC_URL) {
	if(!xxdr_uint(xdrs,&xdrcount1)) return 0;
	if(xdrcount0 != xdrcount1) return 0;
    }
#else
    /* skip the counts */
    expected = expected; /*shut up compiler*/
    if(node->octype == OC_Primitive
       && node->etype != OC_String && node->etype != OC_URL) {
        if(!xxdr_skip(xdrs,2*XDRUNIT)) return 0;
    } else {
	if(!xxdr_skip(xdrs,XDRUNIT)) return 0;
    }
#endif
    return 1;
}

/**************************************************/
/* Moved ocdata.c here */
/**************************************************/

const char StartOfSequence = '\x5A';
const char EndOfSequence = '\xA5';

static int ocerrorstring(XXDR* xdrs);

#define LOCALMEMMAX 1024

/*
Skip arbitrary object based on its octype
and a state 

Cases:

octype		Skip State	actions
-------------------------------------------
Structure
  |Grid
  |DataSet	SKIPINSTANCE	Skip single instance
		  |SKIPFIELDS
		SKIPWHOLE	Skip array of instances
				including leading counts

Sequence	SKIPFIELDS	Skip single record
				(assume leading tag already skipped)
		SKIPINSTANCE	Skip single record
				(including leading tag)
						    
		SKIPWHOLE	Skip all records
				including leading tags
				and trailing end marker

Primitive	<any>		Skip whole primitive array
				including leading counts

Notes:
1. unlisted combinations are not legal/possible.
2. assume that xxdr_getpos is properly positioned.
3. If octype is OC_Sequence, tagp will be set with the
   last tag encountered.
*/

OCerror
ocskipinstance(OCnode* node, XXDR* xdrs, int state, int* tagp)
{
    int i,tag;
    int stat = OC_NOERR;

/* Support switch on combination of octype X state to simply code */
#define CASE(octype,state) ((octype)<<3 | state)

    switch (CASE(node->octype,state)) {

    case CASE(OC_Dataset,SKIPINSTANCE):
    case CASE(OC_Grid,SKIPINSTANCE):
    case CASE(OC_Structure,SKIPINSTANCE):

    case CASE(OC_Dataset,SKIPFIELDS):
    case CASE(OC_Grid,SKIPFIELDS):
    case CASE(OC_Structure,SKIPFIELDS):
    case CASE(OC_Sequence,SKIPFIELDS): /* NOTE this special case */
	if(node->skip.instancesize != OCINDETERMINATE) {
	    if(!xxdr_skip(xdrs,node->skip.instancesize)) OCGOTO(shortxdr);
	} else {/* skip field by field */
	    for(i=0;i<oclistlength(node->subnodes);i++) {
		OCnode* field = (OCnode*)oclistget(node->subnodes,i);
		stat = ocskipinstance(field, xdrs, SKIPWHOLE,NULL);
		if(stat != OC_NOERR) {OCTHROWCHK(stat); goto done;}
	    }
	}
	break;

    case CASE(OC_Dataset,SKIPWHOLE):
    case CASE(OC_Grid,SKIPWHOLE):
    case CASE(OC_Structure,SKIPWHOLE):
	OCASSERT(node->skip.count != OCINDETERMINATE);
	if(node->skip.totalsize != OCINDETERMINATE) {
	    if(!xxdr_skip(xdrs,node->skip.totalsize)) goto badxdr;
	} else {/* skip each instance */
	    if(node->array.rank > 0) {
	        if(!ocskipcounts(xdrs,node,node->skip.count)) goto badxdr;
                for(i=0;i<node->skip.count;i++) {
                    stat = ocskipinstance(node, xdrs, SKIPFIELDS,NULL);
                    if(stat != OC_NOERR) {OCTHROWCHK(stat); goto done;}
                }
            } else { /* scalar */
                stat = ocskipinstance(node, xdrs, SKIPINSTANCE,NULL);
                if(stat != OC_NOERR) {OCTHROWCHK(stat); goto done;}
            }
        }
	break;

    case CASE(OC_Sequence,SKIPINSTANCE): /* Skip record including tag */
	tag = ocgetsequencetag(xdrs); /* always read the tag */
	if(tagp) *tagp = tag;
	if(tag == StartOfSequence) { /* skip record fields */
	    stat = ocskipinstance(node, xdrs, SKIPFIELDS,NULL);
	    if(stat != OC_NOERR) {OCTHROWCHK(stat); break;}
	} /* let caller handle */
	break;
	
    case CASE(OC_Sequence,SKIPWHOLE): /* Skip multiple records including tags */
	for(i=0;;i++) {
	    stat = ocskipinstance(node, xdrs, SKIPINSTANCE, &tag);
	    if(stat != OC_NOERR) {OCTHROWCHK(stat); break;}
	    if(tag == EndOfSequence) break; /* done */
	    if(tag != StartOfSequence) goto badxdr; /* malformed */
	}
	break;

    case CASE(OC_Primitive,SKIPWHOLE):
    case CASE(OC_Primitive,SKIPINSTANCE):
    case CASE(OC_Primitive,SKIPFIELDS):
	OCASSERT(node->skip.count != OCINDETERMINATE);
	if(node->skip.totalsize != OCINDETERMINATE) {
	    /* skip directly past it */
	    if(!xxdr_skip(xdrs,node->skip.totalsize)) goto badxdr;
	} else {/* Walk instance by instance */
	    if(state == SKIPWHOLE) {
		/* read the counts */
		if(!ocskipcounts(xdrs,node,node->skip.count))
		    goto badxdr;
	    }
	    OCASSERT(node->etype == OC_String || node->etype == OC_URL);
	    /* get the count */
	    for(i=0;i<node->skip.count;i++) {
		/* read and skip the string */
	        unsigned int len;
		/* read string size */
		if(!xxdr_uint(xdrs,&len)) OCGOTO(shortxdr);
		/* round up to next XDRUNIT and skip string contents */
		len = RNDUP(len);
		if(!xxdr_skip(xdrs,(size_t)len)) OCGOTO(shortxdr);
            }
	}
	break;

        default:
	    OCPANIC2("ocskipinstance: encountered unexpected node type or state: %d,%d",
			node->octype,state);
	    break;
    }
done:
    return OCTHROW(stat);
shortxdr:
    oc_log(LOGERR,"short xdr packet");
    stat = OC_EXDR;
    goto done;
badxdr:
    oc_log(LOGERR,"malformed xdr packet");
    stat = OC_EXDR;
    goto done;
}

/*
Extract data from the xdr packet into a chunk of memory.
Normally, it is assumed that we are (at least virtually)
"at" a single instance in the xdr packet; which we read.
Virtually because for packed data, we need to point to
the beginning of the packed data and use the index to indicate
which packed element to get. Assume that in any case,
any leading counts have been passed.
*/
OCerror
ocxdrread(OCcontent* content, XXDR* xdrs, char* memory, size_t memsize,
          ocindex_t start, ocindex_t count)
{
    int stat = OC_NOERR;
    unsigned int i;
    size_t elemsize;
    size_t readsize;
    size_t skipsize;
    char localmem[LOCALMEMMAX];
    char* srcmem;    
    unsigned int* p;
    int packed;
    int scalar;
    OCtype octype,etype;
    ocindex_t localstart = start; /* will change if node is cacheing */
    OCnode* node;

    node = content->node;
    octype = node->octype;
    etype = node->etype;

    elemsize = octypesize(etype);

    scalar = (node->array.rank == 0 ? 1 : 0);

    /* check if the data is packed*/
    packed = (octype == OC_Primitive && !scalar
              && (etype == OC_Byte || etype == OC_UByte || etype == OC_Char));
	 
    /* validate memory space*/
    if(memsize < elemsize*count) return OCTHROW(OC_EINVAL);

#ifdef OCIGNORE
    if(!scalar && (!node->cache.cacheable || !node->cache.valid)) {
        unsigned int xdrcount0,xdrcount1;
	/* assume xdr position is correct */
        /* Read leading double count if ! scalar*/
        if(!xxdr_uint(xdrs,&xdrcount0)) OCGOTO(shortxdr);
        if(!xxdr_uint(xdrs,&xdrcount1)) OCGOTO(shortxdr);
        if(xdrcount0 != xdrcount1) return OCTHROW(OC_EXDR);
        if(xdrcount0 < start+count) OCGOTO(shortxdr);
    }
#endif
 
    /* Handle packed data specially*/
    if(packed) {
	readsize = count*1; /* |OC_(Char,Byte,UByte)| == 1 */
	skipsize = start*1; /* |OC_(Char,Byte,UByte)| == 1 */
	/* skip to start of what we want to read */
	if(!xxdr_skip(xdrs,skipsize)) OCGOTO(shortxdr);
	/* read data, keeping xdrs on XDRUNIT boundary */
	if(!xxdr_opaque(xdrs,memory,readsize))
	    OCGOTO(shortxdr);
	return OCTHROW(OC_NOERR);
    }

    /* Not packed */

#ifdef OCIGNORE
    /* If this (primitive) object is cacheable and is valid cache,
       then modify start and set the xdr position accordingly
    */
    if(node->cache.cacheable && node->cache.valid) {
	if(node->cache.index <= start) {
	    localstart -= node->cache.index;
	    if(!xxdr_setpos(xdrs,node->cache.offset)) return xdrerror();
	}
    }
#endif

    /* Compute how much to skip based on the content's cache index */
    localstart = start - content->cache.index;
    if(localstart < 0) localstart = 0;

    /* extract count items; use xxdr_getbytes to speed up*/
    srcmem = memory;
    switch (etype) {
    case OC_Float64: case OC_Int64: case OC_UInt64:
	readsize = count*2*XDRUNIT;
	skipsize = localstart*2*XDRUNIT;
	/* skip to start of what we want to read */
	if(!xxdr_skip(xdrs,skipsize)) OCGOTO(shortxdr);
	if(!xxdr_opaque(xdrs,(char*)srcmem,readsize)) OCGOTO(shortxdr);
	if(etype == OC_Float64) {
	    double* dp;
	    for(dp=(double*)srcmem,i=0;i<count;i++,dp++) {
		double swap;
		xxdrntohdouble((char*)dp,&swap);
		*dp = swap;
	    }
	} else if(!xxdr_network_order) {
	    unsigned long long* llp;
	    for(llp=(unsigned long long*)srcmem,i=0;i<count;i++,p++) {
		swapinline64(llp);
	    }
	}
	break;

    case OC_String: case OC_URL: {
	/* Read string by string */
        char* s = NULL;
	char** pmem = (char**)srcmem;
	/* First skip to the starting string */
	for(i=0;i<localstart;i++) {
	    unsigned int slen;
            if(!xxdr_uint(xdrs,&slen)) OCGOTO(shortxdr);
	    slen = RNDUP(slen);
            if(!xxdr_skip(xdrs,slen)) OCGOTO(shortxdr);
        }
	/* Read count strings */
	for(i=0;i<count;i++) {
	    off_t slen;
	    /* xxdr_string will always alloc the space */	
            if(!xxdr_string(xdrs,&s,&slen)) 
		OCGOTO(shortxdr);
	    pmem[i] = s;
	}
    } break;


    case OC_Char: case OC_Byte: case OC_UByte:
    case OC_Int16: case OC_UInt16:
	/* We need to store the xdr data locally until we can convert it out 
           because  elemsize < sizeof(int) */
	srcmem = localmem;
	if(count*elemsize > sizeof(localmem)) {
	    srcmem = (char*)ocmalloc(count*sizeof(unsigned int));
	    if(srcmem == NULL) {stat = OCTHROW(OC_ENOMEM); goto done;}
	}
	/* fall thru */		
    case OC_Int32: case OC_UInt32:
    case OC_Float32:
        readsize = (count)*XDRUNIT;
        skipsize = (localstart)*XDRUNIT;
	if(!xxdr_skip(xdrs,skipsize)) OCGOTO(shortxdr);
	if(!xxdr_opaque(xdrs,(char*)srcmem,readsize)) OCGOTO(shortxdr);
	if(!xxdr_network_order) {
	    for(p=(unsigned int*)srcmem,i=0;i<count;i++,p++) {
		swapinline32(p);
	    }
	}
	break;

    default: OCPANIC("unexpected etype"); break;
    }

    /* Convert memory to right format */
    switch (etype) {

    case OC_Char: case OC_Byte: case OC_UByte: {
	char* pmem = (char*)memory;
	p = (unsigned int*)srcmem;
	for(i=0;i<count;i++) {
	    unsigned int tmp = *p++;
	    *pmem++ = (unsigned char)tmp;
	}
    } break;

    case OC_Int16: case OC_UInt16: {
	unsigned short* pmem = (unsigned short*)memory;
	p = (unsigned int*)srcmem;
	for(i=0;i<count;i++) {
	    unsigned int tmp = *p++;
	    *pmem++ = (unsigned short)tmp;
	}
    } break;

    default: 
	break; /* already handled above */
    }

    /* set cache */
    content->cache.index = start + count; /* should be our current index */
    content->cache.offset = xxdr_getpos(xdrs); /* should be our current position */

done:
    return OCTHROW(stat);

shortxdr:
    content->cache.valid = 0; /* no longer valid */
    if(!ocerrorstring(xdrs))
        oc_log(LOGERR,"DAP DATADDS packet is apparently too short");
    stat = OCTHROW(OC_EDATADDS);
    goto done;    
}

int
occountrecords(OCnode* node, XXDR* xdrs, size_t* nrecordsp)
{
    int stat = OC_NOERR;
    size_t nrecords = 0;

    if(node->octype != OC_Sequence) return OCTHROW(OC_EINVAL);
    /* checkpoint the xdr position*/
    for(nrecords=0;;nrecords++) {
	int tag = 0;
	stat = ocskipinstance(node,xdrs,SKIPINSTANCE,&tag);
        if(stat != OC_NOERR) break;
        if(tag == EndOfSequence) break;
        if(tag != StartOfSequence) {
            oc_log(LOGERR,"missing/invalid begin/end record marker\n");
            stat = OC_EINVALCOORDS;
            break;
        }
        if(stat != OC_NOERR) break;
    }
    if(nrecordsp != NULL) *nrecordsp = nrecords;
    return OCTHROW(stat);
}



#define tag "Error {\n"

static int
ocerrorstring(XXDR* xdrs)
{
    /* Check to see if the xdrs contains "Error {\n'; assume it is at the beginning of data */
    off_t avail = xxdr_getavail(xdrs);
    char* data = (char*)malloc(avail);
    if(!xxdr_setpos(xdrs,0)) return 0;
    if(!xxdr_opaque(xdrs,data,avail)) return 0;
    /* check for error tag at front */
    if(ocstrncmp(data,tag,sizeof(tag))==0) {
	char* p;
        if((p=strchr(data,'}')) != NULL) *(++p)='\0';
        oc_log(LOGERR,"Server error: %s",data);
        /* Since important, report to stderr as well */
        fprintf(stderr,"Server error: %s",data);
	return 1;
    }
    return 0;
}

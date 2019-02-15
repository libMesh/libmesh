/*********************************************************************
 *   Copyright 2009, UCAR/Unidata
 *   See netcdf/COPYRIGHT file for copying and redistribution conditions.
 *********************************************************************/
/* $Id: data.c,v 1.7 2010/05/24 19:59:56 dmh Exp $ */
/* $Header: /upc/share/CVS/netcdf-3/ncgen/data.c,v 1.7 2010/05/24 19:59:56 dmh Exp $ */

#include        "includes.h"
#include        "ncoffsets.h"
#include        "netcdf_aux.h"
#include        "dump.h"

#undef VERIFY

#define XVSNPRINTF vsnprintf
/*
#define XVSNPRINTF lvsnprintf
extern int lvsnprintf(char*, size_t, const char*, va_list);
*/

#define DATALISTINIT 32

/* Track all known datalist*/
List* alldatalists = NULL;

NCConstant nullconstant;
NCConstant fillconstant;

Datalist nildatalist; /* to support NIL keyword */

Bytebuffer* codebuffer;
Bytebuffer* codetmp;
Bytebuffer* stmt;


/* Forward */
static void setconstlist(NCConstant* con, Datalist* dl);

#ifdef VERIFY
/* index of match */
static int
verify(List* all, Datalist* dl)
{
    int i;
    for(i=0;i<listlength(all);i++) {
	void* pi = listget(all,i);
	if(pi == dl)
	    return i;
    }
    return -1;
}
#endif

/**************************************************/
/**************************************************/

NCConstant*
nullconst(void)
{
    NCConstant* n = ecalloc(sizeof(NCConstant));
    return n;
}

int
isstringable(nc_type nctype)
{
    switch (nctype) {
    case NC_CHAR: case NC_STRING:
    case NC_BYTE: case NC_UBYTE:
    case NC_FILLVALUE:
	return 1;
    default: break;
    }
    return 0;
}

NCConstant*
list2const(Datalist* list)
{
    NCConstant* con = nullconst();
    ASSERT(list != NULL);
    con->nctype = NC_COMPOUND;
    con->lineno = list->data[0]->lineno;
    setconstlist(con,list);
    con->filled = 0;
    return con;
}

Datalist*
const2list(NCConstant* con)
{
    Datalist* list;
    ASSERT(con != NULL);
    list = builddatalist(1);
    if(list != NULL) {
        dlappend(list,con);
    }
    return list;
}

/**************************************************/
#ifdef GENDEBUG
void
report(char* lead, Datalist* list)
{
extern void bufdump(Datalist*,Bytebuffer*);
Bytebuffer* buf = bbNew();
bufdump(list,buf);
fprintf(stderr,"\n%s::%s\n",lead,bbContents(buf));
fflush(stderr);
bbFree(buf);
}

#endif

/**************************************************/

static void
setconstlist(NCConstant* con, Datalist* dl)
{
#ifdef VERIFY
    int pos = verify(alldatalists,dl);
    if(pos >= 0) {
        dumpdatalist(listget(alldatalists,pos),"XXX");
    }
#endif
    con->value.compoundv = dl;
}


/* Deep constant cloning; return struct not pointer to struct*/
NCConstant*
cloneconstant(NCConstant* con)
{
    NCConstant* newcon = NULL;
    Datalist* newdl = NULL;
    char* s = NULL;

    newcon = nullconst();
    if(newcon == NULL) return newcon;
    *newcon = *con;
    switch (newcon->nctype) {
    case NC_STRING:
	if(newcon->value.stringv.len == 0)
	    s = NULL;
	else {
	    s = (char*)ecalloc(newcon->value.stringv.len+1);
	    if(newcon->value.stringv.len > 0)
	        memcpy(s,newcon->value.stringv.stringv,newcon->value.stringv.len);
	    s[newcon->value.stringv.len] = '\0';
	}
	newcon->value.stringv.stringv = s;
	break;
    case NC_OPAQUE:
	s = (char*)ecalloc(newcon->value.opaquev.len+1);
	if(newcon->value.opaquev.len > 0)
	    memcpy(s,newcon->value.opaquev.stringv,newcon->value.opaquev.len);
	s[newcon->value.opaquev.len] = '\0';
	newcon->value.opaquev.stringv = s;
	break;
    case NC_COMPOUND:
	newdl = clonedatalist(con->value.compoundv);
	setconstlist(newcon,newdl);
	break;
    default: break;
    }
    return newcon;
}

/* Deep constant clear*/
void
clearconstant(NCConstant* con)
{
    if(con == NULL) return;
    switch (con->nctype) {
    case NC_STRING:
	if(con->value.stringv.stringv != NULL)
	    efree(con->value.stringv.stringv);
	break;
    case NC_OPAQUE:
	if(con->value.opaquev.stringv != NULL)
	    efree(con->value.opaquev.stringv);
	break;
    case NC_COMPOUND:
	con->value.compoundv = NULL;
	break;
    default: break;
    }
    memset((void*)con,0,sizeof(NCConstant));
}

void
freeconstant(NCConstant* con, int shallow)
{
    if(!shallow) clearconstant(con);
    nullfree(con);
}

/**************************************************/

#if 0
Datalist*
datalistclone(Datalist* dl)
{
    int i;
    Datalist* clone = builddatalist(dl->length);
    for(i=0;i<dl->length;i++) {
	clone->data[i] = cloneconstant(dl->data[i]);
    }
    return clone;
}

Datalist*
datalistappend(Datalist* dl, NCConstant* con)
{
    NCConstant** vector;
    ASSERT(dl != NULL);
    if(con == NULL) return dl;
    vector = (NCConstant**)erealloc(dl->data,sizeof(NCConstant*)*(dl->length+1));
    if(vector == NULL) return NULL;
    vector[dl->length] = cloneconstant(con);
    dl->length++;
    dl->data = vector;
    return dl;
}

Datalist*
datalistreplace(Datalist* dl, unsigned int index, NCConstant* con)
{
    ASSERT(dl != NULL);
    ASSERT(index < dl->length);
    ASSERT(con != NULL);
    dl->data[index] = cloneconstant(con);
    return dl;
}
#endif

int
datalistline(Datalist* ds)
{
    if(ds == NULL || ds->length == 0) return 0;
    return ds->data[0]->lineno;
}


/* Go thru a databuf of possibly nested constants
   and insert commas as needed; ideally, this
   operation should be idempotent so that
   the caller need not worry about it having already
   been applied. Also, handle situation where there may be missing
   matching right braces.
*/

static char* commifyr(char* p, Bytebuffer* buf);
static char* wordstring(char* p, Bytebuffer* buf, int quote);

void
commify(Bytebuffer* buf)
{
    char* list,*p;

    if(bbLength(buf) == 0) return;
    list = bbDup(buf);
    p = list;
    bbClear(buf);
    commifyr(p,buf);
    bbNull(buf);
    efree(list);
}

/* Requires that the string be balanced
   WRT to braces
*/
static char*
commifyr(char* p, Bytebuffer* buf)
{
    int comma = 0;
    int c;
    while((c=*p++)) {
	if(c == ' ') continue;
	if(c == ',') continue;
	else if(c == '}') {
	    break;
	}
	if(comma) bbCat(buf,", "); else comma=1;
	if(c == '{') {
	    bbAppend(buf,'{');
	    p = commifyr(p,buf);
	    bbAppend(buf,'}');
	} else if(c == '\'' || c == '\"') {
	    p = wordstring(p,buf,c);
	} else {
	    bbAppend(buf,c);
	    p=word(p,buf);
	}
    }    
    return p;
}

char*
word(char* p, Bytebuffer* buf)
{
    int c;
    while((c=*p++)) {
	if(c == '}' || c == ' ' || c == ',') break;
	if(c == '\\') {
	    bbAppend(buf,c);
	    c=*p++;
	    if(!c) break;
	}
	bbAppend(buf,(char)c);	
    }	
    p--; /* leave terminator for parent */
    return p;
}

static char*
wordstring(char* p, Bytebuffer* buf, int quote)
{
    int c;
    bbAppend(buf,quote);
    while((c=*p++)) {	    
	if(c == '\\') {
	    bbAppend(buf,c);
	    c = *p++;
	    if(c == '\0') return --p;
	} else if(c == quote) {
	    bbAppend(buf,c);
	    return p;
	}
	bbAppend(buf,c);
    }
    return p;
}


static const char zeros[] =
    "\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0";

void
alignbuffer(NCConstant* prim, Bytebuffer* buf)
{
    int alignment,pad,offset;

    ASSERT(prim->nctype != NC_COMPOUND);

    if(prim->nctype == NC_ECONST)
        alignment = ncaux_class_alignment(prim->value.enumv->typ.typecode);
    else if(usingclassic && prim->nctype == NC_STRING)
        alignment = ncaux_class_alignment(NC_CHAR);
    else if(prim->nctype == NC_CHAR)
        alignment = ncaux_class_alignment(NC_CHAR);
    else
        alignment = ncaux_class_alignment(prim->nctype);
    offset = bbLength(buf);
    pad = getpadding(offset,alignment);
    if(pad > 0) {
	bbAppendn(buf,(void*)zeros,pad);
    }
}



/*
Following routines are in support of language-oriented output
*/

void
codedump(Bytebuffer* buf)
{
   bbCatbuf(codebuffer,buf);
   bbClear(buf);
}

void
codepartial(const char* txt)
{
    bbCat(codebuffer,txt);
}

void
codeline(const char* line)
{
    codepartial(line);
    codepartial("\n");
}

void
codelined(int n, const char* txt)
{
    bbindent(codebuffer,n);
    bbCat(codebuffer,txt);
    codepartial("\n");
}

void
codeflush(void)
{
    if(bbLength(codebuffer) > 0) {
        bbNull(codebuffer);
        fputs(bbContents(codebuffer),stdout);
        fflush(stdout);
        bbClear(codebuffer);
    }
}

void
bbindent(Bytebuffer* buf, const int n)
{
    bbCat(buf,indented(n));
}

/* Provide an restrict snprintf that writes to an expandable buffer */
/* Simulates a simple snprintf because apparently
   the IRIX one is broken wrt return value.
   Supports only %u %d %f %s and %% specifiers
   with optional leading hh or ll.
*/

static void
vbbprintf(Bytebuffer* buf, const char* fmt, va_list argv)
{
    char tmp[128];
    const char* p;
    int c;
    int hcount;
    int lcount;

    char* text;

    for(p=fmt;(c=*p++);) {
	hcount = 0; lcount = 0;
	switch (c) {
	case '%':
retry:	    switch ((c=*p++)) {
	    case '\0': bbAppend(buf,'%'); p--; break;
	    case '%': bbAppend(buf,c); break;
	    case 'h':
		hcount++;
		while((c=*p) && (c == 'h')) {hcount++; p++;}
		if(hcount > 2) hcount = 2;
		goto retry;	        
	    case 'l':
		lcount++;
		while((c=*p) && (c == 'l')) {
		    lcount++;
		    p++;
		}
		if(lcount > 2) lcount = 2;
		goto retry;	        
	    case 'u':
		if(hcount == 2) {
   	            snprintf(tmp,sizeof(tmp),"%hhu",
			(unsigned char)va_arg(argv,unsigned int));
		} else if(hcount == 1) {
   	            snprintf(tmp,sizeof(tmp),"%hu",
			(unsigned short)va_arg(argv,unsigned int));
		} else if(lcount == 2) {
   	            snprintf(tmp,sizeof(tmp),"%llu",
			(unsigned long long)va_arg(argv,unsigned long long));
		} else if(lcount == 1) {
   	            snprintf(tmp,sizeof(tmp),"%lu",
			(unsigned long)va_arg(argv,unsigned long));
		} else {
   	            snprintf(tmp,sizeof(tmp),"%u",
			(unsigned int)va_arg(argv,unsigned int));
		}
		bbCat(buf,tmp);
		break;
	    case 'd':
		if(hcount == 2) {
   	            snprintf(tmp,sizeof(tmp),"%hhd",
			(signed char)va_arg(argv,signed int));
		} else if(hcount == 1) {
   	            snprintf(tmp,sizeof(tmp),"%hd",
			(signed short)va_arg(argv,signed int));
		} else if(lcount == 2) {
   	            snprintf(tmp,sizeof(tmp),"%lld",
			(signed long long)va_arg(argv,signed long long));
		} else if(lcount == 1) {
   	            snprintf(tmp,sizeof(tmp),"%ld",
			(signed long)va_arg(argv,signed long));
		} else {
   	            snprintf(tmp,sizeof(tmp),"%d",
			(signed int)va_arg(argv,signed int));
		}
		bbCat(buf,tmp);
		break;
            case 'f':
		if(lcount > 0) {
   	            snprintf(tmp,sizeof(tmp),"((double)%.16g)",
			(double)va_arg(argv,double));
		} else {
   	            snprintf(tmp,sizeof(tmp),"((float)%.8g)",
			(double)va_arg(argv,double));
		}
		bbCat(buf,tmp);
	        break;
	    case 's':
		text = va_arg(argv,char*);
		bbCat(buf,text);
		break;		
	    case 'c':
		c = va_arg(argv,int);
		bbAppend(buf,(char)c);
		break;		
            default:
		PANIC1("vbbprintf: unknown specifier: %c",(char)c);
	    }
	    break;
	default: 
	    bbAppend(buf,c);
	}
    }
}

void
bbprintf(Bytebuffer* buf, const char *fmt, ...)
{
    va_list argv;
    va_start(argv,fmt);
    vbbprintf(buf,fmt,argv);
    va_end(argv);
}

void
bbprintf0(Bytebuffer* buf, const char *fmt, ...)
{
    va_list argv;
    va_start(argv,fmt);
    bbClear(buf);
    vbbprintf(buf,fmt,argv);
    va_end(argv);
}

void
codeprintf(const char *fmt, ...)
{
    va_list argv;
    va_start(argv,fmt);
    vbbprintf(codebuffer,fmt,argv);
    va_end(argv);
}

NCConstant*
emptycompoundconst(int lineno)
{
    NCConstant* c = nullconst();
    c->lineno = lineno;
    c->nctype = NC_COMPOUND;
    setconstlist(c,builddatalist(0));
    c->filled = 0;
    return c;    
}

/* Make an empty string constant*/
NCConstant*
emptystringconst(int lineno)
{
    NCConstant* c = nullconst();
    ASSERT(c != NULL);
    c->lineno = lineno;
    c->nctype = NC_STRING;
    c->value.stringv.len = 0;
    c->value.stringv.stringv = NULL;
    c->filled = 0;
    return c;    
}

#define INDENTMAX 256
static char* dent = NULL;

char*
indented(int n)
{
    char* indentation;
    if(dent == NULL) {
	dent = (char*)ecalloc(INDENTMAX+1);
	memset((void*)dent,' ',INDENTMAX);
	dent[INDENTMAX] = '\0';	
    }
    if(n*4 >= INDENTMAX) n = INDENTMAX/4;
    indentation = dent+(INDENTMAX - 4*n);
    return indentation;
}

void
dlextend(Datalist* dl)
{
    size_t newalloc;
    NCConstant** newdata = NULL;
    newalloc = (dl->alloc > 0?2*dl->alloc:2);
    newdata = (NCConstant**)ecalloc(newalloc*sizeof(NCConstant*));
    if(dl->length > 0)
        memcpy(newdata,dl->data,sizeof(NCConstant*)*dl->length);
    dl->alloc = newalloc;
    nullfree(dl->data);
    dl->data = newdata;
}


void
capture(Datalist* dl)
{
    if(alldatalists == NULL) alldatalists = listnew();
    listpush(alldatalists,dl);
}

Datalist*
builddatalist(int initial)
{
    Datalist* ci;
    if(initial <= 0) initial = DATALISTINIT;
    initial++; /* for header*/
    ci = (Datalist*)ecalloc(sizeof(Datalist));
    if(ci == NULL) semerror(0,"out of memory\n");
    ci->data = (NCConstant**)ecalloc(sizeof(NCConstant*)*initial);
    ci->alloc = initial;
    ci->length = 0;
    return ci;
}

void
dlappend(Datalist* dl, NCConstant* constant)
{
    if(dl->length >= dl->alloc)
	dlextend(dl);
    dl->data[dl->length++] = (constant);
}

void
dlset(Datalist* dl, size_t pos, NCConstant* constant)
{
    ASSERT(pos < dl->length);
    dl->data[pos] = (constant);
}

/* Convert a datalist to a compound constant */
NCConstant*
builddatasublist(Datalist* dl)
{

  NCConstant* d = nullconst();
  d->nctype = NC_COMPOUND;
  d->lineno = (dl->length > 0?dl->data[0]->lineno:0);
  setconstlist(d,dl);
  d->filled = 0;
  return d;

}

/* Deep copy */
Datalist*
clonedatalist(Datalist* dl)
{
    int i;
    size_t len;
    Datalist* newdl;

    if(dl == NULL) return NULL;
    len = datalistlen(dl);
    newdl = builddatalist(len);
    /* initialize */
    for(i=0;i<len;i++) {
	NCConstant* con = datalistith(dl,i);
	con = cloneconstant(con);
	dlappend(newdl,con);
    }
#if 0
    newdl->vlen = dl->vlen;
#endif
    newdl->readonly = dl->readonly;
    return newdl;
}


/* recursive helpers */

#if 0
static int
isdup(Datalist* dl)
{
    int i;
    size_t limit = listlength(alldatalists);
    for(i=0;i<limit;i++) {
	Datalist* di = listget(alldatalists,i);
	if(di == dl) return 1;
    }
    return 0;
}
#endif

void
reclaimconstant(NCConstant* con)
{
    if(con == NULL) return;
    switch (con->nctype) {
    case NC_STRING:
	if(con->value.stringv.stringv != NULL)
	    efree(con->value.stringv.stringv);
	break;
    case NC_OPAQUE:
	if(con->value.opaquev.stringv != NULL)
	    efree(con->value.opaquev.stringv);
	break;
    case NC_COMPOUND:
#ifdef VERIFY
	{int pos;
	if((pos=verify(alldatalists,con->value.compoundv)) >= 0) {
	    dumpdatalist(listget(alldatalists,pos),"XXX");
	    abort();
	}
	}
#endif
        reclaimdatalist(con->value.compoundv);
	con->value.compoundv = NULL;
	break;
    default: break;
    }
    efree(con);

}

void
reclaimdatalist(Datalist* list)
{
   int i;
   if(list == NULL) return;
   if(list->data != NULL) {
       for(i=0;i<list->length;i++) {
	    NCConstant* con = list->data[i];
	    if(con != NULL) reclaimconstant(con);
       }
       efree(list->data);
       list->data = NULL;
   }
   efree(list);
}

void
reclaimalldatalists(void)
{
    int i;
#if 0
    int j;
    /* Remove duplicates */
    for(i=0;i<listlength(alldatalists);i++) {
	Datalist* di = listget(alldatalists,i);
	if(di == NULL) continue;
        for(j=i;j<listlength(alldatalists);j++) {
	    Datalist* dj = listget(alldatalists,j);
	    if(dj == di) {
	        listset(alldatalists,j,NULL);
fprintf(stderr,"XXX\n");
	    }
        }
    }
#endif
    for(i=0;i<listlength(alldatalists);i++) {
        Datalist* di = listget(alldatalists,i);
	if(di != NULL)
	    reclaimdatalist(di);
    }
    efree(alldatalists);
    alldatalists = NULL;    
}

/* Obsolete */
#if 0
/* return 1 if the next element in the datasrc is compound*/
int
issublist(Datasrc* datasrc) {return istype(datasrc,NC_COMPOUND);}

/* return 1 if the next element in the datasrc is a string*/
int
isstring(Datasrc* datasrc) {return istype(datasrc,NC_STRING);}

/* return 1 if the next element in the datasrc is a fill value*/
int
isfillvalue(Datasrc* datasrc)
{
return srcpeek(datasrc) == NULL || istype(datasrc,NC_FILLVALUE);
}

/* return 1 if the next element in the datasrc is nc_type*/
int
istype(Datasrc* datasrc , nc_type nctype)
{
    NCConstant* ci = srcpeek(datasrc);
    if(ci != NULL && ci->nctype == nctype) return 1;
    return 0;
}

/**************************************************/

void
freedatasrc(Datasrc* src)
{
    efree(src);
}

Datasrc*
allocdatasrc(void)
{
    Datasrc* src;
    src = ecalloc(sizeof(Datasrc));
    src->data = NULL;
    src->index = 0;
    src->length = 0;
    src->prev = NULL;
    return src;
}

Datasrc*
datalist2src(Datalist* list)
{
    Datasrc* src;
    ASSERT(list != NULL);
    src = allocdatasrc();
    src->data = list->data;
    src->index = 0;
    src->length = list->length;
    DUMPSRC(src,"#");
    return src;
}

Datasrc*
const2src(NCConstant* con)
{
    Datasrc* src;
    ASSERT(con != NULL);
    src = allocdatasrc();
    src->data = emalloc(sizeof(NCConstant*));
    src->data[0] = con;
    src->index = 0;
    src->length = 1;
    DUMPSRC(src,"#");
    return src;
}

NCConstant*
srcpeek(Datasrc* ds)
{
    if(ds == NULL) return NULL;
    if(ds->index < ds->length)
	return ds->data[ds->index];
    if(ds->spliced)
	return srcpeek(ds->prev);
    return NULL;
}

void
srcreset(Datasrc* ds)
{
    ds->index = 0;
}

NCConstant*
srcnext(Datasrc* ds)
{
    DUMPSRC(ds,"!");
    if(ds == NULL) return NULL;
    if(ds->index < ds->length)
	return ds->data[ds->index++];
    if(ds->spliced) {
	srcpop(ds);
	return srcnext(ds);
    }
    return NULL;
}

int
srcmore(Datasrc* ds)
{
    if(ds == NULL) return 0;
    if(ds->index < ds->length) return 1;
    if(ds->spliced) return srcmore(ds->prev);
    return 0;
}

int
srcline(Datasrc* ds)
{
    int index = ds->index;
    int len = ds->length;
    /* pick closest available entry*/
    if(len == 0) return 0;
    if(index >= len) index = len-1;
    return ds->data[index]->lineno;
}

void
srcpush(Datasrc* src)
{
    NCConstant* con;
    ASSERT(src != NULL);
    con = srcnext(src);
    ASSERT(con->nctype == NC_COMPOUND);
    srcpushlist(src,con->value.compoundv);
}

void
srcpushlist(Datasrc* src, Datalist* dl)
{
    Datasrc* newsrc;
    ASSERT(src != NULL && dl != NULL);
    newsrc = allocdatasrc();
    *newsrc = *src;
    src->prev = newsrc;
    src->index = 0;
    src->data = dl->data;
    src->length = dl->length;
    DUMPSRC(src,">!");
}

void
srcpop(Datasrc* src)
{
    if(src != NULL) {
        Datasrc* prev = src->prev;
	*src = *prev;
        freedatasrc(prev);
    }
    DUMPSRC(src,"<");
}

void
srcsplice(Datasrc* ds, Datalist* list)
{
    srcpushlist(ds,list);
    ds->spliced = 1;    
}

void
srcsetfill(Datasrc* ds, Datalist* list)
{
    if(ds->index >= ds->length) PANIC("srcsetfill: no space");
    if(ds->data[ds->index]->nctype != NC_FILLVALUE) PANIC("srcsetfill: not fill");
    ds->data[ds->index]->nctype = NC_COMPOUND;
    setconstlist(ds->data[ds->index],list);
}

#endif /*0*/

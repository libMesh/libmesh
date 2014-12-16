/*********************************************************************
 *   Copyright 1993, UCAR/Unidata
 *   See netcdf/COPYRIGHT file for copying and redistribution conditions.
 *********************************************************************/

#include "config.h"

#ifdef HAVE_SYS_TIME_H
#include <sys/time.h>
#endif

#include "oc.h"
extern int oc_dumpnode(OClink, OCddsnode);

#include "ncdap.h"
#include "dapalign.h"

#define LBRACKET '['
#define RBRACKET ']'


static char* repairname(const char* name, const char* badchars);

/**************************************************/
/**
 * Provide a hidden interface to allow utilities
 * to check if a given path name is really an ncdap3 url.
 * If no, return null, else return basename of the url
 * minus any extension.
 */

int
nc__testurl(const char* path, char** basenamep)
{
    NCURI* uri;
    int ok = ncuriparse(path,&uri);
    if(ok) {
	char* slash = (uri->file == NULL ? NULL : strrchr(uri->file, '/'));
	char* dot;
	if(slash == NULL) slash = (char*)path; else slash++;
        slash = nulldup(slash);
	dot = strrchr(slash, '.');
        if(dot != NULL &&  dot != slash) *dot = '\0';
	if(basenamep) *basenamep=slash ; else free(slash);
        ncurifree(uri);
    }
    return ok;
}

/**************************************************/

/*
Given a legal dap name with arbitrary characters,
convert to equivalent legal cdf name.
Currently, the only change is to convert '/'
names to %2f.
*/

char*
cdflegalname(char* name)
{
    return repairname(name,"/");
}

/* Define the type conversion of the DAP variables
   to the external netCDF variable type.
   The proper way is to, for example, convert unsigned short
   to an int to maintain the values.
   Unfortuneately, libnc-dap does not do this:
   it translates the types directly. For example
   libnc-dap upgrades the DAP byte type, which is unsigned char,
   to NC_BYTE, which signed char.
   Oh well.
   For netcdf-4, we can do proper type conversion.
*/
nc_type
nctypeconvert(NCDAPCOMMON* drno, nc_type nctype)
{
    nc_type upgrade = NC_NAT;
    if(drno->controls.flags & NCF_NC3) {
	/* libnc-dap mimic invariant is to maintain type size */
	switch (nctype) {
	case NC_CHAR:    upgrade = NC_CHAR; break;
	case NC_BYTE:    upgrade = NC_BYTE; break;
	case NC_UBYTE:   upgrade = NC_BYTE; break;
	case NC_SHORT:   upgrade = NC_SHORT; break;
	case NC_USHORT:  upgrade = NC_SHORT; break;
	case NC_INT:     upgrade = NC_INT; break;
	case NC_UINT:    upgrade = NC_INT; break;
	case NC_INT64:   upgrade = NC_INT64; break;
	case NC_UINT64:  upgrade = NC_UINT64; break;
	case NC_FLOAT:   upgrade = NC_FLOAT; break;
	case NC_DOUBLE:  upgrade = NC_DOUBLE; break;
	case NC_URL:
	case NC_STRING:  upgrade = NC_CHAR; break;
	default: break;
	}
    } else if(drno->controls.flags & NCF_NC4) {
	/* netcdf-4 conversion is more correct */
	switch (nctype) {
	case NC_CHAR:    upgrade = NC_CHAR; break;
	case NC_BYTE:    upgrade = NC_BYTE; break;
	case NC_UBYTE:   upgrade = NC_UBYTE; break;
	case NC_SHORT:   upgrade = NC_SHORT; break;
	case NC_USHORT:  upgrade = NC_USHORT; break;
	case NC_INT:     upgrade = NC_INT; break;
	case NC_UINT:    upgrade = NC_UINT; break;
	case NC_INT64:   upgrade = NC_INT64; break;
	case NC_UINT64:  upgrade = NC_UINT64; break;
	case NC_FLOAT:   upgrade = NC_FLOAT; break;
	case NC_DOUBLE:  upgrade = NC_DOUBLE; break;
	case NC_URL:
	case NC_STRING:  upgrade = NC_STRING; break;
	default: break;
	}
    }
    return upgrade;
}

nc_type
octypetonc(OCtype etype)
{
    switch (etype) {
    case OC_Char:	return NC_CHAR;
    case OC_Byte:	return NC_UBYTE;
    case OC_UByte:	return NC_UBYTE;
    case OC_Int16:	return NC_SHORT;
    case OC_UInt16:	return NC_USHORT;
    case OC_Int32:	return NC_INT;
    case OC_UInt32:	return NC_UINT;
    case OC_Int64:	return NC_INT64;
    case OC_UInt64:	return NC_UINT64;
    case OC_Float32:	return NC_FLOAT;
    case OC_Float64:	return NC_DOUBLE;
    case OC_String:	return NC_STRING;
    case OC_URL:	return NC_STRING;
    case OC_Dataset:	return NC_Dataset;
    case OC_Sequence:	return NC_Sequence;
    case OC_Structure:	return NC_Structure;
    case OC_Grid:	return NC_Grid;
    case OC_Dimension:	return NC_Dimension;
    case OC_Atomic:	return NC_Atomic;
    default: break;
    }
    return NC_NAT;
}

OCtype
nctypetodap(nc_type nctype)
{
    switch (nctype) {
    case NC_CHAR:	return OC_Char;
    case NC_BYTE:	return OC_Byte;
    case NC_UBYTE:	return OC_UByte;
    case NC_SHORT:	return OC_Int16;
    case NC_USHORT:	return OC_UInt16;
    case NC_INT:	return OC_Int32;
    case NC_UINT:	return OC_UInt32;
    case NC_INT64:	return OC_Int64;
    case NC_UINT64:	return OC_UInt64;
    case NC_FLOAT:	return OC_Float32;
    case NC_DOUBLE:	return OC_Float64;
    case NC_STRING:	return OC_String;
    default : break;
    }
    return OC_NAT;
}

size_t
nctypesizeof(nc_type nctype)
{
    switch (nctype) {
    case NC_CHAR:	return sizeof(char);
    case NC_BYTE:	return sizeof(signed char);
    case NC_UBYTE:	return sizeof(unsigned char);
    case NC_SHORT:	return sizeof(short);
    case NC_USHORT:	return sizeof(unsigned short);
    case NC_INT:	return sizeof(int);
    case NC_UINT:	return sizeof(unsigned int);
    case NC_INT64:	return sizeof(long long);
    case NC_UINT64:	return sizeof(unsigned long long);
    case NC_FLOAT:	return sizeof(float);
    case NC_DOUBLE:	return sizeof(double);
    case NC_STRING:	return sizeof(char*);
    default: PANIC("nctypesizeof");
    }
    return 0;
}

char*
nctypetostring(nc_type nctype)
{
    switch (nctype) {
    case NC_NAT:	return "NC_NAT";
    case NC_BYTE:	return "NC_BYTE";
    case NC_CHAR:	return "NC_CHAR";
    case NC_SHORT:	return "NC_SHORT";
    case NC_INT:	return "NC_INT";
    case NC_FLOAT:	return "NC_FLOAT";
    case NC_DOUBLE:	return "NC_DOUBLE";
    case NC_UBYTE:	return "NC_UBYTE";
    case NC_USHORT:	return "NC_USHORT";
    case NC_UINT:	return "NC_UINT";
    case NC_INT64:	return "NC_INT64";
    case NC_UINT64:	return "NC_UINT64";
    case NC_STRING:	return "NC_STRING";
    case NC_VLEN:	return "NC_VLEN";
    case NC_OPAQUE:	return "NC_OPAQUE";
    case NC_ENUM:	return "NC_ENUM";
    case NC_COMPOUND:	return "NC_COMPOUND";
    case NC_URL:	return "NC_URL";
    case NC_SET:	return "NC_SET";
    case NC_Dataset:	return "NC_Dataset";
    case NC_Sequence:	return "NC_Sequence";
    case NC_Structure:	return "NC_Structure";
    case NC_Grid:	return "NC_Grid";
    case NC_Dimension:	return "NC_Dimension";
    case NC_Atomic:	return "NC_Atomic";
    default: break;
    }
    return NULL;
}

/* Pad a buffer */
int
dapalignbuffer(NCbytes* buf, int alignment)
{
    int pad;
    unsigned long len;
    if(buf == NULL) return 0;
    len = ncbyteslength(buf);
    pad = nccpadding(len,alignment);

#ifdef TEST
    for(;pad > 0;pad--)
        ncbytesappend(buf,0x3a); /* 0x3a was chosen at random */
#else
    ncbytessetlength(buf,len+pad);
#endif
    return 1;
}

size_t
dapdimproduct(NClist* dimensions)
{
    size_t size = 1;
    unsigned int i;
    if(dimensions == NULL) return size;
    for(i=0;i<nclistlength(dimensions);i++) {
	CDFnode* dim = (CDFnode*)nclistget(dimensions,i);
	size *= dim->dim.declsize;
    }
    return size;
}


/* Return value of param or NULL if not found */
const char*
dapparamvalue(NCDAPCOMMON* nccomm, const char* key)
{
    const char* value;

    if(nccomm == NULL || key == NULL) return 0;
    if(!ncurilookup(nccomm->oc.url,key,&value))
	return NULL;
    return value;
}

static const char* checkseps = "+,:;";

/* Search for substring in value of param. If substring == NULL; then just
   check if param is defined.
*/
int
dapparamcheck(NCDAPCOMMON* nccomm, const char* key, const char* subkey)
{
    const char* value;
    char* p;

    if(nccomm == NULL || key == NULL) return 0;
    if(!ncurilookup(nccomm->oc.url,key,&value))
	return 0;
    if(subkey == NULL) return 1;
    p = strstr(value,subkey);
    if(p == NULL) return 0;
    p += strlen(subkey);
    if(*p != '\0' && strchr(checkseps,*p) == NULL) return 0;
    return 1;
}


/* This is NOT UNION */ 
int
nclistconcat(NClist* l1, NClist* l2)
{
    unsigned int i;
    for(i=0;i<nclistlength(l2);i++) nclistpush(l1,nclistget(l2,i));
    return 1;
}

int
nclistminus(NClist* l1, NClist* l2)
{
    unsigned int i,len,found;
    len = nclistlength(l2);
    found = 0;
    for(i=0;i<len;i++) {
	if(nclistdeleteall(l1,nclistget(l2,i))) found = 1;
    }
    return found;
}

int
nclistdeleteall(NClist* l, void* elem)
{
    int i; /* do not make unsigned */
    unsigned int len,found;
    found = 0;
    len = nclistlength(l);
    for(i=len-1;i>=0;i--) {
	void* test = nclistget(l,i);
	if(test==elem) {
	    nclistremove(l,i);
            found=1;
        }
    }
    return found;    
}

/* Collect the set of container nodes ending in "container"*/
void
collectnodepath(CDFnode* node, NClist* path, int withdataset)
{
    if(node == NULL) return;
    nclistpush(path,(void*)node);
    while(node->container != NULL) {
	node = node->container;
	if(!withdataset && node->nctype == NC_Dataset) break;
	nclistinsert(path,0,(void*)node);
    }
}

/* Like collectnodepath3, but in ocspace */
void
collectocpath(OClink conn, OCddsnode node, NClist* path)
{
    OCddsnode container;
    OCtype octype;
    if(node == NULL) return;
    oc_dds_class(conn,node,&octype);
    if(octype != OC_Dataset) {
        oc_dds_container(conn,node,&container);
        if(container != NULL)
            collectocpath(conn,container,path);
    }
    nclistpush(path,(void*)node);
}

char*
makeocpathstring(OClink conn, OCddsnode node, const char* sep)
{
    int i,len,first;
    char* result;
    char* name;
    OCtype octype;
    NClist* ocpath = NULL;
    NCbytes* pathname = NULL;

    /* If we are asking for the dataset path only,
       then nclude it, otherwise elide it
    */
    oc_dds_type(conn,node,&octype);
    if(octype == OC_Dataset) {
        oc_dds_name(conn,node,&name);
	return nulldup(name);
    }

    ocpath = nclistnew();
    collectocpath(conn,node,ocpath);
    len = nclistlength(ocpath);
    assert(len > 0); /* dataset at least */

    pathname = ncbytesnew();
    for(first=1,i=1;i<len;i++) { /* start at 1 to skip dataset name */
	OCddsnode node = (OCddsnode)nclistget(ocpath,i);
	char* name;
        oc_dds_type(conn,node,&octype);
        oc_dds_name(conn,node,&name);
	if(!first) ncbytescat(pathname,sep);
	ncbytescat(pathname,name);
	nullfree(name);
	first = 0;
    }
    result = ncbytesextract(pathname);
    ncbytesfree(pathname);
    nclistfree(ocpath);
    return result;
}

char*
makepathstring(NClist* path, const char* separator, int flags)
{
    int i,len,first;
    NCbytes* pathname = NULL;
    char* result;
    CDFnode* node;

    len = nclistlength(path);
    ASSERT(len > 0); /* dataset at least */

    if(len == 1) {/* dataset only */
        node = (CDFnode*)nclistget(path,0);
	return nulldup(node->ncbasename);
    }

    pathname = ncbytesnew();
    for(first=1,i=0;i<len;i++) {
	CDFnode* node = (CDFnode*)nclistget(path,i);
	char* name;
	if(!node->elided || (flags & PATHELIDE)==0) {
    	    if(node->nctype != NC_Dataset) {
                name = node->ncbasename;
		assert(name != NULL);
	        if(!first) ncbytescat(pathname,separator);
                ncbytescat(pathname,name);
	        first = 0;
	    }
	}
    }
    result = ncbytesextract(pathname);
    ncbytesfree(pathname);
    return result;
}

/* convert path to string using the ncname field */
char*
makecdfpathstring(CDFnode* var, const char* separator)
{
    char* spath;
    NClist* path = nclistnew();
    collectnodepath(var,path,WITHDATASET); /* <= note */
    spath = makepathstring(path,separator,PATHNC);
    nclistfree(path);
    return spath;
}

/* Collect the set names of container nodes ending in "container"*/
void
clonenodenamepath(CDFnode* node, NClist* path, int withdataset)
{
    if(node == NULL) return;
    /* stop at the dataset container as well*/
    if(node->nctype != NC_Dataset)
        clonenodenamepath(node->container,path,withdataset);
    if(node->nctype != NC_Dataset || withdataset)
        nclistpush(path,(void*)nulldup(node->ncbasename));
}

char*
simplepathstring(NClist* names,  char* separator)
{
    int i;
    size_t len;
    char* result;
    if(nclistlength(names) == 0) return nulldup("");
    for(len=0,i=0;i<nclistlength(names);i++) {
	char* name = (char*)nclistget(names,i);
	len += strlen(name);
	len += strlen(separator);
    }
    len++; /* null terminator */
    result = (char*)malloc(len);
    result[0] = '\0';
    for(i=0;i<nclistlength(names);i++) {
	char* segment = (char*)nclistget(names,i);
	if(i > 0) strcat(result,separator);
	strcat(result,segment);
    }
    return result;
}

/* Define a number of location tests */

/* Is node contained (transitively) in a sequence ? */
BOOL
dapinsequence(CDFnode* node)
{
    if(node == NULL || node->container == NULL) return TRUE;
    for(node=node->container;node->nctype != NC_Dataset;node=node->container) {
       if(node->nctype == NC_Sequence) return TRUE;
    }
    return FALSE;
}

/* Is node contained (transitively) in a structure array */
BOOL
dapinstructarray(CDFnode* node)
{
    if(node == NULL) return TRUE;
    for(node=node->container;node->nctype != NC_Dataset;node=node->container) {
       if(node->nctype == NC_Structure
	  && nclistlength(node->array.dimset0) > 0)
	    return TRUE;
    }
    return FALSE;
}

/* Is node a map field of a grid? */
BOOL
dapgridmap(CDFnode* node)
{
    if(node != NULL && node->container != NULL
       && node->container->nctype == NC_Grid) {
	CDFnode* array = (CDFnode*)nclistget(node->container->subnodes,0);
	return (node != array);
    }
    return FALSE;
}

/* Is node an array field of a grid? */
BOOL
dapgridarray(CDFnode* node)
{
    if(node != NULL && node->container != NULL
       && node->container->nctype == NC_Grid) {
	CDFnode* array = (CDFnode*)nclistget(node->container->subnodes,0);
	return (node == array);
    }
    return FALSE;
}

BOOL
dapgridelement(CDFnode* node)
{
    return dapgridarray(node)
           || dapgridmap(node);
}

/* Is node a top-level grid node? */
BOOL
daptopgrid(CDFnode* grid)
{
    if(grid == NULL || grid->nctype != NC_Grid) return FALSE;
    return daptoplevel(grid);
}

/* Is node a top-level sequence node? */
BOOL
daptopseq(CDFnode* seq)
{
    if(seq == NULL || seq->nctype != NC_Sequence) return FALSE;
    return daptoplevel(seq);
}

/* Is node a top-level node? */
BOOL
daptoplevel(CDFnode* node)
{
    if(node->container == NULL 
       || node->container->nctype != NC_Dataset) return FALSE;
    return TRUE;
}

unsigned int
modeldecode(int translation, const char* smodel,
            const struct NCTMODEL* models,
            unsigned int dfalt)
{
    for(;models->translation;models++) {
	if(translation != models->translation) continue;
	if(smodel == models->model
	   || (models->model != NULL && strcasecmp(smodel,models->model)==0)) {
	    /* We have a match */
            return models->flags;
	}
    }
    return dfalt;
}

unsigned long
getlimitnumber(const char* limit)
{
    size_t slen;
    unsigned long multiplier = 1;
    unsigned long lu;

    if(limit == NULL) return 0;
    slen = strlen(limit);
    if(slen == 0) return 0;
    switch (limit[slen-1]) {
    case 'G': case 'g': multiplier = GIGBYTE; break;
    case 'M': case 'm': multiplier = MEGBYTE; break;
    case 'K': case 'k': multiplier = KILOBYTE; break;
    default: break;
    }
    sscanf(limit,"%lu",&lu);
    return (lu*multiplier);
}

void
dapexpandescapes(char *termstring)
{
    char *s, *t, *endp;
    
    /* expand "\" escapes, e.g. "\t" to tab character;
       will only shorten string length, never increase it
    */
    s = termstring;
    t = termstring;
    while(*t) {
	if (*t == '\\') {
	    t++;
	    switch (*t) {
	      case 'a':
		*s++ = '\007'; t++; /* will use '\a' when STDC */
		break;
	      case 'b':
		*s++ = '\b'; t++;
		break;
	      case 'f':
		*s++ = '\f'; t++;
		break;
	      case 'n':
		*s++ = '\n'; t++;
		break;
	      case 'r':
		*s++ = '\r'; t++;
		break;
	      case 't':
		*s++ = '\t'; t++;
		break;
	      case 'v':
		*s++ = '\v'; t++;
		break;
	      case '\\':
		*s++ = '\\'; t++;
		break;
	      case '?':
		*s++ = '\177'; t++;
		break;
	      case 'x':
		t++; /* now t points to one or more hex digits */
		*s++ = (char) strtol(t, &endp, 16);
		t = endp;
		break;
	      case '0':
	      case '1':
	      case '2':
	      case '3':
	      case '4':
	      case '5':
	      case '6':
	      case '7': {
		/* t should now point to 3 octal digits */
		int c;
		c = t[0];
		if(c == 0 || c < '0' || c > '7') goto normal;
		c = t[1];
		if(c == 0 || c < '0' || c > '7') goto normal;
		c = t[2];
		if(c == 0 || c < '0' || c > '7') goto normal;
		c = ((t[0]-'0')<<6)+((t[1]-'0')<<3)+(t[2]-'0');
		*s++ = (char)c;
		t += 3;
		} break;
	      default:
		if(*t == 0)
		    *s++ = '\\';
		else
		    *s++ = *t++;
		break;
	    }
	} else {
normal:	    *s++ = *t++;
	}
    }
    *s = '\0';
    return;
}

#ifdef HAVE_GETTIMEOFDAY
static struct timeval time0;
static struct timeval time1;

static double
deltatime()
{
    double t0, t1;
    t0 = ((double)time0.tv_sec);
    t0 += ((double)time0.tv_usec) / 1000000.0;
    t1 = ((double)time1.tv_sec);
    t1 += ((double)time1.tv_usec) / 1000000.0;
    return (t1 - t0);
}
#endif

/* Provide a wrapper for oc_fetch so we can log what it does */
NCerror
dap_fetch(NCDAPCOMMON* nccomm, OClink conn, const char* ce,
             OCdxd dxd, OCddsnode* rootp)
{
    NCerror ncstat = NC_NOERR;
    OCerror ocstat = OC_NOERR;
    char* ext = NULL;
    OCflags flags = 0;
    int httpcode = 0;

    if(dxd == OCDDS) ext = ".dds";
    else if(dxd == OCDAS) ext = ".das";
    else ext = ".dods";

    if(ce != NULL && strlen(ce) == 0)
	ce = NULL;

    if(FLAGSET(nccomm->controls,NCF_UNCONSTRAINABLE)) {
	ce = NULL;
    }

    if(FLAGSET(nccomm->controls,NCF_ONDISK)) {
	flags |= OCONDISK;
    }

    if(SHOWFETCH) {
	/* Build uri string minus the constraint and #tag */
	char* baseurl = ncuribuild(nccomm->oc.url,NULL,ext,0);
	if(ce == NULL)
            LOG1(NCLOGNOTE,"fetch: %s",baseurl);
	else	
            LOG2(NCLOGNOTE,"fetch: %s?%s",baseurl,ce);
	nullfree(baseurl);
#ifdef HAVE_GETTIMEOFDAY
	gettimeofday(&time0,NULL);
#endif
    }
    ocstat = oc_fetch(conn,ce,dxd,flags,rootp);
    if(FLAGSET(nccomm->controls,NCF_SHOWFETCH)) {
#ifdef HAVE_GETTIMEOFDAY
        double secs;
	gettimeofday(&time1,NULL);
	secs = deltatime();
	nclog(NCLOGNOTE,"fetch complete: %0.3f secs",secs);
#else
	nclog(NCLOGNOTE,"fetch complete.");
#endif
    }
#ifdef DEBUG2
fprintf(stderr,"fetch: dds:\n");
oc_dumpnode(conn,*rootp);
#endif
    
    /* Look at the HTTP return code */
    httpcode = oc_httpcode(conn);
    if(httpcode < 400) {
        ncstat = ocerrtoncerr(ocstat);
    } else if(httpcode >= 500) {
        ncstat = NC_EDAPSVC;
    } else if(httpcode == 401) {
	ncstat = NC_EAUTH;
    } else if(httpcode == 404) {
	ncstat = NC_ENOTFOUND;
    } else {
	ncstat = NC_EACCESS;
    }
    return ncstat;
}

/* Check a name to see if it contains illegal dap characters
*/

static char* baddapchars = "./";

int
dap_badname(char* name)
{
    char* p;
    if(name == NULL) return 0;
    for(p=baddapchars;*p;p++) {
        if(strchr(name,*p) != NULL)
	    return 1;
    }
    return 0;
}

/* Repair a dap name */
char*
dap_repairname(char* name)
{
    /* assume that dap_badname was called on this name and returned 1 */
    return repairname(name,baddapchars);
}

/* Check a name to see if it contains illegal dap characters
   and repair them
*/

static const char* hexdigits = "0123456789abcdef";

static char*
repairname(const char* name, const char* badchars)
{
    char* newname;
    const char *p;
    char *q;
    int c;

    if(name == NULL) return NULL;
    newname = (char*)malloc(1+(3*strlen(name))); /* max needed */
    newname[0] = '\0'; /* so we can use strcat */
    for(p=name,q=newname;(c=*p);p++) {
        if(strchr(badchars,c) != NULL) {
	    int digit;
            char newchar[4];
	    newchar[0] = '%';	
            digit = (c & 0xf0) >> 4;
	    newchar[1] = hexdigits[digit];
            digit = (c & 0x0f);
	    newchar[2] = hexdigits[digit];
	    newchar[3] = '\0';
            strcat(newname,newchar);
            q += 3; /*strlen(newchar)*/
        } else
            *q++ = c;
	*q = '\0'; /* so we can always do strcat */
    }
    *q = '\0'; /* ensure trailing null */
    return newname;
}

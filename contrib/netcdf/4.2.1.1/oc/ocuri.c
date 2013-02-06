/*********************************************************************
 *   Copyright 2010, UCAR/Unidata
 *   See netcdf/COPYRIGHT file for copying and redistribution conditions.
 *   $Header$
 *********************************************************************/

#include "config.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "oc.h"
#include "ocuri.h"

#define OCURIDEBUG

#define LBRACKET '['
#define RBRACKET ']'

#ifndef FIX
#define FIX(s) ((s)==NULL?"":(s))
#endif

#ifndef NILLEN
#define NILLEN(s) ((s)==NULL?0:strlen(s))
#endif

#ifdef HAVE_STRDUP
#  ifndef nulldup
#  define nulldup(s) ((s)==NULL?NULL:strdup(s))
#  endif
#endif

#ifndef HAVE_STRDUP
static char* nulldup(char* s)
{
    char* dup = NULL;
    if(s != NULL) {
	dup = (char*)malloc(strlen(s)+1);
	if(dup != NULL)
	    strcpy(dup,s);
    }
    return dup;
}
#endif

#ifdef OCIGNORE
/* Not all systems have strndup, so provide one*/
static char*
ocstrndup(const char* s, size_t len)
{
    char* dup;
    if(s == NULL) return NULL;
    dup = (char*)ocmalloc(len+1);
    MEMCHECK(dup,NULL);
    memcpy((void*)dup,s,len);
    dup[len] = '\0';
    return dup;
}
#endif

/* Do not trust strncmp */
static int
ocuristrncmp(const char* s1, const char* s2, size_t len)
{
    const char *p,*q;
    if(s1 == s2) return 0;
    if(s1 == NULL) return -1;
    if(s2 == NULL) return +1;
    for(p=s1,q=s2;len > 0;p++,q++,len--) {
	if(*p != *q)
	    return (*p - *q);	
	if(*p == 0) return 0; /* *p == *q == 0 */
    }
    /* 1st len chars are same */
    return 0;
}

static char* legalprotocols[] = {
"file:",
"http:",
"https:",
"ftp:",
NULL /* NULL terminate*/
};

/* Allowable character sets for encode */
static char* fileallow = 
"0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ!#$&'()*+,-./:;=?@_~";

static char* queryallow = 
"0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ!#$&'()*+,-./:;=?@_~";

static void ocparamfree(char** params);
static int ocfind(char** params, const char* key);

/* Do a simple uri parse: return 0 if fail, 1 otherwise*/
int
ocuriparse(const char* uri0, OCURI** ocurip)
{
    OCURI* ocuri = NULL;
    char* uri;
    char** pp;
    char* p;
    char* p1;
    int c;

    /* accumulate parse points*/
    char* protocol = NULL;
    char* params = NULL;
    char* host = NULL;
    char* port = NULL;
    char* constraint = NULL;
    char* user = NULL;
    char* pwd = NULL;
    char* file = NULL;

    if(uri0 == NULL)
	return OC_EBADURL;

    ocuri = (OCURI*)calloc(1,sizeof(OCURI));
    if(ocuri == NULL) return 0;

    /* make local copy of uri */
    uri = nulldup(uri0);

    /* remove all whitespace*/
    p = uri;
    p1 = uri;
    while((c=*p1++)) {if(c != ' ' && c != '\t') *p++ = c;}

    p = uri;

    /* break up the uri string into pieces*/

    /* 1. leading bracketed parameters */
    if(*p == LBRACKET) {
	params = p+1;
	/* find end of the clientparams*/
        for(;*p;p++) {if(p[0] == RBRACKET && p[1] != LBRACKET) break;}
	if(*p == 0) goto fail; /* malformed client params*/
	*p = '\0'; /* leave off the trailing rbracket for now */
	p++; /* move past the params*/
    }

    /* verify that the uri starts with an acceptable protocol*/
    for(pp=legalprotocols;*pp;pp++) {
        if(ocuristrncmp(p,*pp,strlen(*pp))==0) break;
    }
    if(*pp == NULL) goto fail; /* illegal protocol*/
    /* save the protocol */
    protocol = *pp;

    /* 4. skip protocol */
    p += strlen(protocol);

    /* 5. skip // */
    if(*p != '/' && *(p+1) != '/')
	goto fail;
    p += 2;

    /* 6. Mark the end of the host section */
    file = strchr(p,'/');
    if(file) {
	*file++ = '\0'; /* warning: we just overwrote the leading / */
    } else
        goto fail;

    /* 7. extract any user:pwd */
    p1 = strchr(p,'@');
    if(p1) {/* Assume we have user:pwd@ */
	*p1 = '\0';
	user = p;
	pwd = strchr(p,':');
	if(!pwd) goto fail; /* malformed */
	*pwd++ = '\0';
	p = pwd+strlen(pwd)+1;
    }

    /* 8. extract host and port */
    host = p;
    port = strchr(p,':');
    if(port) {
	*port++ = '\0';
    }

    /* 9. Look for '?' */
    constraint = strchr(file,'?');
    if(constraint) {
	*constraint++ = '\0';
    }

    /* assemble the component pieces*/
    if(uri0 && strlen(uri0) > 0)
        ocuri->uri = nulldup(uri0);
    if(protocol && strlen(protocol) > 0) {
        ocuri->protocol = nulldup(protocol);
        /* remove trailing ':' */
        ocuri->protocol[strlen(protocol)-1] = '\0';
    }
    if(user && strlen(user) > 0)
        ocuri->user = nulldup(user);
    if(pwd && strlen(pwd) > 0)
        ocuri->password = nulldup(pwd);
    if(host && strlen(host) > 0)
        ocuri->host = nulldup(host);
    if(port && strlen(port) > 0)
        ocuri->port = nulldup(port);
    if(file && strlen(file) > 0) {
	/* Add back the leading / */
        ocuri->file = malloc(strlen(file)+2);
	strcpy(ocuri->file,"/");
        strcat(ocuri->file,file);
    }
    if(constraint && strlen(constraint) > 0)
        ocuri->constraint = nulldup(constraint);
    ocurisetconstraints(ocuri,constraint);
    if(params != NULL && strlen(params) > 0) {
        ocuri->params = (char*)malloc(1+2+strlen(params));
        strcpy(ocuri->params,"[");
        strcat(ocuri->params,params);
        strcat(ocuri->params,"]");
    }

#ifdef OCXDEBUG
	{
        fprintf(stderr,"ocuri:");
        fprintf(stderr," params=|%s|",FIX(ocuri->params));
        fprintf(stderr," protocol=|%s|",FIX(ocuri->protocol));
        fprintf(stderr," host=|%s|",FIX(ocuri->host));
        fprintf(stderr," port=|%s|",FIX(ocuri->port));
        fprintf(stderr," file=|%s|",FIX(ocuri->file));
        fprintf(stderr," constraint=|%s|",FIX(ocuri->constraint));
        fprintf(stderr,"\n");
    }
#endif
    free(uri);
    if(ocurip != NULL) *ocurip = ocuri;
    return 1;

fail:
    if(ocuri) ocurifree(ocuri);
    if(uri != NULL) free(uri);
    return 0;
}

void
ocurifree(OCURI* ocuri)
{
    if(ocuri == NULL) return;
    if(ocuri->uri != NULL) {free(ocuri->uri);}
    if(ocuri->protocol != NULL) {free(ocuri->protocol);}
    if(ocuri->user != NULL) {free(ocuri->user);}
    if(ocuri->password != NULL) {free(ocuri->password);}
    if(ocuri->host != NULL) {free(ocuri->host);}
    if(ocuri->port != NULL) {free(ocuri->port);}
    if(ocuri->file != NULL) {free(ocuri->file);}
    if(ocuri->constraint != NULL) {free(ocuri->constraint);}
    if(ocuri->projection != NULL) {free(ocuri->projection);}
    if(ocuri->selection != NULL) {free(ocuri->selection);}
    if(ocuri->params != NULL) {free(ocuri->params);}
    if(ocuri->paramlist != NULL) ocparamfree(ocuri->paramlist);
    free(ocuri);
}

/* Replace the constraints */
void
ocurisetconstraints(OCURI* duri,const char* constraints)
{
    char* proj = NULL;
    char* select = NULL;
    const char* p;

    if(duri->constraint == NULL) free(duri->constraint);
    if(duri->projection != NULL) free(duri->projection);
    if(duri->selection != NULL) free(duri->selection);
    duri->constraint = NULL;	
    duri->projection = NULL;	
    duri->selection = NULL;

    if(constraints == NULL || strlen(constraints)==0) return;

    duri->constraint = nulldup(constraints);
    if(*duri->constraint == '?')
	strcpy(duri->constraint,duri->constraint+1);

    p = duri->constraint;
    proj = (char*) p;
    select = strchr(proj,'&');
    if(select != NULL) {
        size_t plen = (select - proj);
	if(plen == 0) {
	    proj = NULL;
	} else {
	    proj = (char*)malloc(plen+1);
	    memcpy((void*)proj,p,plen);
	    proj[plen] = '\0';
	}
	select = nulldup(select);
    } else {
	proj = nulldup(proj);
	select = NULL;
    }
    duri->projection = proj;
    duri->selection = select;
}


/* Construct a complete OC URI without the client params
   and optionally with the constraints;
   caller frees returned string.
   Optionally encode the pieces.
*/

char*
ocuribuild(OCURI* duri, const char* prefix, const char* suffix, int flags)
{
    size_t len = 0;
    char* newuri;
    char* tmpfile;
    char* tmpsuffix;
    char* tmpquery;

    int withparams = ((flags&OCURIPARAMS)
			&& duri->params != NULL);
    int withuserpwd = ((flags&OCURIUSERPWD)
	               && duri->user != NULL && duri->password != NULL);
    int withconstraints = ((flags&OCURICONSTRAINTS)
	                   && duri->constraint != NULL);
#ifdef NEWESCAPE
    int encode = (flags&OCURIENCODE);
#else
    int encode = 0;
#endif

    if(prefix != NULL) len += NILLEN(prefix);
    if(withparams) {
	len += NILLEN("[]");
	len += NILLEN(duri->params);
    }
    len += (NILLEN(duri->protocol)+NILLEN("://"));
    if(withuserpwd) {
	len += (NILLEN(duri->user)+NILLEN(duri->password)+NILLEN(":@"));
    }
    len += (NILLEN(duri->host));
    if(duri->port != NULL) {
	len += (NILLEN(":")+NILLEN(duri->port));
    }
    
    tmpfile = duri->file;
    if(encode)
	tmpfile = ocuriencode(tmpfile,fileallow);
    len += (NILLEN(tmpfile));

    if(suffix != NULL) {
        tmpsuffix = (char*)suffix;
        if(encode)
	    tmpsuffix = ocuriencode(tmpsuffix,fileallow);
        len += (NILLEN(tmpsuffix));
    }

    if(withconstraints) {
	tmpquery = duri->constraint;
        if(encode)
	    tmpquery = ocuriencode(tmpquery,queryallow);
        len += (NILLEN("?")+NILLEN(tmpquery));
    }

    len += 1; /* null terminator */
    
    newuri = (char*)malloc(len);
    if(newuri == NULL) return NULL;

    newuri[0] = '\0';
    if(prefix != NULL) strcat(newuri,prefix);
    if(withparams) {
	strcat(newuri,"[");
	strcat(newuri,duri->params);
	strcat(newuri,"]");
    }
    if(duri->protocol != NULL)
	strcat(newuri,duri->protocol);
    strcat(newuri,"://");
    if(withuserpwd) {
        strcat(newuri,duri->user);
        strcat(newuri,":");
        strcat(newuri,duri->password);	
        strcat(newuri,"@");
    }
    if(duri->host != NULL) { /* may be null if using file: protocol */
        strcat(newuri,duri->host);	
    }
    if(duri->port != NULL) {
        strcat(newuri,":");
        strcat(newuri,duri->port);
    }

    strcat(newuri,tmpfile);
    if(suffix != NULL) strcat(newuri,tmpsuffix);
    if(withconstraints) {
	strcat(newuri,"?");
	strcat(newuri,tmpquery);
    }
    return newuri;
}

/**************************************************/
/* Parameter support */

/*
Client parameters are assumed to be
one or more instances of bracketed pairs:
e.g "[...][...]...".
The bracket content in turn is assumed to be a
comma separated list of <name>=<value> pairs.
e.g. x=y,z=,a=b.
If the same parameter is specifed more than once,
then the first occurrence is used; this is so that
is possible to forcibly override user specified
parameters by prefixing.
IMPORTANT: client parameter string is assumed to
have blanks compress out.
Returns 1 if parse suceeded, 0 otherwise;
*/

int
ocuridecodeparams(OCURI* ocuri)
{
    char* cp;
    char* cq;
    int c;
    int i;
    int nparams;
    char* params0;
    char* params;
    char* params1;
    char** plist;

    if(ocuri == NULL) return 0;
    if(ocuri->params == NULL) return 1;

    params0 = ocuri->params;

    /* Pass 1 to replace beginning '[' and ending ']' */
    if(params0[0] == '[') 
	params = nulldup(params0+1);
    else
	params = nulldup(params0);	

    if(params[strlen(params)-1] == ']')
	params[strlen(params)-1] = '\0';

    /* Pass 2 to replace "][" pairs with ','*/
    params1 = nulldup(params);
    cp=params; cq = params1;
    while((c=*cp++)) {
	if(c == RBRACKET && *cp == LBRACKET) {cp++; c = ',';}
	*cq++ = c;
    }
    *cq = '\0';
    free(params);
    params = params1;

    /* Pass 3 to break string into pieces and count # of pairs */
    nparams=0;
    for(cp=params;(c=*cp);cp++) {
	if(c == ',') {*cp = '\0'; nparams++;}
    }
    nparams++; /* for last one */

    /* plist is an env style list */
    plist = (char**)calloc(1,sizeof(char*)*(2*nparams+1)); /* +1 for null termination */

    /* Pass 4 to break up each pass into a (name,value) pair*/
    /* and insert into the param list */
    /* parameters of the form name name= are converted to name=""*/
    cp = params;
    for(i=0;i<nparams;i++) {
	char* next = cp+strlen(cp)+1; /* save ptr to next pair*/
	char* vp;
	/*break up the ith param*/
	vp = strchr(cp,'=');
	if(vp != NULL) {*vp = '\0'; vp++;} else {vp = "";}
	plist[2*i] = nulldup(cp);	
	plist[2*i+1] = nulldup(vp);
	cp = next;
    }
    plist[2*nparams] = NULL;
    free(params);
    if(ocuri->paramlist != NULL)
	ocparamfree(ocuri->paramlist);
    ocuri->paramlist = plist;
    return 1;
}

const char*
ocurilookup(OCURI* uri, const char* key)
{
    int i;
    if(uri == NULL || key == NULL || uri->params == NULL) return NULL;
    if(uri->paramlist == NULL) {
	i = ocuridecodeparams(uri);
	if(!i) return 0;
    }
    i = ocfind(uri->paramlist,key);
    if(i >= 0)
	return uri->paramlist[(2*i)+1];
    return NULL;
}

int
ocurisetparams(OCURI* uri, const char* newparams)
{
    if(uri == NULL) return 0;
    if(uri->paramlist != NULL) ocparamfree(uri->paramlist);
    uri->paramlist = NULL;
    if(uri->params != NULL) free(uri->params);
    uri->params = nulldup(newparams);
    return 1;
}

/* Internal version of lookup; returns the paired index of the key */
static int
ocfind(char** params, const char* key)
{
    int i;
    char** p;
    for(i=0,p=params;*p;p+=2,i++) {
	if(strcmp(key,*p)==0) return i;
    }
    return -1;
}

static void
ocparamfree(char** params)
{
    char** p;
    if(params == NULL) return;
    for(p=params;*p;p+=2) {
	free(*p);
	if(p[1] != NULL) free(p[1]);
    }
    free(params);
}

#ifdef OCIGNORE
/*
Delete the entry.
return value = 1 => found and deleted;
               0 => param not found
*/
int
ocparamdelete(char** params, const char* key)
{
    int i;
    char** p;
    char** q;
    if(params == NULL || key == NULL) return 0;
    i = ocfind(params,key);
    if(i < 0) return 0;
    p = params+(2*i);
    for(q=p+2;*q;) {	
	*p++ = *q++;
    }
    *p = NULL;
    return 1;
}

static int
oclength(char** params)
{
    int i = 0;
    if(params != NULL) {
	while(*params) {params+=2; i++;}
    }
    return i;
}

/*
Insert new client param (name,value);
return value = 1 => not already defined
               0 => param already defined (no change)
*/
char**
ocparaminsert(char** params, const char* key, const char* value)
{
    int i;
    char** newp;
    size_t len;
    if(params == NULL || key == NULL) return 0;
    i = ocfind(params,key);
    if(i >= 0) return 0;
    /* not found, append */
    i = oclength(params);
    len = sizeof(char*)*((2*i)+1);
    newp = realloc(params,len+2*sizeof(char*));
    memcpy(newp,params,len);
    newp[2*i] = nulldup(key);
    newp[2*i+1] = (value==NULL?NULL:nulldup(value));
    return newp;
}

/*
Replace new client param (name,value);
return value = 1 => replacement performed
               0 => key not found (no change)
*/
int
ocparamreplace(char** params, const char* key, const char* value)
{
    int i;
    if(params == NULL || key == NULL) return 0;
    i = ocfind(params,key);
    if(i < 0) return 0;
    if(params[2*i+1] != NULL) free(params[2*i+1]);
    params[2*i+1] = nulldup(value);
    return 1;
}
#endif


/* Provide % encoders and decoders */


static char* hexchars = "0123456789abcdefABCDEF";

static void
toHex(unsigned int b, char hex[2])
{
    hex[0] = hexchars[(b >> 4) & 0xff];
    hex[1] = hexchars[(b) & 0xff];
}


static unsigned int
fromHex(int c)
{
    if(c >= '0' && c <= '9') return (c - '0');
    if(c >= 'a' && c <= 'f') return (10 + (c - 'a'));
    if(c >= 'A' && c <= 'F') return (10 + (c - 'A'));
    return -1;
}


/* Return a string representing encoding of input; caller must free;
   watch out: will encode whole string, so watch what you give it.
   Allowable argument specifies characters that do not need escaping.
 */

char*
ocuriencode(char* s, char* allowable)
{
    size_t slen;
    char* encoded;
    char* inptr;
    char* outptr;

    if(s == NULL) return NULL;

    slen = strlen(s);
    encoded = (char*)malloc((3*slen) + 1); /* max possible size */

    for(inptr=s,outptr=encoded;*inptr;) {
	int c = *inptr++;
        if(c == ' ') {
	    *outptr++ = '+';
        } else {
            /* search allowable */
            int c2;
	    char* a = allowable;
	    while((c2=*a++)) {
		if(c == c2) break;
	    }
            if(c2) {*outptr++ = c;}
            else {
		char hex[2];
		toHex(c,hex);
		*outptr++ = '%';
		*outptr++ = hex[0];
		*outptr++ = hex[1];
            }
        }
    }
    *outptr = '\0';
    return encoded;
}

/* Return a string representing decoding of input; caller must free;*/
char*
ocuridecode(char* s)
{
    return ocuridecodeonly(s,NULL);
}

/* Return a string representing decoding of input only for specified
   characters;  caller must free
*/
char*
ocuridecodeonly(char* s, char* only)
{
    size_t slen;
    char* decoded;
    char* outptr;
    char* inptr;
    unsigned int c;
    
    if (s == NULL) return NULL;
    if(only == NULL) only = "";

    slen = strlen(s);
    decoded = (char*)malloc(slen+1); /* Should be max we need */

    outptr = decoded;
    inptr = s;
    while((c = *inptr++)) {
	if(c == '+' && strchr(only,'+') != NULL)
	    *outptr++ = ' ';
	else if(c == '%') {
            /* try to pull two hex more characters */
	    if(inptr[0] != '\0' && inptr[1] != '\0'
		&& strchr(hexchars,inptr[0]) != NULL
		&& strchr(hexchars,inptr[1]) != NULL) {
		/* test conversion */
		int xc = (fromHex(inptr[0]) << 4) | (fromHex(inptr[1]));
		if(strchr(only,xc) != NULL) {
		    inptr += 2; /* decode it */
		    c = xc;
                }
            }
        }
        *outptr++ = c;
    }
    *outptr = '\0';
    return decoded;
}


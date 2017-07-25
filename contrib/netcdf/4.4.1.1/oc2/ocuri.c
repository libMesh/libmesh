/*********************************************************************
 *   Copyright 2010, UCAR/Unidata
 *   See netcdf/COPYRIGHT file for copying and redistribuution conditions.
 *   $Header$
 *********************************************************************/

#include "config.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>

#include "oc.h"
#include "ocuri.h"

#undef OCURIDEBUG

#ifdef OCURIDEBUG
static int failpoint = 0;
#define THROW(n) {failpoint=(n); goto fail;}
#else
#define THROW(n)
#endif


#define PADDING 8

#define LBRACKET '['
#define RBRACKET ']'
#define EOFCHAR '\0'

#ifndef FIX
#define FIX(s) ((s)==NULL?"NULL":(s))
#endif

#ifndef NILLEN
#define NILLEN(s) ((s)==NULL?0:strlen(s))
#endif

#ifdef HAVE_STRDUP
#ifndef nulldup
#define nulldup(s) ((s)==NULL?NULL:strdup(s))
#endif
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

#define terminate(p) {*(p) = EOFCHAR;}

#define endof(p) ((p)+strlen(p))

static struct OC_ProtocolInfo {
char* name;
int   filelike; /* 1=>this protocol has no host, user+pwd, or port */
} legalprotocols[] = {
{"file",1},
{"http",0},
{"https",0},
{"ftp",0},
};

/* Allowable character sets for encode */
static char* fileallow =
"0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ!#$&'()*+,-./:;=?@_~";

static char* queryallow =
"0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ!#$&'()*+,-./:;=?@_~";

/* Forward */
static void ocparamfree(char** params);
static int ocfind(char** params, const char* key);
static void oclshift1(char* p);
static void ocrshift1(char* p);
static char* oclocate(char* p, const char* charlist);
static void ocappendparams(char* newuri, char** p);

/* Do a simple uri parse: return 0 if fail, 1 otherwise*/
int
ocuriparse(const char* uri0, OCURI** durip)
{
    OCURI* duri = NULL;
    char* uri = NULL;
    char* p;
    char* q;
    struct OC_ProtocolInfo* proto;
    int i,nprotos;

    /* accumulate parse points*/
    char* protocol = NULL;
    char* host = NULL;
    char* port = NULL;
    char* constraint = NULL;
    char* userpwd = NULL;
    char* file = NULL;
    char* prefixparams = NULL;
    char* suffixparams = NULL;

    if(uri0 == NULL || strlen(uri0) == 0)
	{THROW(1); goto fail;}

    duri = (OCURI*)calloc(1,sizeof(OCURI));
    if(duri == NULL)
	{THROW(2); goto fail;}

    /* save original uri */
    duri->uri = nulldup(uri0);

    /* make local copy of uri */
    uri = (char*)malloc(strlen(uri0)+1+PADDING); /* +1 for trailing null,
                                                    +PADDING for shifting */
    if(uri == NULL)
	{THROW(3); goto fail;}

    /* strings will be broken into pieces with intermixed '\0; characters;
       first char is guaranteed to be '\0' */

    duri->strings = uri;
    uri++;

    /* dup the incoming url */
    strcpy(uri,uri0);

    /* Walk the uri and do the following:
	1. remove all whitespace
	2. remove all '\\' (Temp hack to remove escape characters
                            inserted by Windows or MinGW)
    */
    for(q=uri,p=uri;*p;p++) {
	if(*p != '\\' && *p >= ' ') /* compress out */
	    *q++=*p;
    }
    p = uri;

    /* break up the uri string into big chunks: prefixparams, protocol,
       host section, and the file section (i.e. remainder)
    */

    /* collect any prefix bracketed parameters */
    if(*p == LBRACKET) {
	p++;
	prefixparams = p;
	/* find end of the clientparams; convert LB,RB to '&' */
        for(q=p;*p;p++) {
	    if(p[0] == RBRACKET && p[1] == LBRACKET) {
		*q++ = '&';
		p++;
	    } else if(p[0] == RBRACKET && p[1] != LBRACKET)
		break;
	    else
		*q++=*p;
	}
	if(*p == 0)
	    {THROW(4); goto fail; /* malformed client params*/}
        terminate(q); /* nul term the prefixparams */
	p++; /* move past the final RBRACKET */
    }

    /* Tag the protocol */
    protocol = p;
    p = strchr(p,':');
    if(!p)
	{THROW(5); goto fail;}
    terminate(p); /*overwrite colon*/
    p++; /* skip the colon */

    /* verify that the uri starts with an acceptable protocol*/
    nprotos = (sizeof(legalprotocols)/sizeof(struct OC_ProtocolInfo));
    proto = NULL;
    for(i=0;i<nprotos;i++) {
        if(strcmp(protocol,legalprotocols[i].name)==0) {
	    proto = &legalprotocols[i];
	    break;
	}
    }
    if(proto == NULL)
	{THROW(6); goto fail; /* illegal protocol*/}

    /* skip // */
    if(p[0] != '/' || p[1] != '/')
	{THROW(7); goto fail;}
    p += 2;

    /* If this is all we have (proto://) then fail */
    if(*p == EOFCHAR)
	{THROW(8); goto fail;}

    /* establish the start of the file section */
    if(proto->filelike) {/* everything after proto:// */
	file = p;
	host = NULL; /* and no host section */
    } else { /*!proto->filelike => This means there should be a host section */
        /* locate the end of the host section and therefore the start
           of the file section */
	host = p;
        p  = oclocate(p,"/?#");
	if(p == NULL) {
	    file = endof(host); /* there is no file section */
	} else {
	    ocrshift1(p); /* make room to terminate the host section
                             without overwriting the leading character */
	    terminate(p); /* terminate the host section */
	    file = p+1; /* +1 becauseof the shift */
	}
    }

    /* If you shift in the code below, you must reset file beginning */

    if(host != NULL) {/* Parse the host section */
	/* Check for leading user:pwd@ */
        p = strchr(host,'@');
        if(p) {
	    if(p == host)
		{THROW(9); goto fail; /* we have proto://@ */}
	    userpwd = host;
	    terminate(p); /* overwrite '@' */
	    host = p+1; /* start of host ip name */
	}

        /* extract host and port */
	p = host;
        p = strchr(p,':');
        if(p != NULL) {
	    terminate(p);
	    p++;
	    port = p;
	    if(*port == EOFCHAR)
		{THROW(11); goto fail; /* we have proto://...:/ */}
	    /* The port must look something like a number */
	    for(;*p;p++) {
	        if(strchr("0123456789-",*p) == NULL)
		    {THROW(12); goto fail;  /* probably not a real port, fail */}
	    }
	} /* else *p == NULL */


        /* check for empty host section */
	if(*host == EOFCHAR)
	    {THROW(13); goto fail;}

    }

    assert(file != NULL);
    p = file;

    /* find the end of the file section and the start of the
       constraints and/or suffixparams
    */
    p = oclocate(p,"?#");
    if(p != NULL) { /* we have constraint and/or suffixparams */
	char* fileend = p; /* save the end of the file section */
	char* constraintend = NULL;
	if(*p == '?')
            constraint = p+1;
	else
	    constraint = NULL;
	p = strchr(p,'#'); /* may repeat effect of oclocate above */
	if(p != NULL) {
	    constraintend = p;
	    suffixparams = p+1;
	} else
	    suffixparams = NULL;
	/* Ok, terminate the pieces */
	terminate(fileend); /* terminate file section */
	if(constraint != NULL && constraintend != NULL)
	    terminate(constraintend);
	/* Suffix params are already terminated
           since they should be the last section
           of the original url
        */
    }

    /* check for empty sections */
    if(file != NULL && *file == EOFCHAR)
	file = NULL; /* empty file section */
    if(constraint != NULL && *constraint == EOFCHAR)
	constraint = NULL; /* empty constraint section */
    if(suffixparams != NULL && *suffixparams == EOFCHAR)
	suffixparams = NULL; /* empty suffixparams section */

    if(suffixparams != NULL) {
	if(*suffixparams == EOFCHAR)
	    suffixparams = NULL; /* suffixparams are empty */
    }

    /* do last minute empty check */
    if(protocol != NULL && *protocol == EOFCHAR) protocol = NULL;
    if(userpwd != NULL && *userpwd == EOFCHAR) userpwd = NULL;
    if(host != NULL && *host == EOFCHAR) host = NULL;
    if(port != NULL && *port == EOFCHAR) port = NULL;
    if(file != NULL && *file == EOFCHAR) file = NULL;
    if(constraint != NULL && *constraint == EOFCHAR) constraint = NULL;

    /* assemble the component pieces */
    duri->protocol = protocol;
    duri->userpwd = userpwd;
    duri->host = host;
    duri->port = port;
    duri->file = file;

    ocurisetconstraints(duri,constraint);

    /* concat suffix and prefix params */
    if(prefixparams != NULL || suffixparams != NULL) {
	size_t plen = prefixparams ? strlen(prefixparams) : 0;
	size_t slen = suffixparams ? strlen(suffixparams) : 0;
	size_t space = plen + slen + 1;
	/* add 1 for an extra ampersand if both are defined */
    if(plen > 0 && slen > 0) space++;
    /* Add an extra char for null termination. */
    duri->params = (char*)malloc(space+1);
    if(duri->params == NULL)
      return 0;
    duri->params[0] = EOFCHAR; /* so we can use strcat */
	if(plen > 0) {
      strncat(duri->params,prefixparams,space);
      if(slen > 0)
		strncat(duri->params,"&",space);
	}
	if(slen > 0)
      strncat(duri->params,suffixparams,space);
    }

#ifdef OCURIDEBUG
	{
	int i,nparms;
	char** p;
        fprintf(stderr,"duri:");
        fprintf(stderr," protocol=|%s|",FIX(duri->protocol));
        fprintf(stderr," host=|%s|",FIX(duri->host));
        fprintf(stderr," port=|%s|",FIX(duri->port));
        fprintf(stderr," file=|%s|",FIX(duri->file));
        fprintf(stderr," constraint=|%s|",FIX(duri->constraint));
        fprintf(stderr," params=|%s|",FIX(duri->params));
        fprintf(stderr,"\n");
	if(duri->paramlist == NULL) {
	    if(!ocuridecodeparams(duri)) {
		fprintf(stderr,"DEBUG: param decode failed\n");
		duri->paramlist = NULL;
	    }
	}
	if(duri->paramlist != NULL) {
	    for(p=duri->paramlist,nparms=0;*p;p++,nparms++);
	    nparms = nparms / 2;
	    fprintf(stderr,"params:");
	    for(i=0;i<nparms;i++) {
	        char** pos = duri->paramlist+(i*2);
	        fprintf(stderr," %s=|%s|",pos[0],pos[1]);
	    }
            fprintf(stderr,"\n");
	}
    }
#endif
    if(durip != NULL) *durip = duri; else free(duri);
    return 1;

fail:
    if(duri != NULL) {
	ocurifree(duri);
    }
    return 0;
}

void
ocurifree(OCURI* duri)
{
    if(duri == NULL) return;
    if(duri->uri != NULL) {free(duri->uri);}
    if(duri->params != NULL) {free(duri->params);}
    if(duri->paramlist != NULL) ocparamfree(duri->paramlist);
    if(duri->strings != NULL) {free(duri->strings);}
    if(duri->constraint != NULL) {free(duri->constraint);}
    if(duri->projection != NULL) {free(duri->projection);}
    if(duri->selection != NULL) {free(duri->selection);}
    free(duri);
}

/* Replace the constraints */
void
ocurisetconstraints(OCURI* duri,const char* constraints)
{
    char* proj = NULL;
    char* select = NULL;
    const char* p;

    if(duri->constraint != NULL) free(duri->constraint);
    if(duri->projection != NULL) free(duri->projection);
    if(duri->selection != NULL) free(duri->selection);
    duri->constraint = NULL;
    duri->projection = NULL;
    duri->selection = NULL;

    if(constraints == NULL || strlen(constraints)==0) return;

    duri->constraint = nulldup(constraints);
    if(*duri->constraint == '?')
	oclshift1(duri->constraint);

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
	    proj[plen] = EOFCHAR;
	}
	select = nulldup(select);
    } else {
	proj = nulldup(proj);
	select = NULL;
    }
    duri->projection = proj;
    duri->selection = select;
}


/* Construct a complete OC URI.
   Optionally with the constraints.
   Optionally with the user parameters.
   Caller frees returned string.
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
    int nparams = 0;
    int paramslen = 0;

    /* if both are specified, prefix has priority */
    int withsuffixparams = ((flags&OCURISUFFIXPARAMS)!=0
				&& duri->params != NULL);
    int withprefixparams = ((flags&OCURIPREFIXPARAMS)!=0
				&& duri->params != NULL);
    int withuserpwd = ((flags&OCURIUSERPWD)!=0
	               && duri->userpwd != NULL);
    int withconstraints = ((flags&OCURICONSTRAINTS)!=0
	                   && duri->constraint != NULL);
#ifdef NEWESCAPE
    int encode = (flags&OCURIENCODE);
#else
    int encode = 0;
#endif

    if(prefix != NULL) len += NILLEN(prefix);
    len += (NILLEN(duri->protocol)+NILLEN("://"));
    if(withuserpwd)
	len += (NILLEN(duri->userpwd)+NILLEN("@"));
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

    if(withprefixparams || withsuffixparams) {
	char** p;
	if(duri->paramlist == NULL)
	    if(!ocuridecodeparams(duri))
		return NULL;
	for(paramslen=0,nparams=0,p=duri->paramlist;*p;p++) {
	    nparams++;
	    paramslen += NILLEN(*p);
	}
	if(nparams % 2 == 1)
	    return NULL; /* malformed */
	nparams = (nparams / 2);
	len += paramslen;
	len += 3*nparams; /* for brackets for every param plus possible = */
	if(withsuffixparams)
	    len += strlen("#");
    }

    len += 1; /* null terminator */

    newuri = (char*)malloc(len);
    if(newuri == NULL) return NULL;

    newuri[0] = EOFCHAR;
    if(prefix != NULL) strcat(newuri,prefix);
    if(withprefixparams) {
	ocappendparams(newuri,duri->paramlist);
    }
    if(duri->protocol != NULL)
	strcat(newuri,duri->protocol);
    strcat(newuri,"://");
    if(withuserpwd) {
        strcat(newuri,duri->userpwd);
        strcat(newuri,"@");
    }
    if(duri->host != NULL) { /* may be null if using file: protocol */
        strcat(newuri,duri->host);
    }
    if(duri->port != NULL) {
        strcat(newuri,":");
        strcat(newuri,duri->port);
    }

    if(tmpfile != NULL) {
        strcat(newuri,tmpfile);
        if(suffix != NULL) strcat(newuri,tmpsuffix);
    }
    if(withconstraints) {
	strcat(newuri,"?");
	strcat(newuri,tmpquery);
    }
    if(withsuffixparams & !withprefixparams) {
	strcat(newuri,"#");
	ocappendparams(newuri,duri->paramlist);
    }
    return newuri;
}

static void
ocappendparams(char* newuri, char** p)
{
	while(*p) {
	    strcat(newuri,"[");
	    strcat(newuri,*p++);
	    if(strlen(*p) > 0) {
	        strcat(newuri,"=");
	        strcat(newuri,*p);
	    }
	    p++;
	    strcat(newuri,"]");
	}
}

/**************************************************/
/* Parameter support */

/*
In the original url, client parameters are assumed to be one
or more instances of bracketed pairs: e.g "[...][...]...".
prefixed to the url. This model has been extended to support
specification of the parameters as semicolon separated key=value
pairs in the fragment part of the url.  The fragment part
starts with a '#' and is the last part of the url.

After the url is parsed, the parameter list
is converted to a semicolon separated list with all
whitespace removed.
In any case, each parameter in turn is assumed to be a
of the form <name>=<value> or <name>.
e.g. x=y,z,a=b,w.  If the same parameter is specified more
than once, then the last occurrence is used; this is so
that is possible to forcibly override user specified
parameters by suffixing.  IMPORTANT: client parameter string
is assumed to have blanks compressed out.  Returns 1 if parse
succeeded, 0 otherwise; */

int
ocuridecodeparams(OCURI* ocuri)
{
    char* p;
    int i,c;
    int nparams;
    char* params = NULL;
    char** plist;

    if(ocuri == NULL) return 0;
    if(ocuri->params == NULL) return 1;

    params = strdup(ocuri->params);
    if(params == NULL)
	return 0; /* no memory */

    /* Pass 1:  break string into pieces at the ampersands
       and count # of pairs */
    nparams=0;
    for(p=params;*p;p++) {
	c = *p;
	if(c == '&') {*p = EOFCHAR; nparams++;}
    }
    nparams++; /* for last one */

    /* plist will be an env style list */
    plist = (char**)calloc(1,sizeof(char*)*(2*nparams+1)); /* +1 for null termination */
    if(plist == NULL) {
	free(params);
	return 0;
    }

    /* Break up each param into a (name,value) pair*/
    /* and insert into the param list */
    /* parameters of the form name name= are converted to name=""*/
    for(p=params,i=0;i<nparams;i++) {
      char* next = p+strlen(p)+1; /* save ptr to next pair*/
      char* vp;
      /*break up the ith param*/
      vp = strchr(p,'=');
      if(vp != NULL) {*vp = EOFCHAR; vp++;} else {vp = "";}
      plist[2*i] = nulldup(p);
      plist[2*i+1] = nulldup(vp);
      p = next;
    }
    plist[2*nparams] = NULL;
    free(params);
    if(ocuri->paramlist != NULL)
	ocparamfree(ocuri->paramlist);
    ocuri->paramlist = plist;
    return 1;
}

int
ocurilookup(OCURI* uri, const char* key, const char** resultp)
{
    int i;
    char* value = NULL;
    if(uri == NULL || key == NULL || uri->params == NULL) return 0;
    if(uri->paramlist == NULL) {
	i = ocuridecodeparams(uri);
	if(!i) return 0;
    }
    i = ocfind(uri->paramlist,key);
    if(i < 0)
	return 0;
    value = uri->paramlist[(2*i)+1];
    if(resultp) *resultp = value;
    return 1;
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


/* Return the ptr to the first occurrence of
   any char in the list. Return NULL if no
   occurrences
*/
static char*
oclocate(char* p, const char* charlist)
{
    for(;*p;p++) {
	if(strchr(charlist,*p) != NULL)
	    return p;
    }
    return NULL;
}


/* Shift every char starting at p 1 place to the left */
static void
oclshift1(char* p)
{
    if(p != NULL && *p != EOFCHAR) {
	char* q = p++;
	while((*q++=*p++));
    }
}

/* Shift every char starting at p 1 place to the right */
static void
ocrshift1(char* p)
{
    char cur;
    cur = 0;
    do {
	char next = *p;
	*p++ = cur;
	cur = next;
    } while(cur != 0);
    *p = 0; /* make sure we are still null terminated */
}


/* Provide % encoders and decoders */


static char* hexchars = "0123456789abcdefABCDEF";

static void
toHex(int b, char* hex)
{
    hex[0] = hexchars[(b >> 4) & 0xff];
    hex[1] = hexchars[(b) & 0xff];
}


static int
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
    *outptr = EOFCHAR;
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
	    if(inptr[0] != EOFCHAR && inptr[1] != EOFCHAR
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
    *outptr = EOFCHAR;
    return decoded;
}

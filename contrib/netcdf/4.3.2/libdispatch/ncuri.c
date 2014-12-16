/*********************************************************************
 *   Copyright 2010, UCAR/Unidata
 *   See netcdf/COPYRIGHT file for copying and redistribution conditions.
 *   $Header$
 *********************************************************************/
#include "config.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>

#include "ncuri.h"

#define NCURIDEBUG

#ifdef NCURIDEBUG
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

#ifndef nulldup
#define nulldup(s) ((s)==NULL?NULL:strdup(s))
#endif

#define terminate(p) {*(p) = EOFCHAR;}

#define endof(p) ((p)+strlen(p))

static struct NC_ProtocolInfo {
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
static void ncparamfree(char** params);
static int ncfind(char** params, const char* key);
static void nclshift1(char* p);
static void ncrshift1(char* p);
static char* nclocate(char* p, const char* charlist);
static void ncappendparams(char* newuri, char** p);

/* Do a simple uri parse: return 0 if fail, 1 otherwise*/
int
ncuriparse(const char* uri0, NCURI** durip)
{
    NCURI* duri = NULL;
    char* uri = NULL;
    char* p;
    struct NC_ProtocolInfo* proto;
    int i,nprotos;

    /* accumulate parse points*/
    char* protocol = NULL;
    char* host = NULL;
    char* port = NULL;
    char* constraint = NULL;
    char* user = NULL;
    char* pwd = NULL;
    char* file = NULL;
    char* prefixparams = NULL;
    char* suffixparams = NULL;

    if(uri0 == NULL || strlen(uri0) == 0)
	{THROW(1); goto fail;}

    duri = (NCURI*)calloc(1,sizeof(NCURI));
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
    for(p=uri;*p;p++) {
	if(*p == '\\' || *p < ' ')
	    nclshift1(p); /* compress out */
    }	

    p = uri;

    /* break up the uri string into big chunks: prefixparams, protocol,
       host section, and the file section (i.e. remainder)
    */

    /* collect any prefix bracketed parameters */
    if(*p == LBRACKET) {
	prefixparams = p+1;
	/* find end of the clientparams; convert LB,RB to '&' */
        for(;*p;p++) {
	    if(p[0] == RBRACKET && p[1] == LBRACKET) {
		p[0] = '&';
		nclshift1(p+1);
	    } else if(p[0] == RBRACKET && p[1] != LBRACKET)
		break;
	}
	if(*p == 0)
	    {THROW(4); goto fail; /* malformed client params*/}
        terminate(p); /* nul term the prefixparams (overwrites
                         the final RBRACKET) */
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
    nprotos = (sizeof(legalprotocols)/sizeof(struct NC_ProtocolInfo));
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
    if(p[0] != '/' && p[1] != '/')
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
        p  = nclocate(p,"/?#");
	if(p == NULL) {
	    file = endof(host); /* there is no file section */
	} else {
	    ncrshift1(p); /* make room to terminate the host section
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
	    user = host;
	    terminate(p); /* overwrite '@' */
	    host = p+1; /* start of host ip name */
	    p = strchr(user,':');
 	    if(p == NULL)
		{THROW(10); goto fail; /* malformed */}
	    terminate(p); /*overwrite colon */
	    pwd = p+1;
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
    p = nclocate(p,"?#");
    if(p != NULL) { /* we have constraint and/or suffixparams */
	char* fileend = p; /* save the end of the file section */
	char* constraintend = NULL; 
	if(*p == '?')
            constraint = p+1;
	else
	    constraint = NULL;
	p = strchr(p,'#'); /* may repeat effect of nclocate above */
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
	/* there really are suffix params; so rebuild the suffix params */
	if(*suffixparams == LBRACKET) suffixparams++;
        p = suffixparams;
	/* convert RBRACKET LBRACKET to '&' */
        for(;*p;p++) {
	    if(p[0] == RBRACKET && p[1] == LBRACKET) {
	        p[0] = '&';
		nclshift1(p+1);
	    } else if(p[0] == RBRACKET && p[1] != LBRACKET) {
		/* terminate suffixparams */
		*p = EOFCHAR;
		break;
	    }
	}
	if(*suffixparams == EOFCHAR)
	    suffixparams = NULL; /* suffixparams are empty */
    }

    /* do last minute empty check */
    
    if(protocol != NULL && *protocol == EOFCHAR) protocol = NULL;
    if(user != NULL && *user == EOFCHAR) user = NULL;
    if(pwd != NULL && *pwd == EOFCHAR) pwd = NULL;
    if(host != NULL && *host == EOFCHAR) host = NULL;
    if(port != NULL && *port == EOFCHAR) port = NULL;
    if(file != NULL && *file == EOFCHAR) file = NULL;
    if(constraint != NULL && *constraint == EOFCHAR) constraint = NULL;

    /* assemble the component pieces */
    duri->protocol = protocol;
    duri->user = user;
    duri->password = pwd;
    duri->host = host;
    duri->port = port;
    duri->file = file;

    ncurisetconstraints(duri,constraint);

    /* concat suffix and prefix params */
    if(prefixparams != NULL || suffixparams != NULL) {
	int plen = prefixparams ? strlen(prefixparams) : 0;
	int slen = suffixparams ? strlen(suffixparams) : 0;
	int space = plen + slen + 1;
	/* add 1 for an extra ampersand if both are defined */
        space++;
        duri->params = (char*)malloc(space);
	duri->params[0] = EOFCHAR; /* so we can use strcat */
	if(plen > 0) {
            strcat(duri->params,prefixparams);
	    if(slen > 0)
		strcat(duri->params,"&");
	}
	if(slen > 0)
            strcat(duri->params,suffixparams);
    }

#ifdef NCXDEBUG
	{
        fprintf(stderr,"duri:");
        fprintf(stderr," params=|%s|",FIX(duri->params));
        fprintf(stderr," protocol=|%s|",FIX(duri->protocol));
        fprintf(stderr," host=|%s|",FIX(duri->host));
        fprintf(stderr," port=|%s|",FIX(duri->port));
        fprintf(stderr," file=|%s|",FIX(duri->file));
        fprintf(stderr," constraint=|%s|",FIX(duri->constraint));
        fprintf(stderr,"\n");
    }
#endif
    if(durip != NULL) 
      *durip = duri;
    else
      ncurifree(duri);

    return 1;

fail:
    if(duri != NULL) {
	ncurifree(duri);
    }
    return 0;
}

void
ncurifree(NCURI* duri)
{
    if(duri == NULL) return;
    if(duri->uri != NULL) {free(duri->uri);}
    if(duri->params != NULL) {free(duri->params);}
    if(duri->paramlist != NULL) ncparamfree(duri->paramlist);
    if(duri->strings != NULL) {free(duri->strings);}
    if(duri->constraint != NULL) {free(duri->constraint);}
    if(duri->projection != NULL) {free(duri->projection);}
    if(duri->selection != NULL) {free(duri->selection);}
    free(duri);
}

/* Replace the constraints */
void
ncurisetconstraints(NCURI* duri,const char* constraints)
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
	nclshift1(duri->constraint);

    p = duri->constraint;
    proj = (char*) p;
    select = strchr(proj,'&');
    if(select != NULL) {
        size_t plen = (size_t)(select - proj);
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


/* Construct a complete NC URI.
   Optionally with the constraints.
   Optionally with the user parameters.
   Caller frees returned string.
   Optionally encode the pieces.
*/

char*
ncuribuild(NCURI* duri, const char* prefix, const char* suffix, int flags)
{
    size_t len = 0;
    char* newuri;
    char* tmpfile;
    char* tmpsuffix;
    char* tmpquery;
    size_t nparams = 0;
    size_t paramslen = 0;

    /* if both are specified, prefix has priority */
    int withsuffixparams = ((flags&NCURISUFFIXPARAMS)!=0
				&& duri->params != NULL);
    int withprefixparams = ((flags&NCURIPREFIXPARAMS)!=0
				&& duri->params != NULL);
    int withuserpwd = ((flags&NCURIUSERPWD)!=0
	               && duri->user != NULL && duri->password != NULL);
    int withconstraints = ((flags&NCURICONSTRAINTS)!=0
	                   && duri->constraint != NULL);
#ifdef NEWESCAPE
    int encode = (flags&NCURIENCODE);
#else
    int encode = 0;
#endif

    if(prefix != NULL) len += NILLEN(prefix);
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
	tmpfile = ncuriencode(tmpfile,fileallow);
    len += (NILLEN(tmpfile));

    if(suffix != NULL) {
        tmpsuffix = (char*)suffix;
        if(encode)
	    tmpsuffix = ncuriencode(tmpsuffix,fileallow);
        len += (NILLEN(tmpsuffix));
    }

    if(withconstraints) {
	tmpquery = duri->constraint;
        if(encode)
	    tmpquery = ncuriencode(tmpquery,queryallow);
        len += (NILLEN("?")+NILLEN(tmpquery));
    }

    if(withprefixparams || withsuffixparams) {
	char** p;
	if(duri->paramlist == NULL)
	    if(!ncuridecodeparams(duri))
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
	ncappendparams(newuri,duri->paramlist);
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
	ncappendparams(newuri,duri->paramlist);
    }
    return newuri;
}

static void
ncappendparams(char* newuri, char** p)
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
They may occur either at the front, or suffixed after
a trailing # character. After processing, the list is
converted to an ampersand separated list of the combination
of prefix and suffix parameters.

After the url is parsed, the parameter list
is converted to an ampersand separated list with all
whitespace removed.
In any case, each parameter in turn is assumed to be a
of the form <name>=<value> or [<name>].
e.g. [x=y][z][a=b][w].  If the same parameter is specified more
than once, then the first occurrence is used; this is so
that is possible to forcibly override user specified
parameters by prefixing.  IMPORTANT: client parameter string
is assumed to have blanks compressed out.  Returns 1 if parse
suceeded, 0 otherwise; */

int
ncuridecodeparams(NCURI* ncuri)
{
    char* cp;
    int i,c;
    size_t nparams;
    char* params;
    char** plist;

    if(ncuri == NULL) return 0;
    if(ncuri->params == NULL) return 1;

    params = strdup(ncuri->params); /* so we can modify */

    /* Pass 1 to break string into pieces at the ampersands
       and count # of pairs */
    nparams=0;
    for(cp=params;(c=*cp);cp++) {
	if(c == '&') {*cp = EOFCHAR; nparams++;}
    }
    nparams++; /* for last one */

    /* plist is an env style list */
    plist = (char**)calloc(1,sizeof(char*)*(2*nparams+1)); /* +1 for null termination */
    if(plist == NULL)
	return 0;

    /* Break up each param into a (name,value) pair*/
    /* and insert into the param list */
    /* parameters of the form name name= are converted to name=""*/
    for(cp=params,i=0;i<nparams;i++) {
	char* next = cp+strlen(cp)+1; /* save ptr to next pair*/
	char* vp;
	/*break up the ith param*/
	vp = strchr(cp,'=');
	if(vp != NULL) {*vp = EOFCHAR; vp++;} else {vp = "";}
	plist[2*i] = nulldup(cp);	
	plist[2*i+1] = nulldup(vp);
	cp = next;
    }
    plist[2*nparams] = NULL;
    free(params);
    if(ncuri->paramlist != NULL)
	ncparamfree(ncuri->paramlist);
    ncuri->paramlist = plist;
    return 1;
}

int
ncurilookup(NCURI* uri, const char* key, const char** resultp)
{
    int i;
    char* value = NULL;
    if(uri == NULL || key == NULL || uri->params == NULL) return 0;
    if(uri->paramlist == NULL) {
	i = ncuridecodeparams(uri);
	if(!i) return 0;
    }
    i = ncfind(uri->paramlist,key);
    if(i < 0)
	return 0;
    value = uri->paramlist[(2*i)+1];
    if(resultp) *resultp = value;
    return 1;
}

int
ncurisetparams(NCURI* uri, const char* newparams)
{
    if(uri == NULL) return 0;
    if(uri->paramlist != NULL) ncparamfree(uri->paramlist);
    uri->paramlist = NULL;
    if(uri->params != NULL) free(uri->params);
    uri->params = nulldup(newparams);
    return 1;
}

/* Internal version of lookup; returns the paired index of the key */
static int
ncfind(char** params, const char* key)
{
    int i;
    char** p;
    for(i=0,p=params;*p;p+=2,i++) {
	if(strcmp(key,*p)==0) return i;
    }
    return -1;
}

static void
ncparamfree(char** params)
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
nclocate(char* p, const char* charlist)
{
    for(;*p;p++) {
	if(strchr(charlist,*p) != NULL)
	    return p;
    }
    return NULL;
}


/* Shift every char starting at p 1 place to the left */
static void
nclshift1(char* p)
{
    if(p != NULL && *p != EOFCHAR) {
	char* q = p++;
	while((*q++=*p++));
    }
}

/* Shift every char starting at p 1 place to the right */
static void
ncrshift1(char* p)
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
toHex(unsigned int b, char hex[2])
{
    hex[0] = hexchars[(b >> 4) & 0xff];
    hex[1] = hexchars[(b) & 0xff];
}


static int
fromHex(int c)
{
    if(c >= '0' && c <= '9') return (int) (c - '0');
    if(c >= 'a' && c <= 'f') return (int) (10 + (c - 'a'));
    if(c >= 'A' && c <= 'F') return (int) (10 + (c - 'A'));
    return 0;
}


/* Return a string representing encoding of input; caller must free;
   watch out: will encode whole string, so watch what you give it.
   Allowable argument specifies characters that do not need escaping.
 */

char*
ncuriencode(char* s, char* allowable)
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
            if(c2) {*outptr++ = (char)c;}
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
ncuridecode(char* s)
{
    return ncuridecodeonly(s,NULL);
}

/* Return a string representing decoding of input only for specified
   characters;  caller must free
*/
char*
ncuridecodeonly(char* s, char* only)
{
    size_t slen;
    char* decoded;
    char* outptr;
    char* inptr;
    unsigned int c;
    
    if (s == NULL) return NULL;

    slen = strlen(s);
    decoded = (char*)malloc(slen+1); /* Should be max we need */

    outptr = decoded;
    inptr = s;
    while((c = (unsigned int)*inptr++)) {
	if(c == '+' && only != NULL && strchr(only,'+') != NULL)
	    *outptr++ = ' ';
	else if(c == '%') {
            /* try to pull two hex more characters */
	    if(inptr[0] != EOFCHAR && inptr[1] != EOFCHAR
		&& strchr(hexchars,inptr[0]) != NULL
		&& strchr(hexchars,inptr[1]) != NULL) {
		/* test conversion */
		int xc = (fromHex(inptr[0]) << 4) | (fromHex(inptr[1]));
		if(only == NULL || strchr(only,xc) != NULL) {
		    inptr += 2; /* decode it */
		    c = (unsigned int)xc;
                }
            }
        }
        *outptr++ = (char)c;
    }
    *outptr = EOFCHAR;
    return decoded;
}


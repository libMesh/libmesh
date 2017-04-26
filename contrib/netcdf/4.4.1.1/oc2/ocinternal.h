/* Copyright 2009, UCAR/Unidata and OPeNDAP, Inc.
   See the COPYRIGHT file for more information. */

#ifndef OCINTERNAL_H
#define OCINTERNAL_H

#include "config.h"


#if defined(_WIN32) || defined(_WIN64)
#include <malloc.h>
#endif

/* Required for getcwd, other functions. */
#ifdef _MSC_VER
#include <direct.h>
#define getcwd _getcwd
#endif

#ifdef _AIX
#include <netinet/in.h>
#endif

#include <stdlib.h>
#include <assert.h>
#include <string.h>
#ifdef HAVE_STRINGS_H
#include <strings.h>
#endif
#ifdef HAVE_STDARG_H
#include <stdarg.h>
#endif
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>
#endif

#define CURL_DISABLE_TYPECHECK 1
#include <curl/curl.h>

#include "oclist.h"
#include "ocbytes.h"
#include "ocuri.h"

#ifndef HAVE_STRNDUP
/* Not all systems have strndup, so provide one*/
#define strndup ocstrndup
#endif

#define OCCACHEPOS

#ifndef HAVE_STRNDUP
/* Not all systems have strndup, so provide one*/
#define strndup ocstrndup
#endif

#define OCPATHMAX 8192

#ifndef nullfree
#define nullfree(x) {if((x)!=NULL) free(x);}
#endif

/* Forwards */
typedef struct OCstate OCstate;
typedef struct OCnode OCnode;
typedef struct OCdata OCdata;
struct OCTriplestore;

/* Define the internal node classification values */
#define OC_None  0
#define OC_State 1
#define OC_Node  2
#define OC_Data  3

/* Define a magic number to mark externally visible oc objects */
#define OCMAGIC ((unsigned int)0x0c0c0c0c) /*clever, huh*/

/* Max rc file line size */
#define MAXRCLINESIZE 4096

/* Max number of triples from an rc file */
#define MAXRCLINES 4096

/* Define a struct that all oc objects must start with */
/* OCheader must be defined here to make it available in other headers */
typedef struct OCheader {
    unsigned int magic;
    unsigned int occlass;
} OCheader;

#include "oc.h"
#include "ocx.h"
#include "ocnode.h"
#include "ocdata.h"
#include "occonstraints.h"
#include "ocutil.h"
#include "oclog.h"
#include "xxdr.h"
#include "ocdatatypes.h"
#include "occompile.h"

#ifndef nulldup
#define nulldup(s) (s==NULL?NULL:strdup(s))
#endif

/*
 * Macros for dealing with flag bits.
 */
#define fset(t,f)       ((t) |= (f))
#define fclr(t,f)       ((t) &= ~(f))
#define fisset(t,f)     ((t) & (f))

#define nullstring(s) (s==NULL?"(null)":s)
#define PATHSEPARATOR "."

#define OCDIR "oc"

/* Define infinity for memory size */
#if SIZEOF_SIZE_T == 4 
#define OCINFINITY ((size_t)0xffffffff)
#else
#define OCINFINITY ((size_t)0xffffffffffffffff)
#endif

/* Extend the OCdxd type for internal use */
#define OCVER 3

/* Default initial memory packet size */
#define DFALTPACKETSIZE 0x20000 /*approximately 100k bytes*/

/* Default maximum memory packet size */
#define DFALTMAXPACKETSIZE 0x3000000 /*approximately 50M bytes*/

/* Default user agent string (will have version appended)*/
#define DFALTUSERAGENT "oc"

/* Hold known curl flags */

enum OCCURLFLAGTYPE {CF_UNKNOWN=0,CF_OTHER=1,CF_STRING=2,CF_LONG=3};
struct OCCURLFLAG {
    const char* name;
    int flag;
    int value;
    enum OCCURLFLAGTYPE type;
};

struct OCTriplestore {
    int ntriples;
    struct OCTriple {
        char host[MAXRCLINESIZE]; /* includes port if specified */
        char key[MAXRCLINESIZE];
        char value[MAXRCLINESIZE];
   } triples[MAXRCLINES];
};

/* Collect global state info in one place */
extern struct OCGLOBALSTATE {
    int initialized;
    struct {
        int proto_file;
        int proto_https;
    } curl;
    char* tempdir; /* track a usable temp dir */
    char* home; /* track $HOME for use in creating $HOME/.oc dir */
    struct {
	int ignore; /* if 1, then do not use any rc file */
	int loaded;
        struct OCTriplestore daprc; /* the rc file triple store fields*/
        char* rcfile; /* specified rcfile; overrides anything else */
    } rc;
} ocglobalstate;

/*! Specifies the OCstate = non-opaque version of OClink */
struct OCstate {
    OCheader header; /* class=OC_State */
    OClist* trees; /* list<OCNODE*> ; all root objects */
    OCURI* uri; /* base URI*/
    OCbytes* packet; /* shared by all trees during construction */
    struct OCerrdata {/* Hold info for an error return from server */
	char* code;
	char* message;
	long  httpcode;
	char  curlerrorbuf[CURL_ERROR_SIZE]; /* CURLOPT_ERRORBUFFER*/
    } error;
    CURL* curl; /* curl handle*/
    char curlerror[CURL_ERROR_SIZE];
    struct OCcurlflags {
        int proto_file; /* Is file: supported? */
        int proto_https; /* is https: supported? */
	int compress; /*CURLOPT_ENCODING*/
	int verbose; /*CURLOPT_ENCODING*/
	int timeout; /*CURLOPT_TIMEOUT*/
	int maxredirs; /*CURLOPT_MAXREDIRS*/
	char* useragent; /*CURLOPT_USERAGENT*/
	/* track which of these are created by oc */
#define COOKIECREATED 1
#define NETRCCREATED 2
	int createdflags;
	char* cookiejar; /*CURLOPT_COOKIEJAR,CURLOPT_COOKIEFILE*/
	char* netrc; /*CURLOPT_NETRC,CURLOPT_NETRC_FILE*/
    } curlflags;
    struct OCSSL {
	int   verifypeer; /* CURLOPT_SSL_VERIFYPEER;
                             do not do this when cert might be self-signed
                             or temporarily incorrect */
	int   verifyhost; /* CURLOPT_SSL_VERIFYHOST; for client-side verification */
        char* certificate; /*CURLOPT_SSLCERT*/
	char* key; /*CURLOPT_SSLKEY*/
	char* keypasswd; /*CURLOPT_SSLKEYPASSWD*/
        char* cainfo; /* CURLOPT_CAINFO; certificate authority */
	char* capath;  /*CURLOPT_CAPATH*/
    } ssl;
    struct OCproxy {
	char *host; /*CURLOPT_PROXY*/
	int port; /*CURLOPT_PROXYPORT*/
	char* userpwd; /*CURLOPT_PROXYUSERPWD*/
    } proxy;
    struct OCcredentials {
	char *userpwd; /*CURLOPT_USERPWD*/
    } creds;
    void* usercurldata;
    long ddslastmodified;
    long datalastmodified;
};

/*! OCtree holds extra state info about trees */

typedef struct OCtree
{
    OCdxd  dxdclass;
    char* constraint;
    char* text;
    struct OCnode* root; /* cross link */
    struct OCstate* state; /* cross link */
    OClist* nodes; /* all nodes in tree*/
    /* when dxdclass == OCDATADDS */
    struct {
	char*   memory;   /* allocated memory if OC_ONDISK is not set */
        char*   filename; /* If OC_ONDISK is set */
        FILE*   file;
        off_t   datasize; /* xdr size on disk or in memory */
        off_t   bod;      /* offset of the beginning of packet data */
        off_t   ddslen;   /* length of ddslen (assert(ddslen <= bod)) */
        XXDR*   xdrs;		/* access either memory or file */
        OCdata* data;
    } data;
} OCtree;

/* (Almost) All shared procedure definitions are kept here
   except for: ocdebug.h ocutil.h
   The true external interface is defined in oc.h
*/

#if 0
/* Location: ceparselex.c*/
extern int cedebug;
extern OClist* CEparse(OCstate*,char* input);
#endif

extern OCerror ocopen(OCstate** statep, const char* url);
extern void occlose(OCstate* state);
extern OCerror ocfetch(OCstate*, const char*, OCdxd, OCflags, OCnode**);
extern int oc_network_order;
extern int oc_invert_xdr_double;
extern OCerror ocinternalinitialize(void);

extern OCerror ocupdatelastmodifieddata(OCstate* state);

extern OCerror ocset_useragent(OCstate* state, const char* agent);
extern OCerror ocset_netrc(OCstate* state, const char* path);

/* From ocrc.c */
extern OCerror ocrc_load(); /* find, read, and compile */
extern OCerror ocrc_process(OCstate* state); /* extract relevant triples */
extern char* ocrc_lookup(char* key, char* url);
extern struct OCTriple* ocrc_triple_iterate(char* key, char* url, struct OCTriple* prevp);
extern int ocparseproxy(OCstate* state, char* v);

#endif /*COMMON_H*/

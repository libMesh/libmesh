/* Copyright 2009, UCAR/Unidata and OPeNDAP, Inc.
   See the COPYRIGHT file for more information. */

#ifndef OCINTERNAL_H
#define OCINTERNAL_H

#include "config.h"

#ifdef _AIX
#include <netinet/in.h>
#endif

#include <stdlib.h>
#include <assert.h>
#include <string.h>
#ifndef WIN32
#include <strings.h>
#endif
#include <stdarg.h>
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

#define OCCACHEPOS

#include "oc.h"
#include "ocdatatypes.h"
#include "occonstraints.h"
#include "ocnode.h"
#include "ocutil.h"
#include "oclog.h"
#include "xxdr.h"
#include "ocdata.h"

#ifndef nulldup
#define nulldup(s) (s==NULL?NULL:strdup(s))
#endif

#define nullstring(s) (s==NULL?"(null)":s)
#define PATHSEPARATOR "."

/* Default initial memory packet size */
#define DFALTPACKETSIZE 0x20000 /*approximately 100k bytes*/

/* Default maximum memory packet size */
#define DFALTMAXPACKETSIZE 0x3000000 /*approximately 50M bytes*/

/* Extend the OCdxd type */
#define OCVER 3

/* Define a magic number to mark externally visible oc objects */
#define OCMAGIC ((unsigned int)0x0c0c0c0c) /*clever, huh?*/

/*! Specifies the OCstate. */
typedef struct OCstate
{
    unsigned int magic; /* Mark each structure type */
    CURL* curl; /* curl handle*/
    OClist* trees; /* list<OCnode*> ; all root objects */
    OCURI* uri; /* base URI*/
    OCbytes* packet; /* shared by all trees during construction */
    /* OCContent information */
    struct OCcontent* contentlist;
    struct OCerrdata {/* Hold info for an error return from server */
	char* code;
	char* message;
	long  httpcode;
	char  curlerrorbuf[CURL_ERROR_SIZE]; /* to get curl error message */
    } error;
    /* Store .rc file info */
    struct OCcurlflags {
	int compress;
	int verbose;
	int timeout;
	int followlocation;
	int maxredirs;
	char* useragent;
	char* cookiejar;
	char* cookiefile;
    } curlflags;
    struct OCSSL {
	int   validate;
        char* certificate;
	char* key;
	char* keypasswd;
        char* cainfo; /* certificate authority */
	char* capath; 
	int   verifypeer;
    } ssl;
    struct OCproxy {
	char *host;
	int port;
    } proxy;
    struct OCcredentials {
	char *username;
	char *password;
    } creds;
    long ddslastmodified;
    long datalastmodified;
} OCstate;


/*! Specifies all the info about a particular DAP tree
    i.e. DAS, DDS, or DATADDS as obtained from a fetch response
    This is associated with the root object.
*/
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
	char* memory;           /* allocated memory if OC_INMEMORY is set */
        char* filename;
        FILE* file;
        unsigned long datasize; /* xdr size on disk or in memory */
        unsigned long bod;      /* offset of the beginning of packet data */
        unsigned long ddslen;   /* length of ddslen (assert(ddslen <= bod)) */
        XXDR* xdrs;		/* access either memory or file */
    } data;
} OCtree;

/* (Almost) All shared procedure definitions are kept here
   except for: ocdebug.h ocutil.h
   The true external interface is defined in oc.h
*/

/* Location: ocnode.c */
extern OCnode* ocmakenode(char* name, OCtype ptype, OCnode* root);
extern void occollectpathtonode(OCnode* node, OClist* path);
extern void occomputefullnames(OCnode* root);
extern void occomputesemantics(OClist*);
extern void ocaddattribute(OCattribute* attr, OCnode* parent);
extern OCattribute* ocmakeattribute(char* name, OCtype ptype, OClist* values);
extern size_t ocsetsize(OCnode* node);
extern OCerror occorrelate(OCnode*,OCnode*);
extern OCerror occomputeskipdata(OCstate*, OCnode*);
extern void ocmarkcacheable(OCstate* state, OCnode* ddsroot);

/* Location: dapparselex.c*/
extern int dapdebug;
extern OCerror DAPparse(OCstate*, struct OCtree*, char*);
extern char* dimnameanon(char* basename, unsigned int index);

/* Location: ceparselex.c*/
extern int cedebug;
extern OClist* CEparse(OCstate*,char* input);

/* Location: ocinternal.c*/
extern OCerror ocopen(OCstate** statep, const char* url);
extern void occlose(OCstate* state);

extern OCerror ocfetchf(OCstate*, const char*, OCdxd, OCflags, OCnode**);

/* Location: ocinternal.c */
extern int oc_network_order;
extern int oc_invert_xdr_double;
extern int ocinternalinitialize(void);

/* Location: ocnode.c */
extern void ocfreetree(OCtree* tree);
extern void ocfreeroot(OCnode* root);
extern void ocfreenodes(OClist*);

extern void ocddsclear(struct OCstate*);
extern void ocdasclear(struct OCstate*);
extern void ocdataddsclear(struct OCstate*);
extern void* oclinearize(OCtype etype, unsigned int, char**);

/* Merge DAS with DDS or DATADDS*/
extern int ocddsdasmerge(struct OCstate*, OCnode* das, OCnode* dds);

extern OCerror ocupdatelastmodifieddata(OCstate* state);

extern int ocinternalinitialize(void);


extern OCerror ocsetrcfile(char* rcfile);

/* Global stateflags */
extern int oc_curl_file_supported;
extern int oc_curl_https_supported;

#endif /*COMMON_H*/

#include "ncdispatch.h"
#include "ncuri.h"

#define MAXSERVERURL 4096

/* Define vectors of zeros and ones for use with various nc_get_varX function*/
size_t nc_sizevector0[NC_MAX_VAR_DIMS];
size_t nc_sizevector1[NC_MAX_VAR_DIMS];
ptrdiff_t nc_ptrdiffvector1[NC_MAX_VAR_DIMS];
size_t NC_coord_zero[NC_MAX_VAR_DIMS];
size_t NC_coord_one[NC_MAX_VAR_DIMS];

/* Define the known protocols and their manipulations */
static struct NCPROTOCOLLIST {
    char* protocol;
    char* substitute;
    int   model;
} ncprotolist[] = {
    {"http",NULL,0},
    {"https",NULL,0},
    {"file",NULL,NC_FORMATX_DAP2},
    {"dods","http",NC_FORMATX_DAP2},
    {"dodss","https",NC_FORMATX_DAP2},
    {NULL,NULL,0} /* Terminate search */
};

/* Define the default servers to ping in order;
   make the order attempt to optimize
   against future changes.
*/
static const char* default_servers[] = {
"http://remotetest.unidata.ucar.edu",
NULL
};

/*
static nc_type longtype = (sizeof(long) == sizeof(int)?NC_INT:NC_INT64);
static nc_type ulongtype = (sizeof(unsigned long) == sizeof(unsigned int)?NC_UINT:NC_UINT64);
*/

/* Allow dispatch to do general initialization and finalization */
int
NCDISPATCH_initialize(void)
{
    int status = NC_NOERR;
    int i;
    for(i=0;i<NC_MAX_VAR_DIMS;i++) {
	nc_sizevector0[i] = 0;
        nc_sizevector1[i] = 1;
        nc_ptrdiffvector1[i] = 1;
    }
    for(i=0;i<NC_MAX_VAR_DIMS;i++) {
	NC_coord_one[i] = 1;
	NC_coord_zero[i] = 0;
    }
    return status;
}

int
NCDISPATCH_finalize(void)
{
    int status = NC_NOERR;
    int i;
    return status;
}

/* search list of servers and return first that succeeds when
   concatenated with the specified path part.
   Search list can be prefixed by the second argument.
*/
char*
NC_findtestserver(const char* path, const char** servers)
{
#ifdef USE_DAP
#ifdef ENABLE_DAP_REMOTE_TESTS
    /* NCDAP_ping is defined in libdap2/ncdap.c */
    const char** svc;
    int stat;
    char* url = (char*)malloc(MAXSERVERURL);

    if(path == NULL) path = "";
    if(strlen(path) > 0 && path[0] == '/')
	path++;

    if(servers != NULL) {
        for(svc=servers;*svc != NULL;svc++) {
            snprintf(url,MAXSERVERURL,"%s/%s",*svc,path);
            stat = NCDAP_ping(url);
            if(stat == NC_NOERR)
                return url;
        }
    }
    /* not found in user supplied list; try defaults */
    for(svc=default_servers;*svc != NULL;svc++) {
        snprintf(url,MAXSERVERURL,"%s/%s",*svc,path);
        stat = NCDAP_ping(url);
        if(stat == NC_NOERR)
            return url;
    }
    if(url) free(url);
#endif
#endif
    return NULL;
}


/* return 1 if path looks like a url; 0 otherwise */
int
NC_testurl(const char* path)
{
    int isurl = 0;
    NCURI* tmpurl = NULL;
    char* p;

    if(path == NULL) return 0;

    /* find leading non-blank */
    for(p=(char*)path;*p;p++) {if(*p != ' ') break;}

    /* Do some initial checking to see if this looks like a file path */
    if(*p == '/') return 0; /* probably an absolute file path */

    /* Ok, try to parse as a url */
    if(ncuriparse(path,&tmpurl)) {
	/* Do some extra testing to make sure this really is a url */
        /* Look for a knownprotocol */
        struct NCPROTOCOLLIST* protolist;
        for(protolist=ncprotolist;protolist->protocol;protolist++) {
	    if(strcmp(tmpurl->protocol,protolist->protocol) == 0) {
	        isurl=1;
		break;
	    }
	}
	ncurifree(tmpurl);
	return isurl;
    }
    return 0;
}

/*
Return an NC_FORMATX_... value.
Assumes that the path is known to be a url
*/

int
NC_urlmodel(const char* path)
{
    int model = 0;
    NCURI* tmpurl = NULL;
    struct NCPROTOCOLLIST* protolist;

    model = NC_FORMATX_DAP2;
    return model;
}

#ifdef OBSOLETE
/* Override dispatch table management */
static NC_Dispatch* NC_dispatch_override = NULL;

/* Override dispatch table management */
NC_Dispatch*
NC_get_dispatch_override(void) {
    return NC_dispatch_override;
}

void NC_set_dispatch_override(NC_Dispatch* d)
{
    NC_dispatch_override = d;
}
#endif

/* OBSOLETE
   Overlay by treating the tables as arrays of void*.
   Overlay rules are:
        overlay    base    merge
        -------    ----    -----
          null     null     null
          null      y        y
           x       null      x
           x        y        x
*/

#ifdef OBSOLETE
int
NC_dispatch_overlay(const NC_Dispatch* overlay, const NC_Dispatch* base, NC_Dispatch* merge)
{
    void** voverlay = (void**)overlay;
    void** vmerge;
    int i;
    size_t count = sizeof(NC_Dispatch) / sizeof(void*);
    /* dispatch table must be exact multiple of sizeof(void*) */
    assert(count * sizeof(void*) == sizeof(NC_Dispatch));
    *merge = *base;
    vmerge = (void**)merge;
    for(i=0;i<count;i++) {
        if(voverlay[i] == NULL) continue;
        vmerge[i] = voverlay[i];
    }
    /* Finally, the merge model should always be the overlay model */
    merge->model = overlay->model;
    return NC_NOERR;
}
#endif

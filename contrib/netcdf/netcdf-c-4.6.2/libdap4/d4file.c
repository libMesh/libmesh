/*********************************************************************
 *   Copyright 2016, UCAR/Unidata
 *   See netcdf/COPYRIGHT file for copying and redistribution conditions.
 *********************************************************************/

#include "ncdispatch.h"
#include "ncd4dispatch.h"
#include "d4includes.h"
#include "d4read.h"
#include "d4curlfunctions.h"

#ifdef _MSC_VER
#include <process.h>
#include <direct.h>
#endif

#ifdef HAVE_SYS_STAT_H
#include <sys/stat.h>
#endif

/**************************************************/
/* Forward */

static void applyclientmetacontrols(NCD4meta* meta);
static int constrainable(NCURI*);
static void freeCurl(NCD4curl*);
static void freeInfo(NCD4INFO*);
static int paramcheck(NCD4INFO*, const char* key, const char* subkey);
static const char* getparam(NCD4INFO* info, const char* key);
static int set_curl_properties(NCD4INFO*);

/**************************************************/
/* Constants */

static const char* checkseps = "+,:;";

/**************************************************/
int
NCD4_open(const char * path, int mode,
          int basepe, size_t *chunksizehintp,
          void *mpidata, NC_Dispatch *dispatch, NC *nc)
{
    int ret = NC_NOERR;
    NCD4INFO* d4info = NULL;
    const char* value;
    NCD4meta* meta;

    if(path == NULL)
	return THROW(NC_EDAPURL);

    assert(dispatch != NULL);

    /* Setup our NC and NCDAPCOMMON state*/

    d4info = (NCD4INFO*)calloc(1,sizeof(NCD4INFO));
    if(d4info == NULL) {ret = NC_ENOMEM; goto done;}

    nc->dispatchdata = d4info;
    nc->int_ncid = nc__pseudofd(); /* create a unique id */
    d4info->controller = (NC*)nc;

    /* Parse url and params */
    if(ncuriparse(nc->path,&d4info->uri) != NCU_OK)
	{ret = NC_EDAPURL; goto done;}

    /* Load auth info from rc file */
    if((ret = NC_authsetup(&d4info->auth, d4info->uri)))
	goto done;
    NCD4_curl_protocols(d4info);

    if(!constrainable(d4info->uri))
	SETFLAG(d4info->controls.flags,NCF_UNCONSTRAINABLE);

    /* fail if we are unconstrainable but have constraints */
    if(FLAGSET(d4info->controls.flags,NCF_UNCONSTRAINABLE)) {
	if(d4info->uri->query != NULL) {
	    nclog(NCLOGWARN,"Attempt to constrain an unconstrainable data source: %s",
		   d4info->uri->query);
	    ret = THROW(NC_EDAPCONSTRAINT);
	    goto done;
	}
    }

    /* process control client parameters */
    NCD4_applyclientparamcontrols(d4info);

    /* Use libsrc4 code (netcdf-4) for storing metadata */
    {
	char tmpname[NC_MAX_NAME];

        /* Create fake file name: exact name must be unique,
           but is otherwise irrelevant because we are using NC_DISKLESS
        */
	if(strlen(d4info->controls.substratename) > 0)
            snprintf(tmpname,sizeof(tmpname),"%s",d4info->controls.substratename);
	else
            snprintf(tmpname,sizeof(tmpname),"tmp_%d",nc->int_ncid);

        /* Now, use the file to create the hidden substrate netcdf file.
	   We want this hidden file to always be NC_NETCDF4, so we need to
           force default format temporarily in case user changed it.
	   Since diskless is enabled, create file in-memory.
	*/
	{
	    int new = NC_NETCDF4;
	    int old = 0;
	    int ncid = 0;
	    int ncflags = NC_NETCDF4|NC_CLOBBER;
	    ncflags |= NC_DISKLESS;
	    if(FLAGSET(d4info->controls.debugflags,NCF_DEBUG_COPY)) {
		/* Cause data to be dumped to real file */
		ncflags |= NC_WRITE;
		ncflags &= ~(NC_DISKLESS); /* use real file */
	    }
	    nc_set_default_format(new,&old); /* save and change */
            ret = nc_create(tmpname,ncflags,&ncid);
	    nc_set_default_format(old,&new); /* restore */
	    d4info->substrate.realfile = ((ncflags & NC_DISKLESS) == 0);
	    d4info->substrate.filename = strdup(tmpname);
	    if(d4info->substrate.filename == NULL) ret = NC_ENOMEM;
	    d4info->substrate.nc4id = ncid;
	}
        if(ret != NC_NOERR) goto done;
	/* Avoid fill */
	nc_set_fill(getnc4id(nc),NC_NOFILL,NULL);
    }

    /* Turn on logging; only do this after oc_open*/
    if((value = ncurilookup(d4info->uri,"log")) != NULL) {
	ncloginit();
        if(nclogopen(value))
	    ncsetlogging(1);
	ncloginit();
        if(nclogopen(value))
	    ncsetlogging(1);
    }

    /* Setup a curl connection */
    {
        CURL* curl = NULL; /* curl handle*/
	d4info->curl = (NCD4curl*)calloc(1,sizeof(NCD4curl));
	if(d4info->curl == NULL)
	    {ret = NC_ENOMEM; goto done;}
	/* create the connection */
        if((ret=NCD4_curlopen(&curl))!= NC_NOERR) goto done;
	d4info->curl->curl = curl;
        /* Load misc rc properties */
        NCD4_get_rcproperties(d4info);
        if((ret=set_curl_properties(d4info))!= NC_NOERR) goto done;	
        /* Set the one-time curl flags */
        if((ret=NCD4_set_flags_perlink(d4info))!= NC_NOERR) goto done;
#if 1 /* temporarily make per-link */
        if((ret=NCD4_set_flags_perfetch(d4info))!= NC_NOERR) goto done;
#endif
    }

    d4info->curl->packet = ncbytesnew();
    ncbytessetalloc(d4info->curl->packet,DFALTPACKETSIZE); /*initial reasonable size*/

    /* fetch the dmr + data*/
    {
	int inmem = FLAGSET(d4info->controls.flags,NCF_ONDISK) ? 0 : 1;
        if((ret = NCD4_readDAP(d4info,inmem))) goto done;
    }

    /* if the url goes astray to a random web page, then try to just dump it */
    {
	char* response = ncbytescontents(d4info->curl->packet);
	size_t responselen = ncbyteslength(d4info->curl->packet);

        /* Apply some heuristics to see what we have.
           The leading byte will have the chunk flags, which should
           be less than 0x0f (for now). However, it will not be zero if
           the data was little-endian
	*/
        if(responselen == 0 || response[0] >= ' ') {
	    /* does not look like a chunk, so probable server failure */
	    if(responselen == 0)
	        nclog(NCLOGERR,"Empty DAP4 response");
	    else {/* probable html response */
		nclog(NCLOGERR,"Unexpected DAP response:");
		nclog(NCLOGERR,"==============================");
		nclogtext(NCLOGERR,response);
		nclog(NCLOGERR,"==============================\n");
	    }
	    ret = NC_EDAPSVC;
  	    fflush(stderr);
	    goto done;
	}
    }

    /* Build the meta data */
    if((d4info->substrate.metadata=NCD4_newmeta(ncbyteslength(d4info->curl->packet),
        ncbytescontents(d4info->curl->packet)))==NULL)
	{ret = NC_ENOMEM; goto done;}
    meta = d4info->substrate.metadata;
    meta->controller = d4info;
    meta->ncid = getnc4id(nc); /* Transfer netcdf ncid */

    /* process meta control parameters */
    applyclientmetacontrols(meta);

    /* Infer the mode */
    if((ret=NCD4_infermode(meta))) goto done;

    if((ret=NCD4_dechunk(meta))) goto done;

#ifdef D4DUMPDMR
  {
    fprintf(stderr,"=============\n");
    fputs(d4info->substrate.metadata->serial.dmr,stderr);
    fprintf(stderr,"\n=============\n");
    fflush(stderr);
  }
#endif

    if((ret = NCD4_parse(d4info->substrate.metadata))) goto done;
#ifdef D4DEBUGMETA
  {
    fprintf(stderr,"\n/////////////\n");
    NCbytes* buf = ncbytesnew();
    NCD4_print(d4info->substrate.metadata,buf);
    ncbytesnull(buf);
    fputs(ncbytescontents(buf),stderr);
    ncbytesfree(buf);
    fprintf(stderr,"\n/////////////\n");
    fflush(stderr);
  }
#endif
    if((ret = NCD4_metabuild(d4info->substrate.metadata,d4info->substrate.metadata->ncid))) goto done;
    if(ret != NC_NOERR && ret != NC_EVARSIZE) goto done;
    if((ret = NCD4_processdata(d4info->substrate.metadata))) goto done;

    return THROW(ret);

done:
    if(ret) {
	freeInfo(d4info);
        nc->dispatchdata = NULL;
    }
    return THROW(ret);
}

int
NCD4_close(int ncid, void* ignore)
{
    int ret = NC_NOERR;
    NC* nc;
    NCD4INFO* d4info;
    int substrateid;

    ret = NC_check_id(ncid, (NC**)&nc);
    if(ret != NC_NOERR) goto done;
    d4info = (NCD4INFO*)nc->dispatchdata;
    substrateid = makenc4id(nc,ncid);

    /* We call abort rather than close to avoid trying to write anything,
       except if we are debugging
     */
    if(FLAGSET(d4info->controls.debugflags,NCF_DEBUG_COPY)) {
        /* Dump the data into the substrate */
	if((ret = NCD4_debugcopy(d4info)))
	    goto done;
        ret = nc_close(substrateid);
    } else {
        ret = nc_abort(substrateid);
    }

    freeInfo(d4info);

done:
    return THROW(ret);
}

int
NCD4_abort(int ncid)
{
    return NCD4_close(ncid,NULL);
}

/**************************************************/

/* Reclaim an NCD4INFO instance */
static void
freeInfo(NCD4INFO* d4info)
{
    if(d4info == NULL) return;
    d4info->controller = NULL; /* break link */
    nullfree(d4info->rawurltext);
    nullfree(d4info->urltext);
    ncurifree(d4info->uri);
    freeCurl(d4info->curl);
    nullfree(d4info->data.memory);
    nullfree(d4info->data.ondiskfilename);
    if(d4info->data.ondiskfile != NULL)
	fclose(d4info->data.ondiskfile);
    nullfree(d4info->fileproto.filename);
    if(d4info->substrate.realfile
	&& !FLAGSET(d4info->controls.debugflags,NCF_DEBUG_COPY)) {
	/* We used real file, so we need to delete the temp file
           unless we are debugging.
	   Assume caller has done nc_close|nc_abort on the ncid.
           Note that in theory, this should not be necessary since
           AFAIK the substrate file is still in def mode, and
           when aborted, it should be deleted. But that is not working
           for some reason, so we delete it ourselves.
	*/
#if 0
	if(d4info->substrate.filename != NULL) {
	    unlink(d4info->substrate.filename);
	}
#endif
    }
    nullfree(d4info->substrate.filename); /* always reclaim */
    NCD4_reclaimMeta(d4info->substrate.metadata);
    NC_authclear(&d4info->auth);
    nclistfree(d4info->blobs);
    free(d4info);    
}

static void
freeCurl(NCD4curl* curl)
{
    if(curl == NULL) return;
    NCD4_curlclose(curl->curl);
    ncbytesfree(curl->packet);
    nullfree(curl->errdata.code);
    nullfree(curl->errdata.message);
    free(curl);
}

/* Define the set of protocols known to be constrainable */
static const char* constrainableprotocols[] = {"http", "https",NULL};

static int
constrainable(NCURI* durl)
{
   const char** protocol = constrainableprotocols;
   for(;*protocol;protocol++) {
	if(strcmp(durl->protocol,*protocol)==0)
	    return 1;
   }
   return 0;
}

/*
Set curl properties for link based on rc files etc.
*/
static int
set_curl_properties(NCD4INFO* d4info)
{
    int ret = NC_NOERR;

    if(d4info->auth.curlflags.useragent == NULL) {
	char* agent;
        size_t len = strlen(DFALTUSERAGENT) + strlen(VERSION);
	len++; /*strlcat nul*/
	agent = (char*)malloc(len+1);
	strncpy(agent,DFALTUSERAGENT,len);
	strlcat(agent,VERSION,len);
        d4info->auth.curlflags.useragent = agent;
    }

    /* Some servers (e.g. thredds and columbia) appear to require a place
       to put cookies in order for some security functions to work
    */
    if(d4info->auth.curlflags.cookiejar != NULL
       && strlen(d4info->auth.curlflags.cookiejar) == 0) {
	free(d4info->auth.curlflags.cookiejar);
	d4info->auth.curlflags.cookiejar = NULL;
    }

    if(d4info->auth.curlflags.cookiejar == NULL) {
	/* If no cookie file was defined, define a default */
        char* path = NULL;
        char* newpath = NULL;
        int len;
	errno = 0;
	/* Create the unique cookie file name */
        len =
	  strlen(ncrc_globalstate.tempdir)
	  + 1 /* '/' */
	  + strlen("ncd4cookies");
        path = (char*)malloc(len+1);
        if(path == NULL) return NC_ENOMEM;
	snprintf(path,len,"%s/nc4cookies",ncrc_globalstate.tempdir);
	/* Create the unique cookie file name */
        newpath = NC_mktmp(path);
        free(path);
	if(newpath == NULL) {
	    fprintf(stderr,"Cannot create cookie file\n");
	    goto fail;
	}
	d4info->auth.curlflags.cookiejar = newpath;
	d4info->auth.curlflags.cookiejarcreated = 1;
	errno = 0;
    }
    assert(d4info->auth.curlflags.cookiejar != NULL);

    /* Make sure the cookie jar exists and can be read and written */
    {
	FILE* f = NULL;
	char* fname = d4info->auth.curlflags.cookiejar;
	/* See if the file exists already */
        f = fopen(fname,"r");
	if(f == NULL) {
	    /* Ok, create it */
	    f = fopen(fname,"w+");
	    if(f == NULL) {
	        fprintf(stderr,"Cookie file cannot be read and written: %s\n",fname);
	        {ret= NC_EPERM; goto fail;}
	    }
	} else { /* test if file can be written */
	    fclose(f);
	    f = fopen(fname,"r+");
	    if(f == NULL) {
	        fprintf(stderr,"Cookie file is cannot be written: %s\n",fname);
	        {ret = NC_EPERM; goto fail;}
	    }
	}
	if(f != NULL) fclose(f);
    }

    return THROW(ret);

fail:
    return THROW(ret);
}

void
NCD4_applyclientparamcontrols(NCD4INFO* info)
{
    const char* value;

    /* clear the flags */
    CLRFLAG(info->controls.flags,NCF_CACHE);
    CLRFLAG(info->controls.flags,NCF_SHOWFETCH);
    CLRFLAG(info->controls.flags,NCF_NC4);
    CLRFLAG(info->controls.flags,NCF_NCDAP);
    CLRFLAG(info->controls.flags,NCF_FILLMISMATCH);

    /* Turn on any default on flags */
    SETFLAG(info->controls.flags,DFALT_ON_FLAGS);
    SETFLAG(info->controls.flags,(NCF_NC4|NCF_NCDAP));

    if(paramcheck(info,"show","fetch"))
	SETFLAG(info->controls.flags,NCF_SHOWFETCH);

    if(paramcheck(info,"translate","nc4"))
	info->controls.translation = NCD4_TRANSNC4;

    /* Look at the debug flags */
    if(paramcheck(info,"debug","copy"))
	SETFLAG(info->controls.debugflags,NCF_DEBUG_COPY); /* => close */

    value = getparam(info,"substratename");
    if(value != NULL)
	strncpy(info->controls.substratename,value,NC_MAX_NAME);

    info->controls.opaquesize = DFALTOPAQUESIZE;
    value = getparam(info,"opaquesize");
    if(value != NULL) {
	long long len = 0;
	if(sscanf(value,"%lld",&len) != 1 || len == 0)
	    nclog(NCLOGWARN,"bad [opaquesize] tag: %s",value);
	else
	    info->controls.opaquesize = (size_t)len;	    
    }
    
    value = getparam(info,"fillmismatch");
    if(value != NULL)
	SETFLAG(info->controls.flags,NCF_FILLMISMATCH);

    value = getparam(info,"nofillmismatch");
    if(value != NULL)
	CLRFLAG(info->controls.debugflags,NCF_FILLMISMATCH);
}

static void
applyclientmetacontrols(NCD4meta* meta)    
{
    NCD4INFO* info = meta->controller;
    const char* value = getparam(info,"checksummode");
    if(value != NULL) {
        if(strcmp(value,"ignore")==0)
	    meta->ignorechecksums = 1;
    }
}

/* Search for substring in value of param. If substring == NULL; then just
   check if param is defined.
*/
static int
paramcheck(NCD4INFO* info, const char* key, const char* subkey)
{
    const char* value;
    char* p;

    value = getparam(info, key);
    if(value == NULL)
	return 0;
    if(subkey == NULL) return 1;
    p = strstr(value,subkey);
    if(p == NULL) return 0;
    p += strlen(subkey);
    if(*p != '\0' && strchr(checkseps,*p) == NULL) return 0;
    return 1;
}

/*
Given a parameter key, return its value or NULL if not defined.
*/
static const char*
getparam(NCD4INFO* info, const char* key)
{
    const char* value;

    if(info == NULL || key == NULL) return NULL;
    if((value=ncurilookup(info->uri,key)) == NULL)
	return NULL;
    return value;
}


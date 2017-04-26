/* Copyright 2009, UCAR/Unidata and OPeNDAP, Inc.
   See the COPYRIGHT file for more information. */

#include "config.h"
#include <stdio.h>

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#ifdef HAVE_FCNTL_H
#include <fcntl.h>
#endif
#include <errno.h>

#ifdef _MSC_VER
typedef int pid_t;
#endif

#include "ocinternal.h"
#include "ocdebug.h"
#include "occlientparams.h"
#include "occurlfunctions.h"
#include "ochttp.h"
#include "ocread.h"
#include "dapparselex.h"

#define DATADDSFILE "datadds"

#if 0
/* Note: TMPPATH must end in '/' */
#ifdef __CYGWIN__
#define TMPPATH1 "/cygdrive/c/temp/datadds"
#define TMPPATH2 "./datadds"
#elif defined(_WIN32) || defined(_WIN64)
#define TMPPATH1 "c:\\temp\\datadds"
#define TMPPATH2 ".\\datadds"
#else
#define TMPPATH1 "/tmp/datadds"
#define TMPPATH2 "./datadds"
#endif
#endif

#define CLBRACE '{'
#define CRBRACE '}'

static OCerror ocextractddsinmemory(OCstate*,OCtree*,int);
static OCerror ocextractddsinfile(OCstate*,OCtree*,int);
static char* constraintescape(const char* url);
static OCerror createtempfile(OCstate*,OCtree*);
static int dataError(XXDR* xdrs, OCstate*);

static OCerror ocset_curlproperties(OCstate*);

extern OCnode* makeunlimiteddimension(void);

#if defined(_WIN32) || defined(_WIN64)
#ifdef HAVE_FCNTL_H
#include <fcntl.h>
#endif
#define _S_IREAD 256
#define _S_IWRITE 128
#else
#ifdef HAVE_SYS_STAT_H
#include <sys/stat.h>
#endif
#endif

/* Collect global state info in one place */
struct OCGLOBALSTATE ocglobalstate;

OCerror
ocinternalinitialize(void)
{
    int stat = OC_NOERR;

#if 0
    if(sizeof(off_t) != sizeof(void*)) {
      fprintf(stderr,"OC xxdr depends on the assumption that sizeof(off_t) == sizeof(void*)\n");
      /*
	Commenting out for now, as this does not hold true on 32-bit linux systems.
      OCASSERT(sizeof(off_t) == sizeof(void*));
	*/
    }
#endif

    if(!ocglobalstate.initialized) {
      memset((void*)&ocglobalstate,0,sizeof(ocglobalstate));
      ocglobalstate.initialized = 1;
    }

    /* Capture temp dir*/
    {
	char* tempdir;
	char* p;
	char* q;
	char cwd[4096];
#if defined(_WIN32) || defined(_WIN64)
        tempdir = getenv("TEMP");
#else
	tempdir = "/tmp";
#endif
        if(tempdir == NULL) {
	    fprintf(stderr,"Cannot find a temp dir; using ./\n");
	    tempdir = getcwd(cwd,sizeof(cwd));
	    if(tempdir == NULL || *tempdir == '\0') tempdir = ".";
	}
        ocglobalstate.tempdir= (char*)malloc(strlen(tempdir) + 1);
	for(p=tempdir,q=ocglobalstate.tempdir;*p;p++,q++) {
	    if((*p == '/' && *(p+1) == '/')
	       || (*p == '\\' && *(p+1) == '\\')) {p++;}
	    *q = *p;
	}
	*q = '\0';
#if defined(_WIN32) || defined(_WIN64)
#else
        /* Canonicalize */
	for(p=ocglobalstate.tempdir;*p;p++) {
	    if(*p == '\\') {*p = '/'; };
	}
	*q = '\0';
#endif
    }

    /* Capture $HOME */
    {
	char* p;
	char* q;
        char* home = getenv("HOME");

        if(home == NULL) {
	    /* use tempdir */
	    home = ocglobalstate.tempdir;
	}
        ocglobalstate.home = (char*)malloc(strlen(home) + 1);
	for(p=home,q=ocglobalstate.home;*p;p++,q++) {
	    if((*p == '/' && *(p+1) == '/')
	       || (*p == '\\' && *(p+1) == '\\')) {p++;}
	    *q = *p;
	}
	*q = '\0';
#if defined(_WIN32) || defined(_WIN64)
#else
        /* Canonicalize */
	for(p=home;*p;p++) {
	    if(*p == '\\') {*p = '/'; };
	}
#endif
    }

    /* Compute some xdr related flags */
    xxdr_init();

    ocloginit();

    oc_curl_protocols(&ocglobalstate); /* see what protocols are supported */

    return OCTHROW(stat);
}


/**************************************************/
OCerror
ocopen(OCstate** statep, const char* url)
{
    int stat = OC_NOERR;
    OCstate * state = NULL;
    OCURI* tmpurl = NULL;
    CURL* curl = NULL; /* curl handle*/

    if(!ocuriparse(url,&tmpurl)) {OCTHROWCHK(stat=OC_EBADURL); goto fail;}

    stat = occurlopen(&curl);
    if(stat != OC_NOERR) {OCTHROWCHK(stat); goto fail;}

    state = (OCstate*)ocmalloc(sizeof(OCstate)); /* ocmalloc zeros memory*/
    if(state == NULL) {OCTHROWCHK(stat=OC_ENOMEM); goto fail;}

    /* Setup DAP state*/
    state->header.magic = OCMAGIC;
    state->header.occlass = OC_State;
    state->curl = curl;
    state->trees = oclistnew();
    state->uri = tmpurl;
    if(!ocuridecodeparams(state->uri)) {
	oclog(OCLOGWARN,"Could not parse client parameters");
    }
    state->packet = ocbytesnew();
    ocbytessetalloc(state->packet,DFALTPACKETSIZE); /*initial reasonable size*/

    /* capture curl properties for this link from rc file1*/
    stat = ocset_curlproperties(state);
    if(stat != OC_NOERR) goto fail;

    /* Set the one-time curl flags */
    if((stat=ocset_flags_perlink(state))!= OC_NOERR) goto fail;
#if 1 /* temporarily make per-link */
    if((stat=ocset_flags_perfetch(state))!= OC_NOERR) goto fail;
#endif

    if(statep) *statep = state;
    else {
      if(state != NULL) ocfree(state);
    }
    return OCTHROW(stat);

fail:
    ocurifree(tmpurl);
    if(state != NULL) ocfree(state);
    if(curl != NULL) occurlclose(curl);
    return OCTHROW(stat);
}

OCerror
ocfetch(OCstate* state, const char* constraint, OCdxd kind, OCflags flags,
        OCnode** rootp)
{
    OCtree* tree = NULL;
    OCnode* root = NULL;
    OCerror stat = OC_NOERR;

    tree = (OCtree*)ocmalloc(sizeof(OCtree));
    MEMCHECK(tree,OC_ENOMEM);
    memset((void*)tree,0,sizeof(OCtree));
    tree->dxdclass = kind;
    tree->state = state;
    tree->constraint = constraintescape(constraint);
    if(tree->constraint == NULL)
	tree->constraint = nulldup(constraint);

    /* Set per-fetch curl properties */
#if 0 /* temporarily make per-link */
    if((stat=ocset_flags_perfetch(state))!= OC_NOERR) goto fail;
#endif

    ocbytesclear(state->packet);

    switch (kind) {
    case OCDAS:
        stat = readDAS(state,tree);
	if(stat == OC_NOERR) {
            tree->text = ocbytesdup(state->packet);
	    if(tree->text == NULL) stat = OC_EDAS;
	}
	break;
    case OCDDS:
        stat = readDDS(state,tree);
	if(stat == OC_NOERR) {
            tree->text = ocbytesdup(state->packet);
	    if(tree->text == NULL) stat = OC_EDDS;
	}
	break;
    case OCDATADDS:
	if((flags & OCONDISK) != 0) {/* store in file */
	    /* Create the datadds file immediately
               so that DRNO can reference it*/
            /* Make the tmp file*/
            stat = createtempfile(state,tree);
            if(stat) {OCTHROWCHK(stat); goto fail;}
            stat = readDATADDS(state,tree,flags);
	    if(stat == OC_NOERR) {
                /* Separate the DDS from data and return the dds;
                   will modify packet */
                stat = ocextractddsinfile(state,tree,flags);
	    }
	} else { /*inmemory*/
            stat = readDATADDS(state,tree,flags);
	    if(stat == OC_NOERR) {
                /* Separate the DDS from data and return the dds;
               will modify packet */
            stat = ocextractddsinmemory(state,tree,flags);
	}
	}
	break;
    default:
	break;
    }/*switch*/
    /* Obtain any http code */
    state->error.httpcode = ocfetchhttpcode(state->curl);
    if(stat != OC_NOERR) {
	if(state->error.httpcode >= 400) {
	    oclog(OCLOGWARN,"oc_open: Could not read url; http error = %l",state->error.httpcode);
	} else {
	    oclog(OCLOGWARN,"oc_open: Could not read url");
	}
	goto fail;
    }

    tree->nodes = NULL;
    stat = DAPparse(state,tree,tree->text);
    /* Check and report on an error return from the server */
    if(stat == OC_EDAPSVC  && state->error.code != NULL) {
	oclog(OCLOGERR,"oc_open: server error retrieving url: code=%s message=\"%s\"",
		  state->error.code,
		  (state->error.message?state->error.message:""));
    }
    if(stat) {OCTHROWCHK(stat); goto fail;}
    root = tree->root;
    /* make sure */
    tree->root = root;
    root->tree = tree;

    /* Verify the parse */
    switch (kind) {
    case OCDAS:
        if(root->octype != OC_Attributeset)
	    {OCTHROWCHK(stat=OC_EDAS); goto fail;}
	break;
    case OCDDS:
        if(root->octype != OC_Dataset)
	    {OCTHROWCHK(stat=OC_EDDS); goto fail;}
	break;
    case OCDATADDS:
        if(root->octype != OC_Dataset)
	    {OCTHROWCHK(stat=OC_EDATADDS); goto fail;}
	/* Modify the tree kind */
	tree->dxdclass = OCDATADDS;
	break;
    default: return OC_EINVAL;
    }

    if(kind != OCDAS) {
        /* Process ocnodes to mark those that are cacheable */
        ocmarkcacheable(state,root);
        /* Process ocnodes to handle various semantic issues*/
        occomputesemantics(tree->nodes);
    }

    /* Process ocnodes to compute name info*/
    occomputefullnames(tree->root);

     if(kind == OCDATADDS) {
	if((flags & OCONDISK) != 0) {
            tree->data.xdrs = xxdr_filecreate(tree->data.file,tree->data.bod);
	} else {
#ifdef OCDEBUG
fprintf(stderr,"ocfetch.datadds.memory: datasize=%lu bod=%lu\n",
	(unsigned long)tree->data.datasize,(unsigned long)tree->data.bod);
#endif
	    /* Switch to zero based memory */
            tree->data.xdrs
		= xxdr_memcreate(tree->data.memory,tree->data.datasize,tree->data.bod);
	}
        MEMCHECK(tree->data.xdrs,OC_ENOMEM);
	/* Do a quick check to see if server returned an ERROR {}
           at the beginning of the data
         */
	if(dataError(tree->data.xdrs,state)) {
	    stat = OC_EDATADDS;
	    oclog(OCLOGERR,"oc_open: server error retrieving url: code=%s message=\"%s\"",
		  state->error.code,
		  (state->error.message?state->error.message:""));
	    goto fail;
	}

	/* Compile the data into a more accessible format */
	stat = occompile(state,tree->root);
	if(stat != OC_NOERR)
	    goto fail;
    }

    /* Put root into the state->trees list */
    oclistpush(state->trees,(void*)root);

    if(rootp) *rootp = root;
    return stat;

fail:
    if(root != NULL)
	ocroot_free(root);
    else if(tree != NULL)
	octree_free(tree);
    return OCTHROW(stat);
}

static OCerror
createtempfile(OCstate* state, OCtree* tree)
{
    int stat = OC_NOERR;
    char* path = NULL;
    char* name = NULL;
    int len;

    len =
	  strlen(ocglobalstate.tempdir)
	  + 1 /* '/' */
	  + strlen(DATADDSFILE);
    path = (char*)malloc(len+1);
    if(path == NULL) return OC_ENOMEM;
    occopycat(path,len,3,ocglobalstate.tempdir,"/",DATADDSFILE);
    stat = ocmktmp(path,&name);
    free(path);
    if(stat != OC_NOERR) goto fail;
#ifdef OCDEBUG
    oclog(OCLOGNOTE,"oc_open: creating tmp file: %s",name);
#endif
    tree->data.filename = name; /* remember our tmp file name */
    name = NULL;
    tree->data.file = fopen(tree->data.filename,"w+");
    if(tree->data.file == NULL) return OC_EOPEN;
    /* unlink the temp file so it will automatically be reclaimed */
    if(ocdebug == 0) unlink(tree->data.filename);
    return stat;

fail:
    if(name != NULL) {
        oclog(OCLOGERR,"oc_open: attempt to create tmp file failed: %s",name);
	free(name);
    } else {
        oclog(OCLOGERR,"oc_open: attempt to create tmp file failed: NULL");
    }
    return OCTHROW(stat);
}

void
occlose(OCstate* state)
{
    unsigned int i;
    if(state == NULL) return;

    /* Warning: ocfreeroot will attempt to remove the root from state->trees */
    /* Ok in this case because we are popping the root out of state->trees */
    for(i=0;i<oclistlength(state->trees);i++) {
	OCnode* root = (OCnode*)oclistpop(state->trees);
	ocroot_free(root);
    }
    oclistfree(state->trees);
    ocurifree(state->uri);
    ocbytesfree(state->packet);
    ocfree(state->error.code);
    ocfree(state->error.message);
    ocfree(state->curlflags.useragent);
    if(state->curlflags.cookiejar) {
#if 0
        if(state->curlflags.createdflags & COOKIECREATED)
	    unlink(state->curlflags.cookiejar);
#endif
	ocfree(state->curlflags.cookiejar);
    }
    if(state->curlflags.netrc != NULL) {
#if 0
        if(state->curlflags.createdflags & NETRCCREATED)
	    unlink(state->curlflags.netrc);
#endif
	ocfree(state->curlflags.netrc);
    }
    ocfree(state->ssl.certificate);
    ocfree(state->ssl.key);
    ocfree(state->ssl.keypasswd);
    ocfree(state->ssl.cainfo);
    ocfree(state->ssl.capath);
    ocfree(state->proxy.host);
    ocfree(state->proxy.userpwd);
    ocfree(state->creds.userpwd);
    if(state->curl != NULL) occurlclose(state->curl);
    ocfree(state);
}

static OCerror
ocextractddsinmemory(OCstate* state, OCtree* tree, OCflags flags)
{
    OCerror stat = OC_NOERR;
    size_t ddslen, bod, bodfound;
    /* Read until we find the separator (or EOF)*/
    bodfound = ocfindbod(state->packet,&bod,&ddslen);
    if(!bodfound) {/* No BOD; pretend */
	bod = tree->data.bod;
	ddslen = tree->data.datasize;
    }
    tree->data.bod = bod;
    tree->data.ddslen = ddslen;
    /* copy out the dds */
    if(ddslen > 0) {
        tree->text = (char*)ocmalloc(ddslen+1);
        memcpy((void*)tree->text,(void*)ocbytescontents(state->packet),ddslen);
        tree->text[ddslen] = '\0';
    } else
	tree->text = NULL;
    /* Extract the inmemory contents */
    tree->data.memory = ocbytesextract(state->packet);
#ifdef OCIGNORE
    /* guarantee the data part is on an 8 byte boundary */
    if(tree->data.bod % 8 != 0) {
        unsigned long count = tree->data.datasize - tree->data.bod;
        memcpy(tree->xdrmemory,tree->xdrmemory+tree->data.bod,count);
        tree->data.datasize = count;
	tree->data.bod = 0;
	tree->data.ddslen = 0;
    }
#endif
    if(tree->text == NULL) stat = OC_EDATADDS;
    return OCTHROW(stat);
}

static OCerror
ocextractddsinfile(OCstate* state, OCtree* tree, OCflags flags)
{
    OCerror stat = OC_NOERR;
    size_t ddslen, bod, bodfound;

    /* Read until we find the separator (or EOF)*/
    ocbytesclear(state->packet);
    rewind(tree->data.file);
    bodfound = 0;
    do {
        char chunk[1024];
	size_t count;
	/* read chunks of the file until we find the separator*/
        count = fread(chunk,1,sizeof(chunk),tree->data.file);
	if(count <= 0) break; /* EOF;*/
        ocbytesappendn(state->packet,chunk,count);
	bodfound = ocfindbod(state->packet,&bod,&ddslen);
    } while(!bodfound);
    if(!bodfound) {/* No BOD; pretend */
	bod = tree->data.bod;
	ddslen = tree->data.datasize;
#ifdef OCDEBUG
fprintf(stderr,"missing bod: ddslen=%lu bod=%lu\n",
(unsigned long)ddslen,(unsigned long)bod);
#endif
    }
    tree->data.bod = bod;
    tree->data.ddslen = ddslen;
    /* copy out the dds */
    if(ddslen > 0) {
        tree->text = (char*)ocmalloc(ddslen+1);
        memcpy((void*)tree->text,(void*)ocbytescontents(state->packet),ddslen);
        tree->text[ddslen] = '\0';
    } else
	tree->text = NULL;
    /* reset the position of the tmp file*/
    if(fseek(tree->data.file,(long)tree->data.bod,SEEK_SET) < 0
       || tree->text == NULL)
	stat = OC_EDATADDS;
    return OCTHROW(stat);
}

/* Allow these (non-alpha-numerics) to pass thru */
static char okchars[] = "&/:;,.=?@'\"<>{}!|\\^[]`~";
static char hexdigits[] = "0123456789abcdef";

/* Modify constraint to use %XX escapes */
static char*
constraintescape(const char* url)
{
    size_t len;
    char* p;
    int c;
    char* eurl;

    if(url == NULL) return NULL;
    len = strlen(url);
    eurl = ocmalloc(1+3*len); /* worst case: c -> %xx */
    MEMCHECK(eurl,NULL);
    p = eurl;
    *p = '\0';
    while((c=*url++)) {
	if(c >= '0' && c <= '9') {*p++ = c;}
	else if(c >= 'a' && c <= 'z') {*p++ = c;}
	else if(c >= 'A' && c <= 'Z') {*p++ = c;}
	else if(strchr(okchars,c) != NULL) {*p++ = c;}
	else {
	    *p++ = '%';
	    *p++ = hexdigits[(c & 0xf0)>>4];
	    *p++ = hexdigits[(c & 0xf)];
	}
    }
    *p = '\0';
    return eurl;
}

OCerror
ocupdatelastmodifieddata(OCstate* state)
{
    OCerror status = OC_NOERR;
    long lastmodified;
    char* base = NULL;
    base = ocuribuild(state->uri,NULL,NULL,OCURIENCODE);
    status = ocfetchlastmodified(state->curl, base, &lastmodified);
    free(base);
    if(status == OC_NOERR) {
	state->datalastmodified = lastmodified;
    }
    return OCTHROW(status);
}

/*
    Set curl properties for link based on rc files etc.
*/
static OCerror
ocset_curlproperties(OCstate* state)
{
    OCerror stat = OC_NOERR;

    /* extract the relevant triples int state */
    ocrc_process(state);

    if(state->curlflags.useragent == NULL) {
        size_t len = strlen(DFALTUSERAGENT) + strlen(VERSION) + 1;
	char* agent = (char*)malloc(len+1);
	if(occopycat(agent,len,2,DFALTUSERAGENT,VERSION))
	    state->curlflags.useragent = agent;
	else
	    free(agent);
    }

    /* Some servers (e.g. thredds and columbia) appear to require a place
       to put cookies in order for some security functions to work
    */
    if(state->curlflags.cookiejar != NULL
       && strlen(state->curlflags.cookiejar) == 0) {
	free(state->curlflags.cookiejar);
	state->curlflags.cookiejar = NULL;
    }

    if(state->curlflags.cookiejar == NULL) {
	/* If no cookie file was defined, define a default */
	char tmp[OCPATHMAX+1];
        int stat;
	pid_t pid = getpid();
	snprintf(tmp,sizeof(tmp)-1,"%s/%s.%ld/",ocglobalstate.tempdir,OCDIR,(long)pid);
#ifdef _WIN32
	stat = mkdir(tmp);
#else
	stat = mkdir(tmp,S_IRUSR | S_IWUSR | S_IXUSR);
#endif
	if(stat != 0 && errno != EEXIST) {
	    fprintf(stderr,"Cannot create cookie directory\n");
	    goto fail;
	}
	errno = 0;
	/* Create the unique cookie file name */
	stat = ocmktmp(tmp,&state->curlflags.cookiejar);
	state->curlflags.createdflags |= COOKIECREATED;
	if(stat != OC_NOERR && errno != EEXIST) {
	    fprintf(stderr,"Cannot create cookie file\n");
	    goto fail;
	}
	errno = 0;
    }
    OCASSERT(state->curlflags.cookiejar != NULL);

    /* Make sure the cookie jar exists and can be read and written */
    {
	FILE* f = NULL;
	char* fname = state->curlflags.cookiejar;
	/* See if the file exists already */
        f = fopen(fname,"r");
	if(f == NULL) {
	    /* Ok, create it */
	    f = fopen(fname,"w+");
	    if(f == NULL) {
	        fprintf(stderr,"Cookie file cannot be read and written: %s\n",fname);
	        {stat = OC_EPERM; goto fail;}
	    }
	} else { /* test if file can be written */
	    fclose(f);
	    f = fopen(fname,"r+");
	    if(f == NULL) {
	        fprintf(stderr,"Cookie file is cannot be written: %s\n",fname);
	        {stat = OC_EPERM; goto fail;}
	    }
	}
	if(f != NULL) fclose(f);
    }

#if 0
    /* Create a netrc file if specified  and required,
       where required => >1 NETRC triples exist */
    if(ocrc_netrc_required(state)) {
	/* WARNING: it appears that a user+pwd was specified specifically, then
           the netrc file will be completely disabled. */
	if(state->creds.userpwd != NULL) {
  	    oclog(OCLOGWARN,"The rc file specifies both netrc and user+pwd; this will cause curl to ignore the netrc file");
	}
	stat = oc_build_netrc(state);
    }
#endif

    return stat;

fail:
    return OCTHROW(stat);
}

static char* ERROR_TAG = "Error ";

static int
dataError(XXDR* xdrs, OCstate* state)
{
    int depth=0;
    int errfound = 0;
    off_t ckp=0,avail=0;
    int i=0;
    char* errmsg = NULL;
    char errortext[16]; /* bigger thant |ERROR_TAG|*/
    avail = xxdr_getavail(xdrs);
    if(avail < strlen(ERROR_TAG))
	goto done; /* assume it is ok */
    ckp = xxdr_getpos(xdrs);
    /* Read enough characters to test for 'ERROR ' */
    errortext[0] = '\0';
    xxdr_getbytes(xdrs,errortext,(off_t)strlen(ERROR_TAG));
    if(ocstrncmp(errortext,ERROR_TAG,strlen(ERROR_TAG)) != 0)
	goto done; /* not an immediate error */
    /* Try to locate the whole error body */
    xxdr_setpos(xdrs,ckp);
    for(depth=0,i=0;i<avail;i++) {
	xxdr_getbytes(xdrs,errortext,(off_t)1);
	if(errortext[0] == CLBRACE) depth++;
	else if(errortext[0] == CRBRACE) {
	    depth--;
	    if(depth == 0) {i++; break;}
	}
    }
    errmsg = (char*)malloc((size_t)i+1);
    if(errmsg == NULL) {errfound = 1; goto done;}
    xxdr_setpos(xdrs,ckp);
    xxdr_getbytes(xdrs,errmsg,(off_t)i);
    errmsg[i] = '\0';
    state->error.message = errmsg;
    state->error.code = strdup("?");
    state->error.httpcode = 404;
    xxdr_setpos(xdrs,ckp);
    errfound = 1;
done:
    xxdr_setpos(xdrs,ckp);
    return errfound;
}

/**************************************************/
/* Curl option functions */
/*
Note that if we set the option in curlflags,
then we need to also invoke the ocset_curlopt
to update the curl flags in libcurl.
*/

OCerror
ocset_useragent(OCstate* state, const char* agent)
{
    OCerror stat = OC_NOERR;
    if(state->curlflags.useragent != NULL)
	free(state->curlflags.useragent);
    state->curlflags.useragent = strdup(agent);
    if(state->curlflags.useragent == NULL)
	return OCTHROW(OC_ENOMEM);
    stat = ocset_curlflag(state,CURLOPT_USERAGENT);
    return stat;
}

OCerror
ocset_netrc(OCstate* state, const char* path)
{
    OCerror stat = OC_NOERR;
    if(state->curlflags.netrc != NULL)
	free(state->curlflags.netrc);
    state->curlflags.netrc = strdup(path);
    if(state->curlflags.netrc == NULL)
	return OCTHROW(OC_ENOMEM);
    stat = ocset_curlflag(state,CURLOPT_NETRC);
    return stat;
}

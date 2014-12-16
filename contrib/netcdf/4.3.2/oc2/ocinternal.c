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

#include "ocinternal.h"
#include "ocdebug.h"
#include "occlientparams.h"
#include "ocrc.h"
#include "occurlfunctions.h"
#include "ochttp.h"
#include "ocread.h"
#include "dapparselex.h"

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

#define CLBRACE '{'
#define CRBRACE '}'


/* Define default rc files and aliases*/
static char* rcfilenames[4] = {".daprc",".dodsrc",".ocrc",NULL};

static int ocextractddsinmemory(OCstate*,OCtree*,int);
static int ocextractddsinfile(OCstate*,OCtree*,int);
static char* constraintescape(const char* url);
static OCerror createtempfile(OCstate*,OCtree*);
static int dataError(XXDR* xdrs, OCstate*);

static int ocsetcurlproperties(OCstate*);

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

int
ocinternalinitialize(void)
{
    int stat = OC_NOERR;

    if(!ocglobalstate.initialized) {
        memset((void*)&ocglobalstate,0,sizeof(ocglobalstate));
	ocglobalstate.initialized = 1;
    }

    /* Capture $HOME */
    {
	char* p;
	char* q;
        char* home = getenv("HOME");
	char cwd[4096];
        if(ocglobalstate.home == NULL) {
#if defined(_WIN32) || defined(_WIN64)
	    home = getenv("TEMP");
#else
	    home = "/tmp";
#endif
	}
        if(home == NULL) {
	    home = getcwd(cwd,sizeof(cwd));
	    if(home == NULL || *home == '\0') home = ".";
	}

        /* Convert '\' to '/' */
        ocglobalstate.home = (char*)malloc(strlen(home) + 1);
	for(p=home,q=ocglobalstate.home;*p;p++,q++) {
	    if(*p == '\\') {*q = '/'; } else {*q = *p;}
	}
	*q = '\0';
    }

    /* Compute some xdr related flags */
    xxdr_init();

    ocloginit();

    oc_curl_protocols(&ocglobalstate); /* see what protocols are supported */

    /* compile the .dodsrc, if any */
    {
        char* path = NULL;
	char** alias;
	FILE* f = NULL;
        /* locate the configuration files: . first in '.',  then $HOME */
	for(alias=rcfilenames;*alias;alias++) {
	    size_t pathlen = strlen("./")+strlen(*alias)+1;
            path = (char*)malloc(pathlen);
	    if(path == NULL) return OC_ENOMEM;
	    if(!occopycat(path,pathlen,2,"./",*alias)) {
	        if(path) free(path);
		return OC_EOVERRUN;
	    }
  	    /* see if file is readable */
	    f = fopen(path,"r");
	    if(f != NULL) break;
    	    if(path != NULL) {free(path); path = NULL;} /* cleanup */
	}
	if(f == NULL) { /* try $HOME */
	    OCASSERT(path == NULL);
	    for(alias=rcfilenames;*alias;alias++) {
		size_t pathlen = strlen(ocglobalstate.home)+1+strlen(*alias)+1;
	        path = (char*)malloc(pathlen);
	        if(path == NULL) return OC_ENOMEM;
		if(!occopycat(path,pathlen,3,ocglobalstate.home,"/",*alias)) {
		    if(path) free(path);
		    return OC_EOVERRUN;
		}
		f = fopen(path,"r");
		if(f != NULL) break;
 	        if(path != NULL) {free(path); path=NULL;}
            }
	}
        if(f == NULL) {
            oclog(OCLOGDBG,"Cannot find runtime configuration file");
	} else {
	    OCASSERT(path != NULL);
       	    fclose(f);
            if(ocdebug > 1)
		fprintf(stderr, "DODS RC file: %s\n", path);
            if(ocdodsrc_read(*alias,path) == 0)
	        oclog(OCLOGERR, "Error parsing %s\n",path);
        }
        if(path != NULL) free(path);
    }

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

    /* set curl properties for this link */
    stat = ocsetcurlproperties(state);

    if(statep) *statep = state;
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

    /* Set curl properties: pwd, flags, proxies, ssl */
    if((stat=ocset_user_password(state))!= OC_NOERR) goto fail;
    if((stat=ocset_curl_flags(state)) != OC_NOERR) goto fail;
    if((stat=ocset_proxy(state)) != OC_NOERR) goto fail;
    if((stat=ocset_ssl(state)) != OC_NOERR) goto fail;
    if(state->usercurl)
	state->usercurl((void*)state->curl,state->usercurldata);

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
    int fd = 0;
    char* name = NULL;

    stat = ocmktmp(TMPPATH1,&name, &fd);
    if(stat != OC_NOERR)
        stat = ocmktmp(TMPPATH2,&name,&fd);
    if(stat != OC_NOERR) goto fail;
#ifdef OCDEBUG
    oclog(OCLOGNOTE,"oc_open: using tmp file: %s",name);
#endif
    tree->data.filename = name; /* remember our tmp file name */
    tree->data.file = fdopen(fd,"w+");
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
    
    return stat;
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
	unlink(state->curlflags.cookiejar);
	ocfree(state->curlflags.cookiejar);
    }
    ocfree(state->ssl.certificate);
    ocfree(state->ssl.key);
    ocfree(state->ssl.keypasswd);
    ocfree(state->ssl.cainfo);
    ocfree(state->ssl.capath); 
    ocfree(state->proxy.host);
    ocfree(state->creds.username);
    ocfree(state->creds.password);
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
    fseek(tree->data.file,(long)tree->data.bod,SEEK_SET);
    if(tree->text == NULL) stat = OC_EDATADDS;
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
    return status;
}

/*
    Set curl properties for link based on rc files etc.
*/
static int
ocsetcurlproperties(OCstate* state)
{
    CURLcode cstat = CURLE_OK;

    /* process the triple store wrt to this state */
    if(ocdodsrc_process(state) != OC_NOERR) {
	oclog(OCLOGERR,"Malformed .opendaprc configuration file");
	goto fail;
    }
    if(state->creds.username == NULL && state->creds.password == NULL) {
        if(state->uri->user != NULL && state->uri->password != NULL) {
	    /* this overrides .dodsrc */
            if(state->creds.password) free(state->creds.password);
            state->creds.password = nulldup(state->uri->password);
            if(state->creds.username) free(state->creds.username);
            state->creds.username = nulldup(state->uri->user);
	}
    }
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
    if(state->curlflags.cookiejar == NULL 
       || *state->curlflags.cookiejar) {
#if 1
	/* Apparently anything non-null will work */
	state->curlflags.cookiejar = strdup("");
#else
	/* If no cookie file was defined, define a default */
	char* tmp;
	int fd;
        int stat;
		
        tmp = (char*)malloc(strlen(ocglobalstate.home)
				  +strlen("/")
				  +strlen(OCDIR)
				  +strlen("/")
				  +1);
	if(tmp == NULL)
	    return OC_ENOMEM;
	strcpy(tmp,ocglobalstate.home);
	strcat(tmp,"/");
	strcat(tmp,OCDIR);
	strcat(tmp,"/");
	stat = mkdir(tmp,S_IRUSR | S_IWUSR | S_IXUSR);
	if(stat != 0 && errno != EEXIST) {
	    fprintf(stderr,"Cannot create cookie file\n");
	    return stat;
	}
	errno = 0;
	/* Create the actual cookie file */
	stat = ocmktmp(tmp,&state->curlflags.cookiejar,&fd);
	close(fd);	

#if 0
	fd = creat(tmp,S_IRUSR | S_IWUSR);
	if(fd < 0) {
	    fprintf(stderr,"Cannot create cookie file\n");
	    return OC_EPERM;
	}else
	    close(fd);
#endif
#endif
    }
    return OC_NOERR;

fail:
    if(cstat != CURLE_OK)
	oclog(OCLOGERR, "curl error: %s", curl_easy_strerror(cstat));
    return OC_ECURL;
}

OCerror
ocsetuseragent(OCstate* state, const char* agent)
{
    if(state->curlflags.useragent != NULL)
	free(state->curlflags.useragent);
    state->curlflags.useragent = strdup(agent);
    if(state->curlflags.useragent == NULL)
	return OC_ENOMEM;
    return OC_NOERR;
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

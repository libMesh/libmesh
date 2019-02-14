#include "d4includes.h"
#include "d4curlfunctions.h"
#ifdef HAVE_FCNTL_H
#include <fcntl.h>
#endif
#include "ncwinpath.h"

/* Do conversion if this code was compiled via Vis. Studio or Mingw */

/*Forward*/
static int readpacket(NCD4INFO* state, NCURI*, NCbytes*, NCD4mode, long*);
static int readfile(NCD4INFO* state, const NCURI*, const char* suffix, NCbytes* packet);
static int readfiletofile(NCD4INFO* state, const NCURI*, const char* suffix, FILE* stream, d4size_t*);

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

int
NCD4_readDMR(NCD4INFO* state)
{
    int stat = NC_NOERR;
    long lastmodified = -1;

    stat = readpacket(state,state->uri,state->curl->packet,NCD4_DMR,&lastmodified);
    if(stat == NC_NOERR)
	state->data.dmrlastmodified = lastmodified;
    return THROW(stat);
}

int
NCD4_readDAP(NCD4INFO* state, int flags)
{
    int stat = NC_NOERR;
    long lastmod = -1;

    if((flags & NCF_ONDISK) == 0) {
        stat = readpacket(state,state->uri,state->curl->packet,NCD4_DAP,&lastmod);
        if(stat == NC_NOERR)
            state->data.daplastmodified = lastmod;
    } else { /*((flags & NCF_ONDISK) != 0) */
        NCURI* url = state->uri;
        int fileprotocol = (strcmp(url->protocol,"file")==0);
        if(fileprotocol) {
            stat = readfiletofile(state, url, ".dap", state->data.ondiskfile, &state->data.datasize);
        } else {
	    char* readurl = NULL;
            int flags = 0;
            if(!fileprotocol) flags |= NCURIQUERY;
            flags |= NCURIENCODE;
	    flags |= NCURIPWD;
#ifdef FIX
            ncurisetconstraints(url,state->constraint);
#endif
	    readurl = ncuribuild(url,NULL,".dods",NCURISVC);
	    if(readurl == NULL)
		return THROW(NC_ENOMEM);
            stat = NCD4_fetchurl_file(state->curl, readurl, state->data.ondiskfile,
                                   &state->data.datasize, &lastmod);
            nullfree(readurl);
            if(stat == NC_NOERR)
                state->data.daplastmodified = lastmod;
        }
    }
    return THROW(stat);
}

static const char*
dxxextension(int dxx)
{
    switch(dxx) {
    case NCD4_DMR: return ".dmr";
    case NCD4_DAP: return ".dap";
    default: break;
    }
    return NULL;
}

static int
readpacket(NCD4INFO* state, NCURI* url, NCbytes* packet, NCD4mode dxx, long* lastmodified)
{
    int stat = NC_NOERR;
    int fileprotocol = 0;
    const char* suffix = dxxextension(dxx);
    CURL* curl = state->curl->curl;

    fileprotocol = (strcmp(url->protocol,"file")==0);

    if(fileprotocol) {
	/* Short circuit file://... urls*/
	/* We do this because the test code always needs to read files*/
	stat = readfile(state, url,suffix,packet);
    } else {
        char* fetchurl = NULL;
	int flags = NCURIBASE;
	if(!fileprotocol) flags |= NCURIQUERY;
	flags |= NCURIENCODE;
        fetchurl = ncuribuild(url,NULL,suffix,flags);
	MEMCHECK(fetchurl);
	if(FLAGSET(state->controls.flags,NCF_SHOWFETCH)) {
	    nclog(NCLOGDBG,"fetch url=%s",fetchurl);
#ifdef HAVE_GETTIMEOFDAY
   	    gettimeofday(&time0,NULL);
#endif
	}
        stat = NCD4_fetchurl(curl,fetchurl,packet,lastmodified);
        nullfree(fetchurl);
	if(stat) goto fail;
	if(FLAGSET(state->controls.flags,NCF_SHOWFETCH)) {
            double secs = 0;
#ifdef HAVE_GETTIMEOFDAY
   	    gettimeofday(&time1,NULL);
	    secs = deltatime();
#endif
            nclog(NCLOGDBG,"fetch complete: %0.3f",secs);
	}
    }
#ifdef D4DEBUG
  {
fprintf(stderr,"readpacket: packet.size=%lu\n",
		(unsigned long)ncbyteslength(packet));
  }
#endif
fail:
    return THROW(stat);
}

static int
readfiletofile(NCD4INFO* state, const NCURI* uri, const char* suffix, FILE* stream, d4size_t* sizep)
{
    int stat = NC_NOERR;
    NCbytes* packet = ncbytesnew();
    size_t len;
    stat = readfile(state, uri,suffix,packet);
#ifdef D4DEBUG
fprintf(stderr,"readfiletofile: packet.size=%lu\n",
		(unsigned long)ncbyteslength(packet));
#endif
    if(stat != NC_NOERR) goto unwind;
    len = nclistlength(packet);
    if(stat == NC_NOERR) {
	size_t written;
        fseek(stream,0,SEEK_SET);
	written = fwrite(ncbytescontents(packet),1,len,stream);
	if(written != len) {
#ifdef D4DEBUG
fprintf(stderr,"readfiletofile: written!=length: %lu :: %lu\n",
	(unsigned long)written,(unsigned long)len);
#endif
	    stat = NC_EIO;
	}
    }
    if(sizep != NULL) *sizep = len;
unwind:
    ncbytesfree(packet);
    return THROW(stat);
}

static int
readfile(NCD4INFO* state, const NCURI* uri, const char* suffix, NCbytes* packet)
{
    int stat = NC_NOERR;
    NCbytes* tmp = ncbytesnew();
    char* filename = NULL;

    ncbytescat(tmp,uri->path);
    if(suffix != NULL) ncbytescat(tmp,suffix);
    ncbytesnull(tmp);
    filename = ncbytesextract(tmp);
    ncbytesfree(tmp);

    state->fileproto.filename = filename; /* filename is alloc'd here anyway */

    if(FLAGSET(state->controls.flags,NCF_SHOWFETCH)) {
	char* surl = NULL;
#ifdef HAVE_GETTIMEOFDAY
	gettimeofday(&time0,NULL);
#endif
        surl = ncuribuild((NCURI*)uri,NULL,NULL,NCURIALL);
	nclog(NCLOGDBG,"fetch uri=%s file=%s",surl,filename);
    }
    stat = NC_readfile(filename,packet);
    if(FLAGSET(state->controls.flags,NCF_SHOWFETCH)) {
	double secs;
#ifdef HAVE_GETTIMEOFDAY
   	gettimeofday(&time1,NULL);
	secs = deltatime();
#endif
        nclog(NCLOGDBG,"fetch complete: %0.3f",secs);
    }
    return THROW(stat);
}

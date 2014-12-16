/* Copyright 2009, UCAR/Unidata and OPeNDAP, Inc.
   See the COPYRIGHT file for more information. */

#include "config.h"
#ifdef HAVE_SYS_STAT_H
#include <sys/stat.h>
#endif
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#ifdef HAVE_FCNTL_H
#include <fcntl.h>
#endif
#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>
#endif
#ifdef HAVE_NETINET_IN_H
#include <netinet/in.h>
#endif
#ifdef _WIN32
#include <io.h>
#endif
#include "ocinternal.h"
#include "ocdebug.h"
#include "ochttp.h"
#include "ocread.h"
#include "ocrc.h"
#include "occurlfunctions.h"

/*Forward*/
static int readpacket(OCstate* state, OCURI*, OCbytes*, OCdxd, long*);
static int readfile(const char* path, const char* suffix, OCbytes* packet);
static int readfiletofile(const char* path, const char* suffix, FILE* stream, off_t*);

int
readDDS(OCstate* state, OCtree* tree)
{
    int stat = OC_NOERR;
    long lastmodified = -1;

    ocurisetconstraints(state->uri,tree->constraint);

#ifdef OCDEBUG
fprintf(stderr,"readDDS:\n");
#endif
    stat = readpacket(state,state->uri,state->packet,OCDDS,
			&lastmodified);
    if(stat == OC_NOERR) state->ddslastmodified = lastmodified;

    return stat;
}

int
readDAS(OCstate* state, OCtree* tree)
{
    int stat = OC_NOERR;

    ocurisetconstraints(state->uri,tree->constraint);
#ifdef OCDEBUG
fprintf(stderr,"readDAS:\n");
#endif
    stat = readpacket(state,state->uri,state->packet,OCDAS,NULL);

    return stat;
}

#if 0
int
readversion(OCstate* state, OCURI* url, OCbytes* packet)
{
   return readpacket(state,url,packet,OCVER,NULL);
}
#endif

const char*
ocdxdextension(OCdxd dxd)
{
    switch(dxd) {
    case OCDDS: return ".dds";
    case OCDAS: return ".das";
    case OCDATADDS: return ".dods";
    default: break;
    }
    return NULL;
}

static int
readpacket(OCstate* state, OCURI* url,OCbytes* packet,OCdxd dxd,long* lastmodified)
{
   int stat = OC_NOERR;
   int fileprotocol = 0;
   const char* suffix = ocdxdextension(dxd);
   char* fetchurl = NULL;
   CURL* curl = state->curl;

   fileprotocol = (strcmp(url->protocol,"file")==0);

   if(fileprotocol && !state->curlflags.proto_file) {
        /* Short circuit file://... urls*/
	/* We do this because the test code always needs to read files*/
	fetchurl = ocuribuild(url,NULL,NULL,0);
	stat = readfile(fetchurl,suffix,packet);
    } else {
	int flags = 0;
	if(!fileprotocol) {
	    flags |= OCURICONSTRAINTS;
	}
	flags |= OCURIENCODE;
        fetchurl = ocuribuild(url,NULL,suffix,flags);
	MEMCHECK(fetchurl,OC_ENOMEM);
	if(ocdebug > 0)
            {fprintf(stderr,"fetch url=%s\n",fetchurl); fflush(stderr);}
        stat = ocfetchurl(curl,fetchurl,packet,lastmodified,&state->creds);
	if(stat)
	    oc_curl_printerror(state);
	if(ocdebug > 0)
            {fprintf(stderr,"fetch complete\n"); fflush(stderr);}
    }
    free(fetchurl);
#ifdef OCDEBUG
  {
fprintf(stderr,"readpacket: packet.size=%lu\n",
		(unsigned long)ocbyteslength(packet));
  }
#endif
    return OCTHROW(stat);
}

int
readDATADDS(OCstate* state, OCtree* tree, OCflags flags)
{
    int stat = OC_NOERR;
    long lastmod = -1;

#ifdef OCDEBUG
fprintf(stderr,"readDATADDS:\n");
#endif
    if((flags & OCONDISK) == 0) {
        ocurisetconstraints(state->uri,tree->constraint);
        stat = readpacket(state,state->uri,state->packet,OCDATADDS,&lastmod);
        if(stat == OC_NOERR)
            state->datalastmodified = lastmod;
        tree->data.datasize = ocbyteslength(state->packet);
    } else { /*((flags & OCONDISK) != 0) */
        OCURI* url = state->uri;
        int fileprotocol = 0;
        char* readurl = NULL;

        fileprotocol = (strcmp(url->protocol,"file")==0);

        if(fileprotocol && !state->curlflags.proto_file) {
            readurl = ocuribuild(url,NULL,NULL,0);
            stat = readfiletofile(readurl, ".dods", tree->data.file, &tree->data.datasize);
        } else {
            int flags = 0;
            if(!fileprotocol) flags |= OCURICONSTRAINTS;
            flags |= OCURIENCODE;
	    flags |= OCURIUSERPWD;
            ocurisetconstraints(url,tree->constraint);
            readurl = ocuribuild(url,NULL,".dods",flags);
            MEMCHECK(readurl,OC_ENOMEM);
            if (ocdebug > 0) 
                {fprintf(stderr, "fetch url=%s\n", readurl);fflush(stderr);}
            stat = ocfetchurl_file(state->curl, readurl, tree->data.file,
                                   &tree->data.datasize, &lastmod);
            if(stat == OC_NOERR)
                state->datalastmodified = lastmod;
            if (ocdebug > 0) 
                {fprintf(stderr,"fetch complete\n"); fflush(stderr);}
        }
        free(readurl);
    }
    return OCTHROW(stat);
}

static int
readfiletofile(const char* path, const char* suffix, FILE* stream, off_t* sizep)
{
    int stat = OC_NOERR;
    OCbytes* packet = ocbytesnew();
    size_t len;
    /* check for leading file:/// */
    if(ocstrncmp(path,"file:///",8)==0) path += 7; /* assume absolute path*/
    stat = readfile(path,suffix,packet);
#ifdef OCDEBUG
fprintf(stderr,"readfiletofile: packet.size=%lu\n",
		(unsigned long)ocbyteslength(packet));
#endif
    if(stat != OC_NOERR) goto unwind;
    len = oclistlength(packet);
    if(stat == OC_NOERR) {
	size_t written;
        fseek(stream,0,SEEK_SET);
	written = fwrite(ocbytescontents(packet),1,len,stream);
	if(written != len) {
#ifdef OCDEBUG
fprintf(stderr,"readfiletofile: written!=length: %lu :: %lu\n",
	(unsigned long)written,(unsigned long)len);
#endif
	    stat = OC_EIO;
	}
    }
    if(sizep != NULL) *sizep = len;
unwind:
    ocbytesfree(packet);
    return OCTHROW(stat);
}

static int
readfile(const char* path, const char* suffix, OCbytes* packet)
{
    int stat = OC_NOERR;
    char buf[1024];
    char filename[1024];
    int fd = -1;
    int flags = 0;
    off_t filesize = 0;
    off_t totalread = 0;
    /* check for leading file:/// */
    if(ocstrncmp(path,"file://",7)==0) path += 7; /* assume absolute path*/
    if(!occopycat(filename,sizeof(filename),2,path,(suffix != NULL ? suffix : "")))
	return OCTHROW(OC_EOVERRUN);
    flags = O_RDONLY;
#ifdef O_BINARY
    flags |= O_BINARY;
#endif
    fd = open(filename,flags);
    if(fd < 0) {
	oclog(OCLOGERR,"open failed:%s",filename);
	return OCTHROW(OC_EOPEN);
    }
    /* Get the file size */
    filesize = lseek(fd,(off_t)0,SEEK_END);
    if(filesize < 0) {
	stat = OC_EIO;
	oclog(OCLOGERR,"lseek failed: %s",filename);
	goto done;
    }
    /* Move file pointer back to the beginning of the file */
    (void)lseek(fd,(off_t)0,SEEK_SET);
    stat = OC_NOERR;
    for(totalread=0;;) {
	off_t count = (off_t)read(fd,buf,sizeof(buf));
	if(count == 0)
	    break; /*eof*/
	else if(count <  0) {
	    stat = OC_EIO;
	    oclog(OCLOGERR,"read failed: %s",filename);
	    goto done;
	}
	ocbytesappendn(packet,buf,(unsigned long)count);
	totalread += count;
    }
    if(totalread < filesize) {
	stat = OC_EIO;
	oclog(OCLOGERR,"short read: |%s|=%lu read=%lu\n",
		filename,(unsigned long)filesize,(unsigned long)totalread);
        goto done;
    }

done:
#ifdef OCDEBUG
fprintf(stderr,"readfile: filesize=%lu totalread=%lu\n",
		(unsigned long)filesize,(unsigned long)totalread);
#endif
    if(fd >= 0) close(fd);
    return OCTHROW(stat);
}



/* Copyright 2009, UCAR/Unidata and OPeNDAP, Inc.
   See the COPYRIGHT file for more information. */

#include "config.h"
#include <sys/stat.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#include <fcntl.h>
#include <sys/types.h>
#ifdef HAVE_NETINET_IN_H
#include <netinet/in.h>
#endif
#include "ocinternal.h"
#include "ocdebug.h"
#include "ochttp.h"
#include "ocread.h"
#include "ocrc.h"
#include "occurlfunctions.h"

extern int oc_curl_file_supported;

/*Forward*/
static int readpacket(CURL*, OCURI*, OCbytes*, OCdxd, long*);
static int readfile(char* path, char* suffix, OCbytes* packet);
static int readfiletofile(char* path, char* suffix, FILE* stream, unsigned long*);

int
readDDS(OCstate* state, OCtree* tree)
{
    int stat = OC_NOERR;
    long lastmodified = -1;


    ocurisetconstraints(state->uri,tree->constraint);

ocset_user_password(state);

    stat = readpacket(state->curl,state->uri,state->packet,OCDDS,
			&lastmodified);
    if(stat == OC_NOERR) state->ddslastmodified = lastmodified;

    return stat;
}

int
readDAS(OCstate* state, OCtree* tree)
{
    int stat = OC_NOERR;

    ocurisetconstraints(state->uri,tree->constraint);
    stat = readpacket(state->curl,state->uri,state->packet,OCDAS,NULL);

    return stat;
}

int
readversion(CURL* curl, OCURI* url, OCbytes* packet)
{
   return readpacket(curl,url,packet,OCVER,NULL);
}

static
char* ocdxdextension[] ={
".dds",  /*OCDDS*/
".das",  /*OCDAS*/
".dods", /*OCDATADDS*/
".vers", /*OCVERS*/
};

static int
readpacket(CURL* curl,OCURI* url,OCbytes* packet,OCdxd dxd,long* lastmodified)
{
   int stat = OC_NOERR;
   int fileprotocol = 0;
   char* suffix = ocdxdextension[dxd];
   char* fetchurl = NULL;

   fileprotocol = (strcmp(url->protocol,"file")==0);

   if(fileprotocol && !oc_curl_file_supported) {
        /* Short circuit file://... urls*/
	/* We do this because the test code always needs to read files*/
	fetchurl = ocuribuild(url,NULL,NULL,0);
	stat = readfile(fetchurl,suffix,packet);
    } else {
	int flags = 0;
	if(!fileprotocol) flags |= OCURICONSTRAINTS;
	flags |= OCURIENCODE;
        fetchurl = ocuribuild(url,NULL,suffix,flags);
	MEMCHECK(fetchurl,OC_ENOMEM);
	if(ocdebug > 0)
            {fprintf(stderr,"fetch url=%s\n",fetchurl); fflush(stderr);}
        stat = ocfetchurl(curl,fetchurl,packet,lastmodified);
	if(ocdebug > 0)
            {fprintf(stderr,"fetch complete\n"); fflush(stderr);}
    }
    free(fetchurl);
    return OCTHROW(stat);
}

int
readDATADDS(OCstate* state, OCtree* tree, OCflags flags)
{
    int stat = OC_NOERR;
    long lastmod = -1;

    if((flags & OCONDISK) == 0) {
        ocurisetconstraints(state->uri,tree->constraint);
        stat = readpacket(state->curl,state->uri,state->packet,OCDATADDS,&lastmod);
        if(stat == OC_NOERR)
            state->datalastmodified = lastmod;
        tree->data.datasize = ocbyteslength(state->packet);
    } else {
        OCURI* url = state->uri;
        int fileprotocol = 0;
        char* readurl = NULL;

        fileprotocol = (strcmp(url->protocol,"file")==0);

        if(fileprotocol && !oc_curl_file_supported) {
            readurl = ocuribuild(url,NULL,NULL,0);
            stat = readfiletofile(readurl, ".dods", tree->data.file, &tree->data.datasize);
        } else {
            int flags = 0;
            if(!fileprotocol) flags |= OCURICONSTRAINTS;
            flags |= OCURIENCODE;
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
readfiletofile(char* path, char* suffix, FILE* stream, unsigned long* sizep)
{
    int stat = OC_NOERR;
    OCbytes* packet = ocbytesnew();
    size_t len;
    /* check for leading file:/// */
    if(ocstrncmp(path,"file:///",8)==0) path += 7; /* assume absolute path*/
    stat = readfile(path,suffix,packet);
    if(stat != OC_NOERR) goto unwind;
    len = oclistlength(packet);
    if(stat == OC_NOERR) {
	size_t written;
        fseek(stream,0,SEEK_SET);
	written = fwrite(ocbytescontents(packet),1,len,stream);
	if(written != len) stat = OC_EIO;
    }
    if(sizep != NULL) *sizep = len;
unwind:
    ocbytesfree(packet);
    return OCTHROW(stat);
}

static int
readfile(char* path, char* suffix, OCbytes* packet)
{
    int stat = OC_NOERR;
    char buf[1024];
    char filename[1024];
    int count,size,fd;
    /* check for leading file:/// */
    if(ocstrncmp(path,"file://",7)==0) path += 7; /* assume absolute path*/
    strcpy(filename,path);
    if(suffix != NULL) strcat(filename,suffix);
    fd = open(filename,O_RDONLY);
    if(fd < 0) {
	oc_log(LOGERR,"open failed:%s",filename);
	return OCTHROW(OC_EOPEN);
    }
    size=0;
    stat = OC_NOERR;
    for(;;) {
	count = read(fd,buf,sizeof(buf));
	if(count <= 0)
	    break;
	else if(count <  0) {
	    stat = OC_EIO;
	    oc_log(LOGERR,"read failed: %s",filename);
	    break;
	}
	ocbytesappendn(packet,buf,count);
	size += count;
    }
    close(fd);
    return OCTHROW(stat);
}



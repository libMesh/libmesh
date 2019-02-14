/*********************************************************************
 *   Copyright 2016, UCAR/Unidata
 *   See netcdf/COPYRIGHT file for copying and redistribution conditions.
 *********************************************************************/

#include "d4includes.h"
#include "d4chunk.h"

/**************************************************/

/* Header flags */
#define LAST_CHUNK          (1)
#define ERR_CHUNK           (2)
#define LITTLE_ENDIAN_CHUNK (4)
#define NOCHECKSUM_CHUNK    (8)

#define ALL_CHUNK_FLAGS (LAST_CHUNK|ERR_CHUNK|LITTLE_ENDIAN_CHUNK|NOCHECKSUM_CHUNK)

/**************************************************/

/*
Given a packet as read from the wire via http (or a file), convert in
place from chunked format to a single continguous set of bytes. If an
error packet is recovered, then make that available to the caller and
return an error. Also return whether the data was big endian encoded
and whether it has checksums.
Notes:
*/

/* Define a local struct for convenience */
struct HDR {unsigned int flags;	unsigned int count;};

/* Forward */
static void* getheader(void* p, struct HDR* hdr, int hostlittleendian);
static int processerrchunk(NCD4meta* metadata, void* errchunk, unsigned int count);

/**************************************************/

int
NCD4_dechunk(NCD4meta* metadata)
{
    unsigned char* p;
    unsigned char* q;
    struct HDR hdr;

    if(metadata->mode == NCD4_DSR) 
        return THROW(NC_EDMR);

    metadata->serial.errdata = NULL;
    metadata->serial.dmr = NULL;
    metadata->serial.dap = NULL;
    metadata->serial.hostlittleendian = NCD4_isLittleEndian();
    metadata->serial.remotelittleendian = 0; /* do not actually know yet */
    metadata->serial.remotechecksumming = 0; /* do not actually know yet */
    metadata->localchecksumming = 0; /* do not actually know yet */

    /* Assume proper mode has been inferred already. */

    /* Verify the mode; assume that the <?xml...?> is optional */
    q = metadata->serial.rawdata;
    if(memcmp(q,"<?xml",strlen("<?xml"))==0
       || memcmp(q,"<Dataset",strlen("<Dataset"))==0) {
	if(metadata->mode != NCD4_DMR) 
	    return THROW(NC_EDMR);
	/* setup as dmr only */
	metadata->serial.dmr = (char*)metadata->serial.rawdata; /* temp */
	metadata->serial.dmr[metadata->serial.rawsize-1] = '\0';
	metadata->serial.dmr = strdup((char *)q);
	if(metadata->serial.dmr == NULL)
	    return THROW(NC_ENOMEM);
	return THROW(NC_NOERR); 
    }

    /* We must be processing a DAP mode packet */
    p = metadata->serial.rawdata;
    metadata->serial.dap = p;

#ifdef D4DUMPRAW
    NCD4_tagdump(metadata->serial.rawsize,metadata->serial.rawdata,0,"RAW");
#endif

    /* Get the DMR chunk header*/
    p = getheader(p,&hdr,metadata->serial.hostlittleendian);
    if(hdr.count == 0)
	return THROW(NC_EDMR);
    if(hdr.flags & ERR_CHUNK) {
        return processerrchunk(metadata, (void*)p, hdr.count);
    }

    metadata->serial.remotechecksumming = ((hdr.flags & NOCHECKSUM_CHUNK) ? 0 : 1);
    metadata->localchecksumming = metadata->serial.remotechecksumming;

    metadata->serial.remotelittleendian = ((hdr.flags & LITTLE_ENDIAN_CHUNK) ? 1 : 0);
    metadata->serial.dmr = (char*)p;
    metadata->serial.dmr[hdr.count-1] = '\0';
    metadata->serial.dmr = strdup(metadata->serial.dmr);
    if(metadata->serial.dmr == NULL)
	return THROW(NC_ENOMEM);
    p += hdr.count;

    if(hdr.flags & LAST_CHUNK)
        return THROW(NC_ENODATA);
    /* Read and compress the data chunks */
    q = metadata->serial.dap;
    for(;;) {
	p = getheader(p,&hdr,metadata->serial.hostlittleendian);
	if(hdr.flags & ERR_CHUNK) {
            return processerrchunk(metadata, (void*)p, hdr.count);
	}
	/* data chunk; possibly last; possibly empty */
	if(hdr.count > 0) {
	    d4memmove(q,p,hdr.count); /* will overwrite the header */
	    p += hdr.count;
	    q += hdr.count;
	}
	if(hdr.flags & LAST_CHUNK) break;
    }
    metadata->serial.dapsize = (size_t)DELTA(q,metadata->serial.dap);
#ifdef D4DUMPDMR
    fprintf(stderr,"%s\n",metadata->serial.dmr);
    fflush(stderr);
#endif
#ifdef D4DUMPDAP
    NCD4_tagdump(metadata->serial.dapsize,metadata->serial.dap,0,"DAP");
#endif
    return THROW(NC_NOERR);    
}

static int
processerrchunk(NCD4meta* metadata, void* errchunk, unsigned int count)
{
    metadata->serial.errdata = (char*)d4alloc(count+1);
    if(metadata->serial.errdata == NULL)
        return THROW(NC_ENOMEM);
    memcpy(metadata->serial.errdata,errchunk,count);
    metadata->serial.errdata[count] = '\0';
    return THROW(NC_ENODATA); /* slight lie */
}

/* At the moment, the Hyrax test server
       is serving up the chunk data as little endian.
       So use a heuristic to see which endianness
       makes the more sense.
   This fails for very small dap data sections.
*/
static void*
getheader(void* p, struct HDR* hdr, int hostlittleendian)
{
    unsigned char bytes[4];
#ifdef HYRAXHACK
    struct HDR hyrax;
    unsigned char orig[4];
    memcpy(orig,p,sizeof(bytes));/* save a copy */
#endif
    memcpy(bytes,p,sizeof(bytes));
    p = INCR(p,4); /* on-the-wire hdr is 4 bytes */
    /* assume header is network (big) order */
    hdr->flags = bytes[0]; /* big endian => flags are in byte 0 */
    bytes[0] = 0; /* so we can do byte swap to get count */
    if(hostlittleendian)
        swapinline32(bytes); /* host is little endian */
    hdr->count = *(unsigned int*)bytes; /* get count */
#ifdef HYRAXHACK
    memcpy(bytes,orig,sizeof(bytes)); /* restore */
    hyrax.flags = bytes[3];
    bytes[3] = 0; /* so we can do byte swap to get count */
    if(!hostlittleendian)
        swapinline32(bytes); /* host is big endian */
    hyrax.count = *(unsigned int*)bytes; /* get count */
    /* See which makes more sense */
    if(hyrax.flags <= ALL_CHUNK_FLAGS && hyrax.count >= 0 && hyrax.count < hdr->count) {
	/* Use hyrax version */
	*hdr = hyrax;	
    }
#endif
    return p;
}

/**
Given a raw response, attempt to infer the mode: DMR, DAP, DSR.
Since DSR is not standardizes, it becomes the default.
*/
int
NCD4_infermode(NCD4meta* meta)
{
    d4size_t size = meta->serial.rawsize;
    char* raw = meta->serial.rawdata;

    if(size < 16)
	return THROW(NC_EDAP); /* must have at least this to hold a hdr + partial dmr*/	
    if(memcmp(raw,"<?xml",strlen("<?xml"))==0
       || memcmp(raw,"<Dataset",strlen("<Dataset"))==0) {
	meta->mode = NCD4_DMR;
	goto done;
    }
    raw += 4; /* Pretend we have a DAP hdr */
    if(memcmp(raw,"<?xml",strlen("<?xml"))==0
       || memcmp(raw,"<Dataset",strlen("<Dataset"))==0) {
	meta->mode = NCD4_DAP;
	goto done;
    }
    /* Default to DSR */
    meta->mode = NCD4_DSR;

done:
    return NC_NOERR;
}

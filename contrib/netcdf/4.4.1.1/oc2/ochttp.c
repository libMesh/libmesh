/* Copyright 2009, UCAR/Unidata and OPeNDAP, Inc.
 See the COPYRIGHT file for more information. */

#include "config.h"
#include <sys/stat.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#include <fcntl.h>
#include "ocinternal.h"
#include "ocdebug.h"
#include "ochttp.h"

static size_t WriteFileCallback(void*, size_t, size_t, void*);
static size_t WriteMemoryCallback(void*, size_t, size_t, void*);

struct Fetchdata {
	FILE* stream;
	size_t size;
};

long
ocfetchhttpcode(CURL* curl)
{
    long httpcode = 200;
    CURLcode cstat = CURLE_OK;
    /* Extract the http code */
#ifdef HAVE_CURLINFO_RESPONSE_CODE
    cstat = CURLERR(curl_easy_getinfo(curl,CURLINFO_RESPONSE_CODE,&httpcode));
#else
    cstat = curl_easy_getinfo(curl,CURLINFO_HTTP_CODE,&httpcode);
#endif
    if(cstat != CURLE_OK) httpcode = 0;
    return httpcode;
}

OCerror
ocfetchurl_file(CURL* curl, const char* url, FILE* stream,
		off_t* sizep, long* filetime)
{
	int stat = OC_NOERR;
	CURLcode cstat = CURLE_OK;
	struct Fetchdata fetchdata;

	/* Set the URL */
	cstat = CURLERR(curl_easy_setopt(curl, CURLOPT_URL, (void*)url));
	if (cstat != CURLE_OK)
		goto fail;

	/* send all data to this function  */
	cstat = CURLERR(curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, WriteFileCallback));
	if (cstat != CURLE_OK)
		goto fail;

	/* we pass our file to the callback function */
	cstat = CURLERR(curl_easy_setopt(curl, CURLOPT_WRITEDATA, (void *)&fetchdata));
	if (cstat != CURLE_OK)
		goto fail;

        /* One last thing; always try to get the last modified time */
        cstat = CURLERR(curl_easy_setopt(curl, CURLOPT_FILETIME, (long)1));
	if (cstat != CURLE_OK)
		goto fail;

	fetchdata.stream = stream;
	fetchdata.size = 0;
	cstat = CURLERR(curl_easy_perform(curl));

	if (cstat != CURLE_OK)
	    goto fail;

	if (stat == OC_NOERR) {
	    /* return the file size*/
#ifdef OCDEBUG
	    oclog(OCLOGNOTE,"filesize: %lu bytes",fetchdata.size);
#endif
	    if (sizep != NULL)
		*sizep = fetchdata.size;
	    /* Get the last modified time */
	    if(filetime != NULL)
                cstat = curl_easy_getinfo(curl,CURLINFO_FILETIME,filetime);
            if(cstat != CURLE_OK) goto fail;
	}
	return OCTHROW(stat);

fail:
	oclog(OCLOGERR, "curl error: %s", curl_easy_strerror(cstat));
	return OCTHROW(OC_ECURL);
}

OCerror
ocfetchurl(CURL* curl, const char* url, OCbytes* buf, long* filetime,
           struct OCcredentials* creds)
{
	OCerror stat = OC_NOERR;
	CURLcode cstat = CURLE_OK;
	size_t len;
        long httpcode = 0;

	/* Set the URL */
	cstat = CURLERR(CURLERR(curl_easy_setopt(curl, CURLOPT_URL, (void*)url)));
	if (cstat != CURLE_OK)
		goto fail;
	
#if 0
	if(creds != NULL && creds->password != NULL  && creds->username != NULL) {
	    /* Set user and password */
#if defined (HAVE_CURLOPT_USERNAME) && defined (HAVE_CURLOPT_PASSWORD)
	    cstat = CURLERR(curl_easy_setopt(curl, CURLOPT_USERNAME, creds->username));
	    if (cstat != CURLE_OK)
		goto fail;
	    cstat = CURLERR(curl_easy_setopt(curl, CURLOPT_PASSWORD, creds->password));
	    if (cstat != CURLE_OK)
		goto fail;
#else		
		snprintf(tbuf,1023,"%s:%s",creds->username,creds->password);	
		cstat = CURLERR(curl_easy_setopt(curl, CURLOPT_USERPWD, tbuf));
		if (cstat != CURLE_OK)
			goto fail;
#endif
	}
#endif

	/* send all data to this function  */
	cstat = CURLERR(curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, WriteMemoryCallback));
	if (cstat != CURLE_OK)
		goto fail;

	/* we pass our file to the callback function */
	cstat = CURLERR(curl_easy_setopt(curl, CURLOPT_WRITEDATA, (void*)buf));
	if (cstat != CURLE_OK)
		goto fail;

        /* One last thing; always try to get the last modified time */
	cstat = CURLERR(curl_easy_setopt(curl, CURLOPT_FILETIME, (long)1));

	cstat = CURLERR(curl_easy_perform(curl));

	if(cstat == CURLE_PARTIAL_FILE) {
	    /* Log it but otherwise ignore */
	    oclog(OCLOGWARN, "curl error: %s; ignored",
		   curl_easy_strerror(cstat));
	    cstat = CURLE_OK;
	}
        httpcode = ocfetchhttpcode(curl);

	if(cstat != CURLE_OK) goto fail;

        /* Get the last modified time */
	if(filetime != NULL)
            cstat = CURLERR(curl_easy_getinfo(curl,CURLINFO_FILETIME,filetime));
        if(cstat != CURLE_OK) goto fail;

	/* Null terminate the buffer*/
	len = ocbyteslength(buf);
	ocbytesappend(buf, '\0');
	ocbytessetlength(buf, len); /* dont count null in buffer size*/
#ifdef OCDEBUG
	oclog(OCLOGNOTE,"buffersize: %lu bytes",(off_t)ocbyteslength(buf));
#endif

	return OCTHROW(stat);

fail:
	oclog(OCLOGERR, "curl error: %s", curl_easy_strerror(cstat));
	switch (httpcode) {
	case 401: stat = OC_EAUTH; break;
	case 404: stat = OC_ENOFILE; break;
	case 500: stat = OC_EDAPSVC; break;
	case 200: break;
	default: stat = OC_ECURL; break;
	}
	return OCTHROW(stat);
}

static size_t
WriteFileCallback(void* ptr, size_t size, size_t nmemb,	void* data)
{
	size_t realsize = size * nmemb;
	size_t count;
	struct Fetchdata* fetchdata;
	fetchdata = (struct Fetchdata*) data;
        if(realsize == 0)
	    oclog(OCLOGWARN,"WriteFileCallback: zero sized chunk");
	count = fwrite(ptr, size, nmemb, fetchdata->stream);
	if (count > 0) {
		fetchdata->size += (count * size);
	} else {
	    oclog(OCLOGWARN,"WriteFileCallback: zero sized write");
	}
#ifdef OCPROGRESS
        oclog(OCLOGNOTE,"callback: %lu bytes",(off_t)realsize);
#endif
	return count;
}

static size_t
WriteMemoryCallback(void *ptr, size_t size, size_t nmemb, void *data)
{
	size_t realsize = size * nmemb;
	OCbytes* buf = (OCbytes*) data;
        if(realsize == 0)
	    oclog(OCLOGWARN,"WriteMemoryCallback: zero sized chunk");
	/* Optimize for reading potentially large dods datasets */
	if(!ocbytesavail(buf,realsize)) {
	    /* double the size of the packet */
	    ocbytessetalloc(buf,2*ocbytesalloc(buf));
	}
	ocbytesappendn(buf, ptr, realsize);
#ifdef OCPROGRESS
        oclog(OCLOGNOTE,"callback: %lu bytes",(off_t)realsize);
#endif
	return realsize;
}

#if 0
static void
assembleurl(DAPURL* durl, OCbytes* buf, int what)
{
	encodeurltext(durl->url,buf);
	if(what & WITHPROJ) {
		ocbytescat(buf,"?");
		encodeurltext(durl->projection,buf);
	}
	if(what & WITHSEL) encodeurltext(durl->selection,buf);

}

static char mustencode="";
static char hexchars[16] = {
	'0', '1', '2', '3',
	'4', '5', '6', '7',
	'8', '9', 'a', 'b',
	'c', 'd', 'e', 'f',
};

static void
encodeurltext(char* text, OCbytes* buf)
{
	/* Encode the URL to handle illegal characters */
	len = strlen(url);
	encoded = ocmalloc(len*4+1); /* should never be larger than this*/
	if(encoded==NULL) return;
	p = url; q = encoded;
	while((c=*p++)) {
		if(strchr(mustencode,c) != NULL) {
			char tmp[8];
			int hex1, hex2;
			hex1 = (c & 0x0F);
			hex2 = (c & 0xF0) >> 4;
			tmp[0] = '0'; tmp[1] = 'x';
			tmp[2] = hexchars[hex2]; tmp[3] = hexchars[hex1];
			tmp[4] = '\0';
			ocbytescat(buf,tmp);
		} else *q++ = (char)c;
	}

}

#endif

OCerror
occurlopen(CURL** curlp)
{
	int stat = OC_NOERR;
	CURLcode cstat = CURLE_OK;
	CURL* curl;
	/* initialize curl*/
	curl = curl_easy_init();
	if (curl == NULL)
		stat = OC_ECURL;
	else {
		cstat = CURLERR(curl_easy_setopt(curl, CURLOPT_NOPROGRESS, 1));
		if (cstat != CURLE_OK)
			stat = OC_ECURL;
	}
	if (curlp)
		*curlp = curl;
	return OCTHROW(stat);
}

void
occurlclose(CURL* curl)
{
	if (curl != NULL)
		curl_easy_cleanup(curl);
}

OCerror
ocfetchlastmodified(CURL* curl, char* url, long* filetime)
{
    int stat = OC_NOERR;
    CURLcode cstat = CURLE_OK;

    /* Set the URL */
    cstat = CURLERR(curl_easy_setopt(curl, CURLOPT_URL, (void*)url));
    if (cstat != CURLE_OK)
        goto fail;

    /* Ask for head */
    cstat = CURLERR(curl_easy_setopt(curl, CURLOPT_TIMEOUT, 30)); /* 30sec timeout*/
    cstat = CURLERR(curl_easy_setopt(curl, CURLOPT_CONNECTTIMEOUT, 2));
    cstat = CURLERR(curl_easy_setopt(curl, CURLOPT_HEADER, 1));
    cstat = CURLERR(curl_easy_setopt(curl, CURLOPT_NOBODY, 1));
    cstat = CURLERR(curl_easy_setopt(curl, CURLOPT_NOPROGRESS, 1));
    cstat = CURLERR(curl_easy_setopt(curl, CURLOPT_FILETIME, (long)1));

    cstat = CURLERR(curl_easy_perform(curl));
    if(cstat != CURLE_OK) goto fail;
    if(filetime != NULL)
        cstat = CURLERR(curl_easy_getinfo(curl,CURLINFO_FILETIME,filetime));
    if(cstat != CURLE_OK) goto fail;

    return OCTHROW(stat);

fail:
    oclog(OCLOGERR, "curl error: %s", curl_easy_strerror(cstat));
    return OCTHROW(OC_ECURL);
}

OCerror
ocping(const char* url)
{
    int stat = OC_NOERR;
    CURLcode cstat = CURLE_OK;
    CURL* curl = NULL;
    OCbytes* buf = NULL;

    /* Create a CURL instance */
    stat = occurlopen(&curl);
    if(stat != OC_NOERR) return stat;    

    /* Use redirects */
    cstat = CURLERR(curl_easy_setopt(curl, CURLOPT_MAXREDIRS, 10L));
    if (cstat != CURLE_OK)
        goto done;
    cstat = CURLERR(curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L));
    if (cstat != CURLE_OK)
        goto done;

    /* use a very short timeout: 10 seconds */
    cstat = CURLERR(curl_easy_setopt(curl, CURLOPT_TIMEOUT, (long)10));
    if (cstat != CURLE_OK)
        goto done;

    /* fail on HTTP 400 code errors */
    cstat = CURLERR(curl_easy_setopt(curl, CURLOPT_FAILONERROR, (long)1));
    if (cstat != CURLE_OK)
        goto done;

    /* Try to get the file */
    buf = ocbytesnew();
    stat = ocfetchurl(curl,url,buf,NULL,NULL);
    if(stat == OC_NOERR) {
	/* Don't trust curl to return an error when request gets 404 */
	long http_code = 0;
	cstat = CURLERR(curl_easy_getinfo(curl,CURLINFO_RESPONSE_CODE, &http_code));
        if (cstat != CURLE_OK)
            goto done;
	if(http_code >= 400) {
	    cstat = CURLE_HTTP_RETURNED_ERROR;
	    goto done;
	}
    } else
        goto done;

done:
    ocbytesfree(buf);
    occurlclose(curl);
    if(cstat != CURLE_OK) {
        oclog(OCLOGERR, "curl error: %s", curl_easy_strerror(cstat));
        stat = OC_EDAPSVC;
    }
    return OCTHROW(stat);
}

/* Copyright 2009, UCAR/Unidata and OPeNDAP, Inc.
   See the COPYRIGHT file for more information. */

#include "config.h"
#include "ncrc.h"
#include "ocinternal.h"
#include "ocdebug.h"
#include "occurlfunctions.h"

#define OC_MAX_REDIRECTS 20L

/* Mnemonic */
#define OPTARG void*

/* Define some .rc file entries of interest*/
#define NETRCFILETAG "HTTP.NETRC"

/* Check return value */
#define CHECK(state,flag,value) {if(check(state,flag,(void*)value) != OC_NOERR) {goto done;}}

static OCerror
check(OCstate* state, int flag, void* value)
{
    OCerror stat = ocset_curlopt(state,flag,value);
    return stat;
}

/*
Set a specific curl flag; primary wrapper for curl_easy_setopt
*/
OCerror
ocset_curlopt(OCstate* state, int flag, void* value)
{
    OCerror stat = OC_NOERR;
    CURLcode cstat = CURLE_OK;
    cstat = OCCURLERR(state,curl_easy_setopt(state->curl,flag,value));
    if(cstat != CURLE_OK)
	stat = OC_ECURL;
    return stat;
}

/*
Update a specific flag from state
*/
OCerror
ocset_curlflag(OCstate* state, int flag)
{
    OCerror stat = OC_NOERR;

    switch (flag) {

    case CURLOPT_USERPWD: /* Does both user and pwd */
        if(state->auth.creds.user != NULL && state->auth.creds.pwd != NULL) {
	    CHECK(state, CURLOPT_USERNAME, state->auth.creds.user);
	    CHECK(state, CURLOPT_PASSWORD, state->auth.creds.pwd);
            CHECK(state, CURLOPT_HTTPAUTH, (OPTARG)CURLAUTH_ANY);
	}
	break;

    case CURLOPT_COOKIEJAR: case CURLOPT_COOKIEFILE:
        if(state->auth.curlflags.cookiejar) {
	    /* Assume we will read and write cookies to same place */
	    CHECK(state, CURLOPT_COOKIEJAR, state->auth.curlflags.cookiejar);
	    CHECK(state, CURLOPT_COOKIEFILE, state->auth.curlflags.cookiejar);
        }
	break;

    case CURLOPT_NETRC: case CURLOPT_NETRC_FILE:
	if(state->auth.curlflags.netrc) {
	    CHECK(state, CURLOPT_NETRC, (OPTARG)CURL_NETRC_REQUIRED);
	    CHECK(state, CURLOPT_NETRC_FILE, state->auth.curlflags.netrc);
        }
	break;

    case CURLOPT_VERBOSE:
	if(state->auth.curlflags.verbose)
	    CHECK(state, CURLOPT_VERBOSE, (OPTARG)1L);
	break;

    case CURLOPT_TIMEOUT:
	if(state->auth.curlflags.timeout)
	    CHECK(state, CURLOPT_TIMEOUT, (OPTARG)((long)state->auth.curlflags.timeout));
	break;

    case CURLOPT_USERAGENT:
        if(state->auth.curlflags.useragent)
	    CHECK(state, CURLOPT_USERAGENT, state->auth.curlflags.useragent);
	break;

    case CURLOPT_FOLLOWLOCATION:
        CHECK(state, CURLOPT_FOLLOWLOCATION, (OPTARG)1L);
	break;

    case CURLOPT_MAXREDIRS:
	CHECK(state, CURLOPT_MAXREDIRS, (OPTARG)OC_MAX_REDIRECTS);
	break;

    case CURLOPT_ERRORBUFFER:
	CHECK(state, CURLOPT_ERRORBUFFER, state->error.curlerrorbuf);
	break;

    case CURLOPT_ENCODING:
#ifdef CURLOPT_ENCODING
	if(state->auth.curlflags.compress) {
	    CHECK(state, CURLOPT_ENCODING,"deflate, gzip");
        }
#endif
	break;

    case CURLOPT_PROXY:
	if(state->auth.proxy.host != NULL) {
	    CHECK(state, CURLOPT_PROXY, state->auth.proxy.host);
	    CHECK(state, CURLOPT_PROXYPORT, (OPTARG)(long)state->auth.proxy.port);
	    if(state->auth.proxy.user != NULL && state->auth.proxy.pwd != NULL) {
                CHECK(state, CURLOPT_PROXYUSERNAME, state->auth.proxy.user);
                CHECK(state, CURLOPT_PROXYPASSWORD, state->auth.proxy.pwd);
#ifdef CURLOPT_PROXYAUTH
	        CHECK(state, CURLOPT_PROXYAUTH, (long)CURLAUTH_ANY);
#endif
	    }
	}
	break;

    case CURLOPT_USE_SSL:
    case CURLOPT_SSLCERT: case CURLOPT_SSLKEY:
    case CURLOPT_SSL_VERIFYPEER: case CURLOPT_SSL_VERIFYHOST:
    {
        struct ssl* ssl = &state->auth.ssl;
        CHECK(state, CURLOPT_SSL_VERIFYPEER, (OPTARG)(ssl->verifypeer?1L:0L));
        CHECK(state, CURLOPT_SSL_VERIFYHOST, (OPTARG)(ssl->verifyhost?1L:0L));
        if(ssl->certificate)
            CHECK(state, CURLOPT_SSLCERT, ssl->certificate);
        if(ssl->key)
            CHECK(state, CURLOPT_SSLKEY, ssl->key);
        if(ssl->keypasswd)
            /* libcurl prior to 7.16.4 used 'CURLOPT_SSLKEYPASSWD' */
            CHECK(state, CURLOPT_KEYPASSWD, ssl->keypasswd);
        if(ssl->cainfo)
            CHECK(state, CURLOPT_CAINFO, ssl->cainfo);
        if(ssl->capath)
            CHECK(state, CURLOPT_CAPATH, ssl->capath);
    }
    break;

#ifdef HAVE_CURLOPT_BUFFERSIZE
    case CURLOPT_BUFFERSIZE:
	CHECK(state, CURLOPT_BUFFERSIZE, (OPTARG)state->curlbuffersize);
	break;
#endif

#ifdef HAVE_CURLOPT_KEEPALIVE
    case CURLOPT_TCP_KEEPALIVE:
	if(state->curlkeepalive.active != 0)
	    CHECK(state, CURLOPT_TCP_KEEPALIVE, (OPTARG)1L);
	if(state->curlkeepalive.idle > 0)
	    CHECK(state, CURLOPT_TCP_KEEPIDLE, (OPTARG)state->curlkeepalive.idle);
	if(state->curlkeepalive.interval > 0)
	    CHECK(state, CURLOPT_TCP_KEEPINTVL, (OPTARG)state->curlkeepalive.interval);
	break;
#endif

    default:
        nclog(NCLOGWARN,"Attempt to update unexpected curl flag: %d",flag);
	break;
    }
done:
    return stat;
}


/* Set various general curl flags per fetch  */
OCerror
ocset_flags_perfetch(OCstate* state)
{
    OCerror stat = OC_NOERR;
    /* currently none */
    return stat;
}

/* Set various general curl flags per link */

OCerror
ocset_flags_perlink(OCstate* state)
{
    OCerror stat = OC_NOERR;

    /* Following are always set */
    if(stat == OC_NOERR) stat = ocset_curlflag(state,CURLOPT_ENCODING);
    if(stat == OC_NOERR) stat = ocset_curlflag(state,CURLOPT_NETRC);
    if(stat == OC_NOERR) stat = ocset_curlflag(state,CURLOPT_VERBOSE);
    if(stat == OC_NOERR) stat = ocset_curlflag(state,CURLOPT_TIMEOUT);
    if(stat == OC_NOERR) stat = ocset_curlflag(state,CURLOPT_USERAGENT);
    if(stat == OC_NOERR) stat = ocset_curlflag(state,CURLOPT_COOKIEJAR);
    if(stat == OC_NOERR) stat = ocset_curlflag(state,CURLOPT_USERPWD);
    if(stat == OC_NOERR) stat = ocset_curlflag(state,CURLOPT_PROXY);
    if(stat == OC_NOERR) stat = ocset_curlflag(state,CURLOPT_USE_SSL);
    if(stat == OC_NOERR) stat = ocset_curlflag(state, CURLOPT_FOLLOWLOCATION);
    if(stat == OC_NOERR) stat = ocset_curlflag(state, CURLOPT_MAXREDIRS);
    if(stat == OC_NOERR) stat = ocset_curlflag(state, CURLOPT_ERRORBUFFER);

#ifdef HAVE_CURLOPT_BUFFERSIZE
    /* Optional */
    if(stat == OC_NOERR && state->curlbuffersize > 0)
	stat = ocset_curlflag(state, CURLOPT_BUFFERSIZE);
#endif
#ifdef HAVE_CURLOPT_KEEPALIVE
    if(stat == NC_NOERR && state->curlkeepalive.active != 0)
        stat = ocset_curlflag(state, CURLOPT_TCP_KEEPALIVE);
#endif
    return stat;
}

void
oc_curl_debug(OCstate* state)
{
    state->auth.curlflags.verbose = 1;
    ocset_curlflag(state,CURLOPT_VERBOSE);
    ocset_curlflag(state,CURLOPT_ERRORBUFFER);
}

/* Misc. */

int
ocrc_netrc_required(OCstate* state)
{
    char* netrcfile = NC_rclookup(NETRCFILETAG,state->uri->uri);
    return (netrcfile != NULL || state->auth.curlflags.netrc != NULL ? 0 : 1);
}

void
oc_curl_printerror(OCstate* state)
{
    fprintf(stderr,"curl error details: %s\n",state->curlerror);
}

/* See if http: protocol is supported */
void
oc_curl_protocols(OCstate* state)
{
    const char* const* proto; /*weird*/
    curl_version_info_data* curldata;
    curldata = curl_version_info(CURLVERSION_NOW);
    for(proto=curldata->protocols;*proto;proto++) {
        if(strcmp("http",*proto)==0)
	    state->auth.curlflags.proto_https=1;
    }
}

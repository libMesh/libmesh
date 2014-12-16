/* Copyright 2009, UCAR/Unidata and OPeNDAP, Inc.
   See the COPYRIGHT file for more information. */

#include "config.h"
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ocinternal.h"
#include "ocdebug.h"
#include "oclog.h"

#include "ocrc.h"

#define RTAG ']'
#define LTAG '['

#define TRIMCHARS " \t\r\n"

#define HTTPPREFIXDEPRECATED "CURL."
#define HTTPPREFIX           "HTTP."

static int parseproxy(OCstate* state, char* v);
static int rcreadline(FILE* f, char* more, int morelen);
static void rctrim(char* text);

static void ocdodsrcdump(char* msg, struct OCTriple*, int ntriples);

static char* curllookup(char* suffix,char* url);

/* The Username and password are in the URL if the URL is of the form:
 * http://<name>:<passwd>@<host>/....
 */
static int
occredentials_in_url(const char *url)
{
	char *pos = strstr(url, "http://");
	if (!pos)
		return 0;
	pos += 7;
	if (strchr(pos, '@') && strchr(pos, ':'))
		return 1;

	return 0;
}

static int
ocextract_credentials(const char *url, char **name, char **pw, char **result_url)
{
	char *pos;
	char *end;
	char *middle;
	size_t up_len = 0;
	size_t mid_len = 0;
	size_t midpas_len = 0;
	size_t url_len = 0;

	if (strchr(url, '@')) {
		pos = strstr(url, "http://");
		if (pos)
			pos += 7;
		middle = strchr(pos, ':');
		mid_len = middle - pos;
		*name = malloc(sizeof(char) * (mid_len + 1));
		strncpy(*name, pos, mid_len);
		(*name)[mid_len] = '\0';

		if (middle)
			middle += 1;

		end = strchr(middle, '@');
		midpas_len = end - middle;
		*pw = malloc(sizeof(char) * (midpas_len + 1));
		strncpy(*pw, middle, midpas_len);
		(*pw)[midpas_len] = '\0';

		up_len = end - pos;
		url_len = strlen(url) - up_len;

		*result_url = malloc(sizeof(char) * (url_len + 1));
		if(*result_url == NULL)
			return OC_ENOMEM;

		strncpy(*result_url, url, (size_t)(pos - url));
		strncpy(*result_url + (pos - url), end + 1, url_len - (pos - url));

#if 0
		fprintf(stderr, "URL without username and password: %s:%d\n", sURL, url_len );
		fprintf(stderr, "URL username and password: %s:%d\n", sUP, up_len);
		fprintf(stderr, "URL username: %s:%d\n", sUser, mid_len);
		fprintf(stderr, "URL password: %s:%d\n", sPassword, midpas_len);
#endif
		(*result_url)[url_len] = '\0';

		return OC_NOERR;
	}
	else {
		return OC_EIO;
	}
}

static int
rcreadline(FILE* f, char* more, int morelen)
{
    int i = 0;
    int c = getc(f);
    if(c < 0) return 0;
    for(;;) {
	if(i < morelen)  /* ignore excess characters */
	    more[i++]=c;
	c = getc(f);
	if(c < 0) break; /* eof */
        if(c == '\n') break; /* eol */
    }
    /* null terminate more */
    more[i] = '\0';
    return 1;
}


/* Trim TRIMCHARS from both ends of text; */
static void
rctrim(char* text)
{
    char* p = text;
    size_t len = strlen(text);
    int i;
    /* locate first non-trimchar */
    for(;*p;p++) {
       if(strchr(TRIMCHARS,*p) == NULL) break; /* hit non-trim char */
    }
    strncpy(text,p,len);
    len = strlen(text);
    /* locate last non-trimchar */
    if(len > 0) {
	for(i=(len-1);i>=0;i--) {
            if(strchr(TRIMCHARS,text[i]) == NULL) {
		text[i+1] = '\0'; /* elide trailing trimchars */
		break;
	    }
	}
    }
}

static int
parseproxy(OCstate* state, char* v)
{
    char *host_pos = NULL;
    char *port_pos = NULL;

    if(strlen(v) == 0) return OC_NOERR; /* nothing there*/
    if (occredentials_in_url(v)) {
        char *result_url = NULL;
        ocextract_credentials(v, &state->creds.username,
                            &state->creds.password,
                            &result_url);
        v = result_url;
    }
    /* allocating a bit more than likely needed ... */
    host_pos = strstr(v, "http://");
    if (host_pos)
        host_pos += strlen("http://");
    else
	host_pos = v;
    port_pos = strchr(host_pos, ':');
    if (port_pos) {
        size_t host_len;
        char *port_sep = port_pos;
        port_pos++;
        *port_sep = '\0';
        host_len = strlen(host_pos);
        state->proxy.host = malloc(sizeof(char) * host_len + 1);
        if(state->proxy.host == NULL)
            return OC_ENOMEM;

        strncpy(state->proxy.host, host_pos, host_len);
        state->proxy.host[host_len] = '\0';

        state->proxy.port = atoi(port_pos);
    } else {
        size_t host_len = strlen(host_pos);
        state->proxy.host = malloc(sizeof(char) * host_len + 1);
        if(state->proxy.host == NULL)
            return OC_ENOMEM;

        strncpy(state->proxy.host, host_pos, host_len);
        state->proxy.host[host_len] = '\0';

        state->proxy.port = 80;
    }
#if 0
    state->proxy.host[v_len] = '\0';
    state->proxy.port = atoi(v);
    s_len = strlen(v);
    state->proxy.user = malloc(sizeof(char) * s_len + 1);
    if(state->proxy.user == NULL)
        return OC_ENOMEM;
     strncpy(state->proxy.user, v, s_len);
     state->proxy.user[s_len] = '\0';
     p_len = strlen(v);
     state->proxy.password = malloc(sizeof(char) * p_len + 1);
     if(state->proxy.password == NULL)
         return OC_ENOMEM;
     strncpy(state->proxy.password, v, p_len);
     state->proxy.password[p_len] = '\0';
#endif /*0*/
     if (ocdebug > 1) {
         oclog(OCLOGNOTE,"host name: %s", state->proxy.host);
         oclog(OCLOGNOTE,"user name: %s", state->creds.username);
#ifdef INSECURE
         oclog(OCLOGNOTE,"password: %s", state->creds.password);
#endif
         oclog(OCLOGNOTE,"port number: %d", state->proxy.port);
    }
    if(v) free(v);
    return OC_NOERR;
}

/* insertion sort the triplestore based on url */
static void
sorttriplestore(void)
{
    int i, nsorted;
    struct OCTriple* sorted = NULL;
    struct OCTriplestore* ocdodsrc = ocglobalstate.ocdodsrc;

    if(ocdodsrc == NULL) return; /* nothing to sort */
    if(ocdodsrc->ntriples <= 1) return; /* nothing to sort */
    if(ocdebug > 2)
        ocdodsrcdump("initial:",ocdodsrc->triples,ocdodsrc->ntriples);

    sorted = (struct OCTriple*)malloc(sizeof(struct OCTriple)*ocdodsrc->ntriples);
    if(sorted == NULL) {
        oclog(OCLOGERR,"sorttriplestore: out of memory");
        return;
    }

    nsorted = 0;
    while(nsorted < ocdodsrc->ntriples) {
	int largest;
	/* locate first non killed entry */
	for(largest=0;largest<ocdodsrc->ntriples;largest++) {
            if(ocdodsrc->triples[largest].key[0] != '\0') break;
	}
        OCASSERT(ocdodsrc->triples[largest].key[0] != '\0');
	for(i=0;i<ocdodsrc->ntriples;i++) {
	    if(ocdodsrc->triples[i].key[0] != '\0') { /* avoid empty slots */
	        int lexorder = strcmp(ocdodsrc->triples[i].url,ocdodsrc->triples[largest].url);
   	        int leni = strlen(ocdodsrc->triples[i].url);
 	        int lenlarge = strlen(ocdodsrc->triples[largest].url);
	        /* this defines the ordering */
	        if(leni == 0 && lenlarge == 0) continue; /* if no urls, then leave in order */
	        if(leni != 0 && lenlarge == 0) largest = i;
	        else if(lexorder > 0) largest = i;
	    }
	}
	/* Move the largest entry */
	OCASSERT(ocdodsrc->triples[largest].key[0] != 0);
	sorted[nsorted] = ocdodsrc->triples[largest];
	ocdodsrc->triples[largest].key[0] = '\0'; /* kill entry */
	nsorted++;
      if(ocdebug > 2)
            ocdodsrcdump("pass:",sorted,nsorted);
    }    

    memcpy((void*)ocdodsrc->triples,(void*)sorted,sizeof(struct OCTriple)*nsorted);
    free(sorted);

    if(ocdebug > 0)
	ocdodsrcdump("final .dodsrc order:",ocdodsrc->triples,ocdodsrc->ntriples);
}

/* Create a triple store from a file */
int
ocdodsrc_read(char* basename, char* path)
{
    char line0[MAXRCLINESIZE];
    FILE *in_file = NULL;
    int linecount = 0;
    struct OCTriplestore* ocdodsrc = ocglobalstate.ocdodsrc;

    if(ocdodsrc == NULL) {
        ocdodsrc = (struct OCTriplestore*)malloc(sizeof(struct OCTriplestore));
        if(ocdodsrc == NULL) {
	    oclog(OCLOGERR,"ocdodsrc_read: out of memory");
	    return 0;
	}
        ocglobalstate.ocdodsrc = ocdodsrc;
    }
    ocdodsrc->ntriples = 0;

    in_file = fopen(path, "r"); /* Open the file to read it */
    if (in_file == NULL) {
	oclog(OCLOGERR, "Could not open configuration file: %s",basename);
	return OC_EPERM;
    }

    for(;;) {
	char *line,*key,*value;
	int c;
        if(!rcreadline(in_file,line0,sizeof(line0))) break;
	linecount++;
	if(linecount >= MAXRCLINES) {
	    oclog(OCLOGERR, ".dodsrc has too many lines");
	    return 0;
	}	    	
	line = line0;
	/* check for comment */
	c = line[0];
        if (c == '#') continue;
	rctrim(line);  /* trim leading and trailing blanks */
	if(strlen(line) >= MAXRCLINESIZE) {
	    oclog(OCLOGERR, "%s line too long: %s",basename,line0);
	    return 0;
	}	    	
        /* setup */
	ocdodsrc->triples[ocdodsrc->ntriples].url[0] = '\0';
	ocdodsrc->triples[ocdodsrc->ntriples].key[0] = '\0';
	ocdodsrc->triples[ocdodsrc->ntriples].value[0] = '\0';
	if(line[0] == LTAG) {
	    char* url = ++line;
	    char* rtag = strchr(line,RTAG);
	    if(rtag == NULL) {
		oclog(OCLOGERR, "Malformed [url] in %s entry: %s",basename,line);
		continue;
	    }	    
	    line = rtag + 1;
	    *rtag = '\0';
	    /* save the url */
	    strncpy(ocdodsrc->triples[ocdodsrc->ntriples].url,url,MAXRCLINESIZE);
	    rctrim(ocdodsrc->triples[ocdodsrc->ntriples].url);
	}
	/* split off key and value */
	key=line;
	value = strchr(line, '=');
	if(value == NULL)
	    value = line + strlen(line);
	else {
	    *value = '\0';
	    value++;
	}
	strncpy(ocdodsrc->triples[ocdodsrc->ntriples].key,key,MAXRCLINESIZE);
	if(*value == '\0')
	    strcpy(ocdodsrc->triples[ocdodsrc->ntriples].value,"1");/*dfalt*/
	else
	    strncpy(ocdodsrc->triples[ocdodsrc->ntriples].value,value,MAXRCLINESIZE);
	rctrim(	ocdodsrc->triples[ocdodsrc->ntriples].key);
	rctrim(	ocdodsrc->triples[ocdodsrc->ntriples].value);
	ocdodsrc->ntriples++;
    }
    fclose(in_file);
    sorttriplestore();
    return 1;
}

int
ocdodsrc_process(OCstate* state)
{
    int stat = 0;
    char* value;
    char* url = ocuribuild(state->uri,NULL,NULL,OCURIENCODE);
    struct OCTriplestore* ocdodsrc = ocglobalstate.ocdodsrc;

    if(ocdodsrc == NULL) goto done;
    value = curllookup("DEFLATE",url);
    if(value != NULL) {
        if(atoi(value)) state->curlflags.compress = 1;
        if(ocdebug > 0)
            oclog(OCLOGNOTE,"Compression: %ld", state->curlflags.compress);
    }
    if((value = curllookup("VERBOSE",url)) != NULL) {
        if(atoi(value)) state->curlflags.verbose = 1;
        if(ocdebug > 0)
            oclog(OCLOGNOTE,"curl.verbose: %ld", state->curlflags.verbose);
    }
    if((value = curllookup("TIMEOUT",url)) != NULL) {
        if(atoi(value)) state->curlflags.timeout = atoi(value);
        if(ocdebug > 0)
            oclog(OCLOGNOTE,"curl.timeout: %ld", state->curlflags.timeout);
    }
    if((value = curllookup("USERAGENT",url)) != NULL) {
        if(atoi(value)) state->curlflags.useragent = strdup(value);
        if(state->curlflags.useragent == NULL) {stat = OC_ENOMEM; goto done;}
        if(ocdebug > 0)
            oclog(OCLOGNOTE,"USERAGENT: %s", state->curlflags.useragent);
    }

    if((value = curllookup("COOKIEJAR",url))
       || (value = curllookup("COOKIE_JAR",url))) {
        state->curlflags.cookiejar = strdup(value);
        if(state->curlflags.cookiejar == NULL) {stat = OC_ENOMEM; goto done;}
        if(ocdebug > 0)
            oclog(OCLOGNOTE,"COOKIEJAR: %s", state->curlflags.cookiejar);
    }

    if((value = curllookup("PROXY_SERVER",url)) != NULL) {
        stat = parseproxy(state,value);
        if(stat != OC_NOERR) goto done;
    }

    if((value = curllookup("SSL.VALIDATE",url)) != NULL) {
        if(atoi(value)) state->ssl.validate = 1;
        if(ocdebug > 0)
            oclog(OCLOGNOTE,"CURL.SSL.VALIDATE: %ld", state->ssl.validate);
    }

    if((value = curllookup("SSL.CERTIFICATE",url)) != NULL) {
        state->ssl.certificate = strdup(value);
        if(state->ssl.certificate == NULL) {stat = OC_ENOMEM; goto done;}
        if(ocdebug > 0)
            oclog(OCLOGNOTE,"CREDENTIALS.SSL.CERTIFICATE: %s", state->ssl.certificate);
    }

    if((value = curllookup("SSL.KEY",url)) != NULL) {
        state->ssl.key = strdup(value);
        if(state->ssl.key == NULL) {stat = OC_ENOMEM; goto done;}
        if(ocdebug > 0)
            oclog(OCLOGNOTE,"CREDENTIALS.SSL.KEY: %s", state->ssl.key);
    }

    if((value = curllookup("SSL.KEYPASSWORD",url)) != NULL) {
        state->ssl.keypasswd = strdup(value);
        if(state->ssl.keypasswd == NULL) {stat = OC_ENOMEM; goto done;}
#ifdef INSECURE
        if(ocdebug > 0)
            oclog(OCLOGNOTE,"CREDENTIALS.SSL.KEYPASSWORD: %s", state->ssl.keypasswd);
#endif
    }

    if((value = curllookup("SSL.CAINFO",url)) != NULL) {
        state->ssl.cainfo = strdup(value);
        if(state->ssl.cainfo == NULL) {stat = OC_ENOMEM; goto done;}
        if(ocdebug > 0)
            oclog(OCLOGNOTE,"SSL.CAINFO: %s", state->ssl.cainfo);
    }

    if((value = curllookup("SSL.CAPATH",url)) != NULL) {
        state->ssl.capath = strdup(value);
        if(state->ssl.capath == NULL) {stat = OC_ENOMEM; goto done;}
        if(ocdebug > 0)
            oclog(OCLOGNOTE,"SSL.CAPATH: %s", state->ssl.capath);
    }

    if((value = curllookup("SSL.VERIFYPEER",url)) != NULL) {
	char* s = strdup(value);
	int tf = 0;
	if(s == NULL || strcmp(s,"0")==0 || strcasecmp(s,"false")==0)
	    tf = 0;
	else if(strcmp(s,"1")==0 || strcasecmp(s,"true")==0)
	    tf = 1;
	else
	    tf = 1; /* default if not null */
        state->ssl.verifypeer = tf;
        if(ocdebug > 0)
            oclog(OCLOGNOTE,"SSL.VERIFYPEER: %d", state->ssl.verifypeer);
    }

    if((value = curllookup("CREDENTIALS.USER",url)) != NULL) {
        state->creds.username = strdup(value);
        if(state->creds.username == NULL) {stat = OC_ENOMEM; goto done;}
        if(ocdebug > 0)
            oclog(OCLOGNOTE,"CREDENTIALS.USER: %s", state->creds.username);
    }

    if((value = curllookup("CREDENTIALS.PASSWORD",url)) != NULL) {
        state->creds.password = strdup(value);
        if(state->creds.password == NULL) {stat = OC_ENOMEM; goto done;}
    }

    /* Support combined case */
    if((value = curllookup("CREDENTIALS.USERPASSWORD",url)) != NULL) {
	char* combined  = value;
		char* sep = NULL;
        if(combined == NULL) {stat = OC_ENOMEM; goto done;}
		sep = (char*)strchr(combined,':');
        if(sep == NULL) {
            oclog(OCLOGERR,"CREDENTIALS.USERPASSWORD: no ':' found");
	    stat = OC_EINVAL;
            goto done;
        }
	*sep = '\0';
        state->creds.username = strdup(combined);
        state->creds.password = strdup(sep+1);
    }

    /* else ignore */    

done:
    if(url != NULL) free(url);
    return stat;
}
    
char*
ocdodsrc_lookup(char* key, char* url)
{
    int i,found;
    struct OCTriplestore* ocdodsrc = ocglobalstate.ocdodsrc;
    struct OCTriple* triple = ocdodsrc->triples;

    if(key == NULL || ocdodsrc == NULL) return NULL;
    if(url == NULL) url = "";
    /* Assume that the triple store has been properly sorted */
    for(found=0,i=0;i<ocdodsrc->ntriples;i++,triple++) {
	size_t triplelen = strlen(triple->url);
	int t;
	if(strcmp(key,triple->key) != 0) continue; /* keys do not match */
	/* If the triple entry has no url, then use it (because we have checked all other cases)*/
	if(triplelen == 0) {found=1;break;}
	/* do url prefix comparison */
	t = ocstrncmp(url,triple->url,triplelen);
	if(t ==  0) {found=1; break;}
    }
    if(ocdebug > 2)
    {
	if(found) {
	    fprintf(stderr,"lookup %s: [%s]%s = %s\n",url,triple->url,triple->key,triple->value);
	}
    }    
    return (found ? triple->value : NULL);
}


static void
ocdodsrcdump(char* msg, struct OCTriple* triples, int ntriples)
{
    int i;
    struct OCTriplestore* ocdodsrc = ocglobalstate.ocdodsrc;

    if(msg != NULL) fprintf(stderr,"%s\n",msg);
    if(ocdodsrc == NULL) {
	fprintf(stderr,"<EMPTY>\n");
	return;
    }
    if(triples == NULL) triples= ocdodsrc->triples;
    if(ntriples < 0 ) ntriples= ocdodsrc->ntriples;
    for(i=0;i<ntriples;i++) {
        fprintf(stderr,"\t%s\t%s\t%s\n",
		(strlen(triples[i].url)==0?"--":triples[i].url),
		triples[i].key,
		triples[i].value);
    }
}

/* Isolate the "HTTP." (or "CURL.")
   prefix to allow changing to something else
*/
static char*
curllookup(char* suffix, char* url)
{
    char key[2048];
    char* value = NULL;
    if(!occopycat(key,sizeof(key),2,HTTPPREFIX,suffix))
	return NULL;
    value = ocdodsrc_lookup(key,url);
    if(value == NULL) {
        if(!occopycat(key,sizeof(key),2,HTTPPREFIXDEPRECATED,suffix))
	    return NULL;
        value = ocdodsrc_lookup(key,url);
    }
    return value;
}

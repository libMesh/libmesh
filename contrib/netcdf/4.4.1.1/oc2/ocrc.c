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

#define OCRCFILEENV "DAPRCFILE"

#define RTAG ']'
#define LTAG '['

#define TRIMCHARS " \t\r\n"

static OCerror rc_search(const char* prefix, const char* rcfile, char** pathp);

static int rcreadline(FILE* f, char* more, int morelen);
static void rctrim(char* text);
static char* combinecredentials(const char* user, const char* pwd);

static void storedump(char* msg, struct OCTriple*, int ntriples);

/* Define default rc files and aliases, also defines search order*/
static char* rcfilenames[] = {".daprc",".dodsrc",NULL};

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

static OCerror
ocextract_credentials(const char *url, char **userpwd, char **result_url)
{
    OCURI* parsed = NULL;
    if(!ocuriparse(url,&parsed))
	return OCTHROW(OC_EBADURL);
    if(parsed->userpwd == NULL) {
	ocurifree(parsed);
	return OCTHROW(OC_EBADURL);
    }
    if(userpwd) *userpwd = strdup(parsed->userpwd);
    ocurifree(parsed);
    return OC_NOERR;
}

char*
occombinehostport(const OCURI* uri)
{
    char* hp;
    int len = 0;

    if(uri->host == NULL)
	return NULL;
    else
	len += strlen(uri->host);
    if(uri->port != NULL)
	len += strlen(uri->port);
    hp = (char*)malloc(len+1);
    if(hp == NULL)
	return NULL;
    if(uri->port == NULL)
        occopycat(hp,len+1,1,uri->host);
    else
        occopycat(hp,len+1,3,uri->host,":",uri->port);
    return hp;
}

static char*
combinecredentials(const char* user, const char* pwd)
{
    int userPassSize;
    char *userPassword;

    if(user == NULL) user = "";
    if(pwd == NULL) pwd = "";

    userPassSize = strlen(user) + strlen(pwd) + 2;
    userPassword = malloc(sizeof(char) * userPassSize);
    if (!userPassword) {
        oclog(OCLOGERR,"Out of Memory\n");
	return NULL;
    }
    occopycat(userPassword,userPassSize-1,3,user,":",pwd);
    return userPassword;
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
    size_t len;
    int i;

    len = strlen(text);
    /* locate first non-trimchar */
    for(;*p;p++) {
       if(strchr(TRIMCHARS,*p) == NULL) break; /* hit non-trim char */
    }
    memmove(text,p,strlen(p)+1);
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

int
ocparseproxy(OCstate* state, char* v)
{
    /* Do not free these; they are pointers into v; free v instead */
    char *host_pos = NULL;
    char *port_pos = NULL;
    if(v == NULL || strlen(v) == 0)
	return OC_NOERR; /* nothing there*/
    if (occredentials_in_url(v)) {
        char *result_url = NULL;
        ocextract_credentials(v, &state->proxy.userpwd, &result_url);
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
        if (state->proxy.host == NULL)
            return OCTHROW(OC_ENOMEM);
        strncpy(state->proxy.host, host_pos, host_len);
        state->proxy.host[host_len] = '\0';
        state->proxy.port = atoi(port_pos);
    } else {
        size_t host_len = strlen(host_pos);
        state->proxy.host = malloc(sizeof(char) * host_len + 1);
        if (state->proxy.host == NULL)
            return OCTHROW(OC_ENOMEM);
        strncpy(state->proxy.host, host_pos, host_len);
        state->proxy.host[host_len] = '\0';
        state->proxy.port = 80;
    }
#if 0
    state->proxy.host[v_len] = '\0';
    state->proxy.port = atoi(v);
    s_len = strlen(v);
    state->proxy.user = malloc(sizeof(char) * s_len + 1);
    if (state->proxy.user == NULL)
        return OC_ENOMEM;
     strncpy(state->proxy.user, v, s_len);
     state->proxy.user[s_len] = '\0';
     p_len = strlen(v);
     state->proxy.password = malloc(sizeof(char) * p_len + 1);
     if (state->proxy.password == NULL)
         return OCTHROW(OC_ENOMEM);
     strncpy(state->proxy.password, v, p_len);
     state->proxy.password[p_len] = '\0';
#endif /*0*/
     if (ocdebug > 1) {
         oclog(OCLOGNOTE,"host name: %s", state->proxy.host);
#ifdef INSECURE
         oclog(OCLOGNOTE,"user+pwd: %s", state->proxy.userpwd);
#endif
         oclog(OCLOGNOTE,"port number: %d", state->proxy.port);
    }
    if(v) free(v);
    return OC_NOERR;
}

/* insertion sort the triplestore based on url */
static void
sorttriplestore(struct OCTriplestore* store)
{
    int i, nsorted;
    struct OCTriple* sorted = NULL;

    if(store == NULL) return; /* nothing to sort */
    if(store->ntriples <= 1) return; /* nothing to sort */
    if(ocdebug > 2)
        storedump("initial:",store->triples,store->ntriples);

    sorted = (struct OCTriple*)malloc(sizeof(struct OCTriple)*store->ntriples);
    if(sorted == NULL) {
        oclog(OCLOGERR,"sorttriplestore: out of memory");
        return;
    }

    nsorted = 0;
    while(nsorted < store->ntriples) {
        int largest;
        /* locate first non killed entry */
        for(largest=0;largest<store->ntriples;largest++) {
            if(store->triples[largest].key[0] != '\0') break;
        }
        OCASSERT(store->triples[largest].key[0] != '\0');
        for(i=0;i<store->ntriples;i++) {
            if(store->triples[i].key[0] != '\0') { /* avoid empty slots */
                int lexorder = strcmp(store->triples[i].host,store->triples[largest].host);
                int leni = strlen(store->triples[i].host);
                int lenlarge = strlen(store->triples[largest].host);
                /* this defines the ordering */
                if(leni == 0 && lenlarge == 0) continue; /* if no urls, then leave in order */
                if(leni != 0 && lenlarge == 0) largest = i;
                else if(lexorder > 0) largest = i;
            }
        }
        /* Move the largest entry */
        OCASSERT(store->triples[largest].key[0] != 0);
        sorted[nsorted] = store->triples[largest];
        store->triples[largest].key[0] = '\0'; /* kill entry */
        nsorted++;
      if(ocdebug > 2)
            storedump("pass:",sorted,nsorted);
    }

    memcpy((void*)store->triples,(void*)sorted,sizeof(struct OCTriple)*nsorted);
    free(sorted);

    if(ocdebug > 1)
        storedump("final .rc order:",store->triples,store->ntriples);
}

/* Create a triple store from a file */
static int
ocrc_compile(const char* path)
{
    char line0[MAXRCLINESIZE+1];
    FILE *in_file = NULL;
    int linecount = 0;
    struct OCTriplestore* ocrc = &ocglobalstate.rc.daprc;

    ocrc->ntriples = 0; /* reset; nothing to free */

    in_file = fopen(path, "r"); /* Open the file to read it */
    if (in_file == NULL) {
        oclog(OCLOGERR, "Could not open configuration file: %s",path);
        return OC_EPERM;
    }

    for(;;) {
        char *line,*key,*value;
        int c;
        if(!rcreadline(in_file,line0,sizeof(line0))) break;
        linecount++;
        if(linecount >= MAXRCLINES) {
            oclog(OCLOGERR, ".rc has too many lines");
            return 0;
        }
        line = line0;
        /* check for comment */
        c = line[0];
        if (c == '#') continue;
        rctrim(line);  /* trim leading and trailing blanks */
	if(strlen(line) == 0) continue;
        if(strlen(line) >= MAXRCLINESIZE) {
            oclog(OCLOGERR, "%s line too long: %s",path,line0);
            continue; /* ignore it */
        }
        /* setup */
        ocrc->triples[ocrc->ntriples].host[0] = '\0';
        ocrc->triples[ocrc->ntriples].key[0] = '\0';
        ocrc->triples[ocrc->ntriples].value[0] = '\0';
        if(line[0] == LTAG) {
	    OCURI* uri;
            char* url = ++line;
            char* rtag = strchr(line,RTAG);
            if(rtag == NULL) {
                oclog(OCLOGERR, "Malformed [url] in %s entry: %s",path,line);
                continue;
            }
            line = rtag + 1;
            *rtag = '\0';
            /* compile the url and pull out the host */
	    if(!ocuriparse(url,&uri)) {
                oclog(OCLOGERR, "Malformed [url] in %s entry: %s",path,line);
		continue;
	    }
            strncpy(ocrc->triples[ocrc->ntriples].host,uri->host,MAXRCLINESIZE-1);
	    if(uri->port != NULL) {
                strncat(ocrc->triples[ocrc->ntriples].host,":",MAXRCLINESIZE-1);
                strncat(ocrc->triples[ocrc->ntriples].host,uri->port,MAXRCLINESIZE-1);
	    }
	    ocurifree(uri);
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
        strncpy(ocrc->triples[ocrc->ntriples].key,key,MAXRCLINESIZE-1);
        if(*value == '\0')
            strcpy(ocrc->triples[ocrc->ntriples].value,"1");/*dfalt*/
        else
          strncpy(ocrc->triples[ocrc->ntriples].value,value,(MAXRCLINESIZE-1));
        rctrim( ocrc->triples[ocrc->ntriples].key);
        rctrim( ocrc->triples[ocrc->ntriples].value);
	OCDBG2("rc: key=%s value=%s",
		ocrc->triples[ocrc->ntriples].key,
		ocrc->triples[ocrc->ntriples].value);
        ocrc->ntriples++;
    }
    fclose(in_file);
    sorttriplestore(&ocglobalstate.rc.daprc);
    return 1;
}

/* read and compile the rc file, if any */
OCerror
ocrc_load(void)
{
    OCerror stat = OC_NOERR;
    char* path = NULL;

    if(ocglobalstate.rc.ignore) {
        oclog(OCLOGDBG,"No runtime configuration file specified; continuing");
	return OC_NOERR;
    }
    if(ocglobalstate.rc.loaded) return OC_NOERR;

    /* locate the configuration files in the following order:
       1. specified by set_rcfile
       2. set by DAPRCFILE env variable
       3. '.'
       4. $HOME
    */
    if(ocglobalstate.rc.rcfile != NULL) { /* always use this */
	path = strdup(ocglobalstate.rc.rcfile);
    } else if(getenv(OCRCFILEENV) != NULL && strlen(getenv(OCRCFILEENV)) > 0) {
        path = strdup(getenv(OCRCFILEENV));
    } else {
	char** rcname;
	int found = 0;
	for(rcname=rcfilenames;!found && *rcname;rcname++) {
	    stat = rc_search(".",*rcname,&path);
    	    if(stat == OC_NOERR && path == NULL)  /* try $HOME */
	        stat = rc_search(ocglobalstate.home,*rcname,&path);
	    if(stat != OC_NOERR)
		goto done;
	    if(path != NULL)
		found = 1;
	}
    }
    if(path == NULL) {
        oclog(OCLOGDBG,"Cannot find runtime configuration file; continuing");
    } else {
	if(ocdebug > 0)
	    fprintf(stderr, "RC file: %s\n", path);
        if(ocrc_compile(path) == 0) {
	    oclog(OCLOGERR, "Error parsing %s\n",path);
	    stat = OC_ERCFILE;
	}
    }
done:
    ocglobalstate.rc.loaded = 1; /* even if not exists */
    if(path != NULL)
	free(path);
    return stat;
}

OCerror
ocrc_process(OCstate* state)
{
    OCerror stat = OC_NOERR;
    char* value = NULL;
    OCURI* uri = state->uri;
    char* url_userpwd = NULL;
    char* url_hostport = NULL;

    if(!ocglobalstate.initialized)
	ocinternalinitialize();
    if(!ocglobalstate.rc.loaded)
	ocrc_load();
    /* Note, we still must do this function even if
       ocglobalstate.rc.ignore is set in order
       to getinfo e.g. user:pwd from url
    */

    url_userpwd = uri->userpwd;
    url_hostport = occombinehostport(uri);
    if(url_hostport == NULL)
	return OC_ENOMEM;

    value = ocrc_lookup("HTTP.DEFLATE",url_hostport);
    if(value != NULL) {
        if(atoi(value)) state->curlflags.compress = 1;
        if(ocdebug > 0)
            oclog(OCLOGNOTE,"HTTP.DEFLATE: %ld", state->curlflags.compress);
    }
    if((value = ocrc_lookup("HTTP.VERBOSE",url_hostport)) != NULL) {
        if(atoi(value)) state->curlflags.verbose = 1;
        if(ocdebug > 0)
            oclog(OCLOGNOTE,"HTTP.VERBOSE: %ld", state->curlflags.verbose);
    }
    if((value = ocrc_lookup("HTTP.TIMEOUT",url_hostport)) != NULL) {
        if(atoi(value)) state->curlflags.timeout = atoi(value);
        if(ocdebug > 0)
            oclog(OCLOGNOTE,"HTTP.TIMEOUT: %ld", state->curlflags.timeout);
    }
    if((value = ocrc_lookup("HTTP.USERAGENT",url_hostport)) != NULL) {
        if(atoi(value)) state->curlflags.useragent = strdup(value);
        if(state->curlflags.useragent == NULL) {stat = OC_ENOMEM; goto done;}
        if(ocdebug > 0)
            oclog(OCLOGNOTE,"HTTP.USERAGENT: %s", state->curlflags.useragent);
    }

    if(
          (value = ocrc_lookup("HTTP.COOKIEFILE",url_hostport))
       || (value = ocrc_lookup("HTTP.COOKIE_FILE",url_hostport))
       || (value = ocrc_lookup("HTTP.COOKIEJAR",url_hostport))
       || (value = ocrc_lookup("HTTP.COOKIE_JAR",url_hostport))
      ) {
        state->curlflags.cookiejar = strdup(value);
        if(state->curlflags.cookiejar == NULL) {stat = OC_ENOMEM; goto done;}
        if(ocdebug > 0)
            oclog(OCLOGNOTE,"HTTP.COOKIEJAR: %s", state->curlflags.cookiejar);
    }

    if((value = ocrc_lookup("HTTP.PROXY_SERVER",url_hostport)) != NULL) {
        stat = ocparseproxy(state,value);
        if(stat != OC_NOERR) goto done;
        if(ocdebug > 0)
            oclog(OCLOGNOTE,"HTTP.PROXY_SERVER: %s", value);
    }

    if((value = ocrc_lookup("HTTP.SSL.VALIDATE",url_hostport)) != NULL) {
        if(atoi(value)) {
	    state->ssl.verifypeer = 1;
	    state->ssl.verifyhost = 1;
            if(ocdebug > 0)
                oclog(OCLOGNOTE,"HTTP.SSL.VALIDATE: %ld", 1);
	}
    }

    if((value = ocrc_lookup("HTTP.SSL.CERTIFICATE",url_hostport)) != NULL) {
        state->ssl.certificate = strdup(value);
        if(state->ssl.certificate == NULL) {stat = OC_ENOMEM; goto done;}
        if(ocdebug > 0)
            oclog(OCLOGNOTE,"HTTP.SSL.CERTIFICATE: %s", state->ssl.certificate);
    }

    if((value = ocrc_lookup("HTTP.SSL.KEY",url_hostport)) != NULL) {
        state->ssl.key = strdup(value);
        if(state->ssl.key == NULL) {stat = OC_ENOMEM; goto done;}
        if(ocdebug > 0)
            oclog(OCLOGNOTE,"HTTP.SSL.KEY: %s", state->ssl.key);
    }

    if((value = ocrc_lookup("HTTP.SSL.KEYPASSWORD",url_hostport)) != NULL) {
        state->ssl.keypasswd = strdup(value);
        if(state->ssl.keypasswd == NULL) {stat = OC_ENOMEM; goto done;}
        if(ocdebug > 0)
            oclog(OCLOGNOTE,"HTTP.SSL.KEYPASSWORD: %s", state->ssl.keypasswd);
    }

    if((value = ocrc_lookup("HTTP.SSL.CAINFO",url_hostport)) != NULL) {
        state->ssl.cainfo = strdup(value);
        if(state->ssl.cainfo == NULL) {stat = OC_ENOMEM; goto done;}
        if(ocdebug > 0)
            oclog(OCLOGNOTE,"HTTP.SSL.CAINFO: %s", state->ssl.cainfo);
    }

    if((value = ocrc_lookup("HTTP.SSL.CAPATH",url_hostport)) != NULL) {
        state->ssl.capath = strdup(value);
        if(state->ssl.capath == NULL) {stat = OC_ENOMEM; goto done;}
        if(ocdebug > 0)
            oclog(OCLOGNOTE,"HTTP.SSL.CAPATH: %s", state->ssl.capath);
    }

    if((value = ocrc_lookup("HTTP.SSL.VERIFYPEER",url_hostport)) != NULL) {
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
            oclog(OCLOGNOTE,"HTTP.SSL.VERIFYPEER: %d", state->ssl.verifypeer);
	free(s);
    }

    if((value = ocrc_lookup("HTTP.NETRC",url_hostport)) != NULL) {
        if(state->curlflags.netrc != NULL)
	    free(state->curlflags.netrc);
        state->curlflags.netrc = strdup(value);
        if(state->curlflags.netrc == NULL) {stat = OC_ENOMEM; goto done;}
        if(ocdebug > 0)
            oclog(OCLOGNOTE,"HTTP.NETRC: %s", state->curlflags.netrc);
    }

    { /* Handle various cases for user + password */
	/* First, see if the user+pwd was in the original url */
	char* userpwd = NULL;
	char* user = NULL;
	char* pwd = NULL;
	if(url_userpwd != NULL)
	    userpwd = url_userpwd;
	else {
   	    user = ocrc_lookup("HTTP.CREDENTIALS.USER",url_hostport);
	    pwd = ocrc_lookup("HTTP.CREDENTIALS.PASSWORD",url_hostport);
	    userpwd = ocrc_lookup("HTTP.CREDENTIALS.USERPASSWORD",url_hostport);
	}
	if(userpwd == NULL && user != NULL && pwd != NULL) {
	    userpwd = combinecredentials(user,pwd);
	    state->creds.userpwd = userpwd;
	} else if(userpwd != NULL)
	    state->creds.userpwd = strdup(userpwd);
    }

done:
    if(url_hostport != NULL) free(url_hostport);
    return stat;
}

static struct OCTriple*
ocrc_locate(char* key, char* hostport)
{
    int i,found;
    struct OCTriplestore* ocrc = &ocglobalstate.rc.daprc;
    struct OCTriple* triple;

    if(ocglobalstate.rc.ignore)
	return NULL;
    if(!ocglobalstate.rc.loaded)
	ocrc_load();

    triple = ocrc->triples;

    if(key == NULL || ocrc == NULL) return NULL;
    if(hostport == NULL) hostport = "";
    /* Assume that the triple store has been properly sorted */
    for(found=0,i=0;i<ocrc->ntriples;i++,triple++) {
        size_t hplen = strlen(triple->host);
        int t;
        if(strcmp(key,triple->key) != 0) continue; /* keys do not match */
        /* If the triple entry has no url, then use it
           (because we have checked all other cases)*/
        if(hplen == 0) {found=1;break;}
        /* do hostport match */
        t = strcmp(hostport,triple->host);
        if(t ==  0) {found=1; break;}
    }
    return (found?triple:NULL);
}

char*
ocrc_lookup(char* key, char* hostport)
{
    struct OCTriple* triple = ocrc_locate(key,hostport);
    if(triple != NULL && ocdebug > 2) {
	fprintf(stderr,"lookup %s: [%s]%s = %s\n",hostport,triple->host,triple->key,triple->value);
    }
    return (triple == NULL ? NULL : triple->value);
}


static void
storedump(char* msg, struct OCTriple* triples, int ntriples)
{
    int i;
    struct OCTriplestore* ocrc = &ocglobalstate.rc.daprc;

    if(msg != NULL) fprintf(stderr,"%s\n",msg);
    if(ocrc == NULL) {
        fprintf(stderr,"<EMPTY>\n");
        return;
    }
    if(triples == NULL) triples= ocrc->triples;
    if(ntriples < 0 ) ntriples= ocrc->ntriples;
    for(i=0;i<ntriples;i++) {
        fprintf(stderr,"\t%s\t%s\t%s\n",
                (strlen(triples[i].host)==0?"--":triples[i].host),
                triples[i].key,
                triples[i].value);
    }
}

#if 0
/*
Lookup against all prefixes
*/

static char*
ocrc_lookup(char* suffix, char* url)
{
    char* value = NULL;
    char key[MAXRCLINESIZE+1];
    const char** p = prefixes;
    for(;*p;p++) {
        if(!occopycat(key,sizeof(key),2,*p,suffix))
            return NULL;
	value = ocrc_lookup(key,url);
	if(value != NULL)
	    return value;
    }
    return value;
}

/* compile the rc file, if any */
static OCerror
ocreadrc(void)
{
    OCerror stat = OC_NOERR;
    char* path = NULL;
    /* locate the configuration files: first if specified,
       then '.',  then $HOME */
    if(ocglobalstate.rc.rcfile != NULL) { /* always use this */
	path = ocglobalstate.rc.rcfile;
    } else {
	char** rcname;
	int found = 0;
	for(rcname=rcfilenames;!found && *rcname;rcname++) {
	    stat = rc_search(".",*rcname,&path);
    	    if(stat == OC_NOERR && path == NULL)  /* try $HOME */
	        stat = rc_search(ocglobalstate.home,*rcname,&path);
	    if(stat != OC_NOERR)
		goto done;
	    if(path != NULL)
		found = 1;
	}
    }
    if(path == NULL) {
        oclog(OCLOGDBG,"Cannot find runtime configuration file; continuing");
    } else {
	if(ocdebug > 0)
	    fprintf(stderr, "DODS RC file: %s\n", path);
        if(ocdodsrc_read(path) == 0) {
	    oclog(OCLOGERR, "Error parsing %s\n",path);
	    stat = OC_ERCFILE;
	}
    }
done:
    if(path != NULL)
	free(path);
    return stat;
}
#endif

/**
 * Prefix must end in '/'
 */
static
OCerror
rc_search(const char* prefix, const char* rcname, char** pathp)
{
    char* path = NULL;
    FILE* f = NULL;
    int plen = strlen(prefix);
    int rclen = strlen(rcname);
    OCerror stat = OC_NOERR;

    size_t pathlen = plen+rclen+1+1; /*+1 for '/' +1 for nul*/
    path = (char*)malloc(pathlen);
    if(path == NULL) {
	stat = OC_ENOMEM;
	goto done;
    }
    if(!occopycat(path,pathlen,3,prefix,"/",rcname)) {
        stat = OC_EOVERRUN;
	goto done;
    }
    /* see if file is readable */
    f = fopen(path,"r");
    if(f != NULL)
        oclog(OCLOGDBG, "Found rc file=%s",path);
done:
    if(f == NULL || stat != OC_NOERR) {
      if(path != NULL)
	    free(path);
      path = NULL;
    }

    if(f != NULL)
      fclose(f);
    if(pathp != NULL)
      *pathp = path;
    else {
      free(path);
      path = NULL;
    }

    return OCTHROW(stat);
}

struct OCTriple*
ocrc_triple_iterate(char* key, char* url, struct OCTriple* prev)
{
    struct OCTriple* next;
    if(prev == NULL)
      next = ocrc_locate(key,url);
    else
      next = prev+1;
    if(next == NULL)
      return NULL;
    for(; strlen(next->key) > 0; next++) {
      /* See if key as prefix still matches */
      int cmp = strcmp(key,next->key);
      if(cmp != 0) {next = NULL; break;} /* key mismatch */
      /* compare url */
      cmp = ocstrncmp(url,next->host,strlen(next->host));
      if(cmp ==  0) break;
    }
    return next;
}

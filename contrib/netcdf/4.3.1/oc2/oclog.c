/*********************************************************************
 *   Copyright 2010, UCAR/Unidata
 *   See netcdf/COPYRIGHT file for copying and redistribution conditions.
 *   $Header$
 *********************************************************************/

#include "config.h"
#include <stdlib.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#include <stdio.h>
#include <string.h>
#ifdef HAVE_STDARG_H
#include <stdarg.h>
#endif
#ifdef HAVE_FCNTL_H
#include <fcntl.h>
#endif

#include "oclog.h"

#define PREFIXLEN 8
#define MAXTAGS 256
#define OCTAGDFALT "Log";

static int oclogginginitialized = 0;
static int oclogging = 0;
static int ocsystemfile = 0; /* 1 => we are logging to file we did not open */
static char* oclogfile = NULL;
static FILE* oclogstream = NULL;

static int octagsize = 0;
static char** octagset = NULL;
static char* octagdfalt = NULL;
static char* octagsetdfalt[] = {"Warning","Error","Note","Debug"};
static char* octagname(int tag);

/*!\defgroup OClog OClog Management
@{*/

/*!\internal
*/

void
ocloginit(void)
{
    const char* file;
    if(oclogginginitialized)
	return;
    oclogginginitialized = 1;
    file = getenv(OCENVFLAG);
    ocsetlogging(0);
    oclogfile = NULL;
    oclogstream = NULL;
    /* Use environment variables to preset oclogging state*/
    /* I hope this is portable*/
    if(file != NULL && strlen(file) > 0) {
        if(oclogopen(file)) {
	    ocsetlogging(1);
	}
    }
    octagdfalt = OCTAGDFALT;
    octagset = octagsetdfalt;
}

/*!
Enable/Disable logging.

\param[in] tf If 1, then turn on logging, if 0, then turn off logging.

\return The previous value of the logging flag.
*/

int
ocsetlogging(int tf)
{
    int was;
    if(!oclogginginitialized) ocloginit();
    was = oclogging;
    oclogging = tf;
    return was;
}

/*!
Specify a file into which to place logging output.

\param[in] file The name of the file into which to place logging output.
If the file has the value NULL, then send logging output to
stderr.

\return zero if the open failed, one otherwise.
*/

int
oclogopen(const char* file)
{
    if(!oclogginginitialized) ocloginit();
    oclogclose();
    if(file == NULL || strlen(file) == 0) {
	/* use stderr*/
	oclogstream = stderr;
	oclogfile = NULL;
	ocsystemfile = 1;
    } else if(strcmp(file,"stdout") == 0) {
	/* use stdout*/
	oclogstream = stdout;
	oclogfile = NULL;
	ocsystemfile = 1;
    } else if(strcmp(file,"stderr") == 0) {
	/* use stderr*/
	oclogstream = stderr;
	oclogfile = NULL;
	ocsystemfile = 1;
    } else {
	int fd;
	oclogfile = strdup(file);
	oclogstream = NULL;
	/* We need to deal with this file carefully
	   to avoid unauthorized access*/
	fd = open(oclogfile,O_WRONLY|O_APPEND|O_CREAT,0600);
	if(fd >= 0) {
	    oclogstream = fdopen(fd,"a");
	} else {
	    free(oclogfile);
	    oclogfile = NULL;
	    oclogstream = NULL;
	    ocsetlogging(0);
	    return 0;
	}
	ocsystemfile = 0;
    }
    return 1;
}

void
oclogclose(void)
{
    if(!oclogginginitialized) ocloginit();
    if(oclogstream != NULL && !ocsystemfile) {
	fclose(oclogstream);
    }
    if(oclogfile != NULL) free(oclogfile);
    oclogstream = NULL;
    oclogfile = NULL;
    ocsystemfile = 0;
}

/*!
Send logging messages. This uses a variable
number of arguments and operates like the stdio
printf function.

\param[in] tag Indicate the kind of this log message.
\param[in] fmt Format specification as with printf.
*/

void
oclog(int tag, const char* fmt, ...)
{
    va_list args;
    char* prefix;

    if(!oclogginginitialized) ocloginit();

    if(!oclogging || oclogstream == NULL) return;

    prefix = octagname(tag);
    fprintf(oclogstream,"%s:",prefix);

    if(fmt != NULL) {
      va_start(args, fmt);
      vfprintf(oclogstream, fmt, args);
      va_end( args );
    }
    fprintf(oclogstream, "\n" );
    fflush(oclogstream);
}

void
oclogtext(int tag, const char* text)
{
    oclogtextn(tag,text,strlen(text));
}

/*!
Send arbitrarily long text as a logging message.
Each line will be sent using oclog with the specified tag.

\param[in] tag Indicate the kind of this log message.
\param[in] text Arbitrary text to send as a logging message.
\param[in] count Maximum ength of the text to write.
*/

void
oclogtextn(int tag, const char* text, size_t count)
{
    if(!oclogging || oclogstream == NULL) return;
    fwrite(text,1,count,oclogstream);
    fflush(oclogstream);
}

/* The tagset is null terminated */
void
oclogsettags(char** tagset, char* dfalt)
{
    octagdfalt = dfalt;
    if(tagset == NULL) {
	octagsize = 0;
    } else {
        int i;
	/* Find end of the tagset */
	for(i=0;i<MAXTAGS;i++) {if(tagset[i]==NULL) break;}
	octagsize = i;
    }
    octagset = tagset;
}

static char*
octagname(int tag)
{
    if(tag < 0 || tag >= octagsize) {
	return octagdfalt;
    } else {
	return octagset[tag];
    }
}

/**@}*/

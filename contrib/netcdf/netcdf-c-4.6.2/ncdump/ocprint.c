/* Copyright 2009, UCAR/Unidata and OPeNDAP, Inc.
   See the COPYRIGHT file for more information. */

#define VALIDATE

#define ALLATONCE
#undef TRACK

#include "config.h"

#include <stdlib.h>
#include <stdio.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#include <assert.h>
#include <string.h>
#ifdef HAVE_STRINGS_H
#include <strings.h>
#endif

#include "oc.h"
#include "ocx.h"

/* Utilities */
#include "netcdf.h"
#include "ncuri.h"
#include "ncbytes.h"
#include "nclog.h"

#ifdef WIN32
/*#include <windows.h>*/
#define snprintf _snprintf
#define strcasecmp stricmp
#endif

#ifdef _MSC_VER
#include "XGetopt.h"
int opterr;
int optind;
#endif

#ifndef nulldup
#define nulldup(s) (s==NULL?NULL:strdup(s))
#endif

#define CHECK(x) check_err((ocstat=(x)),0)
#define FAIL(x) check_err((x),1)

/* Define some classifiers */
#define iscontainer(t) ((t) == OC_Dataset || (t) == OC_Structure || (t) == OC_Sequence || (t) == OC_Grid)

#define isatomic(t) ((t) == OC_Atomic)

#define NORC "NONE"

#define LBRACE "{"
#define RBRACE "}"

/*Mnemonic*/
#define TOPLEVEL 1

int ocdebug;

static OCerror ocstat;
static OClink glink;

/* define a large stack of DUMPPATH datanodes */
struct DUMPPATH {
    OCdatanode datanode;
    OCddsnode   node;
    OCtype      octype;
    size_t      rank;
    size_t      dimsizes[OC_MAX_DIMENSIONS];
    int         indexed; /* 1 => indices is valid */
    size_t      indices[OC_MAX_DIMENSIONS];
} stack[2048];

static size_t stacknext;

/* Define the debug options */
struct OCD {
    int debug;          /* -D<integer:1..> */
    int debuglevel;
    int dumpdds;        /* -DN */
    int dumpdatadds;    /* -DX<level 0..> */
    int dumpdatatree;   /* -DD */
    int dumplevel;
    int curl;           /* -DC - make curl be verbose */
    int verbose;        /* -DV - produce more verbose output */
} debug;

/* Define the -X options; currently unused*/
struct OCX {
    int ignore;
} x;

/* Define the other options */
static struct OCOPT {
    char*   surl; /* full url string */
    NCURI*  url;
    struct OCD debug;
    struct OCX x;
    int     showattributes; /* -A */
    int     logging;	    /* -L */
    char*   netrc ;	    /* -N */
    char*   rcfile ;	    /* -R */
    int     selfsigned ;    /* -S */
    int     octest; 	    /* -T */ /* match original octest output */
    int     generate;	    /* -g|-G */
    int     optdas;	    /* -p */
    int     optdatadds;	    /* -p */
    int     optdds;	    /* -p */
    FILE*   output;         /* -o */
    /* Deprecated */
    char*   constraint;	    /* -C */
    NCbytes* userparams;    /* -U */
} ocopt;

static char blanks[2048];
#define BLANKSPERDENT 2

/* Forward*/
static void usage(char*);
static int fail(char*);
static void check_err(OCerror stat, int dofail);
static void dumpflags(void);

struct OCOPT;
static OCerror processdata(int);
static OCerror printdata(OClink, OCdatanode);

static OCerror printdata_indices(OClink, OCdatanode, NCbytes*,int);
static OCerror printdata_container(OClink, OCdatanode, NCbytes*,int);
static OCerror printdata_leaf(OClink, OCdatanode, NCbytes*,int);

static off_t odom_init(size_t rank, size_t* indices, size_t* dimsizes);
static int odom_more(size_t rank, size_t* indices, size_t* dimsizes);
static void odom_next(size_t rank, size_t* indices, size_t* dimsizes);

static OCerror dumpdatanode(OClink, OCdatanode, size_t count, void* memory, NCbytes*);
static OCerror generatedds(OClink, OCddsnode, NCbytes*, int depth);
static char* generatedas(OClink,OCddsnode);
static OCerror generatedasr(OClink, OCddsnode, NCbytes*, int depth);
static OCerror generateddsattributes(OClink, OCddsnode node, NCbytes*, int);

static OCerror printdims(OClink link, OCddsnode node, NCbytes* buffer);
static char* stringescape(char*);
static char* idescape(char*, char*, size_t);
static int needsescapes(const char* s);
static size_t totaldimsize(size_t,size_t*);
static char* indent(int n);
static void pushstack(OCdatanode datanode);
static void popstack() {stacknext--;}
#ifdef TRACK
static void printstack(char* msg);
#endif

static char* optionmsg =
" [-A]"
" [-D <debugarg>]"
" [-G]"
" [-L]"
" [-N <netrc file>]"
" [-S]"
" [-R <rcfile>]"
" [-T]"
" [-h]"
" [-o <file|->]"
" [-p das|dds|datadds]"
" <url>";

static OCflags ocflags;

EXTERNL int nc_initialize(void);

static void
init()
{
    memset(&ocopt,0,sizeof(ocopt));
    ocopt.generate = 1;             /* -G|-g */
    ocopt.userparams = ncbytesnew(); /* -U */
    nc_initialize();
}

int
main(int argc, char **argv)
{
    int c;
    char* suffix;

#ifdef OCDEBUG
    { int i;
	fprintf(stderr,"argv =");
	for(i=0;i<argc;i++)
	    fprintf(stderr," %s",argv[i]);
	fprintf(stderr,"\n");
    }
#endif

    init();

    opterr = 1;

    while ((c = getopt(argc, argv, "AC:D:GLN:R:STU:X:gho:u:f:p:")) != EOF) {
        switch (c) {
	case 'A': ocopt.showattributes = 1; break;
	case 'C': ocopt.constraint = nulldup(optarg); break;
        case 'G': case 'g': ocopt.generate = 1; break;
        case 'L': ocopt.logging = 1; break;
	case 'N': ocopt.netrc = nulldup(optarg); break;
	case 'R': ocopt.rcfile = nulldup(optarg); break;
        case 'S': ocopt.selfsigned = 1; break;
        case 'T': ocopt.octest = 1; break;
	case 'U':
	    if(optarg != NULL) {
		ncbytesappend(ocopt.userparams,';');
		ncbytescat(ocopt.userparams,optarg);
	    }
	    break;
        case 'D': {
	    int c0;
	    if(optarg == NULL || strlen(optarg) == 0) usage("missing -D argument");
	    c0 = optarg[0];
	    if(c0 >= '0' && c0 <= '9') {/* debug level */
		ocopt.debug.debuglevel = atoi(optarg); break;
	    } else switch (c0) {
	           case 'C': ocopt.debug.curl = 1; break;
	           case 'D': ocopt.debug.dumpdatatree = 1; break;
	           case 'N': ocopt.debug.dumpdds = 1; break;
	           case 'X': ocopt.debug.dumpdatadds = 1;
			     ocopt.debug.dumplevel = atoi(optarg+1);
			     break;
	           case 'V': ocopt.debug.verbose = 1; break;
		   default: usage("unknown -D option");
		   }
	    } break;
        case 'X': {
	    int c0;
	    int so = (optarg == NULL ? 0 : strlen(optarg));
	    if(so == 0) usage("missing -X argument");
	    c0 = optarg[0];
	    switch (c0) {
	    default:
		usage("unknown -X option");
	    }
	} break;

	case 'o':
            if(ocopt.output != NULL) fclose(ocopt.output);
	    if(optarg == NULL)
		usage("-o does not specify a file name");
	    ocopt.output = fopen(optarg,"w");
            if(ocopt.output == NULL)
		usage("-o file not writeable");
	    break;

	case 'u': case 'f':
	    ocopt.surl = optarg;
	    break;

	case 'p':
	    if(optarg == NULL)
		usage("-p does not specify an argument");
	    if(strcasecmp(optarg,"das")==0) ocopt.optdas=1;
	    else if(strcasecmp(optarg,"dds")==0) ocopt.optdds=1;
	    else if(strcasecmp(optarg,"data")==0) ocopt.optdatadds=1;
	    else if(strcasecmp(optarg,"datadds")==0) ocopt.optdatadds=1;
	    else if(strcasecmp(optarg,"all")==0) {
		ocopt.optdas = 1; ocopt.optdds = 1; ocopt.optdatadds = 1;
	    } else usage("unknown -p option");
	    break;

	case 'h': usage(""); exit(0);

        default: usage("unknown option");
        }
    }

    if(ocopt.output == NULL)
	ocopt.output = stdout;

    if (ocopt.debug.debuglevel > 0) {
        ocdebug = ocopt.debug.debuglevel;
    }

    if(ocopt.logging) {
	ncloginit();
	ncsetlogging(1);
	if(!nclogopen(NULL))
	    fprintf(stderr,"Failed to open logging output\n");
    }

    argc -= optind;
    argv += optind;

    if (argc > 0 && ocopt.surl== NULL) {
        ocopt.surl = nulldup(argv[argc - 1]);
    } else
        usage("Multiple urls specified");

    if (ocopt.surl == NULL)
        ocopt.surl = getenv("URLSRC");

    if (ocopt.surl == NULL) {
        usage("no source url specified\n");
    }

    /* Compile the url */
    if(ncuriparse(ocopt.surl,&ocopt.url) != NCU_OK) {
	fprintf(stderr,"malformed source url: %s\n",ocopt.surl);
	exit(1);
    }

    /* For convenience, see if the url has a trailing .dds, .das, or .dods
       and if so, use it
    */
    suffix = strrchr(ocopt.url->path,'.');
    if(suffix != NULL) {
	int match = 0;
	if(strcmp(suffix,".das")==0) {
	    ocopt.optdas = 1;
	    ocopt.optdds = 0;
	    ocopt.optdatadds = 0;
	    match = 1;
	} else if(strcmp(suffix,".dds")==0) {
	    ocopt.optdas = 0;
	    ocopt.optdds = 1;
	    ocopt.optdatadds = 0;
	    match = 1;
	} else if(strcmp(suffix,".dods")==0) {
	    ocopt.optdas = 0;
	    ocopt.optdds = 0;
	    ocopt.optdatadds = 1;
	    match = 1;
	}
	/* Remove the suffix */
	if(match)
	    *suffix = '\0';
    }

    /* If -C was specified, then it has precedence */
    if(ocopt.constraint != NULL) {
	ncurisetquery(ocopt.url,ocopt.constraint);
	nullfree(ocopt.constraint);
	ocopt.constraint = NULL;
    }
    /* Rebuild the url string */
    if(ocopt.surl != NULL) free(ocopt.surl);
    ocopt.surl = ncuribuild(ocopt.url,NULL,NULL,NCURIALL);

    /* Reparse */
    if(ncuriparse(ocopt.surl,&ocopt.url) != NCU_OK) {
	fprintf(stderr,"malformed source url: %s\n",ocopt.surl);
	exit(1);
    }

    if(ocopt.rcfile != NULL) {
    }

    if (ocopt.debug.verbose)
        dumpflags();

    processdata(ocflags);

    return 0;
}

static void
dumpflags(void)
{
    char* tmp;
    if(ocopt.showattributes) fprintf(stderr," -A");
    if(ocopt.debug.debug) fprintf(stderr," -D%d",ocopt.debug.debuglevel);
    if(ocopt.debug.dumpdds) fprintf(stderr," -DN");
    if(ocopt.debug.dumpdatatree) fprintf(stderr," -DD");
    if(ocopt.debug.dumpdatadds) fprintf(stderr," -DX%d",ocopt.debug.dumplevel);
    if(ocopt.debug.verbose) fprintf(stderr," -DV");
    if(ocopt.generate) fprintf(stderr," -G");
    if(ocopt.logging) fprintf(stderr," -L");
    if(ocopt.logging) fprintf(stderr," -N %s",ocopt.netrc);
    if(ocopt.logging) fprintf(stderr," -R %s",ocopt.rcfile);
    if(ocopt.optdas || ocopt.optdds || ocopt.optdatadds) {
	fprintf(stderr," -p");
	if(ocopt.optdas) fprintf(stderr," das");
	if(ocopt.optdds) fprintf(stderr," dds");
	if(ocopt.optdatadds) fprintf(stderr," datadds");
    }
    tmp = ncuribuild(ocopt.url,NULL,NULL,NCURIALL);
    fprintf(stderr,"%s\n",tmp);
    free(tmp);
}

static void
usage(char* msg)
{
    if(msg) fprintf(stderr,"error: %s\n",msg);
    fprintf(stderr,"usage: ocprint %s\n",optionmsg);
    fail(NULL);
}

static int
fail(char* msg)
{
    if(msg) fprintf(stderr,"fatalerror: %s\n",msg);
    fflush(ocopt.output); fflush(stderr);
    exit(1);
}

/******************************************/
static void
check_err(OCerror ocstat, int dofail)
{
    if(ocstat == OC_NOERR) return;
    fprintf(stderr,"error status returned: (%d) %s\n",ocstat,oc_errstring(ocstat));
    if(dofail) fail(NULL);
}

static OCerror
processdata(OCflags flags)
{
    char* totalurl;
    OClink link;
    OCddsnode dasroot, ddsroot, dataddsroot;
    OCdatanode rootdatanode;

    totalurl = ncuribuild(ocopt.url,NULL,NULL,NCURIALL);
    FAIL(oc_open(totalurl,&link));
    free(totalurl);
    glink = link;

    if(ocopt.debug.curl)
	oc_trace_curl(link);

    if(ocopt.netrc)
	oc_set_netrc(link,ocopt.netrc);

#if 0
    if(ocopt.selfsigned)
	oc_set_curlopt(link,"CURLOPT_VERIFYPEER", (void*)0L);
#endif

    if(ocopt.optdas) {
        ocstat = oc_fetch(link,ocopt.url->query,OCDAS,0,&dasroot);
        if(ocstat != OC_NOERR) {
            fprintf(stderr,"error status returned: (%d) %s\n",ocstat,oc_errstring(ocstat));
            fprintf(stderr,"Could not read DAS; continuing.\n");
            ocopt.optdas = 0;
            ocopt.showattributes = 0;
        } else if(ocopt.generate) {
            char* das = generatedas(link,dasroot);
            fprintf(ocopt.output,"%s",das);
            free(das);
        } else {
	    const char* text = oc_tree_text(link,dasroot);
            fprintf(ocopt.output,"%s",(text?text:"null"));
        }
    }
    fflush(ocopt.output);

    if(ocopt.optdds) {
        ocstat = oc_fetch(link,ocopt.url->query,OCDDS,flags,&ddsroot);
        if(ocstat != OC_NOERR) {
            fprintf(stderr,"error status returned: (%d) %s\n",ocstat,oc_errstring(ocstat));
            fprintf(stderr,"Could not read DDS; continuing.\n");
            ocopt.optdds = 0;
        } else {
            if(ocopt.showattributes && !ocopt.optdas) {
                FAIL(oc_fetch(link,ocopt.url->query,OCDAS,flags,&dasroot));
            }
            if(ocopt.showattributes || ocopt.optdas) {
                FAIL(oc_merge_das(link,dasroot,ddsroot));
            }
            if(ocopt.generate) {
	        NCbytes* buffer = ncbytesnew();
                FAIL(generatedds(link,ddsroot,buffer,0));
                fprintf(ocopt.output,"%s",ncbytescontents(buffer));
		ncbytesfree(buffer);
            } else {
                const char* text = oc_tree_text(link,ddsroot);
                fprintf(ocopt.output,"%s",(text?text:"null"));
            }
        }
        if(ocopt.debug.dumpdds)
            oc_dds_ddnode(link,ddsroot);
    }
    fflush(ocopt.output);

    if(ocopt.optdatadds) {
        ocstat = oc_fetch(link,ocopt.url->query,OCDATADDS,flags,&dataddsroot);
        if(ocstat) {
            fprintf(stderr,"Cannot read DATADDS: %s\n",ocopt.surl);
            exit(1);
        }
        if(ocopt.debug.dumpdds)
            oc_dds_ddnode(link,dataddsroot);
        if(ocopt.debug.dumpdatadds)
            oc_dds_dd(link,dataddsroot,ocopt.debug.dumplevel);

        FAIL(oc_dds_getdataroot(link,dataddsroot,&rootdatanode));
        if(ocopt.debug.dumpdatatree)
	    oc_data_ddtree(link,rootdatanode);
        stacknext = 0;
        printdata(link,rootdatanode);
    }
    fflush(ocopt.output);

    oc_close(link);
    return OC_NOERR;
}

/**
This is the main procedure for
printing a data tree.
*/
static OCerror
printdata(OClink link, OCdatanode datanode)
{
    NCbytes* buffer;
    OCtype octype;

    buffer = ncbytesnew();

    /* verify that datanode is a Dataset datanode */
    FAIL(oc_data_octype(link,datanode,&octype));
    assert(octype == OC_Dataset);

    printdata_container(link,datanode,buffer,TOPLEVEL);

    fprintf(ocopt.output,"%s",ncbytescontents(buffer));

    ncbytesfree(buffer);
    return OC_NOERR;
}


/* Print a single container datanode */
static OCerror
printdata_container(OClink link, OCdatanode datanode, NCbytes* buffer, int istoplevel)
{
    OCerror stat = OC_NOERR;
    size_t i;
    OCddsnode node;
    OCtype octype;
    size_t nsubnodes;

    /* Obtain some information about the node */
    FAIL(oc_data_ddsnode(link,datanode,&node));
    FAIL(oc_dds_nsubnodes(link,node,&nsubnodes));
    FAIL(oc_data_octype(link,datanode,&octype));

    /* If this is not a single instance container, then
       defer to the appropriate function */
    if(isatomic(octype))
	return printdata_leaf(link,datanode,buffer,istoplevel);
    if(oc_data_indexable(link,datanode))
	return printdata_indices(link,datanode,buffer,!TOPLEVEL);

    /* Must be a container instance or a record */
    for(i=0;i<nsubnodes;i++) {
        OCdatanode field;
	FAIL(oc_data_ithfield(link,datanode,i,&field));
	pushstack(field);
        FAIL(printdata_indices(link,field,buffer,istoplevel));
	popstack();
	if(oc_data_free(link,field) != OC_NOERR)
	    break;
    }
    return stat;
}

static OCerror
printdata_indices(OClink link, OCdatanode datanode, NCbytes* buffer, int istoplevel)
{
    OCerror stat = OC_NOERR;
    size_t i;
    OCddsnode node;
    size_t rank;
    OCtype octype;
    size_t dimsizes[OC_MAX_DIMENSIONS];
    size_t indices[OC_MAX_DIMENSIONS];

    /* Obtain some information about the node */
    FAIL(oc_data_ddsnode(link,datanode,&node));
    FAIL(oc_dds_octype(link,node,&octype));
    FAIL(oc_dds_rank(link,node,&rank));

    /* If this is not an indexed structure or a sequence then
       defer to the appropriate function */
    if(isatomic(octype))
	return printdata_leaf(link,datanode,buffer,istoplevel);
    if(iscontainer(octype) && !oc_data_indexable(link,datanode))
	return printdata_container(link,datanode,buffer,istoplevel);

    /* Iterate over the datanodes */
    if(octype == OC_Structure) {
	/* Get dimension sizes */
        FAIL(oc_dds_dimensionsizes(link,node,dimsizes));

	/* init odometer and get cross-product */
	odom_init(rank,indices,dimsizes);
        while(odom_more(rank,indices,dimsizes)) {
	    OCdatanode element;
	    FAIL(oc_data_ithelement(link,datanode,indices,&element));
            pushstack(element);
	    /* walk the container */
            printdata_container(link,element,buffer,!TOPLEVEL);
            popstack();
	    oc_data_free(link,element);
	    odom_next(rank,indices,dimsizes);
        }
    } else if(octype == OC_Sequence) {
        /* Dump each element */
        for(i=0;;i++) {
	    OCdatanode record;
	    stat = oc_data_ithrecord(link,datanode,i,&record);
	    if(stat != OC_NOERR) {
	        if(stat == OC_EINDEX) break;
	        return stat;
	    }
            pushstack(record);
            printdata_container(link,record,buffer,!TOPLEVEL); /* print current record */
            popstack();
	    oc_data_free(link,record);
        }
#ifdef VALIDATE
	{
	    size_t count;
	    FAIL(oc_data_recordcount(link,datanode,&count));
	    assert(count == i);
        }
#endif
    } else
        abort();

    return OC_NOERR;
}

static OCerror
printdata_leaf(OClink link, OCdatanode datanode, NCbytes* buffer, int istoplevel)
{
    OCddsnode node;
    OCtype octype,atomtype;
    size_t elemsize;
    size_t memsize;
    char* memory;
    size_t count,rank;

    /* Obtain some information about the node */
    FAIL(oc_data_ddsnode(link,datanode,&node));
    FAIL(oc_dds_octype(link,node,&octype));
    FAIL(oc_dds_atomictype(link,node,&atomtype));
    FAIL(oc_dds_rank(link,node,&rank));

    assert(octype == OC_Atomic);

    /* If this variable is top-level then
       use the oc_dds_read functions
       in order to test them
    */

    elemsize = oc_typesize(atomtype);

    if(rank == 0) {/* Scalar case */
	memory = calloc(elemsize,1); /* reading only one value */
        /* read the scalar */
	if(istoplevel) {
	    FAIL(oc_dds_read(link,node,NULL,NULL,elemsize,memory));
	} else {
	    FAIL(oc_data_read(link,datanode,NULL,NULL,elemsize,memory));
	}
        count = 1;
    } else {
	size_t dimsizes[OC_MAX_DIMENSIONS];
	size_t indices[OC_MAX_DIMENSIONS];
        FAIL(oc_dds_dimensionsizes(link,node,dimsizes));
	/* init odometer and get cross-product */
	count = odom_init(rank,indices,dimsizes);
        memsize = elemsize*count;
        memory = calloc(memsize,1);

#ifdef ALLATONCE /* read all at once */
        /* indices should be all zeros at this point */
	if(istoplevel) {
	    FAIL(oc_dds_read(link,node,indices,dimsizes,memsize,memory));
	} else {
	    FAIL(oc_data_read(link,datanode,indices,dimsizes,memsize,memory));
	}
#else /* BYITEM */
        {
  	    size_t offset;
	    size_t one[OC_MAX_DIMENSIONS];
            /* Initialize the read-by-one counts */
	    for(i=0;i<rank;i++) one[i]=0;
	    one[rank-1] = 1;
            /* Read whole atomic array item by item using an odometer */
	    for(offset=0,i=0;i<count;i++,offset+=elemsize) {
		if(!odom_more(rank,indices,dimsizes))
		    abort();
		if(istoplevel) {
		    FAIL(oc_dds_read(link,node,
                                      indices,one,
				      elemsize,memory+offset));
		} else {
	            FAIL(oc_data_read(link,datanode,
				       indices,one,
                                       elemsize,memory+offset));
		}
		odom_next(rank,indices,dimsizes);
	    }
        }
#endif
    }
    dumpdatanode(link,datanode,count,memory,buffer);
    if(atomtype == OC_String || atomtype == OC_URL)
	oc_reclaim_strings(count,(char**)memory);
    free(memory);
    return OC_NOERR;
}

static OCerror
generatedds(OClink link, OCddsnode node, NCbytes* buffer, int depth)
{
    size_t i,rank,nattr,nsubnodes;
    OCtype octype, atomtype;
    OCddsnode container,field;
    char id1[1024];
    char* name;

    ncbytescat(buffer,indent(depth));

    /* get all info about the node */
    FAIL(oc_dds_properties(link,node,&name,&octype,&atomtype,&container,
                             &rank,&nsubnodes,&nattr));

    if(octype == OC_Atomic) {
        ncbytescat(buffer,oc_typetostring(atomtype));
	ncbytescat(buffer," ");
	ncbytescat(buffer,idescape(name,id1,sizeof(id1)));
        /* dump dim info (if any)*/
	printdims(link,node,buffer);
        ncbytescat(buffer,";\n");
        generateddsattributes(link,node,buffer,depth+1);
    } else { /*must be container*/
	const char* typename = oc_typetostring(octype);
        ncbytescat(buffer,typename);
        ncbytescat(buffer," ");
        ncbytescat(buffer,LBRACE);
        ncbytescat(buffer,"\n");
        for(i=0;i<nsubnodes;i++) {
	    FAIL(oc_dds_ithfield(link,node,i,&field));
	    if(octype == OC_Grid) {
                ncbytescat(buffer,indent(depth));
		switch (i) {
		case 0: ncbytescat(buffer,"Array:\n"); break;
		case 1: ncbytescat(buffer,"Maps:\n"); break;
		default: break;
		}
 	    }
	    generatedds(link,field,buffer,depth+1);
        }
        ncbytescat(buffer,indent(depth));
        ncbytescat(buffer,RBRACE);
        ncbytescat(buffer," ");
        ncbytescat(buffer,idescape(name,id1,sizeof(id1)));
	printdims(link,node,buffer);
        ncbytescat(buffer,";\n");
        generateddsattributes(link,node,buffer,depth+1);
    }
    if(name) free(name);
    return OC_NOERR;
}

static OCerror
printdims(OClink link, OCddsnode node, NCbytes* buffer)
{
    int i;
    size_t rank,size;
    OCddsnode dimids[OC_MAX_DIMENSIONS];
    char tmp[1024];
    char id1[1024];

    FAIL(oc_dds_rank(link,node,&rank));
    if(rank == 0) return OC_NOERR;

    FAIL(oc_dds_dimensions(link,node,dimids));
    for(i=0;i<rank;i++) {
	OCddsnode dim = dimids[i];
	char* dimname = NULL;
	FAIL(oc_dimension_properties(link,dim,&size,&dimname));
	if(dimname == NULL)
	    snprintf(tmp,sizeof(tmp),"[%lu]",(unsigned long)size);
        else
	    snprintf(tmp,sizeof(tmp),"[%s=%lu]",idescape(dimname,id1,sizeof(id1)),(unsigned long)size);
	ncbytescat(buffer,tmp);
        if(dimname) free(dimname);
    }
    return OC_NOERR;
}

static OCerror
generateddsattributes(OClink link, OCddsnode node, NCbytes* buffer, int depth)
{
    size_t i,j;
    char tmp[128];
    size_t nattrs;
    char* aname = NULL;
    char* name = NULL;
    OCtype atomtype;
    size_t nvalues;
    char** values = NULL;
    char id1[1024];

    FAIL(oc_dds_attr_count(link,node,&nattrs));
    FAIL(oc_dds_name(link,node,&name));

    if(ocopt.showattributes && nattrs > 0) {
        for(i=0;i<nattrs;i++) {
            FAIL(oc_dds_attr(link,node,i,NULL,NULL,&nvalues,NULL));
   	    values = (char**)malloc(nvalues*sizeof(char*));
            FAIL(oc_dds_attr(link,node,i,&aname,&atomtype,NULL,values));
            snprintf(tmp,sizeof(tmp),"%s%s %s:%s = ",indent(depth),
                        oc_typetostring(atomtype),idescape(name,id1,sizeof(id1)),aname);
            ncbytescat(buffer,tmp);
            for(j=0;j<nvalues;j++) {
                if(j > 0) ncbytescat(buffer,", ");
		if(needsescapes(values[j])) {
		    char* escaped = stringescape(values[j]);
                    ncbytescat(buffer,"\"");
                    ncbytescat(buffer,escaped);
                    ncbytescat(buffer,"\"");
		    if(escaped) free(escaped);
		} else
                    ncbytescat(buffer,values[j]);
            }
            ncbytescat(buffer,";\n");
	    oc_reclaim_strings(nvalues,values);
	    if(values) free(values);
            if(aname) free(aname);
	    values = NULL;
	    aname = NULL;
        }
    }
    if(name) free(name);
    return OC_NOERR;
}

static char*
generatedas(OClink link, OCddsnode root)
{
    size_t i, nsubnodes;
    char* result;
    NCbytes* buffer = ncbytesnew();

    FAIL(oc_dds_nsubnodes(link,root,&nsubnodes));
    ncbytescat(buffer,"Attributes {\n");
    for(i=0;i<nsubnodes;i++) {
	OCddsnode attr;
	FAIL(oc_dds_ithfield(link,root,i,&attr));
        generatedasr(link,attr,buffer,1);
    }
    ncbytescat(buffer,"}\n");
    result = ncbytesdup(buffer);
    ncbytesfree(buffer);
    return result;
}

static OCerror
generatedasr(OClink link, OCddsnode ddsnode, NCbytes* buffer, int depth)
{
    size_t i,nsubnodes;
    char tmp[256];
    OCtype octype, atomtype;
    char* name = NULL;
    char id1[1024];

    /* get some node info */
    FAIL(oc_dds_name(link,ddsnode,&name));
    FAIL(oc_dds_octype(link,ddsnode,&octype));
    FAIL(oc_dds_atomictype(link,ddsnode,&atomtype));

    if(octype == OC_Attributeset) {
        /* get node subcount */
        FAIL(oc_dds_nsubnodes(link,ddsnode,&nsubnodes));
        snprintf(tmp,sizeof(tmp),"%s%s {\n",indent(depth),idescape(name,id1,sizeof(id1)));
        ncbytescat(buffer,tmp);
        for(i=0;i<nsubnodes;i++) {
	    OCddsnode attr;
	    FAIL(oc_dds_ithfield(link,ddsnode,i,&attr));
            generatedasr(link,attr,buffer,depth+1);
        }
        ncbytescat(buffer,indent(depth));
        ncbytescat(buffer,"}\n");
    } else if(octype == OC_Attribute) {
        /* get some info about the node */
	size_t nvalues;
	FAIL(oc_das_attr_count(link,ddsnode,&nvalues));
        snprintf(tmp,sizeof(tmp),"%s%s %s",indent(depth),
                oc_typetostring(atomtype),idescape(name,id1,sizeof(id1)));
        ncbytescat(buffer,tmp);
        for(i=0;i<nvalues;i++) {
            char* value;
            OCtype ptype;
            FAIL(oc_das_attr(link,ddsnode,i,&ptype,&value));
            if(i > 0) ncbytescat(buffer,",");
            if(ptype == OC_String || ptype == OC_URL) {
                char* se = stringescape(value);
                snprintf(tmp,sizeof(tmp)," \"%s\"",se);
                free(se);
            } else
                snprintf(tmp,sizeof(tmp)," %s",value);
            ncbytescat(buffer,tmp);
            free(value);
        }
        ncbytescat(buffer,";\n");
    } else {
        snprintf(tmp,sizeof(tmp),"ocget DAS: unexpected type: %d",(int)octype);
        ncbytescat(buffer,tmp);
    }
    if(name) free(name);
    return OC_NOERR;
}

static char hexdigits[] = "0123456789abcdef";

/* Add escape characters to a string */
static char*
stringescape(char* s)
{
    size_t len;
    char* p;
    int c;
    char* escapedstring;

    if(s == NULL) return NULL;
    len = strlen(s);
    escapedstring = (char*)malloc(4*len);
    p = escapedstring;
    while((c=*s++)) {
        if(c == '"' || c == '\\') {*p++ = '\\'; *p++ = c;}
        else if (c < ' ' || c >= 127) {
            *p++ = '\\'; *p++ = 'x';
            *p++ = hexdigits[(c & 0xf0)>>4];
            *p++ = hexdigits[(c & 0xf)];
        } else
            *p++ = c;
    }
    *p = '\0';
    return escapedstring;
}

static char idchars[] = "_%";

/* Add escape characters to an identifier */
static char*
idescape(char* id, char* escapeid, size_t esize)
{
    char* p;
    int c;

    if(id == NULL) return NULL;
    p = escapeid;
    *p = '\0';
    esize--; /* leave room for null */
    while(esize-- > 0 && (c=*id++)) {
        if(c >= '0' && c <= '9') {*p++ = c;}
        else if(c >= 'a' && c <= 'z') {*p++ = c;}
        else if(c >= 'A' && c <= 'Z') {*p++ = c;}
        else if(strchr(idchars,c) != NULL) {*p++ = c;}
        else {
            *p++ = '%';
            *p++ = hexdigits[(c & 0xf0)>>4];
            *p++ = hexdigits[(c & 0xf)];
        }
    }
    *p = '\0';
    return escapeid;
}

static char valuechars[] = " \\\"";

/**
Return 1 if the given string, used as a value, should be escaped.
*/
static int
needsescapes(const char* s)
{
    const char* p = s;
    int c;
    while((c=*p++)) {
	if(strchr(valuechars,c) != NULL)
	    return 1; /* needs to be escaped */
    }
    return 0;
}


static OCerror
dumpdatanode(OClink link, OCdatanode datanode, size_t count, void* memory, NCbytes* buffer)
{
    size_t i;
    size_t delta;
    OCddsnode node;
    OCtype atomtype;
    OCtype octype;
    NCbytes* path = NULL;
    char* name = NULL;
    char id[1024];
    char tmp[1024];
    struct DUMPPATH* entry = NULL;

    FAIL(oc_data_ddsnode(link,datanode,&node));
    FAIL(oc_dds_octype(link,node,&octype));
    FAIL(oc_dds_atomictype(link,node,&atomtype));
    delta = oc_typesize(atomtype);

#ifdef TRACK
    printstack("dumpdatanode");
#endif

    /* construct the fully qualified name from the stack; watch out for duplicates
       from e.g. sequence versus record */
    path = ncbytesnew();
    for(i=0;i<stacknext;i++) {
	entry = stack + i;
	/* check for duplicate */
	if(i<(stacknext-1) && entry->node == stack[i+1].node) continue;

	/* Get various pieces of additional node information */
        FAIL(oc_dds_name(glink,entry->node,&name));
        (void)idescape(name,id,sizeof(id));
	if(name) free(name);

        switch (entry->octype) {

        case OC_Dataset:
            break;

        case OC_Structure:
            ncbytescat(path,"/");
            ncbytescat(path,id);
            if(entry->rank > 0) {
		for(i=0;i<entry->rank;i++) {
                    sprintf(tmp,"[%lu]",(unsigned long)entry->indices[i]);
                    ncbytescat(path,tmp);
		}
            }
	    break;

        case OC_Grid:
            ncbytescat(path,"/");
            ncbytescat(path,id);
	    break;

        case OC_Sequence:
            ncbytescat(path,"/");
            ncbytescat(path,id);
            sprintf(tmp,"[%lu]",(unsigned long)entry->indices[0]);
            ncbytescat(path,tmp);
            break;

        case OC_Atomic:
            ncbytescat(path,"/");
            ncbytescat(path,id);
            break; /* deal with below */

        default: ncbytescat(path,"?"); break;
        }
    }
    /* work with the final entry */
    assert(entry == (stack + (stacknext - 1)));
    assert(datanode == entry->datanode);
    snprintf(tmp,sizeof(tmp),"%s %s",
                oc_typetostring(atomtype),
                ncbytescontents(path));
    ncbytescat(buffer,tmp);
    if(entry->rank > 0) {
	if(ocopt.octest) { /* Match the octest output */
	    off_t xproduct;
	    xproduct = totaldimsize(entry->rank,entry->dimsizes);
            snprintf(tmp,sizeof(tmp),"[0..%lu]",(unsigned long)xproduct-1);
            ncbytescat(buffer,tmp);
	} else {
	    for(i=0;i<entry->rank;i++) {
                snprintf(tmp,sizeof(tmp),"[0..%lu]",(unsigned long)entry->dimsizes[i]-1);
                ncbytescat(buffer,tmp);
	    }
	}
    }
    count = totaldimsize(entry->rank,entry->dimsizes);
    for(i=0;i<count;i++) {
        char *memory_local = memory;
        ncbytescat(buffer," ");
        oc_typeprint(atomtype,memory_local+(i*delta),sizeof(tmp),tmp);
        ncbytescat(buffer,tmp);
    }
    ncbytescat(buffer,"\n");
    ncbytesfree(path);
    return OC_NOERR;
}

static off_t
odom_init(size_t rank, size_t* indices, size_t* dimsizes)
{
    int i;
    off_t count;
    for(count=1,i=0;i<rank;i++) {
        indices[i] = 0;
	count *= dimsizes[i];
    }
    return count;
}

static void
odom_next(size_t rank, size_t* indices, size_t* dimsizes)
{
    int i;
    for(i=rank-1;i>=0;i--) {
	indices[i]++;
	if(indices[i] < dimsizes[i]) break;
	if(i > 0) indices[i] = 0;
    }
}

/* Return 0 if we have exhausted the indices, 1 otherwise */
static int
odom_more(size_t rank, size_t* indices, size_t* dimsizes)
{
    if(indices[0] >= dimsizes[0]) return 0;
    return 1;
}

/* Compute total # of elements if dimensioned */
static size_t
totaldimsize(size_t rank, size_t* sizes)
{
    unsigned int i;
    size_t count = 1;
    for(i=0;i<rank;i++) {
        count *= sizes[i];
    }
    return count;
}

static char*
indent(int n)
{
    size_t nblanks = BLANKSPERDENT * n;
    memset(blanks,' ',nblanks);
    blanks[nblanks] = '\0';
    return blanks;
}

static void
pushstack(OCdatanode datanode)
{
    struct DUMPPATH* entry = stack+stacknext;
    entry->datanode = datanode;
    FAIL(oc_data_ddsnode(glink,entry->datanode,&entry->node));
    FAIL(oc_dds_octype(glink,entry->node,&entry->octype));
    FAIL(oc_dds_rank(glink,entry->node,&entry->rank));
    if(entry->rank > 0) {
        FAIL(oc_dds_dimensionsizes(glink,entry->node,entry->dimsizes));
    }
    entry->indexed = oc_data_indexed(glink,entry->datanode);
    if(entry->indexed) {
	FAIL(oc_data_position(glink,entry->datanode,entry->indices));
    }
    stacknext++;
}


#ifdef TRACK
static void printstack(char* msg)
{
    size_t i,j;
    struct DUMPPATH* entry;
    fprintf(stderr,"\n%s@stack: %u\n",msg,stacknext);
    for(entry=stack,i=0;i<stacknext;i++,entry++) {
	OCdatanode datanode = entry->datanode;
	OCddsnode node;
	size_t rank;
	size_t edges[OC_MAX_DIMENSIONS];
	char* name;
        FAIL(oc_dds_rank(glink,entry->node,&rank));
        if(entry->rank > 0)
            FAIL(oc_dds_dimensionsizes(glink,entry->node,edges));
        FAIL(oc_dds_name(glink,node,&name));
        fprintf(stderr,"    [%d] (%s)",(int)i,name)
	for(j=0;j<rank;j++)
            fprintf(stderr,"[%u]",(unsigned int)edges[j]);
        fprintf(stderr,"\n");
    }
}
#endif

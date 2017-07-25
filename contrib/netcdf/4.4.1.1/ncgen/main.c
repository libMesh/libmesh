/*********************************************************************
 *   Copyright 1993, UCAR/Unidata
 *   See netcdf/COPYRIGHT file for copying and redistribution conditions.
 *********************************************************************/
/* $Id: main.c,v 1.33 2010/05/26 21:43:36 dmh Exp $ */
/* $Header: /upc/share/CVS/netcdf-3/ncgen/main.c,v 1.33 2010/05/26 21:43:36 dmh Exp $ */

#include "includes.h"
#include "offsets.h"
#ifdef HAVE_GETOPT_H
#include <getopt.h>
#endif

#ifdef _MSC_VER
#include "XGetopt.h"
#define snprintf _snprintf
int opterr;
int optind;
#endif

/* Default is netcdf-3 mode 1 */
#define DFALTCMODE 0

extern void init_netcdf(void);
extern void parse_init(void);
extern int ncgparse(void);

/* For error messages */
char* progname;
char* cdlname;

/* option flags */
int nofill_flag;
char* mainname; /* name to use for main function; defaults to "main"*/
Language l_flag;
int syntax_only;
int header_only;

/* flags for tracking what output format to use */
int k_flag;    /* > 0  => -k was specified on command line*/
int format_attribute; /* 1=>format came from format attribute */
int enhanced_flag; /* 1 => netcdf-4 */
int cdf5_flag; /* 1 => cdf5 | maybe netcdf-4 */
int specials_flag; /* 1=> special attributes are present */
int usingclassic;
int cmode_modifier;
int diskless;
int ncloglevel;

GlobalSpecialData globalspecials;

char* binary_ext = ".nc";

size_t nciterbuffersize;

struct Vlendata* vlendata;

char *netcdf_name; /* command line -o file name */
char *datasetname; /* name from the netcdf <name> {} */

extern FILE *ncgin;

/* Forward */
static char* ubasename(char*);
void usage( void );
int main( int argc, char** argv );

/* Define tables vs modes for legal -k values*/
struct Kvalues legalkinds[NKVALUES] = {
    /* NetCDF-3 classic format (32-bit offsets) */
    {"classic", NC_FORMAT_CLASSIC}, /* canonical format name */
    {"nc3", NC_FORMAT_CLASSIC},	    /* short format name */
    {"1", NC_FORMAT_CLASSIC},	/* deprecated, use "-3" or "-k nc3" instead */

    /* NetCDF-3 64-bit offset format */
    {"64-bit offset", NC_FORMAT_64BIT_OFFSET}, /* canonical format name */
    {"nc6", NC_FORMAT_64BIT_OFFSET},		/* short format name */
    {"2", NC_FORMAT_64BIT_OFFSET},     /* deprecated, use "-6" or "-k nc6" instead */
    {"64-bit-offset", NC_FORMAT_64BIT_OFFSET}, /* aliases */

    /* NetCDF-4 HDF5-based format */
    {"netCDF-4", NC_FORMAT_NETCDF4}, /* canonical format name */
    {"nc4", NC_FORMAT_NETCDF4},	     /* short format name */
    {"3", NC_FORMAT_NETCDF4},   /* deprecated, use "-4" or "-k nc4" instead */
    {"netCDF4", NC_FORMAT_NETCDF4},  /* aliases */
    {"hdf5", NC_FORMAT_NETCDF4},
    {"enhanced", NC_FORMAT_NETCDF4},
    {"netcdf-4", NC_FORMAT_NETCDF4},
    {"netcdf4", NC_FORMAT_NETCDF4},

    /* NetCDF-4 HDF5-based format, restricted to classic data model */
    {"netCDF-4 classic model", NC_FORMAT_NETCDF4_CLASSIC}, /* canonical format name */
    {"nc7", NC_FORMAT_NETCDF4_CLASSIC}, /* short format name */
    {"4", NC_FORMAT_NETCDF4_CLASSIC}, /* deprecated, use "-7" or -k nc7" instead */
    {"netCDF-4-classic", NC_FORMAT_NETCDF4_CLASSIC}, /* aliases */
    {"netCDF-4_classic", NC_FORMAT_NETCDF4_CLASSIC},
    {"netCDF4_classic", NC_FORMAT_NETCDF4_CLASSIC},
    {"hdf5-nc3", NC_FORMAT_NETCDF4_CLASSIC},
    {"enhanced-nc3", NC_FORMAT_NETCDF4_CLASSIC},

    /* CDF-5 format */
    {"5", NC_FORMAT_64BIT_DATA},
    {"64-bit-data", NC_FORMAT_64BIT_DATA},
    {"64-bit data", NC_FORMAT_64BIT_DATA},
    {"nc5", NC_FORMAT_64BIT_DATA},
    {"cdf5", NC_FORMAT_64BIT_DATA},
    {"cdf-5", NC_FORMAT_64BIT_DATA},

    /* null terminate*/
    {NULL,0}
};

#ifndef _MSC_VER
struct Languages {
    char* name;
    Language flag;
} legallanguages[] = {
{"b", L_BINARY},
{"c", L_C},
{"C", L_C},
{"f77", L_F77},
{"fortran77", L_F77},
{"Fortran77", L_F77},
{"j", L_JAVA},
{"java", L_JAVA},
{NULL,L_UNDEFINED}
};
#else
typedef struct Languages {
		char* name;
		Language flag;
} Languages;

struct Languages legallanguages[] = {
{"b", L_BINARY},
{"c", L_C},
{"C", L_C},
{"f77", L_F77},
{"fortran77", L_F77},
{"Fortran77", L_F77},
{"j", L_JAVA},
{"java", L_JAVA},
{NULL,L_UNDEFINED}
};
#endif

#if 0 /*not used*/
/* BOM Sequences */
static char* U8   = "\xEF\xBB\xBF";    /* UTF-8 */
static char* BE32 = "\x00\x00\xFE\xFF"; /* UTF-32; big-endian */
static char* LE32 = "\xFF\xFE";       /* UTF-32; little-endian */
static char* BE16 = "\xFE\xFF";       /* UTF-16; big-endian */
static char* LE16 = "\xFF\xFE";       /* UTF-16; little-endian */
#endif

/* The default minimum iterator size depends
   on whether we are doing binary or language
   based output.
*/
#define DFALTBINNCITERBUFFERSIZE  0x40000 /* about 250k bytes */
#define DFALTLANGNCITERBUFFERSIZE  0x4000 /* about 15k bytes */

/* strip off leading path */
/* result is malloc'd */

static char *
ubasename(char *logident)
{
    char* sep;

    sep = strrchr(logident,'/');
#ifdef MSDOS
    if(sep == NULL) sep = strrchr(logident,'\\');
#endif
    if(sep == NULL) return logident;
    sep++; /* skip past the separator */
    return sep;
}

void
usage(void)
{
    derror("Usage: %s"
" [-1]"
" [-3]"
" [-4]"
" [-5]"
" [-6]"
" [-7]"
" [-b]"
" [-B buffersize]"
" [-d]"
" [-D debuglevel]"
" [-h]"
" [-k kind ]"
" [-l language=b|c|f77|java]"
" [-M <name>]"
" [-n]"
" [-o outfile]"
" [-P]"
" [-x]"
" [file ... ]",
	   progname);
    derror("netcdf library version %s", nc_inq_libvers());
}

int
main(
	int argc,
	char *argv[])
{
    int c;
    FILE *fp;
	struct Languages* langs;
    char* lang_name;
#ifdef __hpux
    setlocale(LC_CTYPE,"");
#endif

    init_netcdf();

    opterr = 1;			/* print error message if bad option */
    progname = ubasename(argv[0]);
    cdlname = "-";
    netcdf_name = NULL;
    datasetname = NULL;
    l_flag = 0;
    nofill_flag = 0;
    syntax_only = 0;
    header_only = 0;
    mainname = "main";
    nciterbuffersize = 0;

    k_flag = 0;
    format_attribute = 0;
    enhanced_flag = 0;
    cdf5_flag = 0;
    specials_flag = 0;
    diskless = 0;
#ifdef LOGGING
    ncloglevel = NC_TURN_OFF_LOGGING;
#else
    ncloglevel = -1;
#endif
    memset(&globalspecials,0,sizeof(GlobalSpecialData));

#if _CRAYMPP && 0
    /* initialize CRAY MPP parallel-I/O library */
    (void) par_io_init(32, 32);
#endif

    while ((c = getopt(argc, argv, "134567bB:cdD:fhHk:l:M:no:Pv:xL:")) != EOF)
      switch(c) {
	case 'd':
	  debug = 1;
	  break;
	case 'D':
	  debug = atoi(optarg);
	  break;
	case 'c': /* for c output, old version of "-lc" */
	  if(l_flag != 0) {
	    fprintf(stderr,"Please specify only one language\n");
	    return 1;
	  }
	  l_flag = L_C;
	  fprintf(stderr,"-c is deprecated: please use -lc\n");
	  break;
	case 'f': /* for f77 output, old version of "-lf" */
	  if(l_flag != 0) {
	    fprintf(stderr,"Please specify only one language\n");
	    return 1;
	  }
	  l_flag = L_F77;
	  fprintf(stderr,"-f is deprecated: please use -lf77\n");
	  break;
	case 'b': /* for binary netcdf output, ".nc" extension */
	  if(l_flag != 0) {
	    fprintf(stderr,"Please specify only one language\n");
	    return 1;
	  }
	  l_flag = L_BINARY;
	  break;
	case 'h':
	  header_only = 1;
	  break;
	case 'H':
	  usage();
	  exit(0);
        case 'l': /* specify language, instead of using -c or -f or -b */
	{
	    if(l_flag != 0) {
              fprintf(stderr,"Please specify only one language\n");
              return 1;
	    }
            if(!optarg) {
              derror("%s: output language is null", progname);
              return(1);
            }
            lang_name = (char*) emalloc(strlen(optarg)+1);
	    (void)strcpy(lang_name, optarg);
	    for(langs=legallanguages;langs->name != NULL;langs++) {
              if(strcmp(lang_name,langs->name)==0) {
	  	l_flag = langs->flag;
                break;
              }
	    }
	    if(langs->name == NULL) {
              derror("%s: output language %s not implemented",progname, lang_name);
              return(1);
	    }
	}; break;
	case 'L':
	    ncloglevel = atoi(optarg);
	    break;
	case 'n':		/* old version of -b, uses ".cdf" extension */
	  if(l_flag != 0) {
	    fprintf(stderr,"Please specify only one language\n");
	    return 1;
	  }
	  l_flag = L_BINARY;
          binary_ext = ".cdf";
	  break;
	case 'o':		/* to explicitly specify output name */
	  netcdf_name = nulldup(optarg);
	  break;
	case 'x': /* set nofill mode to speed up creation of large files */
	  nofill_flag = 1;
	  break;
        case 'v': /* a deprecated alias for "kind" option */
	    /*FALLTHRU*/
	case 'k': { /* for specifying variant of netCDF format to be generated
		     Possible values are:
		     Format names:
		       "classic" or "nc3"
		       "64-bit offset" or "nc6"
		       "64-bit data" or "nc5" or "cdf-5"
		       "netCDF-4" or "nc4"
		       "netCDF-4 classic model" or "nc7"
		       "netCDF-5" or "nc5" or "cdf5"
		     Format version numbers (deprecated):
		       1 (=> classic)
		       2 (=> 64-bit offset)
		       3 (=> netCDF-4)
		       4 (=> netCDF-4 classic model)
                       5 (=> classic 64 bit data aka CDF-5)
		   */
	    struct Kvalues* kvalue;
	    char *kind_name = (optarg != NULL
				? (char *) emalloc(strlen(optarg)+1)
				: emalloc(1));
	    if (! kind_name) {
		derror ("%s: out of memory", progname);
		return(1);
	    }
            if(optarg != NULL)
              (void)strcpy(kind_name, optarg);
            for(kvalue=legalkinds;kvalue->name;kvalue++) {
              if(strcmp(kind_name,kvalue->name) == 0) {
                k_flag = kvalue->k_flag;
                break;
              }
            }
            if(kvalue->name == NULL) {
                derror("Invalid format: %s",kind_name);
                return 2;
            }
	} break;
	case '3':		/* output format is classic (netCDF-3) */
	    k_flag = NC_FORMAT_CLASSIC;
	    break;
	case '6':		/* output format is 64-bit-offset (netCDF-3 version 2) */
	    k_flag = NC_FORMAT_64BIT_OFFSET;
	    break;
	case '4':		/* output format is netCDF-4 (variant of HDF5) */
	    k_flag = NC_FORMAT_NETCDF4;
	    break;
	case '5':		/* output format is CDF5 */
	    k_flag = NC_FORMAT_CDF5;
	    break;
	case '7':		/* output format is netCDF-4 (restricted to classic model)*/
	    k_flag = NC_FORMAT_NETCDF4_CLASSIC;
	    break;
	case 'M': /* Determine the name for the main function */
	    mainname = nulldup(optarg);
	    break;
	case 'B':
	  nciterbuffersize = atoi(optarg);
	  break;
	case 'P': /* diskless with persistence */
	  diskless = 1;
	  break;
	case '?':
	  usage();
	  return(8);
      }

    if(l_flag == 0) {
	l_flag = L_BINARY; /* default */
	/* Treat -k or -o as an implicit -lb assuming no other -l flags */
        if(k_flag == 0 && netcdf_name == NULL)
	    syntax_only = 1;
    }

    /* Compute/default the iterator buffer size */
    if(l_flag == L_BINARY) {
	if(nciterbuffersize == 0 )
	    nciterbuffersize = DFALTBINNCITERBUFFERSIZE;
    } else {
	if(nciterbuffersize == 0)
	    nciterbuffersize = DFALTLANGNCITERBUFFERSIZE;
    }

#ifndef ENABLE_C
    if(c_flag) {
	  fprintf(stderr,"C not currently supported\n");
	  exit(1);
    }
#endif
#ifndef ENABLE_BINARY
    if(l_flag == L_BINARY) {
	  fprintf(stderr,"Binary netcdf not currently supported\n");
	  exit(1);
    }
#endif
#ifndef ENABLE_JAVA
    if(l_flag == L_JAVA) {
	  fprintf(stderr,"Java not currently supported\n");
	  exit(1);
    }
#else
    if(l_flag == L_JAVA && mainname != NULL && strcmp(mainname,"main")==0)
      mainname = "Main";
#endif
#ifndef ENABLE_F77
    if(l_flag == L_F77) {
	  fprintf(stderr,"F77 not currently supported\n");
	  exit(1);
    }
#endif

    if(l_flag != L_BINARY)
	diskless = 0;

    argc -= optind;
    argv += optind;

    if (argc > 1) {
	derror ("%s: only one input file argument permitted",progname);
	return(6);
    }

    fp = stdin;
    if (argc > 0 && strcmp(argv[0], "-") != 0) {
	char bom[4];
	size_t count;
	if ((fp = fopen(argv[0], "r")) == NULL) {
	    derror ("can't open file %s for reading: ", argv[0]);
	    perror("");
	    return(7);
	}
   	/* Check the leading bytes for an occurrence of a BOM */
        /* re: http://www.unicode.org/faq/utf_bom.html#BOM */
	/* Attempt to read the first four bytes */
	memset(bom,0,sizeof(bom));
	count = fread(bom,1,2,fp);
	if(count == 2) {
	    switch (bom[0]) {
	    case '\x00':
	    case '\xFF':
	    case '\xFE':
	        /* Only UTF-* is allowed; complain and exit */
		fprintf(stderr,"Input file contains a BOM indicating a non-UTF8 encoding\n");
		return 1;
	    case '\xEF':
		/* skip the BOM */
	        (void)fread(bom,1,1,fp);
	        break;
	    default: /* legal printable char, presumably; rewind */
	        rewind(fp);
		break;
	    }
	}

	cdlname = (char*)emalloc(NC_MAX_NAME);
	cdlname = nulldup(argv[0]);
	if(cdlname != NULL) {
	  if(strlen(cdlname) > NC_MAX_NAME)
	    cdlname[NC_MAX_NAME] = '\0';
	}
    }

    parse_init();
    ncgin = fp;
    if(debug >= 2) {ncgdebug=1;}
    if(ncgparse() != 0)
        return 1;

    /* Compute the k_flag (1st pass) using rules in the man page (ncgen.1).*/

#ifndef USE_NETCDF4
    if(enhanced_flag) {
	derror("CDL input is enhanced mode, but --disable-netcdf4 was specified during build");
	return 0;
    }
#endif

    if(l_flag == L_JAVA || l_flag == L_F77) {
        k_flag = 1;
	if(enhanced_flag) {
	    derror("Java or Fortran requires classic model CDL input");
	    return 0;
	}
    }

    if(k_flag == 0)
      k_flag = globalspecials._Format;

    if(cdf5_flag && !enhanced_flag && k_flag == 0)
      k_flag = 5;
    if(enhanced_flag && k_flag == 0)
      k_flag = 3;

    if(enhanced_flag && k_flag != 3) {
      if(enhanced_flag && k_flag != 3 && k_flag != 5) {
        derror("-k or _Format conflicts with enhanced CDL input");
        return 0;
      }
    }

    if(specials_flag > 0 && k_flag == 0)
#ifdef USE_NETCDF4
	k_flag = 3;
#else
	k_flag = 1;
#endif

    if(k_flag == 0)
	k_flag = 1;

    /* Figure out usingclassic */
    switch (k_flag) {
    case NC_FORMAT_64BIT_DATA:
    case NC_FORMAT_CLASSIC:
    case NC_FORMAT_64BIT_OFFSET:
    case NC_FORMAT_NETCDF4_CLASSIC:
	usingclassic = 1;
	break;
    case NC_FORMAT_NETCDF4:
    default:
	usingclassic = 0;
	break;
    }

    /* compute cmode_modifier */
    switch (k_flag) {
    case NC_FORMAT_CLASSIC:
	cmode_modifier = 0; break;
    case NC_FORMAT_64BIT_OFFSET:
	cmode_modifier = NC_64BIT_OFFSET; break;
    case NC_FORMAT_NETCDF4:
	cmode_modifier = NC_NETCDF4; break;
    case NC_FORMAT_NETCDF4_CLASSIC:
	cmode_modifier = NC_NETCDF4 | NC_CLASSIC_MODEL; break;
    case NC_FORMAT_64BIT_DATA:
	cmode_modifier = NC_CDF5; break;
    default: ASSERT(0); /* cannot happen */
    }

    if(diskless)
	cmode_modifier |= (NC_DISKLESS|NC_NOCLOBBER);

    processsemantics();
    if(!syntax_only && error_count == 0)
        define_netcdf();

    return 0;
}
END_OF_MAIN()

void
init_netcdf(void) /* initialize global counts, flags */
{
    compute_alignments();
    memset((void*)&nullconstant,0,sizeof(NCConstant));
    fillconstant = nullconstant;
    fillconstant.nctype = NC_FILLVALUE;

    codebuffer = bbNew();
    stmt = bbNew();
    error_count = 0; /* Track # of errors */
}

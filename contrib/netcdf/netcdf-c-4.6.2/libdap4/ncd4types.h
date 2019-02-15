/*********************************************************************
 *   Copyright 2016, UCAR/Unidata
 *   See netcdf/COPYRIGHT file for copying and redistribution conditions.
 *********************************************************************/

/*
Type declarations and associated constants
are defined here.	
*/

#ifndef D4TYPES_H
#define D4TYPES_H 1

#undef COMPILEBYDEFAULT

#include "ncrc.h"
#include "ncauth.h"

/*
Control if struct fields can be map targets.
Currently turned off because semantics are unclear.
*/
#undef ALLOWFIELDMAPS

#define long64 long long
#define ncerror int

/* Misc. code controls */
#define FILLCONSTRAINT TRUE

#define DEFAULTSTRINGLENGTH 64

/* Size of the checksum */
#define CHECKSUMSIZE 4

/**************************************************/
/* sigh, do the forwards */

typedef struct NCD4INFO NCD4INFO;
typedef enum NCD4CSUM NCD4CSUM;
typedef enum NCD4mode NCD4mode;
typedef enum NCD4translation NCD4translation;
typedef struct NCD4curl NCD4curl;
typedef struct NCD4meta NCD4meta;
typedef struct NCD4node NCD4node;
typedef struct NCD4params NCD4params;

/**************************************************/
/* DMR Tree node sorts */

/* Concepts from netcdf-4 data model */
/* Define as powers of 2 so we can create a set */
typedef enum NCD4sort {
    NCD4_NULL=0,
    NCD4_ATTR=1,
    NCD4_ATTRSET=2,
    NCD4_XML=4, /* OtherXML */
    NCD4_DIM=8,
    NCD4_GROUP=16, /* Including Dataset */
    NCD4_TYPE=32, /* atom types, opaque, enum, struct or seq */
    NCD4_VAR=64, /* Variable or field */
    NCD4_ECONST=128,
} NCD4sort;

#define ISA(sort,flags) ((sort) & (flags))

#define ISATTR(sort) ISA((sort),(NCD4_ATTR))
#define ISDIM(sort) ISA((sort),(NCD4_DIM))
#define ISGROUP(sort) ISA((sort),(NCD4_GROUP))
#define ISTYPE(sort) ISA((sort),(NCD4_TYPE))
#define ISVAR(sort) ISA((sort),(NCD4_VAR))
#define ISECONST(sort) ISA((sort),(NCD4_ECONST))

/* Define some nc_type aliases */
#define NC_NULL NC_NAT
#define NC_SEQ NC_VLEN
#define NC_STRUCT NC_COMPOUND

#define ISCMPD(subsort) ((subsort) == NC_SEQ || (subsort) == NC_STRUCT)

/**************************************************/
/* Special attributes; When defining netcdf-4,
   these are suppressed, except for UCARTAGMAPS
*/

#define RESERVECHAR '_'

#define UCARTAG         "_edu.ucar."
#define UCARTAGVLEN     "_edu.ucar.isvlen"
#define UCARTAGOPAQUE   "_edu.ucar.opaque.size"
#define UCARTAGORIGTYPE "_edu.ucar.orig.type"
#define UCARTAGUNLIM    "_edu.ucar.isunlimited"

/* These are attributes inserted into the netcdf-4 file */
#define NC4TAGMAPS      "_edu.ucar.maps"

/**************************************************/
/* Misc.*/

/* Define possible translation modes */
enum NCD4translation {
NCD4_NOTRANS = 0, /* Apply straight DAP4->NetCDF4 translation */
NCD4_TRANSNC4 = 1, /* Use _edu.ucar flags to achieve better translation */
};

/* Define possible debug flags */
#define NCF_DEBUG_NONE  0
#define NCF_DEBUG_COPY  1 /* Dump data into the substrate and close it rather than abortiing it */

/* Define possible retrieval modes */
enum NCD4mode {
NCD4_DMR = 1,
NCD4_DAP = 2,
NCD4_DSR = 4
};


/* Define storage for all the primitive types (plus vlen) */
union ATOMICS {
    char i8[8];
    unsigned char u8[8];
    short i16[4];
    unsigned short u16[4];
    int i32[2];
    unsigned int u32[2];
    long long i64[1];
    unsigned long long u64[1];
    float f32[2];
    double f64[1];
#if SIZEOF_VOIDP == 4
    char* s[2];
#elif SIZEOF_VOIDP == 8
    char* s[1];
#endif
#if (SIZEOF_VOIDP + SIZEOF_SIZE_T) == 8
    nc_vlen_t vl[1];
#endif
};

/**************************************************/
/* !Node type for the NetCDF-4 metadata produced from
   parsing the DMR tree.
   We only use a single node type tagged with the sort.
   Some information is not kept (e.g. attributes).
*/
struct NCD4node {
    NCD4sort sort; /* gross discriminator */
    nc_type subsort; /* subcases of sort */
    char* name; /* Raw unescaped */
    NCD4node*  container; /* parent object: e.g. group, enum, compound...*/
    int visited; /* for recursive walks of all nodes */
    /* Sort specific fields */
    NClist* groups;   /* NClist<NCD4node*> groups in group */
    NClist* vars;   /* NClist<NCD4node*> vars in group, fields in struct/seq */
    NClist* types;   /* NClist<NCD4node*> types in group */
    NClist* dims;    /* NClist<NCD4node*>; dimdefs in group, dimrefs in vars */
    NClist* attributes; /* NClist<NCD4node*> */
    NClist* maps;       /* NClist<NCD4node*> */
    NClist* xmlattributes; /* NClist<String> */
    NCD4node* basetype;
    struct { /* sort == NCD4_ATTRIBUTE */
        NClist* values;
    } attr;
    struct { /* sort == NCD4_OPAQUE */
	long long size; /* 0 => var length */
    } opaque;
    struct { /* sort == NCD4_DIMENSION */
	long long size;
	int isunlimited;
	int isanonymous;
    } dim;
    struct { /* sort == NCD4_ECONST || sort == NCD4_ENUM */    
        union ATOMICS ecvalue;
        NClist* econsts; /* NClist<NCD4node*> */
    } en;
    struct { /* sort == NCD4_GROUP */
	NClist* elements;   /* NClist<NCD4node*> everything at the top level in a group */
        int isdataset;
        char* dapversion;
        char* dmrversion;
        char* datasetname;
        NClist* varbyid; /* NClist<NCD4node*> indexed by varid */
    } group;
    struct { /* Meta Info */
        int id; /* Relevant netcdf id; interpretation depends on sort; */
	int isfixedsize; /* sort == NCD4_TYPE; Is this a fixed size (recursively) type? */
	d4size_t dapsize; /* size of the type as stored in the dap data; will, as a rule,
                             be same as memsize only for types <= NC_UINT64 */
        nc_type cmpdid; /*netcdf id for the compound type created for seq type */
	size_t memsize; /* size of a memory instance without taking dimproduct into account,
                           but taking compound alignment into account  */
        d4size_t offset; /* computed structure field offset in memory */
        size_t alignment; /* computed structure field alignment in memory */
    } meta;
    struct { /* Data compilation info */
        int flags; /* See d4data for actual flags */
	D4blob dap4data; /* offset and start pos for this var's data in serialization */
        unsigned int remotechecksum; /* toplevel variable checksum as sent by server*/    
        unsigned int localchecksum; /* toplevel variable checksum as computed by client */    
    } data;
    struct { /* Track netcdf-4 conversion info */
	int isvlen;	/*  _edu.ucar.isvlen */
	/* Split UCARTAGORIGTYPE into group plus name */
	struct {
  	    NCD4node* group;
	    char* name;
	} orig;
	/* Represented elsewhare: _edu.ucar.opaque.size */ 
    } nc4;
};

/* Tracking info about the serialized input before and after de-chunking */
typedef struct NCD4serial {
    size_t rawsize; /* |rawdata| */ 
    void* rawdata;
    size_t dapsize; /* |dapdata|; this is transient */
    void* dap; /* pointer into rawdata where dap data starts */ 
    char* dmr;/* copy of dmr */ 
    char* errdata; /* null || error chunk (null terminated) */
    int hostlittleendian; /* 1 if the host is little endian */
    int remotelittleendian; /* 1 if the packet says data is little endian */
    int remotechecksumming; /* 1 if the packet says checksums are included */
} NCD4serial;

/* This will be passed out of the parse */
struct NCD4meta {
    NCD4INFO* controller;
    int ncid; /* root ncid of the substrate netcdf-4 file;
		 warning: copy of NCD4Info.substrate.nc4id */
    NCD4node* root;
    NCD4mode  mode; /* Are we reading DMR (only) or DAP (includes DMR) */
    NClist* allnodes; /*list<NCD4node>*/
    struct Error { /* Content of any error response */
	char* parseerror;
	int   httpcode;
	char* message;
	char* context;
	char* otherinfo;
    } error;
    int debuglevel;
    NCD4serial serial;
    int ignorechecksums; /* 1=> compute but ignore */
    int localchecksumming; /* 1=>compute local checksum */
    int swap; /* 1 => swap data */
    /* Define some "global" (to a DMR) data */
    NClist* groupbyid; /* NClist<NCD4node*> indexed by groupid >> 16; this is global */
    NCD4node* _bytestring; /* If needed */
};

typedef struct NCD4parser {
    char* input;
    int debuglevel;
    NCD4meta* metadata;
    /* Capture useful subsets of dataset->allnodes */
    NClist* types; /*list<NCD4node>*/
    NClist* dims; /*list<NCD4node>*/
    NClist* vars; /*list<NCD4node>*/
    NClist* groups; /*list<NCD4node>*/
    /* Convenience for short cut fqn detection */
    NClist* atomictypes; /*list<NCD4node>*/
    char* used; /* mark indices in atomictypes that have been used */
    NCD4node* dapopaque; /* Single non-fixed-size opaque type */
} NCD4parser;

/**************************************************/

/* Curl info */
struct NCD4curl {
    CURL* curl; /* curl handle*/
    NCbytes* packet; 
    struct errdata {/* Hold info for an error return from server */
	char* code;
	char* message;
	long  httpcode;
	char  errorbuf[CURL_ERROR_SIZE]; /* CURLOPT_ERRORBUFFER*/
    } errdata;
    struct {
	int active; /* Activate keepalive? */
	long idle; /* KEEPIDLE value */
	long interval; /* KEEPINTVL value */
    } keepalive; /* keepalive info */
    long buffersize; /* read buffer size */    
};

/**************************************************/
/* Define a structure holding common info */

struct NCD4INFO {
    NC*   controller; /* Parent instance of NCD4INFO */
    char* rawurltext; /* as given to ncd4_open */
    char* urltext;    /* as modified by ncd4_open */
    NCURI* uri;      /* parse of rawuritext */
    NCD4curl* curl;
    int inmemory; /* store fetched data in memory? */
    struct {
	char*   memory;   /* allocated memory if ONDISK is not set */
        char*   ondiskfilename; /* If ONDISK is set */
        FILE*   ondiskfile;     /* ditto */
        d4size_t   datasize; /* size on disk or in memory */
        long dmrlastmodified;
        long daplastmodified;
    } data;
    struct {
	int realfile; /* 1 => we created actual temp file */
	char* filename; /* of the substrate file */
        int nc4id; /* substrate nc4 file ncid used to hold metadata; not same as external id  */
	NCD4meta* metadata;
    } substrate;
    struct {
        NCCONTROLS  flags;
        NCCONTROLS  debugflags;
	NCD4translation translation;
	char substratename[NC_MAX_NAME];
	size_t opaquesize; /* default opaque size */
    } controls;
    NCauth auth;
    struct {
	char* filename;
    } fileproto;
    NClist* blobs;
};

#endif /*D4TYPES_H*/

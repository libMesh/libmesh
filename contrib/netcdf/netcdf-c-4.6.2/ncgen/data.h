/*********************************************************************
 *   Copyright 2009, UCAR/Unidata
 *   See netcdf/COPYRIGHT file for copying and redistribution conditions.
 *********************************************************************/

#ifndef DATA_H
#define DATA_H 1

#ifndef NO_STDARG
#  include <stdarg.h>
#else
#  include <varargs.h>
#endif

/* nmemonics*/
#define TOPLEVEL 1

/* Forward types */
struct Datalist;
struct Symbol;
struct Dimset;
typedef struct Generator Generator;

/* any one possible value*/
typedef union Constvalue {
    struct Datalist* compoundv; /* NC_COMPOUND*/
    char charv;                 /* NC_CHAR*/
    signed char int8v;          /* NC_BYTE*/
    unsigned char uint8v;       /* NC_UBYTE*/
    short int16v;               /* NC_SHORT*/
    unsigned short uint16v;     /* NC_USHORT*/
    int int32v;                 /* NC_INT*/
    unsigned int uint32v;       /* NC_UINT*/
    long long int64v;           /* NC_INT64*/
    unsigned long long uint64v; /* NC_UINT64*/
    float floatv;               /* NC_FLOAT*/
    double doublev;             /* NC_DOUBLE*/
    struct Stringv {		/* NC_STRING*/
	int len;
	char* stringv; 
    } stringv;
    struct Opaquev {     /* NC_OPAQUE*/
	int len; /* length as originally written (rounded to even number)*/
	char* stringv; /*as  constant was written*/
		      /* (padded to even # chars >= 16)*/
		      /* without leading 0x*/
    } opaquev;
    struct Symbol* enumv;   /* NC_ECONST*/
} Constvalue;

typedef struct NCConstant {
    nc_type 	  nctype; /* NC_INT,... */
    int		  lineno;
    Constvalue    value;
    int           filled; /* was this originally NC_FILLVALUE? */
} NCConstant;

typedef struct Datalist {
    int           readonly; /* data field is shared with another Datalist*/
    size_t  length; /* |data| */
    size_t  alloc;  /* track total allocated space for data field*/
    NCConstant**     data; /* actual list of constants constituting the datalist*/
    /* Track various values associated with the datalist*/
    /* (used to be in Constvalue.compoundv)*/
#if 0
    struct Vlen {
        struct Symbol* schema; /* type/var that defines structure of this*/
        unsigned int count; /* # of vlen basetype instances*/
	unsigned int uid;       /* unique id for NC_VLEN*/
    } vlen;
#endif
} Datalist;

/* Define a structure to track
   location of current read point in the Datalist sequence
   In effect, we are parsing the data sequence.
   Push and pop of data sources is supported (see srcpush() below).*/
typedef struct Datasrc {
    NCConstant**   data;     /* duplicate pointer; so do not free.*/
    int index;        
    int length;
    int spliced;           /* Was this list spliced into our parent ? */
    struct Datasrc* prev; /* linked list for debugging */
} Datasrc;

#if 0
/* Define a holder for passing a start/count array */
struct Vlendata {
    char* data;
    unsigned long count;
};
extern struct Vlendata* vlendata;
#endif

extern List* alldatalists;

/* from: data.c */
extern Datalist* builddatalist(int initialize);
extern void capture(Datalist* dl);
extern void dlappend(Datalist*, NCConstant*);
extern void dlset(Datalist*, size_t, NCConstant*);
extern NCConstant* builddatasublist(Datalist* dl);
extern void dlextend(Datalist* dl);
extern void dlsetalloc(Datalist* dl, size_t newalloc);
extern Datalist* clonedatalist(Datalist* dl);
extern void reclaimalldatalists(void);
extern void reclaimdatalist(Datalist*);
extern void reclaimconstant(NCConstant*);

int       datalistline(Datalist*);
#define   datalistith(dl,i) ((dl)==NULL?NULL:((i) >= (dl)->length?NULL:(dl)->data[i]))
#define   datalistlen(dl) ((dl)==NULL?0:(dl)->length)

NCConstant* list2const(Datalist*);
Datalist* const2list(NCConstant* con);

int isstringable(nc_type nctype);

#define islistconst(con) ((con)!=NULL && (con)->nctype == NC_COMPOUND)
#define isfillconst(con) ((con)!=NULL && (con)->nctype == NC_FILLVALUE)
#define constline(con) (con==NULL?0:(con)->lineno)
#define consttype(con) (con==NULL?NC_NAT:(con)->nctype)

#define isnilconst(con) ((con)!=NULL && (con)->nctype == NC_NIL)
#define   compoundfor(con) ((con)==NULL?NULL:(con)->value.compoundv)

NCConstant* emptycompoundconst(int lineno);
NCConstant* emptystringconst(int);

NCConstant* cloneconstant(NCConstant* con); /* deep clone*/
void clearconstant(NCConstant* con); /* deep clear*/
#define freeconst(con) freeconstant(con,DEEP);
void freeconstant(NCConstant* con, int shallow);

void alignbuffer(struct NCConstant* prim, Bytebuffer* buf);

/* Code dump support procedures */
void bbindent(Bytebuffer*,const int);
void bbprintf(Bytebuffer*,const char *fmt, ...); /* append */
void bbprintf0(Bytebuffer*,const char *fmt, ...); /* clear, then append*/
/* Following dump to codebuffer */
void codeprintf(const char *fmt, ...);
void codedump(Bytebuffer*);
void codepartial(const char*);
void codeline(const char*);
void codelined(int n,const char*);
void codeflush(void); /* flush codebuffer to stdout */

void commify(Bytebuffer* buf);
char* word(char* p, Bytebuffer* buf);

/* Provide buffers for language based generators */
extern Bytebuffer* codebuffer; /* buffer over the std output */
extern Bytebuffer* stmt; /* single stmt text generation */

/* Aliases */
#define srcincr(src) srcnext(src)
#define srcget(src) srcpeek(src)

extern NCConstant* nullconst(void);
extern NCConstant nullconstant;
extern NCConstant fillconstant;
extern NCConstant nilconstant;

/* From genchar.c */
void gen_charattr(Datalist*, Bytebuffer*);
void gen_charseq(Datalist*, Bytebuffer*);
void gen_chararray(struct Dimset*, int, Datalist*, Bytebuffer*, Datalist* fillsrc);

typedef enum ListClass {
    LISTDATA, LISTATTR, LISTVLEN, LISTCOMPOUND, LISTFIELDARRAY
} ListClass;

struct Generator {
    void* globalstate; /* per-generator; per list state is in the method args where needed */
    int (*charconstant)(Generator*,struct Symbol*,Bytebuffer*,...);
    int (*constant)(Generator*,struct Symbol*,NCConstant*,Bytebuffer*,...);
    int (*listbegin)(Generator*,struct Symbol*,void*,ListClass,size_t,Bytebuffer*,int*,...);
    int (*list)(Generator*,struct Symbol*,void*,ListClass,int,size_t,Bytebuffer*,...);
    int (*listend)(Generator*,struct Symbol*,void*,ListClass,int,size_t,Bytebuffer*,...);
    int (*vlendecl)(Generator*,struct Symbol*,Bytebuffer*,int,size_t,...);
    int (*vlenstring)(Generator*,struct Symbol*,Bytebuffer*,int*,size_t*,...);
};

extern int generator_getstate(Generator*,void**);
extern int generator_reset(Generator*,void*);

typedef int (*Writer)(Generator*,struct Symbol*,Bytebuffer*,int,const size_t*,const size_t*);

extern void generate_attrdata(struct Symbol*, Generator*, Writer writer, Bytebuffer*);
extern void generate_vardata(struct Symbol*, Generator*, Writer writer,Bytebuffer*);
extern void generate_basetype(struct Symbol*,NCConstant*,Bytebuffer*,Datalist*,Generator*);


/* Obsolete */
#if 0
Datasrc* datalist2src(Datalist* list);
Datasrc* const2src(NCConstant*);
void freedatasrc(Datasrc* src);
int issublist(Datasrc* src);
int isstring(Datasrc* src);
int isfillvalue(Datasrc* src);
int istype(Datasrc* src, nc_type);
void srcpush(Datasrc*);
void srcpushlist(Datasrc* src, Datalist* cmpd);
void srcpop(Datasrc*);
void srcsetfill(Datasrc* ds, Datalist* list);
NCConstant* srcnext(Datasrc*);
int srcmore(Datasrc*);
int srcline(Datasrc* ds);
void srcreset(Datasrc* ds);
#define srclen(s) ((s)==NULL?0:(s)->length)
#ifdef FASTDATASRC
#define srcpeek(ds) ((ds)==NULL || (ds)->index >= (ds)->max?NULL:(ds)->data+(ds)->index)
#else
NCConstant* srcpeek(Datasrc*);
#endif
#endif /*0*/

#endif /*DATA_H*/


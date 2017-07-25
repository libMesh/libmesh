#ifndef NC_GENLIB_H
#define NC_GENLIB_H
/*********************************************************************
 *   Copyright 1993, UCAR/Unidata
 *   See netcdf/COPYRIGHT file for copying and redistribution conditions.
 *   $Header: /upc/share/CVS/netcdf-3/ncgen/genlib.h,v 1.20 2010/05/17 23:26:45 dmh Exp $
 *********************************************************************/
#include "config.h"
#include <stdlib.h>
#include <limits.h>
#include "generr.h"

/* break if C Line length exceeds this*/
#define C_MAX_STMT	72
/* break if FORTRAN Line length exceeds this*/
#define F77_MAX_STMT	66

#define PATHSEPARATOR   "/"

/* Convenience*/
#define REQUIRED 1
#define DONTCARE -1
#define NOTFLAT 0

#define nulllen(s) ((s)==NULL?0:strlen(s))

#define PRIMNO (NC_STRING - NC_NAT + 1)
extern struct Symbol* primsymbols[PRIMNO];

extern void derror ( const char *fmt, ... )
#ifdef _GNUC_
       __attribute__ ((format (printf, 1, 2)))
#endif
;

extern void	verror ( const char *fmt, ... )
#ifdef _GNUC_
       __attribute__ ((format (printf, 1, 2)))
#endif
;

extern void markcdf4(const char *msg);
extern char* getmarkcdf4(void);

extern void markcdf5(const char *msg);
extern char* getmarkcdf4(void);


/*
All external procedures in ncgen.h have been moved to this file.
*/

/* from: genlib.c */
extern void define_netcdf(void);/* generates all define mode stuff */
extern void close_netcdf ( void ); /* generates close */
extern char* cprefixed(List* prefix, char* suffix, char* separator);
extern void topfqn(Symbol* sym);
extern void nestedfqn(Symbol* sym);
extern void attfqn(Symbol* sym);

/* from: escapes.c */
extern int unescape(char*, const char*, int, int);
extern int unescapeoct(const char* s);
extern int unescapehex(const char* s);
extern char* cescapifychar(unsigned int c, int quote);
extern char* codify(const char *name);
extern char* escapifychar(unsigned int c, char* s0, int quote);
extern char* escapify(char*,int,size_t);
extern char* escapifyname(char* s0);
extern void cquotestring(Bytebuffer*,char quote);
extern void f77quotestring(Bytebuffer*);
extern char* xescapify(char* s0, int quote, size_t len);
extern char* jescapify(char* s0, int quote, size_t len);
extern char* jescapifyname(char* s0);
extern char* fqnescape(const char* s);

/* from: getfill.c */
extern void nc_getfill(NCConstant*);
extern char* nc_dfaltfillname(nc_type);
extern struct Datalist* getfiller(Symbol*); /* symbol isa variable|type */

/* from: ncgen.y */
extern Symbol* install(const char *sname);
extern Symbol* basetypefor(nc_type nctype);/* Convert nctype to a Symbol*/
extern Symbol* makearraytype(Symbol*, Dimset*);

/* from: cvt.c */
extern void convert1(NCConstant*,NCConstant*); /* Convert an arbitrary value to another */
extern void setprimlength(NCConstant* prim, unsigned long len);
extern struct Datalist* convertstringtochars(NCConstant* str);


/* from: semantic.c */
extern  void processsemantics(void);
extern  size_t nctypesize(nc_type);
extern  Symbol* locate(Symbol* refsym);
extern  Symbol* lookup(nc_class objectclass, Symbol* pattern);
extern  Symbol* lookupingroup(nc_class objectclass, char* name, Symbol* grp);
extern  Symbol* lookupgroup(List* prefix);
extern int nounlimited(Dimset* dimset, int from);
extern int lastunlimited(Dimset* dimset);
extern void padstring(NCConstant* con, size_t desiredlength, int fillchar);

extern Datalist* explodestrings(Datalist*,char*);
extern Datalist* implodestrings(Datalist*,char*);
extern int explodestringconst(NCConstant* con, char* tag, NCConstant*);

extern char* indented(int n);

/* Generators for cdf, c, and fortran */

#ifdef ENABLE_BINARY
/* from: genbin.c */
extern Generator* bin_generator;
extern void gen_netcdf(const char *filename);
extern void cl_netcdf(void);
#endif

#ifdef ENABLE_C
/* from: genc.c */
extern Generator* c_generator;
extern void gen_ncc(const char *filename);
extern void cl_c(void);
extern const char* ctypename(Symbol*);
extern const char* nctype(nc_type type);
extern const char* ncctype(nc_type type);
extern const char* ncstype(nc_type type);
extern const char* cname(Symbol* sym);

#endif

#ifdef ENABLE_F77
/* from: genf77.c */
extern Generator* f77_generator;
extern void gen_ncf77(const char *filename);
extern void cl_f77(void);
extern const char* f77name(Symbol*);
extern const char* f77typename(Symbol*);
#endif

#ifdef ENABLE_JAVA
/* from: genj.c */
extern Generator* j_generator;
extern void gen_ncjava(const char *filename);
extern void cl_java(void);
extern void jpartial(char*);
extern void jline(char*);
extern void jlined(int,char*);
extern void jflush(void);
#endif


/* from: main.c */
extern int k_flag;    /* -k value from command line*/
extern int format_attribute; /* 1 if format came from _FORMAT attribute */
extern int enhanced_flag; /* 1 => netcdf-4 constructs appear in the parse */
extern int cdf5_flag; /* 1 => cdf-5 unsigned types in the parse */
extern int specials_flag; /* 1 => special attributes are present */
extern int usingclassic;   /* 1 => k_flag == 1|2|5 */
extern int k_flag;
extern int ncloglevel;
extern GlobalSpecialData globalspecials;

/* Global data */

extern Symbol* symlist;      /* all symbol objects created */
extern Symbol* rootgroup;

/* Track definitions of dims, types, attributes, and vars*/
extern List* grpdefs;
extern List* dimdefs;
extern List* attdefs;
extern List* gattdefs;
extern List* xattdefs;
extern List* typdefs;
extern List* vardefs;
extern List* condefs;

extern int CDFmodel;

extern int lineno;
extern int derror_count;
extern int kflag_flag;
extern int cmode_modifier;
extern Language l_flag;
extern char* binary_ext;
extern int nofill_flag;
extern int header_only;
extern char* mainname;
extern size_t nciterbuffersize;

extern char* progname; /* for error messages*/
extern char *netcdf_name; /* command line -o file name */
extern char *datasetname; /* name from the netcdf <name> {} */
extern char *cdlname; /* name from the command line */

/* from: util.c */
extern void* emalloc (size_t);
extern void* ecalloc (size_t);
extern void* erealloc(void*,size_t);

extern const char* specialname(int tag);

#endif /*!NC_GENLIB_H*/

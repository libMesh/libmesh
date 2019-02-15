/*********************************************************************
  *   Copyright 2016, UCAR/Unidata
  *   See netcdf/COPYRIGHT file for copying and redistribution conditions.
  *********************************************************************/

/*
External functions and state variables are
defined here, including function-like #defines.
*/

#ifndef NCD4_H
#define NCD4_H 1

#include "d4includes.h"
#include "d4util.h"
#include "d4debug.h"
#include "nc4internal.h"

/**************************************************/
/* Constants */

#define RCFILEENV "DAPRCFILE"

/* Figure out a usable max path name max */
#ifdef PATH_MAX /* *nix* */
#define NC_MAX_PATH PATH_MAX
#else
#  ifdef MAX_PATH /*windows*/
#    define NC_MAX_PATH MAX_PATH
#  else
#    define NC_MAX_PATH 4096
#  endif
#endif

#define COUNTERTYPE unsigned long long
#define COUNTERSIZE sizeof(COUNTERTYPE)

/* Clear allocated memory (see e.g. d4data.c); Potentially costly*/
#undef CLEARMEM

/* Clear empty structure alignment space; Potentially costly, but probably less than  CLEARMEM */
#define CLEARSTRUCT

/* Always use fixed size opaques */
#define FIXEDOPAQUE
#define DFALTOPAQUESIZE 16

/**************************************************/

#undef nullfree
#ifndef nullfree
#define nullfree(m) ((m)==NULL?NULL:(free(m),NULL))
#endif
#define nulldup(s) ((s)==NULL?NULL:strdup(s))

/**************************************************/
/* DSP API wrappers */

#ifdef FIX
extern int dsp_getDMR(ND4dsp* dsp, DCR** dcrp);
extern int dsp_getDAP(ND4dsp* dsp, DCR** dcrp);
extern int dsp_close(ND4dsp* dsp);

/* DSP API */
extern int dsp_open(const char* path, ND4dsp** dspp);

/*
extern NCerror dapbuildvaraprojection(CDFnode*,
		     const size_t* startp, const size_t* countp, const ptrdiff_t* stridep,
		     struct DCEprojection** projectionlist);

extern NCerror NCD4_getvarx(int ncid, int varid,
	    const size_t *startp,
	    const size_t *countp,
	    const ptrdiff_t *stridep,
	    void *data,
	    nc_type dsttype0);
*/
#endif

/**************************************************/

/* From d4http.c */
extern long NCD4_fetchhttpcode(CURL* curl);
extern int NCD4_fetchurl_file(CURL* curl, const char* url, FILE* stream, d4size_t* sizep, long* filetime);
extern int NCD4_fetchurl(CURL* curl, const char* url, NCbytes* buf, long* filetime);
extern int NCD4_curlopen(CURL** curlp);
extern void NCD4_curlclose(CURL* curl);
extern int NCD4_fetchlastmodified(CURL* curl, char* url, long* filetime);
extern int NCD4_ping(const char* url);

/* From d4read.c */
extern int NCD4_readDMR(NCD4INFO* state);
extern int NCD4_readDAP(NCD4INFO* state, int flags);

/* From d4parser.c */
extern int NCD4_parse(NCD4meta*);
extern NCD4node* NCD4_findAttr(NCD4node* container, const char* attrname);
extern NCD4node* NCD4_groupFor(NCD4node* node);

/* From d4printer.c */
extern int NCD4_print(NCD4meta*, NCbytes* output);

/* From d4meta.c */
extern NCD4meta* NCD4_newmeta(size_t size, void* rawdata);
extern void NCD4_reclaimMeta(NCD4meta*);
extern void reclaimNode(NCD4node* node);
extern void NCD4_setdebuglevel(NCD4meta*,int);
extern int NCD4_metabuild(NCD4meta*, int ncid);
extern size_t NCD4_computeTypeSize(NCD4meta*, NCD4node* type);

/* From d4chunk.c */
extern int NCD4_dechunk(NCD4meta*);
extern int NCD4_infermode(NCD4meta* meta);

/* From d4swap.c */
extern int NCD4_swapdata(NCD4meta*, NClist* topvars);

/* From d4fix.c */
extern int NCD4_delimit(NCD4meta*, NCD4node* var, void** offsetp);
extern int NCD4_moveto(NCD4meta*, NCD4node* var, d4size_t count, void** offsetp);
extern int NCD4_toposort(NCD4meta*);

/* From d4data.c */
extern int NCD4_processdata(NCD4meta*);
extern int NCD4_fillinstance(NCD4meta*, NCD4node* type, void** offsetp, void** dstp, NClist* blobs);
extern int NCD4_getToplevelVars(NCD4meta* meta, NCD4node* group, NClist* toplevel);

/* From d4util.c */
extern d4size_t NCD4_dimproduct(NCD4node* node);
extern size_t NCD4_typesize(nc_type tid);
extern int NCD4_isLittleEndian(void);/* Return 1 if this machine is little endian */
extern int NCD4_errorNC(int code, const int line, const char* file);
extern int NCD4_error(int code, const int line, const char* file, const char* fmt, ...);
extern char* NCD4_makeFQN(NCD4node* node);
extern char* NCD4_makeName(NCD4node*,const char* sep);
extern int NCD4_parseFQN(const char* fqn0, NClist* pieces);
extern char* NCD4_deescape(const char* esc);
extern char* NCD4_entityescape(const char* s);

/* From d4dump.c */
extern void NCD4_dumpbytes(size_t size, const void* data0, int swap);
extern void NCD4_tagdump(size_t size, const void* data0, int swap, const char* tag);
extern void NCD4_dumpvars(NCD4node* group);
extern union ATOMICS* NCD4_dumpatomic(NCD4node* var, void* data);

/* From d4rc.c */
extern int NCD4_rcload(void);
extern int NCD4_rcprocess(NCD4INFO* info);
extern void NCD4_rcfree(NClist* rc);
extern char* NCD4_rclookup(char* key, char* hostport);
extern int NCD4_parseproxy(NCD4INFO* info, const char* surl);
extern int NCD4_rcdefault(NCD4INFO*);

/* From d4cvt.c */
extern int NCD4_convert(nc_type srctype, nc_type dsttype, char* memory0, char* value0, size_t count);

/* From d4crc32.c */
extern unsigned int NCD4_crc32(unsigned int crc, const void *buf, size_t size);

/* d4file.c */
extern void NCD4_applyclientparamcontrols(NCD4INFO*);

/* Add an extra function whose sole purpose is to allow
   configure(.ac) to test for the presence of this code.
*/
extern int nc__dap4(void);

/**************************************************/
/* Macro defined functions */

#undef NCCHECK
#undef FAIL
#define NCCHECK(expr) if((ret=(expr))) {ret = NCD4_errorNC(ret,__LINE__,__FILE__); goto done;}else{}
#define FAIL(code,fmt,...) do{ret=NCD4_error(code,__LINE__,__FILE__,fmt , ##__VA_ARGS__); goto done;}while(0)

#undef INCR
#undef DECR
#undef DELTA
#define INCR(offset,size) ((void*)(((char*)(offset))+(size)))
#define DECR(offset,size) ((void*)(((char*)(offset))-(size)))
#define DELTA(p1,p2) ((ptrdiff_t)(((char*)(p1))-((char*)(p2))))

#undef GETCOUNTER
#undef SKIPCOUNTER
#define GETCOUNTER(p) ((d4size_t)*((COUNTERTYPE*)(p)))
#define SKIPCOUNTER(p) {p=INCR(p,COUNTERSIZE);}

#undef PUSH
#define PUSH(list,value) do{if((list)==NULL) {(list)=nclistnew();} else{}; nclistpush((list),(value));}while(0)

#define getnc3id(d4) (getdap(d4)->nc4id)

#define ISTOPLEVEL(var) ((var)->container == NULL || (var)->container->sort == NCD4_GROUP)

#define FILEIDPART(NCID) (((unsigned int) (NCID)) >> ID_SHIFT)
#define GROUPIDPART(NCID) (((unsigned int) (NCID)) & GRP_ID_MASK)
#define MAKENCID(grp,file) ((((unsigned int)(file)) << ID_SHIFT) | (grp))

#define getdap(ncp) ((NCD4INFO*)((NC*)ncp)->dispatchdata)
#define getnc4id(ncp) (getdap(ncp)->substrate.nc4id)

/* Convert a dap4 grpid to a substrate id */
#define makenc4id(ncp,dap4id) (((dap4id) & GRP_ID_MASK) | getdap(ncp)->substrate.nc4id)
/* and the inverse */
#define makedap4id(ncp,nc4id) (((nc4id) & GRP_ID_MASK) | (ncp)->ext_ncid)

#ifdef CLEARMEM
#define d4alloc(n) (calloc(1,(size_t)(n)))
#else
#define d4alloc(n) (malloc((size_t)(n)))
#endif

/* A number of hacks have been inserted
   to deal with issues in accessing hyrax
   using DAP4.
*/
#undef HYRAXHACK


#endif /*NCD4_H*/


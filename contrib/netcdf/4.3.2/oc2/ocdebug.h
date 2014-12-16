/* Copyright 2009, UCAR/Unidata and OPeNDAP, Inc.
   See the COPYRIGHT file for more information. */

#ifndef OCOCDBG_H
#define OCOCDBG_H

#ifdef HAVE_STDARG_H
#include <stdarg.h>
#endif

#ifndef OCDEBUG
#undef OCDEBUG
#endif

/* OCCATCHERROR is used to detect errors as close
   to their point of origin as possible. When
   enabled, one can set a breakpoint in ocbreakpoint()
   to catch the failure. Turing it on incurs a significant
   performance penalty, so it is off by default.*/

#define OCCATCHERROR

#define OCPANIC(msg) assert(ocpanic(msg))
#define OCPANIC1(msg,arg) assert(ocpanic(msg,arg))
#define OCPANIC2(msg,arg1,arg2) assert(ocpanic(msg,arg1,arg2))

/* Make it possible to catch assertion failures by breakpointing ocpanic*/
#define OCASSERT(expr) if(!(expr)) {OCPANIC((#expr));} else {}

/* Need some syntactic trickery to make these macros work*/
#ifdef OCDEBUG
#define OCDBG(l,msg) {oclog(OCLOGDBG,msg);}
#define OCDBG1(l,msg,arg) {oclog(OCLOGDBG,msg,arg);}
#define OCDBG2(l,msg,arg1,arg2) {oclog(OCLOGDBG,msg,arg1,arg2);}
#define OCDBGTEXT(l,text) {oclogtext(OCLOGNOTE,text);} else {}
#define OCDBGCODE(l,code) {code;}

#else
#define OCDBG(l,msg)
#define OCDBG1(l,msg,arg)
#define OCDBG2(l,msg,arg1,arg2)
#define OCDBGTEXT(l,text)
#define OCDBGCODE(l,code)
#endif


/*
OCPROGRESS attempts to provide some info
about how IO is getting along.
*/
#undef OCPROGRESS

extern int ocdebug;
extern int cedebug;

/*extern char* dent2(int n);*/
/*/extern char* dent(int n);*/
extern int ocpanic(const char* fmt, ...);

extern int xdrerror(void);

/*
Provide wrapped versions of calloc and malloc.
The wrapped version panics if memory
is exhausted.  It also guarantees that the
memory has been zero'd.
*/

extern void* occalloc(size_t size, size_t nelems);
extern void* ocmalloc(size_t size);
extern void  ocfree(void*);

#define MEMCHECK(var,throw) {if((var)==NULL) return (throw);}
#define MEMFAIL(var) MEMCHECK(var,OCTHROW(OC_ENOMEM))
#define MEMGOTO(var,label) {if((var)==NULL) goto label;}

#ifdef OCCATCHERROR
extern int ocbreakpoint(int err);
extern int octhrow(int err);
/* Place breakpoint on ocbreakpoint to catch errors close to where they occur*/
#define OCTHROW(e) octhrow(e)
#define OCTHROWCHK(e) (void)octhrow(e)
#define OCGOTO(label) {ocbreakpoint(-1); goto label;}
#else
#define OCTHROW(e) (e)
#define OCTHROWCHK(e)
#define OCGOTO(label) goto label
#endif


#endif /*OCOCDBG_H*/


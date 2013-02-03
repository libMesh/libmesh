/* Copyright 2009, UCAR/Unidata and OPeNDAP, Inc.
   See the COPYRIGHT file for more information. */

#ifndef OCUTIL_H
#define OCUTIL_H 1

/* Forward */
struct OCstate;

#define ocmax(x,y) ((x) > (y) ? (x) : (y))

extern char* ocstrndup(const char* s, size_t len);
extern int ocstrncmp(const char* s1, const char* s2, size_t len);

extern size_t octypesize(OCtype etype);
extern char*  octypetostring(OCtype octype);
extern char*  octypetoddsstring(OCtype octype);
extern char* ocerrstring(int err);
extern OCerror ocsvcerrordata(struct OCstate*,char**,char**,long*);
extern OCerror octypeprint(OCtype etype, char* buf, size_t bufsize, void* value);
extern size_t xxdrsize(OCtype etype);

extern size_t totaldimsize(OCnode*);

extern void makedimlist(OClist* path, OClist* dims);

extern int findbod(OCbytes* buffer, size_t*, size_t*);

/* Reclaimers*/
extern void freeOCnode(OCnode*,int);
extern void ocfreeprojectionclause(OCprojectionclause* clause);

/* Misc. */
extern void ocdataddsmsg(struct OCstate*, struct OCtree*);

#endif /*UTIL_H*/

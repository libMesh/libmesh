/* Copyright 2009, UCAR/Unidata and OPeNDAP, Inc.
   See the COPYRIGHT file for more information. */

#ifndef OCBYTES_H
#define OCBYTES_H 1

typedef struct OCbytes {
  int nonextendible; /* 1 => fail if an attempt is made to extend this buffer*/
  unsigned long alloc;
  unsigned long length;
  char* content;
} OCbytes;

#if defined(_CPLUSPLUS_) || defined(__CPLUSPLUS__) || defined(__CPLUSPLUS)
#define EXTERNC extern "C"
#else
#define EXTERNC extern
#endif

EXTERNC OCbytes* ocbytesnew(void);
EXTERNC void ocbytesfree(OCbytes*);
EXTERNC int ocbytessetalloc(OCbytes*,unsigned long);
EXTERNC int ocbytessetlength(OCbytes*,unsigned long);
EXTERNC int ocbytesfill(OCbytes*, char fill);

/* Produce a duplicate of the contents*/
EXTERNC char* ocbytesdup(OCbytes*);
/* Extract the contents and leave buffer empty */
EXTERNC char* ocbytesextract(OCbytes*);

/* Return the ith byte; -1 if no such index */
EXTERNC int ocbytesget(OCbytes*,unsigned long);
/* Set the ith byte */
EXTERNC int ocbytesset(OCbytes*,unsigned long,char);

/* Append one byte */
EXTERNC int ocbytesappend(OCbytes*,char); /* Add at Tail */
/* Append n bytes */
EXTERNC int ocbytesappendn(OCbytes*,const void*,unsigned long); /* Add at Tail */

/* Null terminate the byte string without extending its length (for debugging) */
EXTERNC int ocbytesnull(OCbytes*);

/* Concatenate a null-terminated string to the end of the buffer */
EXTERNC int ocbytescat(OCbytes*,const char*);

/* Set the contents of the buffer; mark the buffer as non-extendible */
EXTERNC int ocbytessetcontents(OCbytes*, char*, unsigned long);

/* Following are always "in-lined"*/
#define ocbyteslength(bb) ((bb)!=NULL?(bb)->length:0)
#define ocbytesalloc(bb) ((bb)!=NULL?(bb)->alloc:0)
#define ocbytescontents(bb) (((bb)!=NULL && (bb)->content!=NULL)?(bb)->content:(char*)"")
#define ocbytesextend(bb,len) ocbytessetalloc((bb),(len)+(bb->alloc))
#define ocbytesclear(bb) ((bb)!=NULL?(bb)->length=0:0)
#define ocbytesavail(bb,n) ((bb)!=NULL?((bb)->alloc - (bb)->length) >= (n):0)

#endif /*OCBYTES_H*/

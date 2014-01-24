/* Copyright 2009, UCAR/Unidata and OPeNDAP, Inc.
   See the COPYRIGHT file for more information. */
#ifndef OCLIST_H
#define OCLIST_H 1

/* Define the type of the elements in the list*/

#if defined(_CPLUSPLUS_) || defined(__CPLUSPLUS__)
#define EXTERNC extern "C"
#else
#define EXTERNC extern
#endif

EXTERNC int oclistnull(void*);

typedef struct OClist {
  unsigned long alloc;
  unsigned long length;
  void** content;
} OClist;

EXTERNC OClist* oclistnew(void);
EXTERNC int oclistfree(OClist*);
EXTERNC int oclistsetalloc(OClist*,unsigned long);
EXTERNC int oclistsetlength(OClist*,unsigned long);

/* Set the ith element */
EXTERNC int oclistset(OClist*,unsigned long,void*);
/* Get value at position i */
EXTERNC void* oclistget(OClist*,unsigned long);/* Return the ith element of l */
/* Insert at position i; will push up elements i..|seq|. */
EXTERNC int oclistinsert(OClist*,unsigned long,void*);
/* Remove element at position i; will move higher elements down */
EXTERNC void* oclistremove(OClist* l, unsigned long i);

/* Tail operations */
EXTERNC int oclistpush(OClist*,void*); /* Add at Tail */
EXTERNC void* oclistpop(OClist*);
EXTERNC void* oclisttop(OClist*);

/* Duplicate and return the content (null terminate) */
EXTERNC void** oclistdup(OClist*);

/* Look for value match */
EXTERNC int oclistcontains(OClist*, void*);

/* Remove element by value; only removes first encountered */
EXTERNC int oclistelemremove(OClist* l, void* elem);

/* remove duplicates */
EXTERNC int oclistunique(OClist*);

/* Create a clone of a list */
EXTERNC OClist* oclistclone(OClist*);

/* Following are always "in-lined"*/
#define oclistclear(l) oclistsetlength((l),0)
#define oclistextend(l,len) oclistsetalloc((l),(len)+(l->alloc))
#define oclistcontents(l)  ((l)==NULL?NULL:(l)->content)
#define oclistlength(l)  ((l)==NULL?0:(l)->length)

#endif /*OCLIST_H*/

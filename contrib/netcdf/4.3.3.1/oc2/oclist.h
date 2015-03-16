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
  size_t alloc;
  size_t length;
  void** content;
} OClist;

EXTERNC OClist* oclistnew(void);
EXTERNC int oclistfree(OClist*);
EXTERNC int oclistsetalloc(OClist*,size_t);
EXTERNC int oclistsetlength(OClist*,size_t);

/* Set the ith element */
EXTERNC int oclistset(OClist*,size_t,void*);
/* Get value at position i */
EXTERNC void* oclistget(OClist*,size_t);/* Return the ith element of l */
/* Insert at position i; will push up elements i..|seq|. */
EXTERNC int oclistinsert(OClist*,size_t,void*);
/* Remove element at position i; will move higher elements down */
EXTERNC void* oclistremove(OClist* l, size_t i);

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

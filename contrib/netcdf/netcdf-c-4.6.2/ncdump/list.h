/* Copyright 2009, UCAR/Unidata and OPeNDAP, Inc.
   See the COPYRIGHT file for more information. */
#ifndef LIST_H
#define LIST_H 1

/* Define the type of the elements in the list*/

#if defined(_CPLUSPLUS_) || defined(__CPLUSPLUS__)
#define EXTERNC extern "C"
#else
#define EXTERNC extern
#endif

EXTERNC int listnull(void*);

typedef struct List {
  unsigned long alloc;
  unsigned long length;
  void** content;
} List;

EXTERNC List* listnew(void);
EXTERNC int listfree(List*);
EXTERNC int listsetalloc(List*,unsigned long);
EXTERNC int listsetlength(List*,unsigned long);

/* Set the ith element */
EXTERNC int listset(List*,unsigned long,void*);
/* Get value at position i */
EXTERNC void* listget(List*,unsigned long);/* Return the ith element of l */
/* Insert at position i; will push up elements i..|seq|. */
EXTERNC int listinsert(List*,unsigned long,void*);
/* Remove element at position i; will move higher elements down */
EXTERNC void* listremove(List* l, unsigned long i);

/* Tail operations */
EXTERNC int listpush(List*,void*); /* Add at Tail */
EXTERNC void* listpop(List*);
EXTERNC void* listtop(List*);

/* Duplicate and return the content (null terminate) */
EXTERNC void** listdup(List*);

/* Look for value match */
EXTERNC int listcontains(List*, void*);

/* Remove element by value; only removes first encountered */
EXTERNC int listelemremove(List* l, void* elem);

/* remove duplicates */
EXTERNC int listunique(List*);

/* Create a clone of a list */
EXTERNC List* listclone(List*);

/* Following are always "in-lined"*/
#define listclear(l) listsetlength((l),0)
#define listextend(l,len) listsetalloc((l),(len)+(l->alloc))
#define listcontents(l)  ((l)==NULL?NULL:(l)->content)
#define listlength(l)  ((l)==NULL?0:(int)(l)->length)

#endif /*LIST_H*/

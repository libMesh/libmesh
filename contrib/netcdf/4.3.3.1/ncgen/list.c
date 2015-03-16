/* Copyright 2009, UCAR/Unidata and OPeNDAP, Inc.
   See the COPYRIGHT file for more information. */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "list.h"

int listnull(void* e) {return e == NULL;}

#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif

#define DEFAULTALLOC 16
#define ALLOCINCR 16

List* listnew(void)
{
  List* l;
/*
  if(!initialized) {
    memset((void*)&DATANULL,0,sizeof(void*));
    initialized = 1;
  }
*/
  l = (List*)malloc(sizeof(List));
  if(l) {
    l->alloc=0;
    l->length=0;
    l->content=NULL;
  }
  return l;
}

int
listfree(List* l)
{
  if(l) {
    l->alloc = 0;
    if(l->content != NULL) {free(l->content); l->content = NULL;}
    free(l);
  }
  return TRUE;
}

int
listsetalloc(List* l, unsigned long sz)
{
  void** newcontent = NULL;
  if(l == NULL) return FALSE;
  if(sz <= 0) {sz = (l->length?2*l->length:DEFAULTALLOC);}
  if(l->alloc >= sz) {return TRUE;}
  newcontent=(void**)calloc(sz,sizeof(void*));
  if(newcontent != NULL && l->alloc > 0 && l->length > 0 && l->content != NULL) {
    memcpy((void*)newcontent,(void*)l->content,sizeof(void*)*l->length);
  }
  if(l->content != NULL) free(l->content);
  l->content=newcontent;
  l->alloc=sz;
  return TRUE;
}

int
listsetlength(List* l, unsigned long sz)
{
  if(l == NULL) return FALSE;
  if(sz > l->alloc && !listsetalloc(l,sz)) return FALSE;
  l->length = sz;
  return TRUE;
}

void*
listget(List* l, unsigned long index)
{
  if(l == NULL || l->length == 0) return NULL;
  if(index >= l->length) return NULL;
  return l->content[index];
}

int
listset(List* l, unsigned long index, void* elem)
{
  if(l == NULL) return FALSE;
  if(index >= l->length) return FALSE;
  l->content[index] = elem;
  return TRUE;
}

/* Insert at position i of l; will push up elements i..|seq|. */
int
listinsert(List* l, unsigned long index, void* elem)
{
  int i; /* do not make unsigned */
  if(l == NULL) return FALSE;
  if(index > l->length) return FALSE;
  listsetalloc(l,0);
  for(i=(int)l->length;i>index;i--) l->content[i] = l->content[i-1];
  l->content[index] = elem;
  l->length++;
  return TRUE;
}

int
listpush(List* l, void* elem)
{
  if(l == NULL) return FALSE;
  if(l->length >= l->alloc) listsetalloc(l,0);
  l->content[l->length] = elem;
  l->length++;
  return TRUE;
}

void*
listpop(List* l)
{
  if(l == NULL || l->length == 0) return NULL;
  l->length--;  
  return l->content[l->length];
}

void*
listtop(List* l)
{
  if(l == NULL || l->length == 0) return NULL;
  return l->content[l->length - 1];
}

void*
listremove(List* l, unsigned long i)
{
  unsigned long len;
  void* elem;
  if(l == NULL || (len=l->length) == 0) return NULL;
  if(i >= len) return NULL;
  elem = l->content[i];
  for(i+=1;i<len;i++) l->content[i-1] = l->content[i];
  l->length--;
  return elem;  
}

/* Duplicate and return the content (null terminate) */
void**
listdup(List* l)
{
    void** result = (void**)malloc(sizeof(void*)*(l->length+1));
    memcpy((void*)result,(void*)l->content,sizeof(void*)*l->length);
    result[l->length] = (void*)0;
    return result;
}

int
listcontains(List* l, void* elem)
{
    unsigned long i;
    for(i=0;i<listlength(l);i++) {
	if(elem == listget(l,i)) return 1;
    }
    return 0;
}

/* Remove element by value; only removes first encountered */
int
listelemremove(List* l, void* elem)
{
  unsigned long len;
  unsigned long i;
  int found = 0;
  if(l == NULL || (len=l->length) == 0) return 0;
  for(i=0;i<listlength(l);i++) {
    void* candidate = l->content[i];
    if(elem == candidate) {
      for(i+=1;i<len;i++) l->content[i-1] = l->content[i];
      l->length--;
      found = 1;
      break;
    }
  }
  return found;
}




/* Extends list to include a unique operator 
   which remove duplicate values; NULL values removed
   return value is always 1.
*/

int
listunique(List* l)
{
    unsigned long i,j,k,len;
    void** content;
    if(l == NULL || l->length == 0) return 1;
    len = l->length;
    content = l->content;
    for(i=0;i<len;i++) {
        for(j=i+1;j<len;j++) {
	    if(content[i] == content[j]) {
		/* compress out jth element */
                for(k=j+1;k<len;k++) content[k-1] = content[k];	
		len--;
	    }
	}
    }
    l->length = len;
    return 1;
}

List*
listclone(List* l)
{
    List* clone = listnew();
    *clone = *l;
    clone->content = listdup(l);
    return clone;
}

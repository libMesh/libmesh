/*********************************************************************
  *   Copyright 1993, UCAR/Unidata
  *   See netcdf/COPYRIGHT file for copying and redistribution conditions.
  *********************************************************************/
/* $Header$ */

#ifndef DCECONSTRAINTS_H
#define DCECONSTRAINTS_H 1

#include "ceconstraints.h"

/* Provide a universal cast type containing common fields */

/* Define the common "supertype */
typedef struct DCEnode {
    CEsort sort;    
} DCEnode;

/* The slice structure is assumed common to all uses */
typedef struct DCEslice {
    DCEnode node;
    size_t first;
    size_t count;
    size_t length; /* count*stride */
    size_t stride;
    size_t stop; /* == first + count*/
    size_t declsize;  /* from defining dimension, if any.*/
} DCEslice;

typedef struct DCEsegment {
    DCEnode node;
    char* name; 
    int slicesdefined; /*1=>slice counts defined, except declsize*/
    int slicesdeclized; /*1=>slice declsize defined */
    size_t rank;
    DCEslice slices[NC_MAX_VAR_DIMS];    
    void* annotation;
} DCEsegment;

typedef struct DCEfcn {
    DCEnode node;
    char* name; 
    NClist* args;
} DCEfcn;

typedef struct DCEvar {
    DCEnode node;
    NClist* segments;
    void* annotation;
} DCEvar;

typedef struct DCEconstant {
    DCEnode node;
    CEsort discrim;
    char* text;
    long long intvalue;
    double floatvalue;
} DCEconstant;

typedef struct DCEvalue {
    DCEnode node;
    CEsort discrim;
    /* Do not bother with a union */
    DCEconstant* constant;
    DCEvar* var;
    DCEfcn* fcn;
} DCEvalue;

typedef struct DCEselection {
    DCEnode node;
    CEsort operator;
    DCEvalue* lhs;
    NClist* rhs;
} DCEselection;

typedef struct DCEprojection {
    DCEnode node;
    CEsort discrim;
    /* Do not bother with a union */
    DCEvar* var;
    DCEfcn* fcn;
} DCEprojection;

typedef struct DCEconstraint {
    DCEnode node;
    NClist* projections;
    NClist* selections;
} DCEconstraint;


extern int dceparseconstraints(char* constraints, DCEconstraint* DCEonstraint);
extern int dceslicemerge(DCEslice* dst, DCEslice* src);
extern int dcemergeprojectionlists(NClist* dst, NClist* src);

extern DCEnode* dceclone(DCEnode* node);
extern NClist* dceclonelist(NClist* list);

extern void dcefree(DCEnode* node);
extern void dcefreelist(NClist* list);

extern char* dcetostring(DCEnode* node);
extern char* dcelisttostring(NClist* list,char*);
extern void dcetobuffer(DCEnode* node, NCbytes* buf);
extern void dcelisttobuffer(NClist* list, NCbytes* buf,char*);

extern NClist* dceallnodes(DCEnode* node, CEsort which);

extern DCEnode* dcecreate(CEsort sort);

extern void dcemakewholeslice(DCEslice* slice, size_t declsize);
extern void dcemakewholeprojection(DCEprojection*);

extern int dceiswholesegment(DCEsegment*);
extern int dceiswholeslice(DCEslice*);
extern int dceiswholeseglist(NClist*);
extern int dceiswholeprojection(DCEprojection*);
extern int dcesamepath(NClist* list1, NClist* list2);
extern int dcemergeprojections(DCEprojection* dst, DCEprojection* src);

extern int dceverbose;

#endif /*DCECONSTRAINTS_H*/

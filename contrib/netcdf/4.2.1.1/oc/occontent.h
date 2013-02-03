/* Copyright 2009, UCAR/Unidata and OPeNDAP, Inc.
   See the COPYRIGHT file for more information. */

#ifndef OCCONTENT_H
#define OCCONTENT_H

/*! Specifies the OCcontent*/
typedef struct OCcontent {
    unsigned int magic;
    OCmode mode;
    struct OCstate* state;
    struct OCnode* node;
    struct OCtree* tree;
    int    packed; /* track OC_Char and OC_Byte specially*/
    struct OCCACHE {
	int valid;
        ocindex_t index;      /* index corresponding to offset */
        ocindex_t maxindex;   /* max allowable index, if known0 => max unknown */
        ocoffset_t offset;    /* location of this content in the xdr data */
    } cache;  /* track last xdr position and index of this content */
    struct OCcontent* next;
} OCcontent;

extern OCcontent* ocnewcontent(OCstate* state);
extern void ocfreecontent(OCstate* state, OCcontent* content);
extern OCmode ocgetmode(OCcontent* content);

extern int ocdataith(struct OCstate*, OCcontent*, size_t, OCcontent*);
extern int ocdatacount(struct OCstate*, OCcontent*, size_t*);

extern int ocrootdata(struct OCstate*, struct OCnode*, struct OCcontent*);
extern int ocgetcontent(struct OCstate*, struct OCcontent*, void* memory,
                        size_t memsize, size_t start, size_t count);

#endif /*OCCONTENT_H*/

/* Copyright 2009, UCAR/Unidata and OPeNDAP, Inc.
   See the COPYRIGHT file for more information. */

#ifndef OCDATA_H
#define OCDATA_H

typedef struct OCdimcounter {
    int rank;
    size_t index[OC_MAX_DIMS];
    size_t size[OC_MAX_DIMS];
} OCdimcounter;

extern const char StartOfSequence;
extern const char EndOfSequence;

/*Forward */
struct OCcontent;

/* Skip arbitrary dimensioned instance; Handles dimensioning.*/
extern int ocskip(OCnode* node, XXDR* xdrs);

extern int occountrecords(OCnode* node, XXDR* xdrs, size_t* nrecordsp);

extern int ocxdrread(struct OCcontent*, XXDR*, char* memory, size_t, ocindex_t index, ocindex_t count);

extern int ocskipinstance(OCnode* node, XXDR* xdrs, int state, int* tagp);

#endif /*OCDATA_H*/

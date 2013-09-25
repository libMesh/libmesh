/* Copyright 2009, UCAR/Unidata and OPeNDAP, Inc.
   See the COPYRIGHT file for more information. */

#ifndef READ_H
#define READ_H


extern int readDDS(OCstate*, OCtree*);
extern int readDAS(OCstate*, OCtree*);

extern int readDATADDS(OCstate*, OCtree*, int inmemory);

#endif /*READ_H*/

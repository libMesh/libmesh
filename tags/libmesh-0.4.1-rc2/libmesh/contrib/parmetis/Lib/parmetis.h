/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * par_metis.h
 *
 * This file includes all necessary header files
 *
 * Started 8/27/94
 * George
 *
 * $Id: parmetis.h,v 1.2 2003-07-15 12:40:11 benkirk Exp $
 */

#define PARMETIS_MAJOR_VERSION        3
#define PARMETIS_MINOR_VERSION        0 

/*
#define DEBUG			1
#define DMALLOC			1
*/

#include <stdheaders.h>

#ifdef DMALLOC
#include <dmalloc.h>
#endif

#include <rename.h>
#include <defs.h>
#include <struct.h>
#include <macros.h>
#include <proto.h>

/**
 * Hmm...  These are copied from proto.h, but for some reason they
 * don't get included?
 */

/* kmetis.c */
void ParMETIS_V3_PartKway(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, int *, float *, float *, int *, int *, idxtype *, MPI_Comm *);

/* ametis.c */
void ParMETIS_V3_AdaptiveRepart(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, int *, float *, float *, float *, int *, int *, idxtype *, MPI_Comm *);

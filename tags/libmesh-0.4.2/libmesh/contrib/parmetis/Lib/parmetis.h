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
 * $Id: parmetis.h,v 1.3 2004-02-10 13:28:06 benkirk Exp $
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

#include "rename.h"
#include "defs.h"
#include "struct.h"
#include "macros.h"
#include "proto.h"

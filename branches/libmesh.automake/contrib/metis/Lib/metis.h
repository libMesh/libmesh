/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * metis.h
 *
 * This file includes all necessary header files
 *
 * Started 8/27/94
 * George
 *
 * $Id$
 */

/*
#define	DEBUG		1
#define	DMALLOC		1
*/
#ifdef LIBMESH_IS_COMPILING_METIS
#  if defined(__GNUC__) && !defined(__INTEL_COMPILER) /* intel can masquerade as GNUC... */
#    pragma GCC diagnostic ignored "-Wimplicit-function-declaration"
#  endif
#endif

#include "./stdheaders.h"

#ifdef DMALLOC
#include <dmalloc.h>
#endif

#include "./defs.h"
#include "./struct.h"
#include "./macros.h"
#include "./rename.h"
#include "./proto.h"


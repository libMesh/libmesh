/*
 *	Copyright 1996, University Corporation for Atmospheric Research
 *      See netcdf/COPYRIGHT file for copying and redistribution conditions.
 */
/* $Id: rnd.h 2501 2007-11-20 02:33:29Z benkirk $ */
#ifndef _RNDUP

/* useful for aligning memory */
#define	_RNDUP(x, unit)  ((((x) + (unit) - 1) / (unit)) \
	* (unit))
#define	_RNDDOWN(x, unit)  ((x) - ((x)%(unit)))

#define M_RND_UNIT	(sizeof(double))
#define	M_RNDUP(x) _RNDUP(x, M_RND_UNIT)
#define	M_RNDDOWN(x)  __RNDDOWN(x, M_RND_UNIT)

#endif

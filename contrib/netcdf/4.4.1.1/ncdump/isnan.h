/*********************************************************************
 *   Copyright 2008, University Corporation for Atmospheric Research
 *   See netcdf/README file for copying and redistribution conditions.
 *   $Id: isnan.h,v 1.3 2008/04/23 22:05:00 russ Exp $
 *********************************************************************/

#include <config.h>

#ifndef NO_FLOAT_H
#include <float.h>		/* for DBL_MAX */
#endif /* NO_FLOAT_H */
#include <math.h>

#if ! (defined(isinf) || HAVE_DECL_ISINF)
#define isinf(x) (DBL_MAX/((double)(x))==0.0)
#endif /* !HAVE_DECL_ISINF */
#if ! (defined(isnan) || HAVE_DECL_ISNAN)
#define isnan(x) ((x)!=(x))
#endif /* !HAVE_DECL_ISNAN */
#if ! (defined(isfinite) || HAVE_DECL_ISFINITE)
#define isfinite(x) (!(isinf(x)||isnan(x)))
#endif /* !HAVE_DECL_ISFINITE */

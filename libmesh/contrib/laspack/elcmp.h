/****************************************************************************/
/*                                 elcmp.h                                  */
/****************************************************************************/
/*                                                                          */
/* type of matrix ELements and vector CoMPonents                            */
/*                                                                          */
/* Copyright (C) 1992-1996 Tomas Skalicky. All rights reserved.             */
/*                                                                          */
/****************************************************************************/
/*                                                                          */
/*        ANY USE OF THIS CODE CONSTITUTES ACCEPTANCE OF THE TERMS          */
/*              OF THE COPYRIGHT NOTICE (SEE FILE COPYRGHT.H)               */
/*                                                                          */
/****************************************************************************/

#ifndef ELCMP_H
#define ELCMP_H

#include <float.h>
#include <math.h>

#include "copyrght.h"

typedef double Real;

#ifdef __BORLANDC__
/* BC 2.0 does not handle IEEE arithmetic correctly */
#define IsZero(a) (fabs(a) < 1.0e20 * DBL_MIN)
#define IsOne(a)  (fabs(a - 1.0) < 10.0 * DBL_EPSILON)
#else
#define IsZero(a) (fabs(a) < 10.0 * DBL_MIN)
#define IsOne(a)  (fabs(a - 1.0) < 10.0 * DBL_EPSILON)
#endif /* __BORLANDC__ */

typedef struct {
    size_t Pos;
    Real Val;
} ElType;

#endif /* ELCMP_H */

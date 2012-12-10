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

/* configurations specific for use with libMesh */
#include "laspack_config.h"


#ifdef _LP_INCLUDED_FROM_CPLUSPLUS
extern "C" {
#endif

  
#include <float.h>
#include <math.h>

#include "copyrght.h"


#if defined(_LP_LIBMESH_USE_COMPLEX_NUMBERS)

   /* complex arithmetic, hua */
#  include<complex>
   typedef std::complex<double> _LPNumber;
   typedef std::complex<double> _LPDouble;
   typedef double               _LPReal;

#  define _LPfabs(c) std::abs(c)
   /* when comparing two complex values, use the real-signed absolute (possibly wrong?) */
#  define _LPmax(x,y) (   ((x).real()/fabs((x).real())*std::abs(x))  \
                        > ((y).real()/fabs((y).real())*std::abs(y))  \
                        ? (x) : (y) )
#  define _LPmin(x,y) (   ((x).real()/fabs((x).real())*std::abs(x))  \
                        < ((y).real()/fabs((y).real())*std::abs(y))  \
                        ? (x) : (y) )
#  define _LPPrintFormat(c) std::abs(c)

   /* To be consistent with PETSc, for complex we only compare the real part.
    * See also: $PETSC_DIR/vec//impls/seq/dvec2.c */
   /* _LPAbsRealPart() returns a complex with fabs(real part) */
#  define _LPAbsRealPart(x) (_LPNumber(fabs((x).real()), (x).imag()))
#  define _LPIsGreater(x,y) ((x).real() > (y).real())
#  define _LPRealPart(x)    (x).real()

#  define _LPNormxNorm(x)   ((x).real() * (x).real() + (x).imag() * (x).imag())

#  define _LPIsZeroNumber(c) (  ( fabs((c).real())     + fabs((c).imag()) ) < 10.0 * DBL_MIN  )
#  define _LPIsOneNumber(c)  (  ( fabs((c).real()-1.0) + fabs((c).imag()) ) < 10.0 * DBL_EPSILON  )
#  define _LPIsZeroReal(a) (fabs(a) < 10.0 * DBL_MIN)
#  define _LPIsOneReal(a)  (fabs(a - 1.0) < 10.0 * DBL_EPSILON)

#else

   /* real arithmetic */
   typedef double _LPNumber;
   typedef double _LPDouble;
   typedef double _LPReal;

#  define _LPfabs(c) (fabs(c))
#  define _LPmax(x, y) ((x) > (y) ? (x) : (y))
#  define _LPmin(x, y) ((x) < (y) ? (x) : (y))
#  define _LPPrintFormat(c) (c)

#  define _LPAbsRealPart(x) (fabs(x))
#  define _LPIsGreater(x,y) ((x) > (y))
#  define _LPRealPart(x)    (x)

#  define _LPNormxNorm(x)   ((x) * (x))

#  define _LPIsZeroNumber(a) (fabs(a) < 10.0 * DBL_MIN)
#  define _LPIsOneNumber(a)  (fabs(a - 1.0) < 10.0 * DBL_EPSILON)
#  define _LPIsZeroReal(a)   (fabs(a) < 10.0 * DBL_MIN)
#  define _LPIsOneReal(a)    (fabs(a - 1.0) < 10.0 * DBL_EPSILON)

#endif



#ifdef __BORLANDC__
   /* BC 2.0 does not handle IEEE arithmetic correctly */
#  if defined(_LP_LIBMESH_USE_COMPLEX_NUMBERS)
      Choke this: not tested.
#  else
      /* deliberately override above-given macros */
#     define _LPIsZero(a) (fabs(a) < 1.0e20 * DBL_MIN)
#     define _LPIsOne(a)  (fabs(a - 1.0) < 10.0 * DBL_EPSILON)
#  endif
#endif /* __BORLANDC__ */



typedef struct {
    size_t Pos;
    _LPNumber Val;
} ElType;

#ifdef _LP_INCLUDED_FROM_CPLUSPLUS
}
#endif

#endif /* ELCMP_H */

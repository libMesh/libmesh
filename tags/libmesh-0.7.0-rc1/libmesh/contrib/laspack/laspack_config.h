/****************************************************************************/
/*                                 elcmp.h                                  */
/****************************************************************************/
/*                                                                          */
/* includes definitions of libMesh                                          */
/*                                                                          */
/* Copyright (C) 1992-1996 Tomas Skalicky. All rights reserved.             */
/*                                                                          */
/****************************************************************************/
/*                                                                          */
/*        ANY USE OF THIS CODE CONSTITUTES ACCEPTANCE OF THE TERMS          */
/*              OF THE COPYRIGHT NOTICE (SEE FILE COPYRGHT.H)               */
/*                                                                          */
/****************************************************************************/

#ifndef LASPACK_CONFIG_H
#define LASPACK_CONFIG_H

#include "libmesh_config.h"

#if defined(LIBMESH_USE_REAL_NUMBERS)
#  undef _LP_LIBMESH_USE_COMPLEX_NUMBERS

#  ifdef __cplusplus
      /* someone included us from C++ */
#     define _LP_INCLUDED_FROM_CPLUSPLUS 1
#  else
      /* compile LASPACK with real arithmetic in C */
#     undef _LP_INCLUDED_FROM_CPLUSPLUS
#  endif /* __cplusplus */

#elif defined(LIBMESH_USE_COMPLEX_NUMBERS)

   /* either way, whether someone included us or LASPACK
    * is compiled, this is a C++ compiler which does not
    * want the extern "C" */
#  define _LP_LIBMESH_USE_COMPLEX_NUMBERS 1
#  undef _LP_INCLUDED_FROM_CPLUSPLUS

#else
   Choke this: something wrong.
#endif


#endif /* LASPACK_CONFIG_H */

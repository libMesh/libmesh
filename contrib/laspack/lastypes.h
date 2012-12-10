/****************************************************************************/
/*                                lastypes.h                                 */
/****************************************************************************/
/*                                                                          */
/* basic LASpack TYPES                                                      */
/*                                                                          */
/* Copyright (C) 1992-1996 Tomas Skalicky. All rights reserved.             */
/*                                                                          */
/****************************************************************************/
/*                                                                          */
/*        ANY USE OF THIS CODE CONSTITUTES ACCEPTANCE OF THE TERMS          */
/*              OF THE COPYRIGHT NOTICE (SEE FILE COPYRGHT.H)               */
/*                                                                          */
/****************************************************************************/

#ifndef LASTYPES_H
#define LASTYPES_H

#include "laspack_config.h"
#ifdef _LP_INCLUDED_FROM_CPLUSPLUS
extern "C" {
#endif

  
#include "copyrght.h"

#if defined(_LP_LIBMESH_USE_COMPLEX_NUMBERS)
   /* complex arithmetic, use bool for internal boolean */
   typedef bool    _LPBoolean;
#  ifndef _LP_DEFINED_BOOLEAN
#     define _LP_DEFINED_BOOLEAN
#     define _LPFalse false
#     define _LPTrue  true
#  endif /* _LP_DEFINED_BOOLEAN */
#else
    /* use LASPACK enum */
#  ifndef _LP_DEFINED_BOOLEAN
#     define _LP_DEFINED_BOOLEAN
         typedef enum {
             _LPFalse = 0,
             _LPTrue  = 1
         } _LPBoolean; /* boolean type */
#  endif /* _LP_DEFINED_BOOLEAN */
#endif /* _LP_LIBMESH_USE_COMPLEX_NUMBERS */

typedef enum {
    Rowws,
    Clmws
} ElOrderType;

typedef enum {
    Normal,
    Tempor
} InstanceType;

#ifdef _LP_INCLUDED_FROM_CPLUSPLUS
}
#endif

#endif /* LASTYPES_H */

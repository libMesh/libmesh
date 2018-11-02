/****************************************************************************/
/*                                 qvector.h                                 */
/****************************************************************************/
/*                                                                          */
/* type VECTOR                                                              */
/*                                                                          */
/* Copyright (C) 1992-1996 Tomas Skalicky. All rights reserved.             */
/*                                                                          */
/****************************************************************************/
/*                                                                          */
/*        ANY USE OF THIS CODE CONSTITUTES ACCEPTANCE OF THE TERMS          */
/*              OF THE COPYRIGHT NOTICE (SEE FILE COPYRGHT.H)               */
/*                                                                          */
/****************************************************************************/

#ifndef QVECTOR_H
#define QVECTOR_H

#include "laspack_config.h"
#ifdef _LP_INCLUDED_FROM_CPLUSPLUS
extern "C" {
#endif

  
#include <stdlib.h>

#include "lastypes.h"
#include "elcmp.h"
#include "copyrght.h"

typedef struct {
    char *Name;
    size_t Dim;
    InstanceType Instance;
    int LockLevel;
    _LPNumber Multipl;   /* "stretch" factor */
    _LPBoolean OwnData;
    _LPNumber *Cmp;
} QVector;

void V_Constr(QVector *V, const char *Name, size_t Dim, InstanceType Instance,
	      _LPBoolean OwnData);
void V_Destr(QVector *V);
void V_SetName(QVector *V, const char *Name);
const char *V_GetName(const QVector *V);
size_t V_GetDim(const QVector *V);
void V_SetCmp(QVector *V, size_t Ind, _LPNumber Val);
void V_SetAllCmp(QVector *V, _LPNumber Val);
void V_SetRndCmp(QVector *V);
_LPNumber V_GetCmp(QVector *V, size_t Ind);
void V_AddCmp(QVector *V, size_t Ind, _LPNumber Val);

/* macros for fast access */
#define     V__SetCmp(PtrV, Ind, Val)       (PtrV)->Cmp[Ind] = (Val)
#define     V__GetCmp(PtrV, Ind)            (PtrV)->Cmp[Ind]
#define     V__AddCmp(PtrV, Ind, Val)       (PtrV)->Cmp[Ind] += (Val)

void V_Lock(QVector *V);
void V_Unlock(QVector *V);

#ifdef _LP_INCLUDED_FROM_CPLUSPLUS
}
#endif

#endif /* QVECTOR_H */

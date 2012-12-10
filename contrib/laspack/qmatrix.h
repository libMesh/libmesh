/****************************************************************************/
/*                                qmatrix.h                                 */
/****************************************************************************/
/*                                                                          */
/* type QMATRIX                                                             */
/*                                                                          */
/* Copyright (C) 1992-1996 Tomas Skalicky. All rights reserved.             */
/*                                                                          */
/****************************************************************************/
/*                                                                          */
/*        ANY USE OF THIS CODE CONSTITUTES ACCEPTANCE OF THE TERMS          */
/*              OF THE COPYRIGHT NOTICE (SEE FILE COPYRGHT.H)               */
/*                                                                          */
/****************************************************************************/

#ifndef QMATRIX_H
#define QMATRIX_H

#include "laspack_config.h"
#ifdef _LP_INCLUDED_FROM_CPLUSPLUS
extern "C" {
#endif

  
#include <stdlib.h>

#include "lastypes.h"
#include "elcmp.h"
#include "qvector.h"
#include "copyrght.h"

typedef struct QMatrixType {
    char *Name;
    size_t Dim;
    _LPBoolean Symmetry;
    ElOrderType ElOrder;
    InstanceType Instance;
    int LockLevel;
    _LPDouble MultiplD;
    _LPDouble MultiplU;
    _LPDouble MultiplL;
    _LPBoolean OwnData;
    size_t *Len;
    ElType **El;
    _LPBoolean *ElSorted;
    _LPBoolean *DiagElAlloc;
    ElType **DiagEl;
    _LPBoolean *ZeroInDiag;
    _LPNumber *InvDiagEl;
    _LPBoolean UnitRightKer;
    _LPNumber *RightKerCmp;
    _LPBoolean UnitLeftKer;
    _LPNumber *LeftKerCmp;
    void *EigenvalInfo;
    _LPBoolean *ILUExists;
    struct QMatrixType *ILU;
} QMatrix;

void Q_Constr(QMatrix *Q, const char *Name, size_t Dim, _LPBoolean Symmetry,
              ElOrderType ElOrder, InstanceType Instance, _LPBoolean OwnData);
void Q_Destr(QMatrix *Q);
void Q_SetName(QMatrix *Q, const char *Name);
const char *Q_GetName(QMatrix *Q);
size_t Q_GetDim(QMatrix *Q);
_LPBoolean Q_GetSymmetry(QMatrix *Q);
ElOrderType Q_GetElOrder(QMatrix *Q);
void Q_SetLen(QMatrix *Q, size_t RoC, size_t Len);
size_t Q_GetLen(QMatrix *Q, size_t RoC);
void Q_SetEntry(QMatrix *Q, size_t RoC, size_t Entry, size_t Pos, _LPNumber Val);
size_t Q_GetPos(QMatrix *Q, size_t RoC, size_t Entry);
_LPNumber Q_GetVal(QMatrix *Q, size_t RoC, size_t Entry);
void Q_AddVal(QMatrix *Q, size_t RoC, size_t Entry, _LPNumber Val);

/* macros for fast access */
#define     Q__GetLen(PtrQ, RoC)               (PtrQ)->Len[RoC]
#define     Q__SetEntry(PtrQ, RoC, Entry, Pos_, Val_) { \
                (PtrQ)->El[RoC][Entry].Pos = (Pos_); \
                (PtrQ)->El[RoC][Entry].Val = (Val_); \
            }
#define     Q__GetPos(PtrQ, RoC, Entry)        (PtrQ)->El[RoC][Entry].Pos
#define     Q__GetVal(PtrQ, RoC, Entry)        (PtrQ)->El[RoC][Entry].Val
#define     Q__AddVal(PtrQ, RoC, Entry, Val_) \
                (PtrQ)->El[RoC][Entry].Val += (Val_)

_LPNumber Q_GetEl(QMatrix *Q, size_t Row, size_t Clm);

void Q_SortEl(QMatrix *Q);
void Q_AllocInvDiagEl(QMatrix *Q);

void Q_SetKer(QMatrix *Q, QVector *RightKer, QVector *LeftKer);
_LPBoolean Q_KerDefined(QMatrix *Q);

void **Q_EigenvalInfo(QMatrix *Q);

void Q_Lock(QMatrix *Q);
void Q_Unlock(QMatrix *Q);

#ifdef _LP_INCLUDED_FROM_CPLUSPLUS
}
#endif

#endif /* QMATRIX_H */

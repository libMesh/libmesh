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

#ifdef __cplusplus
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
    double Multipl;
    Boolean OwnData;
    Real *Cmp;
} QVector;

void V_Constr(QVector *V, char *Name, size_t Dim, InstanceType Instance,
	      Boolean OwnData);
void V_Destr(QVector *V);
void V_SetName(QVector *V, char *Name);
char *V_GetName(QVector *V);
size_t V_GetDim(QVector *V);
void V_SetCmp(QVector *V, size_t Ind, Real Val);
void V_SetAllCmp(QVector *V, Real Val);
void V_SetRndCmp(QVector *V);
Real V_GetCmp(QVector *V, size_t Ind);
void V_AddCmp(QVector *V, size_t Ind, Real Val);

/* macros for fast access */
#define     V__SetCmp(PtrV, Ind, Val)       (PtrV)->Cmp[Ind] = (Val)
#define     V__GetCmp(PtrV, Ind)            (PtrV)->Cmp[Ind]
#define     V__AddCmp(PtrV, Ind, Val)       (PtrV)->Cmp[Ind] += (Val)

void V_Lock(QVector *V);
void V_Unlock(QVector *V);

#ifdef __cplusplus
}
#endif

#endif /* QVECTOR_H */

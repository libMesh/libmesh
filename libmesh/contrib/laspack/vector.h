/****************************************************************************/
/*                                 vector.h                                 */
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

#ifndef VECTOR_H
#define VECTOR_H

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
} Vector;

void V_Constr(Vector *V, char *Name, size_t Dim, InstanceType Instance,
	      Boolean OwnData);
void V_Destr(Vector *V);
void V_SetName(Vector *V, char *Name);
char *V_GetName(Vector *V);
size_t V_GetDim(Vector *V);
void V_SetCmp(Vector *V, size_t Ind, Real Val);
void V_SetAllCmp(Vector *V, Real Val);
void V_SetRndCmp(Vector *V);
Real V_GetCmp(Vector *V, size_t Ind);
void V_AddCmp(Vector *V, size_t Ind, Real Val);

/* macros for fast access */
#define     V__SetCmp(PtrV, Ind, Val)       (PtrV)->Cmp[Ind] = (Val)
#define     V__GetCmp(PtrV, Ind)            (PtrV)->Cmp[Ind]
#define     V__AddCmp(PtrV, Ind, Val)       (PtrV)->Cmp[Ind] += (Val)

void V_Lock(Vector *V);
void V_Unlock(Vector *V);

#endif /* VECTOR_H */

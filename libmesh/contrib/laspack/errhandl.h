/****************************************************************************/
/*                                errhandl.h                                */
/****************************************************************************/
/*                                                                          */
/* ERRor HANDLing routines                                                  */
/*                                                                          */
/* Copyright (C) 1992-1996 Tomas Skalicky. All rights reserved.             */
/*                                                                          */
/****************************************************************************/
/*                                                                          */
/*        ANY USE OF THIS CODE CONSTITUTES ACCEPTANCE OF THE TERMS          */
/*              OF THE COPYRIGHT NOTICE (SEE FILE COPYRGHT.H)               */
/*                                                                          */
/****************************************************************************/

#ifndef ERRHANDL_H
#define ERRHANDL_H

#include <stdio.h>

#include "copyrght.h"

typedef enum {
    LASOK,
    LASMemAllocErr,
    LASLValErr,
    LASDimErr,
    LASRangeErr,
    LASSymStorErr,
    LASMatrCombErr,
    LASMulInvErr,
    LASElNotSortedErr,
    LASZeroInDiagErr,
    LASZeroPivotErr,
    LASILUStructErr,
    LASBreakdownErr,
    LASUserBreak
} LASErrIdType;

void LASError(LASErrIdType ErrId, char *ProcName, char *Object1Name,
              char *Object2Name, char *Object3Name);
void LASBreak(void);
LASErrIdType LASResult(void);
void WriteLASErrDescr(FILE *File);

#endif /* ERRHANDL_H */

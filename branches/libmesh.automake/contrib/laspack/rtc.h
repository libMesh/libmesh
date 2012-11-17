/****************************************************************************/
/*                                  rtc.h                                   */
/****************************************************************************/
/*                                                                          */
/* Residual Termination Control                                             */
/*                                                                          */
/* Copyright (C) 1992-1996 Tomas Skalicky. All rights reserved.             */
/*                                                                          */
/****************************************************************************/
/*                                                                          */
/*        ANY USE OF THIS CODE CONSTITUTES ACCEPTANCE OF THE TERMS          */
/*              OF THE COPYRIGHT NOTICE (SEE FILE COPYRGHT.H)               */
/*                                                                          */
/****************************************************************************/

#ifndef RTC_H
#define RTC_H

#include "laspack_config.h"
#ifdef _LP_INCLUDED_FROM_CPLUSPLUS
extern "C" {
#endif

  
#include "lastypes.h"
#include "qvector.h"
#include "itersolv.h"
#include "copyrght.h"

/* identifiers for iteration methods */

typedef enum {
    /* classical iterative methods */
    JacobiIterId,
    SORForwIterId,
    SORBackwIterId,
    SSORIterId,

    /* semi-iterative methods */
    ChebyshevIterId,

    /* CG and CG-like methods */
    CGIterId,
    CGNIterId,
    GMRESIterId,
    BiCGIterId,
    QMRIterId,
    CGSIterId,
    BiCGSTABIterId,

    /* multigrid and multigrid based methods */
    MGIterId,
    NestedMGIterId,
    MGPCGIterId,
    BPXPCGIterId
} IterIdType;

typedef void (*RTCAuxProcType)(int, _LPReal, _LPReal, IterIdType);

void SetRTCAccuracy(_LPReal Eps);
void SetRTCAuxProc(RTCAuxProcType AuxProc);
_LPBoolean RTCResult(int Iter, _LPReal rNorm, _LPReal bNorm, IterIdType IterId);
int GetLastNoIter(void);
_LPReal GetLastAccuracy(void);


#ifdef _LP_INCLUDED_FROM_CPLUSPLUS
}
#endif

#endif /* RTC_H */

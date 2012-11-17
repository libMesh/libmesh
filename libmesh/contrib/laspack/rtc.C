/****************************************************************************/
/*                                  rtc.c                                   */
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

#include <stddef.h>

#include "rtc.h"
#include "errhandl.h"
#include "elcmp.h"
#include "operats.h"
#include "copyrght.h"

/* accuracy for Residual Termination Control */
static _LPReal RTCEps = 1e-8;

/* auxiliary procedure to be performed by Residual Termination Control */
static RTCAuxProcType RTCAuxProc = NULL;

/* number of iterations performed during last call of a iteration method */
static int LastNoIter = 0;

/* accuracy reached during last call of a iteration method */
static _LPReal LastAcc = 0.0;

void SetRTCAccuracy(_LPReal Eps)
/* set accuracy for the RTC */
{
    RTCEps = Eps;
}

void SetRTCAuxProc(RTCAuxProcType AuxProc)
/* set auxiliary procedure of RTC */
{
    RTCAuxProc = AuxProc;
}

_LPBoolean RTCResult(int Iter, _LPReal rNorm, _LPReal bNorm, IterIdType IterId)
/* get result of RTC */
{
    _LPBoolean Result;

    if (LASResult() == LASOK) {
        if (rNorm < RTCEps * bNorm || (_LPIsZeroReal(bNorm) && _LPIsOneReal(1.0 + rNorm)))
            Result = _LPTrue;
        else
            Result = _LPFalse;

        LastNoIter = Iter;
        if (!_LPIsZeroReal(bNorm))
            LastAcc = rNorm / bNorm;
        else
            LastAcc = 1.0;

        if (RTCAuxProc != NULL)
            (*RTCAuxProc)(Iter, rNorm, bNorm, IterId);
    } else {
        Result = _LPTrue;
    }

    return(Result);
}

int GetLastNoIter()
/* get number of iterations performed during last call of a iteration method */
{
    return(LastNoIter);
}

_LPReal GetLastAccuracy()
/* get accuracy reached during last call of a iteration method */
{
    return(LastAcc);
}

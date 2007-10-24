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
static double RTCEps = 1e-8;

/* auxiliary procedure to be performed by Residual Termination Control */
static RTCAuxProcType RTCAuxProc = NULL;

/* number of iterations performed during last call of a iteration method */
static int LastNoIter = 0;

/* accuracy reached during last call of a iteration method */
static double LastAcc = 0.0;

void SetRTCAccuracy(double Eps)
/* set accuracy for the RTC */
{
    RTCEps = Eps;
}

void SetRTCAuxProc(RTCAuxProcType AuxProc)
/* set auxiliary procedure of RTC */
{
    RTCAuxProc = AuxProc;
}

Boolean RTCResult(int Iter, double rNorm, double bNorm, IterIdType IterId)
/* get result of RTC */
{
    Boolean Result;

    if (LASResult() == LASOK) {
        if (rNorm < RTCEps * bNorm || (IsZero(bNorm) && IsOne(1.0 + rNorm)))
            Result = True;
        else
            Result = False;

        LastNoIter = Iter;
        if (!IsZero(bNorm))
            LastAcc = rNorm / bNorm;
        else
            LastAcc = 1.0;

        if (RTCAuxProc != NULL)
            (*RTCAuxProc)(Iter, rNorm, bNorm, IterId);
    } else {
        Result = True;
    }

    return(Result);
}

int GetLastNoIter()
/* get number of iterations performed during last call of a iteration method */
{
    return(LastNoIter);
}

double GetLastAccuracy()
/* get accuracy reached during last call of a iteration method */
{
    return(LastAcc);
}

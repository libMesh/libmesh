/****************************************************************************/
/*                                 mlsolv.h                                 */
/****************************************************************************/
/*                                                                          */
/* Multi-Level SOLVers                                                      */
/*                                                                          */
/* Copyright (C) 1992-1996 Tomas Skalicky. All rights reserved.             */
/*                                                                          */
/****************************************************************************/
/*                                                                          */
/*        ANY USE OF THIS CODE CONSTITUTES ACCEPTANCE OF THE TERMS          */
/*              OF THE COPYRIGHT NOTICE (SEE FILE COPYRGHT.H)               */
/*                                                                          */
/****************************************************************************/

#ifndef MLSOLV_H
#define MLSOLV_H

#include "laspack_config.h"
#ifdef _LP_INCLUDED_FROM_CPLUSPLUS
extern "C" {
#endif

  
#include "qvector.h"
#include "matrix.h"
#include "qmatrix.h"
#include "itersolv.h"
#include "copyrght.h"

QVector *MGStep(int NoLevels, QMatrix *A, QVector *x, QVector *b,
	    Matrix *R, Matrix *P, int Level, int Gamma,
            IterProcType SmoothProc, int Nu1, int Nu2, 
     	    PrecondProcType PrecondProc, _LPDouble Omega,
            IterProcType SolvProc, int NuC,
	    PrecondProcType PrecondProcC, _LPDouble OmegaC);
QVector *MGIter(int NoLevels, QMatrix *A, QVector *x, QVector *b,
	    Matrix *R, Matrix *P, int MaxIter, int Gamma,
            IterProcType SmoothProc, int Nu1, int Nu2, 
	    PrecondProcType PrecondProc, _LPDouble Omega,
            IterProcType SolvProc, int NuC,
	    PrecondProcType PrecondProcC, _LPDouble OmegaC);
QVector *NestedMGIter(int NoLevels, QMatrix *A, QVector *x, QVector *b,
	    Matrix *R, Matrix *P, int Gamma,
            IterProcType SmoothProc, int Nu1, int Nu2, 
	    PrecondProcType PrecondProc, _LPDouble Omega,
            IterProcType SolvProc, int NuC,
	    PrecondProcType PrecondProcC, _LPDouble OmegaC);
QVector *MGPCGIter(int NoLevels, QMatrix *A, QVector *x, QVector *b,
	    Matrix *R, Matrix *P, int MaxIter, int NoMGIter, int Gamma,
            IterProcType SmoothProc, int Nu1, int Nu2, 
	    PrecondProcType PrecondProc, _LPDouble Omega,
            IterProcType SolvProc, int NuC,
	    PrecondProcType PrecondProcC, _LPDouble OmegaC);
QVector *BPXPrecond(int NoLevels, QMatrix *A, QVector *y, QVector *c,
            Matrix *R, Matrix *P, int Level,
            IterProcType SmoothProc, int Nu, 
	    PrecondProcType PrecondProc, _LPDouble Omega,
            IterProcType SmoothProcC, int NuC,
	    PrecondProcType PrecondProcC, _LPDouble OmegaC);
QVector *BPXPCGIter(int NoLevels, QMatrix *A, QVector *x, QVector *b,
	    Matrix *R, Matrix *P, int MaxIter,
            IterProcType SmoothProc, int Nu, 
	    PrecondProcType PrecondProc, _LPDouble Omega,
            IterProcType SmoothProcC, int NuC,
	    PrecondProcType PrecondProcC, _LPDouble OmegaC);

#ifdef _LP_INCLUDED_FROM_CPLUSPLUS
}
#endif

#endif /* MLSOLV_H */

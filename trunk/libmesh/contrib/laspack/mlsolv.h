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

#ifdef __cplusplus
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
     	    PrecondProcType PrecondProc, double Omega,
            IterProcType SolvProc, int NuC,
	    PrecondProcType PrecondProcC, double OmegaC);
QVector *MGIter(int NoLevels, QMatrix *A, QVector *x, QVector *b,
	    Matrix *R, Matrix *P, int MaxIter, int Gamma,
            IterProcType SmoothProc, int Nu1, int Nu2, 
	    PrecondProcType PrecondProc, double Omega,
            IterProcType SolvProc, int NuC,
	    PrecondProcType PrecondProcC, double OmegaC);
QVector *NestedMGIter(int NoLevels, QMatrix *A, QVector *x, QVector *b,
	    Matrix *R, Matrix *P, int Gamma,
            IterProcType SmoothProc, int Nu1, int Nu2, 
	    PrecondProcType PrecondProc, double Omega,
            IterProcType SolvProc, int NuC,
	    PrecondProcType PrecondProcC, double OmegaC);
QVector *MGPCGIter(int NoLevels, QMatrix *A, QVector *x, QVector *b,
	    Matrix *R, Matrix *P, int MaxIter, int NoMGIter, int Gamma,
            IterProcType SmoothProc, int Nu1, int Nu2, 
	    PrecondProcType PrecondProc, double Omega,
            IterProcType SolvProc, int NuC,
	    PrecondProcType PrecondProcC, double OmegaC);
QVector *BPXPrecond(int NoLevels, QMatrix *A, QVector *y, QVector *c,
            Matrix *R, Matrix *P, int Level,
            IterProcType SmoothProc, int Nu, 
	    PrecondProcType PrecondProc, double Omega,
            IterProcType SmoothProcC, int NuC,
	    PrecondProcType PrecondProcC, double OmegaC);
QVector *BPXPCGIter(int NoLevels, QMatrix *A, QVector *x, QVector *b,
	    Matrix *R, Matrix *P, int MaxIter,
            IterProcType SmoothProc, int Nu, 
	    PrecondProcType PrecondProc, double Omega,
            IterProcType SmoothProcC, int NuC,
	    PrecondProcType PrecondProcC, double OmegaC);

#ifdef __cplusplus
}
#endif

#endif /* MLSOLV_H */

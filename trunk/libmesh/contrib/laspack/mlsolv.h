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

#include "vector.h"
#include "matrix.h"
#include "qmatrix.h"
#include "itersolv.h"
#include "copyrght.h"

Vector *MGStep(int NoLevels, QMatrix *A, Vector *x, Vector *b,
	    Matrix *R, Matrix *P, int Level, int Gamma,
            IterProcType SmoothProc, int Nu1, int Nu2, 
     	    PrecondProcType PrecondProc, double Omega,
            IterProcType SolvProc, int NuC,
	    PrecondProcType PrecondProcC, double OmegaC);
Vector *MGIter(int NoLevels, QMatrix *A, Vector *x, Vector *b,
	    Matrix *R, Matrix *P, int MaxIter, int Gamma,
            IterProcType SmoothProc, int Nu1, int Nu2, 
	    PrecondProcType PrecondProc, double Omega,
            IterProcType SolvProc, int NuC,
	    PrecondProcType PrecondProcC, double OmegaC);
Vector *NestedMGIter(int NoLevels, QMatrix *A, Vector *x, Vector *b,
	    Matrix *R, Matrix *P, int Gamma,
            IterProcType SmoothProc, int Nu1, int Nu2, 
	    PrecondProcType PrecondProc, double Omega,
            IterProcType SolvProc, int NuC,
	    PrecondProcType PrecondProcC, double OmegaC);
Vector *MGPCGIter(int NoLevels, QMatrix *A, Vector *x, Vector *b,
	    Matrix *R, Matrix *P, int MaxIter, int NoMGIter, int Gamma,
            IterProcType SmoothProc, int Nu1, int Nu2, 
	    PrecondProcType PrecondProc, double Omega,
            IterProcType SolvProc, int NuC,
	    PrecondProcType PrecondProcC, double OmegaC);
Vector *BPXPrecond(int NoLevels, QMatrix *A, Vector *y, Vector *c,
            Matrix *R, Matrix *P, int Level,
            IterProcType SmoothProc, int Nu, 
	    PrecondProcType PrecondProc, double Omega,
            IterProcType SmoothProcC, int NuC,
	    PrecondProcType PrecondProcC, double OmegaC);
Vector *BPXPCGIter(int NoLevels, QMatrix *A, Vector *x, Vector *b,
	    Matrix *R, Matrix *P, int MaxIter,
            IterProcType SmoothProc, int Nu, 
	    PrecondProcType PrecondProc, double Omega,
            IterProcType SmoothProcC, int NuC,
	    PrecondProcType PrecondProcC, double OmegaC);

#endif /* MLSOLV_H */

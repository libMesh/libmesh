/****************************************************************************/
/*                                itersolv.h                                */
/****************************************************************************/
/*                                                                          */
/* ITERative SOLVers for systems of linear equations                        */
/*                                                                          */
/* Copyright (C) 1992-1996 Tomas Skalicky. All rights reserved.             */
/*                                                                          */
/****************************************************************************/
/*                                                                          */
/*        ANY USE OF THIS CODE CONSTITUTES ACCEPTANCE OF THE TERMS          */
/*              OF THE COPYRIGHT NOTICE (SEE FILE COPYRGHT.H)               */
/*                                                                          */
/****************************************************************************/

#ifndef ITERSOLV_H
#define ITERSOLV_H

#ifdef __cplusplus
extern "C" {
#endif

  
#include "qvector.h"
#include "qmatrix.h"
#include "precond.h"
#include "eigenval.h"
#include "copyrght.h"

typedef QVector *(*IterProcType)(QMatrix *, QVector *, QVector *, int,
				PrecondProcType, double);

/* classical iterative methods */

QVector *JacobiIter(QMatrix *A, QVector *x, QVector *b, int NoIter,
	    PrecondProcType Dummy, double Omega);
QVector *SORForwIter(QMatrix *A, QVector *x, QVector *b, int NoIter,
            PrecondProcType Dummy, double Omega);
QVector *SORBackwIter(QMatrix *A, QVector *x, QVector *b, int NoIter,
            PrecondProcType Dummy, double Omega);
QVector *SSORIter(QMatrix *A, QVector *x, QVector *b, int NoIter,
            PrecondProcType Dummy, double Omega);

/* semi-iterative methods */

QVector *ChebyshevIter(QMatrix *A, QVector *x, QVector *b, int MaxIter,
            PrecondProcType PrecondProc, double OmegaPrecond);

/* CG and CG-like methods */

QVector *CGIter(QMatrix *A, QVector *x, QVector *b, int MaxIter, 
            PrecondProcType PrecondProc, double OmegaPrecond);
QVector *CGNIter(QMatrix *A, QVector *x, QVector *b, int MaxIter, 
            PrecondProcType PrecondProc, double OmegaPrecond);
QVector *GMRESIter(QMatrix *A, QVector *x, QVector *b, int MaxIter, 
            PrecondProcType PrecondProc, double OmegaPrecond);
QVector *BiCGIter(QMatrix *A, QVector *x, QVector *b, int MaxIter,
            PrecondProcType PrecondProc, double OmegaPrecond);
QVector *QMRIter(QMatrix *A, QVector *x, QVector *b, int MaxIter,
            PrecondProcType PrecondProc, double OmegaPrecond);
QVector *CGSIter(QMatrix *A, QVector *x, QVector *b, int MaxIter,
            PrecondProcType PrecondProc, double OmegaPrecond);
QVector *BiCGSTABIter(QMatrix *A, QVector *x, QVector *b, int MaxIter,
            PrecondProcType PrecondProc, double OmegaPrecond);
void SetGMRESRestart(int MaxSteps);

#ifdef __cplusplus
}
#endif

#endif /* ITERSOLV_H */

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

#include "laspack_config.h"
#ifdef _LP_INCLUDED_FROM_CPLUSPLUS
extern "C" {
#endif

  
#include "qvector.h"
#include "qmatrix.h"
#include "precond.h"
#include "eigenval.h"
#include "copyrght.h"

typedef QVector *(*IterProcType)(QMatrix *, QVector *, QVector *, int,
				PrecondProcType, _LPDouble);

/* classical iterative methods */

QVector *JacobiIter(QMatrix *A, QVector *x, QVector *b, int NoIter,
		    PrecondProcType Dummy, _LPDouble Omega);
QVector *SORForwIter(QMatrix *A, QVector *x, QVector *b, int NoIter,
            PrecondProcType Dummy, _LPDouble Omega);
QVector *SORBackwIter(QMatrix *A, QVector *x, QVector *b, int NoIter,
            PrecondProcType Dummy, _LPDouble Omega);
QVector *SSORIter(QMatrix *A, QVector *x, QVector *b, int NoIter,
            PrecondProcType Dummy, _LPDouble Omega);

/* semi-iterative methods */
#ifndef _LP_LIBMESH_USE_COMPLEX_NUMBERS 
QVector *ChebyshevIter(QMatrix *A, QVector *x, QVector *b, int MaxIter,
            PrecondProcType PrecondProc, _LPDouble OmegaPrecond);
#endif

/* CG and CG-like methods */

QVector *CGIter(QMatrix *A, QVector *x, QVector *b, int MaxIter, 
            PrecondProcType PrecondProc, _LPDouble OmegaPrecond);
QVector *CGNIter(QMatrix *A, QVector *x, QVector *b, int MaxIter, 
            PrecondProcType PrecondProc, _LPDouble OmegaPrecond);
QVector *GMRESIter(QMatrix *A, QVector *x, QVector *b, int MaxIter, 
            PrecondProcType PrecondProc, _LPDouble OmegaPrecond);
QVector *BiCGIter(QMatrix *A, QVector *x, QVector *b, int MaxIter,
            PrecondProcType PrecondProc, _LPDouble OmegaPrecond);
QVector *QMRIter(QMatrix *A, QVector *x, QVector *b, int MaxIter,
            PrecondProcType PrecondProc, _LPDouble OmegaPrecond);
QVector *CGSIter(QMatrix *A, QVector *x, QVector *b, int MaxIter,
            PrecondProcType PrecondProc, _LPDouble OmegaPrecond);
QVector *BiCGSTABIter(QMatrix *A, QVector *x, QVector *b, int MaxIter,
            PrecondProcType PrecondProc, _LPDouble OmegaPrecond);
void SetGMRESRestart(int MaxSteps);

#ifdef _LP_INCLUDED_FROM_CPLUSPLUS
}
#endif

#endif /* ITERSOLV_H */

/****************************************************************************/
/*                                precond.h                                 */
/****************************************************************************/
/*                                                                          */
/* PRECONDitioners for iterative solvers of systems of linear equations     */
/*                                                                          */
/* Copyright (C) 1992-1996 Tomas Skalicky. All rights reserved.             */
/*                                                                          */
/****************************************************************************/
/*                                                                          */
/*        ANY USE OF THIS CODE CONSTITUTES ACCEPTANCE OF THE TERMS          */
/*              OF THE COPYRIGHT NOTICE (SEE FILE COPYRGHT.H)               */
/*                                                                          */
/****************************************************************************/

#ifndef PRECOND_H
#define PRECOND_H

#ifdef __cplusplus
extern "C" {
#endif

  
#include "lastypes.h"
#include "qvector.h"
#include "qmatrix.h"
#include "copyrght.h"

typedef QVector *(*PrecondProcType)(QMatrix *, QVector *, QVector *, double);

/* declaration of preconditioners */

QVector *JacobiPrecond(QMatrix *A, QVector *y, QVector *c, double Omega);
QVector *SSORPrecond(QMatrix *A, QVector *y, QVector *c, double Omega);
QVector *ILUPrecond(QMatrix *A, QVector *y, QVector *c, double Omega);

#ifdef __cplusplus
}
#endif

#endif /* PRECOND_H */

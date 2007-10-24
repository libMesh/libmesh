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

#include "laspack_config.h"
#ifdef _LP_INCLUDED_FROM_CPLUSPLUS
extern "C" {
#endif

  
#include "lastypes.h"
#include "qvector.h"
#include "qmatrix.h"
#include "copyrght.h"

typedef QVector *(*PrecondProcType)(QMatrix *, QVector *, QVector *, _LPDouble);

/* declaration of preconditioners */

QVector *JacobiPrecond(QMatrix *A, QVector *y, QVector *c, _LPDouble Omega);
QVector *SSORPrecond(QMatrix *A, QVector *y, QVector *c, _LPDouble Omega);
QVector *ILUPrecond(QMatrix *A, QVector *y, QVector *c, _LPDouble Omega);

#ifdef _LP_INCLUDED_FROM_CPLUSPLUS
}
#endif

#endif /* PRECOND_H */

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

#include "lastypes.h"
#include "vector.h"
#include "qmatrix.h"
#include "copyrght.h"

typedef Vector *(*PrecondProcType)(QMatrix *, Vector *, Vector *, double);

/* declaration of preconditioners */

Vector *JacobiPrecond(QMatrix *A, Vector *y, Vector *c, double Omega);
Vector *SSORPrecond(QMatrix *A, Vector *y, Vector *c, double Omega);
Vector *ILUPrecond(QMatrix *A, Vector *y, Vector *c, double Omega);

#endif /* PRECOND_H */

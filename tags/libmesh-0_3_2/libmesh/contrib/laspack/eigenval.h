/****************************************************************************/
/*                                eigenval.h                                */
/****************************************************************************/
/*                                                                          */
/* estimation of extremal EIGENVALues                                       */
/*                                                                          */
/* Copyright (C) 1992-1996 Tomas Skalicky. All rights reserved.             */
/*                                                                          */
/****************************************************************************/
/*                                                                          */
/*        ANY USE OF THIS CODE CONSTITUTES ACCEPTANCE OF THE TERMS          */
/*              OF THE COPYRIGHT NOTICE (SEE FILE COPYRGHT.H)               */
/*                                                                          */
/****************************************************************************/

#ifndef EIGENVAL_H
#define EIGENVAL_H

#ifdef __cplusplus
extern "C" {
#endif

  
#include "qvector.h"
#include "qmatrix.h"
#include "precond.h"
#include "copyrght.h"

/* estimation of extremal eigenvalues */

void SetEigenvalAccuracy(double Eps);
double GetMinEigenval(QMatrix *Q, PrecondProcType PrecondProc, double OmegaPrecond);
double GetMaxEigenval(QMatrix *Q, PrecondProcType PrecondProc, double OmegaPrecond);

#ifdef __cplusplus
}
#endif

#endif /* EIGENVAL_H */

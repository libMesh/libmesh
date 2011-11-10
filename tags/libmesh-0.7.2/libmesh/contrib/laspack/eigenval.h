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

#include "laspack_config.h"
#ifdef _LP_INCLUDED_FROM_CPLUSPLUS
extern "C" {
#endif

  
#include "qvector.h"
#include "qmatrix.h"
#include "precond.h"
#include "copyrght.h"

/* estimation of extremal eigenvalues */

void SetEigenvalAccuracy(_LPReal Eps);
_LPDouble GetMinEigenval(QMatrix *Q, PrecondProcType PrecondProc, _LPDouble OmegaPrecond);
_LPDouble GetMaxEigenval(QMatrix *Q, PrecondProcType PrecondProc, _LPDouble OmegaPrecond);

#ifdef _LP_INCLUDED_FROM_CPLUSPLUS
}
#endif


#endif /* EIGENVAL_H */

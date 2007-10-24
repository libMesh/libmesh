/****************************************************************************/
/*                                operats.h                                 */
/****************************************************************************/
/*                                                                          */
/* basic OPERATionS for the types vector, matrix and qmatrix                */
/*                                                                          */
/* Copyright (C) 1992-1996 Tomas Skalicky. All rights reserved.             */
/*                                                                          */
/****************************************************************************/
/*                                                                          */
/*        ANY USE OF THIS CODE CONSTITUTES ACCEPTANCE OF THE TERMS          */
/*              OF THE COPYRIGHT NOTICE (SEE FILE COPYRGHT.H)               */
/*                                                                          */
/****************************************************************************/

#ifndef OPERATS_H
#define OPERATS_H

#include "laspack_config.h"
#ifdef _LP_INCLUDED_FROM_CPLUSPLUS
extern "C" {
#endif

  
#include <stdlib.h>

#include "lastypes.h"
#include "qvector.h"
#include "matrix.h"
#include "qmatrix.h"
#include "copyrght.h"

QVector *Asgn_VV(QVector *V1, QVector *V2);
QVector *AddAsgn_VV(QVector *V1, QVector *V2);
QVector *SubAsgn_VV(QVector *V1, QVector *V2);
QVector *MulAsgn_VS(QVector *V, _LPDouble S);
QVector *Add_VV(QVector *V1, QVector *V2);
QMatrix *Add_QQ(QMatrix *Q1, QMatrix *Q2);
QVector *Sub_VV(QVector *V1, QVector *V2);
QMatrix *Sub_QQ(QMatrix *Q1, QMatrix *Q2);
QVector *Mul_SV(_LPDouble S, QVector *V);
Matrix *Mul_SM(_LPDouble S, Matrix *M);
QMatrix *Mul_SQ(_LPDouble S, QMatrix *Q);
_LPDouble Mul_VV(QVector *V1, QVector *V2);
_LPDouble InnerProd_VV(QVector *V1, QVector *V2);
QVector *Mul_MV(Matrix *M, QVector *V);
QVector *Mul_QV(QMatrix *Q, QVector *V);
QVector *MulInv_QV(QMatrix *Q, QVector *V);
Matrix *Transp_M(Matrix *M);
QMatrix *Transp_Q(QMatrix *Q);
QMatrix *Diag_Q(QMatrix *Q);
QMatrix *Upper_Q(QMatrix *Q);
QMatrix *Lower_Q(QMatrix *Q);
_LPReal l1Norm_V(QVector *V);
_LPReal l2Norm_V(QVector *V);
_LPReal MaxNorm_V(QVector *V);
QVector *OrthoRightKer_VQ(QVector *V, QMatrix *Q);
QVector *OrthoLeftKer_VQ(QVector *V, QMatrix *Q);

#ifdef _LP_INCLUDED_FROM_CPLUSPLUS
}
#endif

#endif /* OPERATS_H */

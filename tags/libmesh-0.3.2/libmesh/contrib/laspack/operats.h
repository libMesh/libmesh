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

#ifdef __cplusplus
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
QVector *MulAsgn_VS(QVector *V, double S);
QVector *Add_VV(QVector *V1, QVector *V2);
QMatrix *Add_QQ(QMatrix *Q1, QMatrix *Q2);
QVector *Sub_VV(QVector *V1, QVector *V2);
QMatrix *Sub_QQ(QMatrix *Q1, QMatrix *Q2);
QVector *Mul_SV(double S, QVector *V);
Matrix *Mul_SM(double S, Matrix *M);
QMatrix *Mul_SQ(double S, QMatrix *Q);
double Mul_VV(QVector *V1, QVector *V2);
QVector *Mul_MV(Matrix *M, QVector *V);
QVector *Mul_QV(QMatrix *Q, QVector *V);
QVector *MulInv_QV(QMatrix *Q, QVector *V);
Matrix *Transp_M(Matrix *M);
QMatrix *Transp_Q(QMatrix *Q);
QMatrix *Diag_Q(QMatrix *Q);
QMatrix *Upper_Q(QMatrix *Q);
QMatrix *Lower_Q(QMatrix *Q);
double l1Norm_V(QVector *V);
double l2Norm_V(QVector *V);
double MaxNorm_V(QVector *V);
QVector *OrthoRightKer_VQ(QVector *V, QMatrix *Q);
QVector *OrthoLeftKer_VQ(QVector *V, QMatrix *Q);

#ifdef __cplusplus
}
#endif

#endif /* OPERATS_H */

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

#include <stdlib.h>

#include "lastypes.h"
#include "vector.h"
#include "matrix.h"
#include "qmatrix.h"
#include "copyrght.h"

Vector *Asgn_VV(Vector *V1, Vector *V2);
Vector *AddAsgn_VV(Vector *V1, Vector *V2);
Vector *SubAsgn_VV(Vector *V1, Vector *V2);
Vector *MulAsgn_VS(Vector *V, double S);
Vector *Add_VV(Vector *V1, Vector *V2);
QMatrix *Add_QQ(QMatrix *Q1, QMatrix *Q2);
Vector *Sub_VV(Vector *V1, Vector *V2);
QMatrix *Sub_QQ(QMatrix *Q1, QMatrix *Q2);
Vector *Mul_SV(double S, Vector *V);
Matrix *Mul_SM(double S, Matrix *M);
QMatrix *Mul_SQ(double S, QMatrix *Q);
double Mul_VV(Vector *V1, Vector *V2);
Vector *Mul_MV(Matrix *M, Vector *V);
Vector *Mul_QV(QMatrix *Q, Vector *V);
Vector *MulInv_QV(QMatrix *Q, Vector *V);
Matrix *Transp_M(Matrix *M);
QMatrix *Transp_Q(QMatrix *Q);
QMatrix *Diag_Q(QMatrix *Q);
QMatrix *Upper_Q(QMatrix *Q);
QMatrix *Lower_Q(QMatrix *Q);
double l1Norm_V(Vector *V);
double l2Norm_V(Vector *V);
double MaxNorm_V(Vector *V);
Vector *OrthoRightKer_VQ(Vector *V, QMatrix *Q);
Vector *OrthoLeftKer_VQ(Vector *V, QMatrix *Q);

#endif /* OPERATS_H */

/****************************************************************************/
/*                                factor.h                                  */
/****************************************************************************/
/*                                                                          */
/* incomplete FACTORization for the type qmatrix                            */
/*                                                                          */
/* Copyright (C) 1992-1996 Tomas Skalicky. All rights reserved.             */
/*                                                                          */
/****************************************************************************/
/*                                                                          */
/*        ANY USE OF THIS CODE CONSTITUTES ACCEPTANCE OF THE TERMS          */
/*              OF THE COPYRIGHT NOTICE (SEE FILE COPYRGHT.H)               */
/*                                                                          */
/****************************************************************************/

#ifndef FACTOR_H
#define FACTOR_H

#include "vector.h"
#include "qmatrix.h"
#include "copyrght.h"

QMatrix *ILUFactor(QMatrix *Q);

#endif /* FACTOR_H */


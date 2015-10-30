/****************************************************************************/
/*                                eigenval.c                                */
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

#include <stddef.h>
#include <math.h>

#include "eigenval.h"
#include "elcmp.h"
#include "errhandl.h"
#include "operats.h"
#include "rtc.h"
#include "copyrght.h"
#ifndef _LP_LIBMESH_USE_COMPLEX_NUMBERS

typedef struct {
    _LPDouble MinEigenval;
    _LPDouble MaxEigenval;
    PrecondProcType PrecondProcUsed;
    _LPDouble OmegaPrecondUsed;
} EigenvalInfoType;

/* accuracy for the estimation of extremal eigenvalues */
static _LPReal EigenvalEps = 1e-4;

static void EstimEigenvals(QMatrix *A, PrecondProcType PrecondProc, _LPDouble OmegaPrecond);
static void SearchEigenval(size_t n, _LPDouble *Alpha, _LPDouble *Beta, size_t k,
	        _LPDouble BoundMin, _LPDouble BoundMax, _LPBoolean *Found, _LPDouble *Lambda);
static size_t NoSmallerEigenvals(size_t n, _LPDouble *Alpha, _LPDouble *Beta, _LPDouble Lambda);

void SetEigenvalAccuracy(_LPReal Eps)
/* set accuracy for the estimation of extremal eigenvalues */
{
    EigenvalEps = Eps;
}

_LPDouble GetMinEigenval(QMatrix *A, PrecondProcType PrecondProc, _LPDouble OmegaPrecond)
/* returns estimate for minimum eigenvalue of the matrix A */
{
    _LPDouble MinEigenval;
    
    EigenvalInfoType *EigenvalInfo;
    
    Q_Lock(A);
    
    if (LASResult() == LASOK) {
        EigenvalInfo = (EigenvalInfoType *)*(Q_EigenvalInfo(A));
        /* if eigenvalues not estimated yet, ... */
        if (EigenvalInfo == NULL) {
            EigenvalInfo = (EigenvalInfoType *)malloc(sizeof(EigenvalInfoType));
            if (EigenvalInfo != NULL) {
	        *(Q_EigenvalInfo(A)) = (void *)EigenvalInfo;
                EstimEigenvals(A, PrecondProc, OmegaPrecond);
            } else {
                LASError(LASMemAllocErr, "GetMinEigenval", Q_GetName(A), NULL, NULL);
            }
        }
	
        /* if eigenvalues estimated with an other preconditioner, ... */
        if (EigenvalInfo->PrecondProcUsed != PrecondProc
	    || EigenvalInfo->OmegaPrecondUsed != OmegaPrecond) {
            EstimEigenvals(A, PrecondProc, OmegaPrecond);
        }

        if (LASResult() == LASOK)
            MinEigenval = EigenvalInfo->MinEigenval;
        else 
            MinEigenval = 1.0;
    } else {
        MinEigenval = 1.0;
    }
    
    return(MinEigenval);             
}

_LPDouble GetMaxEigenval(QMatrix *A, PrecondProcType PrecondProc, _LPDouble OmegaPrecond)
/* returns estimate for maximum eigenvalue of the matrix A */
{
    _LPDouble MaxEigenval;

    EigenvalInfoType *EigenvalInfo;
    
    Q_Lock(A);
    
    if (LASResult() == LASOK) {
        EigenvalInfo = (EigenvalInfoType *)*(Q_EigenvalInfo(A));
        /* if eigenvalues not estimated yet, ... */
        if (EigenvalInfo == NULL) {
            EigenvalInfo = (EigenvalInfoType *)malloc(sizeof(EigenvalInfoType));
            if (EigenvalInfo != NULL) {
	        *(Q_EigenvalInfo(A)) = (void *)EigenvalInfo;
                EstimEigenvals(A, PrecondProc, OmegaPrecond);
            } else {
                LASError(LASMemAllocErr, "GetMaxEigenval", Q_GetName(A), NULL, NULL);
            }
        }
	
        /* if eigenvalues estimated with an other preconditioner, ... */
        if (EigenvalInfo->PrecondProcUsed != PrecondProc
	    || EigenvalInfo->OmegaPrecondUsed != OmegaPrecond) {
            EstimEigenvals(A, PrecondProc, OmegaPrecond);
        }

        if (LASResult() == LASOK)
            MaxEigenval = EigenvalInfo->MaxEigenval;
        else 
            MaxEigenval = 1.0;
    } else {
        MaxEigenval = 1.0;
    }
    
    return(MaxEigenval);             
}

static void EstimEigenvals(QMatrix *A, PrecondProcType PrecondProc, _LPDouble OmegaPrecond)
/* estimates extremal eigenvalues of the matrix A by means of the Lanczos method */
{
    /*
     *  for details to the Lanczos algorithm see
     *
     *  G. H. Golub, Ch. F. van Loan:
     *  Matrix Computations;
     *  North Oxford Academic, Oxford, 1986
     *
     *  (for modification for preconditioned matrices compare with sec. 10.3) 
     *
     */
   
    _LPDouble LambdaMin = 0.0, LambdaMax = 0.0;
    _LPDouble LambdaMinOld, LambdaMaxOld;
    _LPDouble GershBoundMin = 0.0, GershBoundMax = 0.0;
    _LPDouble *Alpha, *Beta;
    size_t Dim, j;
    _LPBoolean Found;
    QVector q, qOld, h, p;

    Q_Lock(A);
    
    Dim = Q_GetDim(A);
    V_Constr(&q, "q", Dim, Normal, _LPTrue);
    V_Constr(&qOld, "qOld", Dim, Normal, _LPTrue);
    V_Constr(&h, "h", Dim, Normal, _LPTrue);
    if (PrecondProc != NULL)
        V_Constr(&p, "p", Dim, Normal, _LPTrue);
   
    if (LASResult() == LASOK) {
        Alpha = (_LPDouble *)malloc((Dim + 1) * sizeof(_LPDouble));
        Beta = (_LPDouble *)malloc((Dim + 1) * sizeof(_LPDouble));
        if (Alpha != NULL && Beta != NULL) {
	    j = 0;
            
            V_SetAllCmp(&qOld, 0.0);
            V_SetRndCmp(&q);
	    if (Q_KerDefined(A))
	        OrthoRightKer_VQ(&q, A);
            if (Q_GetSymmetry(A) && PrecondProc != NULL) {
	        (*PrecondProc)(A, &p, &q, OmegaPrecond);
                MulAsgn_VS(&q, 1.0 / sqrt(Mul_VV(&q, &p)));
	    } else {
                MulAsgn_VS(&q, 1.0 / l2Norm_V(&q));
	    }
            
            Beta[0] = 1.0;
            do {
	        j++;
                if (Q_GetSymmetry(A) && PrecondProc != NULL) {
		    /* p = M^(-1) q */
		    (*PrecondProc)(A, &p, &q, OmegaPrecond);
		    /* h = A p */
                    Asgn_VV(&h, Mul_QV(A, &p));
	            if (Q_KerDefined(A))
	                OrthoRightKer_VQ(&h, A);
		    /* Alpha = p . h */
                    Alpha[j] = Mul_VV(&p, &h);
		    /* r = h - Alpha q - Beta qOld */
                    SubAsgn_VV(&h, Add_VV(Mul_SV(Alpha[j], &q), Mul_SV(Beta[j-1], &qOld)));
                    /* z = M^(-1) r */
		    (*PrecondProc)(A, &p, &h, OmegaPrecond);
		    /* Beta = sqrt(r . z) */
                    Beta[j] = sqrt(Mul_VV(&h, &p));
                    Asgn_VV(&qOld, &q);
		    /* q = r / Beta */
                    Asgn_VV(&q, Mul_SV(1.0 / Beta[j], &h));
		} else {
		    /* h = A p */
  		    if (Q_GetSymmetry(A)) {
                        Asgn_VV(&h, Mul_QV(A, &q));
		    } else {
                        if (PrecondProc != NULL) {
			    (*PrecondProc)(A, &h, Mul_QV(A, &q), OmegaPrecond);
			    (*PrecondProc)(Transp_Q(A), &h, &h, OmegaPrecond);
                            Asgn_VV(&h, Mul_QV(Transp_Q(A), &h));
                        } else {
                            Asgn_VV(&h, Mul_QV(Transp_Q(A), Mul_QV(A, &q)));
                        }
                    }
	            if (Q_KerDefined(A))
	                OrthoRightKer_VQ(&h, A); 
		    /* Alpha = q . h */
                    Alpha[j] = Mul_VV(&q, &h);
		    /* r = h - Alpha q - Beta qOld */
                    SubAsgn_VV(&h, Add_VV(Mul_SV(Alpha[j], &q), Mul_SV(Beta[j-1], &qOld)));
                    /* Beta = || r || */
		    Beta[j] = l2Norm_V(&h);
                    Asgn_VV(&qOld, &q);
		    /* q = r / Beta */
                    Asgn_VV(&q, Mul_SV(1.0 / Beta[j], &h));
		}
		
		LambdaMaxOld = LambdaMax;
                LambdaMinOld = LambdaMin;
		
                /* determination of extremal eigenvalues of the tridiagonal matrix
                   (Beta[i-1] Alpha[i] Beta[i]) (where 1 <= i <= j) 
		   by means of the method of bisection; bounds for eigenvalues
		   are determined after Gershgorin circle theorem */
                if (j == 1) {
		    GershBoundMin = Alpha[1] - _LPfabs(Beta[1]);
	  	    GershBoundMax = Alpha[1] + _LPfabs(Beta[1]);
		    
                    LambdaMin = Alpha[1];
                    LambdaMax = Alpha[1];
		} else {
		    GershBoundMin = _LPmin(Alpha[j] - _LPfabs(Beta[j]) - _LPfabs(Beta[j - 1]),
					GershBoundMin);
		    GershBoundMax = _LPmax(Alpha[j] + _LPfabs(Beta[j]) + _LPfabs(Beta[j - 1]),
				        GershBoundMax);

                    SearchEigenval(j, Alpha, Beta, 1, GershBoundMin, LambdaMin,
		        &Found, &LambdaMin);
		    if (!Found)
                        SearchEigenval(j, Alpha, Beta, 1, GershBoundMin, GershBoundMax,
		            &Found, &LambdaMin);
		    
	            SearchEigenval(j, Alpha, Beta, j, LambdaMax, GershBoundMax,
		        &Found, &LambdaMax);
		    if (!Found)
                        SearchEigenval(j, Alpha, Beta, j, GershBoundMin, GershBoundMax,
		            &Found, &LambdaMax);
                }
            } while (!_LPIsZeroReal(Beta[j]) && j < Dim
		&& (_LPfabs(LambdaMin - LambdaMinOld) > EigenvalEps * LambdaMin
                || _LPfabs(LambdaMax - LambdaMaxOld) > EigenvalEps * LambdaMax)
                && LASResult() == LASOK);
                
	    if (Q_GetSymmetry(A)) {
	        LambdaMin = (1.0 - j * EigenvalEps) * LambdaMin;
	    } else {
	      LambdaMin = (1.0 - sqrt((_LPDouble) j) * EigenvalEps) * sqrt(LambdaMin);
            }
            if (Alpha != NULL)
                free(Alpha);
            if (Beta != NULL)
                free(Beta);
        } else {
            LASError(LASMemAllocErr, "EstimEigenvals", Q_GetName(A), NULL, NULL);
	}

    }
    
    V_Destr(&q);
    V_Destr(&qOld);
    V_Destr(&h);
    if (PrecondProc != NULL)
        V_Destr(&p);
    
    if (LASResult() == LASOK) {
        ((EigenvalInfoType *)*(Q_EigenvalInfo(A)))->MinEigenval = LambdaMin;
        ((EigenvalInfoType *)*(Q_EigenvalInfo(A)))->MaxEigenval = LambdaMax;
        ((EigenvalInfoType *)*(Q_EigenvalInfo(A)))->PrecondProcUsed = PrecondProc;
        ((EigenvalInfoType *)*(Q_EigenvalInfo(A)))->OmegaPrecondUsed = OmegaPrecond;
    } else {
        ((EigenvalInfoType *)*(Q_EigenvalInfo(A)))->MinEigenval = 1.0;
        ((EigenvalInfoType *)*(Q_EigenvalInfo(A)))->MaxEigenval = 1.0;
        ((EigenvalInfoType *)*(Q_EigenvalInfo(A)))->PrecondProcUsed = NULL;
        ((EigenvalInfoType *)*(Q_EigenvalInfo(A)))->OmegaPrecondUsed = 1.0;
    }

    Q_Unlock(A);
}

static void SearchEigenval(size_t n, _LPDouble *Alpha, _LPDouble *Beta, size_t k,
         _LPDouble BoundMin, _LPDouble BoundMax, _LPBoolean *Found, _LPDouble *Lambda)
/* search the k-th eigenvalue of the tridiagonal matrix
   (Beta[i-1] Alpha[i] Beta[i]) (where 1 <= i <= n) 
   by means of the method of bisection */
{
    /*
     *  for details to the method of bisection see
     *
     *  G. H. Golub, Ch. F. van Loan:
     *  Matrix Computations;
     *  North Oxford Academic, Oxford, 1986
     *
     */
   
    if (NoSmallerEigenvals(n, Alpha, Beta, BoundMin) < k
	&& NoSmallerEigenvals(n, Alpha, Beta, BoundMax) >= k) {
        while (_LPfabs(BoundMax - BoundMin) > 0.01 * EigenvalEps 
	    * (_LPfabs(BoundMin) + _LPfabs(BoundMax))) {
            *Lambda = 0.5 * (BoundMin + BoundMax);
	    if (NoSmallerEigenvals(n, Alpha, Beta, *Lambda) >= k) 
	        BoundMax = *Lambda;
   	    else
	        BoundMin = *Lambda;
        }
	*Lambda = BoundMax;

	*Found = _LPTrue;
    } else {
	*Found = _LPFalse;
    }
}

static size_t NoSmallerEigenvals(size_t n, _LPDouble *Alpha, _LPDouble *Beta, _LPDouble Lambda)
/* returns number of eigenvalues of the tridiagonal matrix
   (Beta[i-1] Alpha[i] Beta[i]) (where 1 <= i <= n) 
   which are less then Lambda */
{
    size_t No;
    
    _LPDouble p, pNew, pOld, Sign;
    size_t i;
    
    No = 0;
    
    pOld = 1.0;
    p = (Alpha[1] - Lambda) / _LPfabs(Beta[1]);
    /* check for change of sign */
    if (_LPIsZeroReal(p) || p * pOld < 0)
        No++;
    
    for (i = 2; i <= n; i++) {
        Sign = Beta[i-1] / _LPfabs(Beta[i-1]);
        pNew = ((Alpha[i] - Lambda) * p - Beta[i-1] * Sign * pOld) / _LPfabs(Beta[i]);
        pOld = p;
        p = pNew;
	
	/* check for change of sign */
	if (p * pOld < 0 || (_LPIsZeroReal(p) && !_LPIsZeroReal(pOld)))
	    No++;
    }
   
    return(No);
}



#else

void SetEigenvalAccuracy(_LPReal /* Eps */)
{
    printf("ERROR: Eigenvalue estimation not implemented for complex arithmetic.\n");
    abort();
}


_LPDouble GetMinEigenval(QMatrix * /* Q */, 
			 PrecondProcType /* PrecondProc */, 
			 _LPDouble /* OmegaPrecond */)
{
    printf("ERROR: Eigenvalue estimation not implemented for complex arithmetic.\n");
    abort();
}


_LPDouble GetMaxEigenval(QMatrix * /* Q */, 
			 PrecondProcType /* PrecondProc */, 
			 _LPDouble /* OmegaPrecond */)
{
    printf("ERROR: Eigenvalue estimation not implemented for complex arithmetic.\n");
    abort();
}

#endif /* _LP_LIBMESH_USE_COMPLEX_NUMBERS */


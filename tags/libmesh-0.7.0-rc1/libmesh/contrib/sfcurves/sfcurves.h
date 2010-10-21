#ifndef __sfcurves_h__
#define __sfcurves_h__
 
#ifdef __cplusplus
extern "C" {
#endif


/*--------------------- Prototypes ---------------------------*/
void morton(double *x, double *y, double *z, int *N, int *table);

void hilbert(double *x, double *y, double *z, int *N, int *table);

#ifdef __cplusplus
}
#endif

#endif /* #define __sfcurves_h__ */

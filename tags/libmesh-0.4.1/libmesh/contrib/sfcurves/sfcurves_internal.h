#ifndef __sfcurves_internal_h__
#define __sfcurves_internal_h__

#include "sfcurves.h"
 
#ifdef __cplusplus
extern "C" {
#endif


#define MIN(a,b) ((a)<(b)?(a):(b))
#define MAX(a,b) ((a)>(b)?(a):(b))
#define INT unsigned long int

struct m_str {
  double x;
  double y;
  double z;
  INT index[3];
  int table;
};


/*--------------------- Prototypes ---------------------------*/
void hsfc3d(unsigned coord[], unsigned * nkey, unsigned   key[]);   

void interleave(double x, double y, double z, double deg, INT index[3]);

int cmp_indx(const void *a, const void *b);
 
#ifdef __cplusplus
}
#endif

#endif /* #define __sfcurves_internal_h__ */

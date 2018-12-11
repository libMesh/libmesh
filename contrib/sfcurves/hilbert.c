#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include "sfcurves_internal.h"



/* Bits per unsigned word */

#define MaxBits ( sizeof(unsigned) * CHAR_BIT )

/*--------------------------------------------------------------------*/
/* 3D Hilbert Space-filling curve */

void hsfc3d(
  unsigned   coord[] , /* IN: Normalized integer coordinates */
  unsigned * nkey ,    /* IN: Word length of 'key' */
  unsigned   key[] )   /* OUT: space-filling curve key */
{
  static int init = 0 ;
  static unsigned char gray_inv[ 2*2*2 ] ;

  const unsigned NKey  = ( 3 < *nkey ) ? 3 : (*nkey) ;
  const unsigned NBits = ( MaxBits * NKey ) / 3 ;

  unsigned i ;
  unsigned char axis[3+3] ;

  /* GRAY coding */

  if ( ! init ) {
    unsigned char gray[ 2*2*2 ] ;
    register unsigned k ;
    register unsigned j ;

    gray[0] = 0 ;
    for ( k = 1 ; k < sizeof(gray) ; k <<= 1 ) {
      for ( j = 0 ; j < k ; j++ ) gray[k+j] = k | gray[k-(j+1)] ;
    }
    for ( k = 0 ; k < sizeof(gray) ; k++ ) gray_inv[ gray[k] ] = k ;
    init = 1 ;
  }

  /* Zero out the key */

  for ( i = 0 ; i < NKey ; ++i ) key[i] = 0 ;

  axis[0] = 0 << 1 ;
  axis[1] = 1 << 1 ;
  axis[2] = 2 << 1 ;

  for ( i = 1 ; i <= NBits ; i++ ) {
    const unsigned s = MaxBits - i ;
    const unsigned c = gray_inv[
      (((( coord[ axis[0] >> 1 ] >> s ) ^ axis[0] ) & 01 ) << 0 ) |
      (((( coord[ axis[1] >> 1 ] >> s ) ^ axis[1] ) & 01 ) << 1 ) |
      (((( coord[ axis[2] >> 1 ] >> s ) ^ axis[2] ) & 01 ) << 2 ) ];
    unsigned n ;

    /* Set the 3bits */

    for ( n = 0 ; n < 3 ; ++n ) {
      const unsigned bit   = 01 & ( c >> ( 2 - n ) );  /* Bit value  */
      const unsigned off   = 3 * i + n ;               /* Bit offset */
      const unsigned which = off / MaxBits ;           /* Which word */
      const unsigned shift = MaxBits - off % MaxBits ; /* Which bits */

      /**
       * The previous test was: "if (MaxBits == shift)" but this
       * caused an off-by-one error for some meshes so we've changed
       * it to what you see below. This seems to fix the illegal
       * access problem, but I'm not sure if it's correct wrt the
       * algorithm itself.
       */
      if (which > 2) {
        key[ which - 1 ] |= bit ;
      }
      else {
        key[ which ] |= bit << shift ;
      }
    }

    /* Determine the recursive quadrant */

    axis[3+0] = axis[0] ;
    axis[3+1] = axis[1] ;
    axis[3+2] = axis[2] ;

    switch( c ) {
    case 0:
      axis[0] = axis[3+2];
      axis[1] = axis[3+1];
      axis[2] = axis[3+0];
      break ;
    case 1:
      axis[0] = axis[3+0];
      axis[1] = axis[3+2];
      axis[2] = axis[3+1];
      break ;
    case 2:
      axis[0] = axis[3+0];
      axis[1] = axis[3+1];
      axis[2] = axis[3+2];
      break ;
    case 3:
      axis[0] = axis[3+2] ^ 01 ;
      axis[1] = axis[3+0] ^ 01 ;
      axis[2] = axis[3+1];
      break ;
    case 4:
      axis[0] = axis[3+2];
      axis[1] = axis[3+0] ^ 01 ;
      axis[2] = axis[3+1] ^ 01 ;
      break ;
    case 5:
      axis[0] = axis[3+0];
      axis[1] = axis[3+1];
      axis[2] = axis[3+2];
      break ;
    case 6:
      axis[0] = axis[3+0];
      axis[1] = axis[3+2] ^ 01 ;
      axis[2] = axis[3+1] ^ 01 ;
      break ;
    case 7:
      axis[0] = axis[3+2] ^ 01 ;
      axis[1] = axis[3+1];
      axis[2] = axis[3+0] ^ 01 ;
      break ;
    default:
      exit(-1);
    }
  }
}

/*--------------------------------------------------------------------*/


void fhsfc3d(
  double     coord[] , /* IN: Normalized floating point coordinates */
  unsigned * nkey ,    /* IN: Word length of key */
  unsigned   key[] )   /* OUT: space-filling curve key */
{
  const double imax = UINT_MAX;
  unsigned c[3] ;
  c[0] = (unsigned) (coord[0] * imax) ;
  c[1] = (unsigned) (coord[1] * imax) ;
  c[2] = (unsigned) (coord[2] * imax) ;
  hsfc3d( c , nkey , key );
}



void hilbert(double *x, double *y, double *z, int *N, int *table)
{
  unsigned int index[3];

  double extrx[2];
  double extry[2];
  double extrz[2];
  struct m_str *s;

  int i;

  double temp[3]={0.0, 0.0, 0.0};
  unsigned leng=3;

  s=malloc((*N)*sizeof(struct m_str ));

  extrx[0]=extrx[1]=x[0];
  extry[0]=extry[1]=y[0];
  extrz[0]=extrz[1]=z[0];
  table[0]=1;

  for(i=1;i<*N;i++) {
    extrx[0]=MIN(extrx[0],x[i]);
    extrx[1]=MAX(extrx[1],x[i]);
    extry[0]=MIN(extry[0],y[i]);
    extry[1]=MAX(extry[1],y[i]);
    extrz[0]=MIN(extrz[0],z[i]);
    extrz[1]=MAX(extrz[1],z[i]);
    table[i]=i+1;
  }


  for(i=0;i<*N;i++) {
    temp[0]=(x[i]-extrx[0])/(extrx[1]-extrx[0]);
    temp[1]=(y[i]-extry[0])/(extry[1]-extry[0]);
    temp[2]=(z[i]-extrz[0])/(extrz[1]-extrz[0]);

    fhsfc3d(temp,&leng,index);

    s[i].x=x[i];
    s[i].y=y[i];
    s[i].z=z[i];
    s[i].index[0]=index[0];
    s[i].index[1]=index[1];
    s[i].index[2]=index[2];
    s[i].table=table[i];

  }


  qsort(s,*N,sizeof(struct m_str), cmp_indx);

  for(i=0;i<*N;i++) {

    x[i]=s[i].x;
    y[i]=s[i].y;
    z[i]=s[i].z;
    table[i]=s[i].table;

  }

  free(s);
}

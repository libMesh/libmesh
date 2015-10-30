#include <stdlib.h>
#include <math.h>
#include "sfcurves_internal.h"



void interleave(double x, double y, double z, double deg, unsigned int index[3])
{
  unsigned int tempx,tempy,tempz;
  
  int i,j=2,cnt=0;

  tempx=(unsigned int )(x*deg);
  tempy=(unsigned int )(y*deg);
  tempz=(unsigned int )(z*deg);

  index[0]=0;
  index[1]=0;
  index[2]=0;

  for(i=sizeof(unsigned int)*8-1;i>=0;i--) {
    index[j] += (tempx >> i) & 01;
    index[j]  = index[j] << 01;

    if (( cnt % (sizeof(unsigned int )*8) == 0) && (cnt !=0 ) ) {
      cnt = 0;
      j--;
    }
    else cnt++;
    
    index[j] += (tempy >> i) & 01;
    index[j]  = index[j] << 01;

    if (( cnt % (sizeof(unsigned int )*8) == 0) && (cnt !=0 ) ) {
      cnt = 0;
      j--;
    }
    else cnt++;
    
    index[j] += (tempz >> i) & 01;
    index[j]  = index[j] << 01;

    if (( cnt % (sizeof(unsigned int )*8) == 0) && (cnt !=0 ) ) {
      cnt = 0;
      j--;
    }
    else cnt++;

  }
}



void morton(double *x, double *y, double *z, int *N, int *table)
{  
  unsigned int index[3];

  double extrx[2];
  double extry[2];
  double extrz[2];
  struct m_str *s;
    
  int i;
  
  double big=1.0e9;
  
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



  for(i=0;i<*N;i++){
    interleave((x[i]-extrx[0])/(extrx[1]-extrx[0]),
               (y[i]-extry[0])/(extry[1]-extry[0]),
               (z[i]-extrz[0])/(extrz[1]-extrz[0]),
               big, index);
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




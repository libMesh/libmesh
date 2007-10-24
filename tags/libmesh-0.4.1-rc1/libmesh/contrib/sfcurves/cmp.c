#include "sfcurves_internal.h"

int cmp_indx(const void *a, const void *b)
{    
  const struct m_str *tmpa=a,*tmpb=b;
  
  if      (tmpa[0].index[2] > tmpb[0].index[2]) return(1);
  else if (tmpa[0].index[2] < tmpb[0].index[2]) return(-1);
  else if (tmpa[0].index[2] == tmpb[0].index[2]) {
    if      (tmpa[0].index[1] > tmpb[0].index[1]) return(1);
    else if (tmpa[0].index[1] < tmpb[0].index[1]) return(-1);
    else if (tmpa[0].index[1] == tmpb[0].index[1]) {
      if      (tmpa[0].index[0] > tmpb[0].index[0]) return(1);
      else if (tmpa[0].index[0] < tmpb[0].index[0]) return(-1);
    }
  }
  
  else return(0);      
}

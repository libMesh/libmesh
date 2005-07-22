#include "sfcurves_internal.h"

int cmp_indx(const void *a, const void *b)
{    
  const struct m_str *tmpa=a,*tmpb=b;

  if (tmpa->index[0] > tmpb->index[0]) return  1;
  if (tmpa->index[0] < tmpb->index[0]) return -1;
  
  if (tmpa->index[1] > tmpb->index[1]) return  1;
  if (tmpa->index[1] < tmpb->index[1]) return -1;
  
  if (tmpa->index[2] > tmpb->index[2]) return  1;
  if (tmpa->index[2] < tmpb->index[2]) return -1;
  
  return 0;      
}

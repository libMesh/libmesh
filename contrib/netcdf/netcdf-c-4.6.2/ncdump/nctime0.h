/*********************************************************************
 *   Copyright 2008, University Corporation for Atmospheric Research
 *   See netcdf/COPYRIGHT file for copying and redistribution conditions.
 *   Russ Rew
 *********************************************************************/
#ifndef _NCTIME0_H
#define _NCTIME0_H

#include "nctime.h"


extern void insert_bounds_info(int ncid, int varid, ncatt_t *attp);
extern int is_valid_time_unit(const char *units);
extern int is_bounds_att(ncatt_t *attp);
extern void get_timeinfo(int ncid, int varid, ncvar_t *vp);
extern void print_att_times(int ncid, int varid, const ncatt_t *att);

#endif /* _NCTIME0_H */


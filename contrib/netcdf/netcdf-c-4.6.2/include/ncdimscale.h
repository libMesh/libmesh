/* This is part of the netCDF package.
   Copyright 2011 University Corporation for Atmospheric Research/Unidata
   See COPYRIGHT file for conditions of use.

   Includes for some HDF5 stuff needed by tests.
*/

#ifndef _NCDIMSCALE_H_
#define _NCDIMSCALE_H_

#include <hdf5.h>

typedef struct hdf5_objid 
{
   unsigned long fileno[2]; /* file number */
   haddr_t objno[2]; /* object number */
} HDF5_OBJID_T; 

#endif

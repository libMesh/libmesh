/* Copyright 2018, UCAR/Unidata See netcdf/COPYRIGHT file for copying
 * and redistribution conditions.*/
/**
 * @file @internal This file handles the variable functions for the
 * HDF4 dispatch layer.
 *
 * @author Ed Hartnett
 */

#include "config.h"
#include <nc4internal.h>
#include "hdf4dispatch.h"
#include "nc4dispatch.h"
#include <mfhdf.h>

/**
 * Read an array of values. This is called by nc_get_vara() for
 * netCDF-4 files, as well as all the other nc_get_vara_*
 * functions. HDF4 files are handled as a special case.
 *
 * @param ncid File ID.
 * @param varid Variable ID.
 * @param startp Array of start indices.
 * @param countp Array of counts.
 * @param ip pointer that gets the data.
 * @param memtype The type of these data after it is read into memory.

 * @return ::NC_NOERR for success.
 * @return ::NC_EBADID Bad ncid.
 * @author Ed Hartnett, Dennis Heimbigner
 */
int
NC_HDF4_get_vara(int ncid, int varid, const size_t *startp,
                 const size_t *countp, void *ip, int memtype)
{
   NC_VAR_HDF4_INFO_T *hdf4_var;
   NC_VAR_INFO_T *var;
   int32 start32[NC_MAX_VAR_DIMS], edge32[NC_MAX_VAR_DIMS];
   size_t nelem = 1;
   void *data;
   int retval, d;
   int range_error;

   LOG((2, "%s: ncid 0x%x varid %d memtype %d", __func__, ncid, varid,
        memtype));

   /* No scalars in HDF4 SD API. Caller must also provide place to put
    * data. */
   if (!startp || !countp || !ip)
      return NC_EINVAL;

   /* Find our metadata for this file, group, and var. */
   if ((retval = nc4_find_grp_h5_var(ncid, varid, NULL, NULL, &var)))
      return retval;
   assert(var && var->hdr.name && var->format_var_info);

   /* Get the HDF4 specific var metadata. */
   hdf4_var = (NC_VAR_HDF4_INFO_T *)var->format_var_info;

   /* Convert starts/edges to the int32 type HDF4 wants. Also learn
    * how many elements of data are being read. */
   for (d = 0; d < var->ndims; d++)
   {
      start32[d] = startp[d];
      edge32[d] = countp[d];
      nelem *= countp[d];
   }

   /* If memtype was not give, use variable type. */
   if (memtype == NC_NAT)
      memtype = var->type_info->hdr.id;

   /* If we need to convert data, allocate temp storage. */
   if (var->type_info->hdr.id == memtype)
      data = ip;
   else
      if (!(data = malloc(var->type_info->size * nelem)))
         return NC_ENOMEM;

   /* Read the data with HDF4. */
   if (SDreaddata(hdf4_var->sdsid, start32, NULL, edge32, data))
      return NC_EHDFERR;

   /* Do we need to convert data? */
   if (var->type_info->hdr.id != memtype)
   {
      if ((retval = nc4_convert_type(data, ip, var->type_info->hdr.id, memtype, nelem,
                                     &range_error, NULL, 0)))
         return retval;
      free(data);
      if (range_error)
         return range_error;
   }

   return NC_NOERR;
}

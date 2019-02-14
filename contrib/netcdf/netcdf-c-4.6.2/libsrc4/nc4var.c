/* Copyright 2003-2006, University Corporation for Atmospheric
 * Research. See COPYRIGHT file for copying and redistribution
 * conditions.*/
/**
 * @file
 * @internal This file is part of netcdf-4, a netCDF-like interface
 * for HDF5, or a HDF5 backend for netCDF, depending on your point of
 * view. This file handles the NetCDF-4 variable functions.
 *
 * @author Ed Hartnett, Dennis Heimbigner, Ward Fisher
 */

#include "config.h"
#include <nc4internal.h>
#include "nc4dispatch.h"
#include "hdf5internal.h"
#include <math.h>

/**
 * @internal Set chunk cache size for a variable. This is the internal
 * function called by nc_set_var_chunk_cache().
 *
 * @param ncid File ID.
 * @param varid Variable ID.
 * @param size Size in bytes to set cache.
 * @param nelems Number of elements in cache.
 * @param preemption Controls cache swapping.
 *
 * @returns ::NC_NOERR No error.
 * @returns ::NC_EBADID Bad ncid.
 * @returns ::NC_ENOTVAR Invalid variable ID.
 * @returns ::NC_ESTRICTNC3 Attempting netcdf-4 operation on strict nc3 netcdf-4 file.
 * @returns ::NC_EINVAL Invalid input.
 * @returns ::NC_EHDFERR HDF5 error.
 * @author Ed Hartnett
 */
int
NC4_set_var_chunk_cache(int ncid, int varid, size_t size, size_t nelems,
                        float preemption)
{
   NC *nc;
   NC_GRP_INFO_T *grp;
   NC_FILE_INFO_T *h5;
   NC_VAR_INFO_T *var;
   int retval;

   /* Check input for validity. */
   if (preemption < 0 || preemption > 1)
      return NC_EINVAL;

   /* Find info for this file and group, and set pointer to each. */
   if ((retval = nc4_find_nc_grp_h5(ncid, &nc, &grp, &h5)))
      return retval;
   assert(nc && grp && h5);

   /* Find the var. */
   var = (NC_VAR_INFO_T*)ncindexith(grp->vars,varid);
   if(!var)
      return NC_ENOTVAR;
   assert(var && var->hdr.id == varid);

   /* Set the values. */
   var->chunk_cache_size = size;
   var->chunk_cache_nelems = nelems;
   var->chunk_cache_preemption = preemption;

   if ((retval = nc4_reopen_dataset(grp, var)))
      return retval;

   return NC_NOERR;
}

/**
 * @internal A wrapper for NC4_set_var_chunk_cache(), we need this
 * version for fortran. Negative values leave settings as they are.
 *
 * @param ncid File ID.
 * @param varid Variable ID.
 * @param size Size in bytes to set cache.
 * @param nelems Number of elements in cache.
 * @param preemption Controls cache swapping.
 *
 * @returns ::NC_NOERR for success
 * @author Ed Hartnett
 */
int
nc_set_var_chunk_cache_ints(int ncid, int varid, int size, int nelems,
                            int preemption)
{
   size_t real_size = H5D_CHUNK_CACHE_NBYTES_DEFAULT;
   size_t real_nelems = H5D_CHUNK_CACHE_NSLOTS_DEFAULT;
   float real_preemption = CHUNK_CACHE_PREEMPTION;

   if (size >= 0)
      real_size = ((size_t) size) * MEGABYTE;

   if (nelems >= 0)
      real_nelems = nelems;

   if (preemption >= 0)
      real_preemption = preemption / 100.;

   return NC4_set_var_chunk_cache(ncid, varid, real_size, real_nelems,
                                  real_preemption);
}

/**
 * @internal This is called by nc_get_var_chunk_cache(). Get chunk
 * cache size for a variable.
 *
 * @param ncid File ID.
 * @param varid Variable ID.
 * @param sizep Gets size in bytes of cache.
 * @param nelemsp Gets number of element slots in cache.
 * @param preemptionp Gets cache swapping setting.
 *
 * @returns ::NC_NOERR No error.
 * @returns ::NC_EBADID Bad ncid.
 * @returns ::NC_ENOTVAR Invalid variable ID.
 * @returns ::NC_ENOTNC4 Not a netCDF-4 file.
 * @author Ed Hartnett
 */
int
NC4_get_var_chunk_cache(int ncid, int varid, size_t *sizep,
                        size_t *nelemsp, float *preemptionp)
{
   NC *nc;
   NC_GRP_INFO_T *grp;
   NC_FILE_INFO_T *h5;
   NC_VAR_INFO_T *var;
   int retval;

   /* Find info for this file and group, and set pointer to each. */
   if ((retval = nc4_find_nc_grp_h5(ncid, &nc, &grp, &h5)))
      return retval;
   assert(nc && grp && h5);

   /* Find the var. */
   var = (NC_VAR_INFO_T*)ncindexith(grp->vars,varid);
   if(!var)
      return NC_ENOTVAR;
   assert(var && var->hdr.id == varid);

   /* Give the user what they want. */
   if (sizep)
      *sizep = var->chunk_cache_size;
   if (nelemsp)
      *nelemsp = var->chunk_cache_nelems;
   if (preemptionp)
      *preemptionp = var->chunk_cache_preemption;

   return NC_NOERR;
}

/**
 * @internal A wrapper for NC4_get_var_chunk_cache(), we need this
 * version for fortran.
 *
 * @param ncid File ID.
 * @param varid Variable ID.
 * @param sizep Gets size in bytes of cache.
 * @param nelemsp Gets number of element slots in cache.
 * @param preemptionp Gets cache swapping setting.
 *
 * @returns ::NC_NOERR No error.
 * @returns ::NC_EBADID Bad ncid.
 * @returns ::NC_ENOTVAR Invalid variable ID.
 * @returns ::NC_ENOTNC4 Not a netCDF-4 file.
 * @author Ed Hartnett
 */
int
nc_get_var_chunk_cache_ints(int ncid, int varid, int *sizep,
                            int *nelemsp, int *preemptionp)
{
   size_t real_size, real_nelems;
   float real_preemption;
   int ret;

   if ((ret = NC4_get_var_chunk_cache(ncid, varid, &real_size,
                                      &real_nelems, &real_preemption)))
      return ret;

   if (sizep)
      *sizep = real_size / MEGABYTE;
   if (nelemsp)
      *nelemsp = (int)real_nelems;
   if(preemptionp)
      *preemptionp = (int)(real_preemption * 100);

   return NC_NOERR;
}

/**
 * @internal Get all the information about a variable. Pass NULL for
 * whatever you don't care about. This is the internal function called
 * by nc_inq_var(), nc_inq_var_deflate(), nc_inq_var_fletcher32(),
 * nc_inq_var_chunking(), nc_inq_var_chunking_ints(),
 * nc_inq_var_fill(), nc_inq_var_endian(), nc_inq_var_filter(), and
 * nc_inq_var_szip().
 *
 * @param ncid File ID.
 * @param varid Variable ID.
 * @param name Gets name.
 * @param xtypep Gets type.
 * @param ndimsp Gets number of dims.
 * @param dimidsp Gets array of dim IDs.
 * @param nattsp Gets number of attributes.
 * @param shufflep Gets shuffle setting.
 * @param deflatep Gets deflate setting.
 * @param deflate_levelp Gets deflate level.
 * @param fletcher32p Gets fletcher32 setting.
 * @param contiguousp Gets contiguous setting.
 * @param chunksizesp Gets chunksizes.
 * @param no_fill Gets fill mode.
 * @param fill_valuep Gets fill value.
 * @param endiannessp Gets one of ::NC_ENDIAN_BIG ::NC_ENDIAN_LITTLE
 * ::NC_ENDIAN_NATIVE
 * @param idp Pointer to memory to store filter id.
 * @param nparamsp Pointer to memory to store filter parameter count.
 * @param params Pointer to vector of unsigned integers into which
 * to store filter parameters.
 *
 * @returns ::NC_NOERR No error.
 * @returns ::NC_EBADID Bad ncid.
 * @returns ::NC_ENOTVAR Bad varid.
 * @returns ::NC_ENOMEM Out of memory.
 * @returns ::NC_EINVAL Invalid input.
 * @author Ed Hartnett, Dennis Heimbigner
 */
int
NC4_inq_var_all(int ncid, int varid, char *name, nc_type *xtypep,
                int *ndimsp, int *dimidsp, int *nattsp,
                int *shufflep, int *deflatep, int *deflate_levelp,
                int *fletcher32p, int *contiguousp, size_t *chunksizesp,
                int *no_fill, void *fill_valuep, int *endiannessp,
                unsigned int* idp, size_t* nparamsp, unsigned int* params
   )
{
   NC *nc;
   NC_GRP_INFO_T *grp;
   NC_FILE_INFO_T *h5;
   NC_VAR_INFO_T *var;
   int d;
   int retval;

   LOG((2, "%s: ncid 0x%x varid %d", __func__, ncid, varid));

   /* Find info for this file and group, and set pointer to each. */
   if ((retval = nc4_find_nc_grp_h5(ncid, &nc, &grp, &h5)))
      return retval;

   assert(nc);
   assert(grp && h5);

   /* Walk through the list of vars, and return the info about the one
      with a matching varid. If the varid is -1, find the global
      atts and call it a day. */
   if (varid == NC_GLOBAL)
   {
      if (nattsp)
      {
         /* Do we need to read the atts? */
         if (grp->atts_not_read)
            if ((retval = nc4_read_atts(grp, NULL)))
               return retval;

         *nattsp = ncindexcount(grp->att);
      }
      return NC_NOERR;
   }

   /* Find the var. */
   var = (NC_VAR_INFO_T*)ncindexith(grp->vars,varid);
   if(!var)
      return NC_ENOTVAR;
   assert(var && var->hdr.id == varid);

   /* Copy the data to the user's data buffers. */
   if (name)
      strcpy(name, var->hdr.name);
   if (xtypep)
      *xtypep = var->type_info->hdr.id;
   if (ndimsp)
      *ndimsp = var->ndims;
   if (dimidsp)
      for (d = 0; d < var->ndims; d++)
         dimidsp[d] = var->dimids[d];
   if (nattsp)
   {
      if (var->atts_not_read)
         if ((retval = nc4_read_atts(grp, var)))
            return retval;
      *nattsp = ncindexcount(var->att);
   }

   /* Chunking stuff. */
   if (!var->contiguous && chunksizesp)
      for (d = 0; d < var->ndims; d++)
      {
         chunksizesp[d] = var->chunksizes[d];
         LOG((4, "chunksizesp[%d]=%d", d, chunksizesp[d]));
      }

   if (contiguousp)
      *contiguousp = var->contiguous ? NC_CONTIGUOUS : NC_CHUNKED;

   /* Filter stuff. */
   if (deflatep)
      *deflatep = (int)var->deflate;
   if (deflate_levelp)
      *deflate_levelp = var->deflate_level;
   if (shufflep)
      *shufflep = (int)var->shuffle;
   if (fletcher32p)
      *fletcher32p = (int)var->fletcher32;

   if (idp)
      *idp = var->filterid;
   if (nparamsp)
      *nparamsp = (var->params == NULL ? 0 : var->nparams);
   if (params && var->params != NULL)
      memcpy(params,var->params,var->nparams*sizeof(unsigned int));

   /* Fill value stuff. */
   if (no_fill)
      *no_fill = (int)var->no_fill;

   /* Don't do a thing with fill_valuep if no_fill mode is set for
    * this var, or if fill_valuep is NULL. */
   if (!var->no_fill && fill_valuep)
   {
      /* Do we have a fill value for this var? */
      if (var->fill_value)
      {
         if (var->type_info->nc_type_class == NC_STRING)
         {
            assert(*(char **)var->fill_value);
            /* This will allocate memeory and copy the string. */
            if (!(*(char **)fill_valuep = strdup(*(char **)var->fill_value)))
            {
               free(*(char **)fill_valuep);
               return NC_ENOMEM;
            }
         }
         else
         {
            assert(var->type_info->size);
            memcpy(fill_valuep, var->fill_value, var->type_info->size);
         }
      }
      else
      {
         if (var->type_info->nc_type_class == NC_STRING)
         {
            if (!(*(char **)fill_valuep = calloc(1, sizeof(char *))))
               return NC_ENOMEM;

            if ((retval = nc4_get_default_fill_value(var->type_info, (char **)fill_valuep)))
            {
               free(*(char **)fill_valuep);
               return retval;
            }
         }
         else
         {
            if ((retval = nc4_get_default_fill_value(var->type_info, fill_valuep)))
               return retval;
         }
      }
   }

   /* Does the user want the endianness of this variable? */
   if (endiannessp)
      *endiannessp = var->type_info->endianness;

   return NC_NOERR;
}

/**
 * @internal Inquire about chunking settings for a var. This is used
 * by the fortran API.
 *
 * @param ncid File ID.
 * @param varid Variable ID.
 * @param contiguousp Gets contiguous setting.
 * @param chunksizesp Gets chunksizes.
 *
 * @returns ::NC_NOERR No error.
 * @returns ::NC_EBADID Bad ncid.
 * @returns ::NC_ENOTVAR Invalid variable ID.
 * @returns ::NC_EINVAL Invalid input
 * @returns ::NC_ENOMEM Out of memory.
 * @author Ed Hartnett
 */
int
nc_inq_var_chunking_ints(int ncid, int varid, int *contiguousp, int *chunksizesp)
{
   NC_VAR_INFO_T *var;
   size_t *cs = NULL;
   int i, retval;

   /* Get pointer to the var. */
   if ((retval = nc4_find_grp_h5_var(ncid, varid, NULL, NULL, &var)))
      return retval;
   assert(var);

   /* Allocate space for the size_t copy of the chunksizes array. */
   if (var->ndims)
      if (!(cs = malloc(var->ndims * sizeof(size_t))))
         return NC_ENOMEM;

   /* Call the netcdf-4 version directly. */
   retval = NC4_inq_var_all(ncid, varid, NULL, NULL, NULL, NULL, NULL,
                            NULL, NULL, NULL, NULL, contiguousp, cs, NULL,
                            NULL, NULL, NULL, NULL, NULL);

   /* Copy from size_t array. */
   if (!retval && chunksizesp && var->contiguous == NC_CHUNKED)
   {
      for (i = 0; i < var->ndims; i++)
      {
         chunksizesp[i] = (int)cs[i];
         if (cs[i] > NC_MAX_INT)
            retval = NC_ERANGE;
      }
   }

   if (var->ndims)
      free(cs);
   return retval;
}

/**
 * @internal Find the ID of a variable, from the name. This function
 * is called by nc_inq_varid().
 *
 * @param ncid File ID.
 * @param name Name of the variable.
 * @param varidp Gets variable ID.

 * @returns ::NC_NOERR No error.
 * @returns ::NC_EBADID Bad ncid.
 * @returns ::NC_ENOTVAR Bad variable ID.
 */
int
NC4_inq_varid(int ncid, const char *name, int *varidp)
{
   NC *nc;
   NC_GRP_INFO_T *grp;
   NC_VAR_INFO_T *var;
   char norm_name[NC_MAX_NAME + 1];
   int retval;

   if (!name)
      return NC_EINVAL;
   if (!varidp)
      return NC_NOERR;

   LOG((2, "%s: ncid 0x%x name %s", __func__, ncid, name));

   /* Find info for this file and group, and set pointer to each. */
   if ((retval = nc4_find_nc_grp_h5(ncid, &nc, &grp, NULL)))
      return retval;

   /* Normalize name. */
   if ((retval = nc4_normalize_name(name, norm_name)))
      return retval;

   /* Find var of this name. */
   var = (NC_VAR_INFO_T*)ncindexlookup(grp->vars,norm_name);
   if(var)
   {
      *varidp = var->hdr.id;
      return NC_NOERR;
   }
   return NC_ENOTVAR;
}

/**
 * @internal
 *
 * This function will change the parallel access of a variable from
 * independent to collective.
 *
 * @param ncid File ID.
 * @param varid Variable ID.
 * @param par_access NC_COLLECTIVE or NC_INDEPENDENT.
 *
 * @returns ::NC_NOERR No error.
 * @returns ::NC_EBADID Invalid ncid passed.
 * @returns ::NC_ENOTVAR Invalid varid passed.
 * @returns ::NC_ENOPAR LFile was not opened with nc_open_par/nc_create_var.
 * @returns ::NC_EINVAL Invalid par_access specified.
 * @returns ::NC_NOERR for success
 * @author Ed Hartnett, Dennis Heimbigner
 */
int
NC4_var_par_access(int ncid, int varid, int par_access)
{
#ifndef USE_PARALLEL4
   return NC_ENOPAR;
#else
   NC *nc;
   NC_GRP_INFO_T *grp;
   NC_FILE_INFO_T *h5;
   NC_VAR_INFO_T *var;
   int retval;

   LOG((1, "%s: ncid 0x%x varid %d par_access %d", __func__, ncid,
        varid, par_access));

   if (par_access != NC_INDEPENDENT && par_access != NC_COLLECTIVE)
      return NC_EINVAL;

   /* Find info for this file and group, and set pointer to each. */
   if ((retval = nc4_find_nc_grp_h5(ncid, &nc, &grp, &h5)))
      return retval;

   /* This function only for files opened with nc_open_par or nc_create_par. */
   if (!h5->parallel)
      return NC_ENOPAR;

   /* Find the var, and set its preference. */
   var = (NC_VAR_INFO_T*)ncindexith(grp->vars,varid);
   if (!var) return NC_ENOTVAR;
   assert(var->hdr.id == varid);

   if (par_access)
      var->parallel_access = NC_COLLECTIVE;
   else
      var->parallel_access = NC_INDEPENDENT;
   return NC_NOERR;
#endif /* USE_PARALLEL4 */
}

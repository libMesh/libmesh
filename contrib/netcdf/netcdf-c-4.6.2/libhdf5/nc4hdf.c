/* Copyright 2003, University Corporation for Atmospheric
 * Research. See the COPYRIGHT file for copying and redistribution
 * conditions. */
/**
 * @file
 * @internal This file is part of netcdf-4, a netCDF-like interface
 * for HDF5, or a HDF5 backend for netCDF, depending on your point of
 * view.
 *
 * This file contains functions internal to the netcdf4 library. None of
 * the functions in this file are exposed in the exetnal API. These
 * functions handle the HDF interface.
 *
 * @author Ed Hartnett, Dennis Heimbigner, Ward Fisher
 */

#include "config.h"
#include "hdf5internal.h"
#include <math.h>

#ifdef HAVE_INTTYPES_H
#define __STDC_FORMAT_MACROS
#include <inttypes.h>
#endif

#define NC_HDF5_MAX_NAME 1024 /**< @internal Max size of HDF5 name. */

/**
 * @internal Flag attributes in a linked list as dirty.
 *
 * @param attlist List of attributes, may be NULL.
 *
 * @return NC_NOERR No error.
 * @author Dennis Heimbigner
 */
static int
flag_atts_dirty(NCindex *attlist) {

   NC_ATT_INFO_T *att = NULL;
   int i;

   if(attlist == NULL) {
      return NC_NOERR;
   }

   for(i=0;i<ncindexsize(attlist);i++) {
      att = (NC_ATT_INFO_T*)ncindexith(attlist,i);
      if(att == NULL) continue;
      att->dirty = NC_TRUE;
   }

   return NC_NOERR;
}

/**
 * @internal This function is needed to handle one special case: what
 * if the user defines a dim, writes metadata, then goes back into
 * define mode and adds a coordinate var for the already existing
 * dim. In that case, I need to recreate the dim's dimension scale
 * dataset, and then I need to go to every var in the file which uses
 * that dimension, and attach the new dimension scale.
 *
 * @param grp Pointer to group info struct.
 * @param dimid Dimension ID.
 * @param dimscaleid HDF5 dimension scale ID.
 *
 * @returns NC_NOERR No error.
 * @returns NC_EHDFERR HDF5 returned an error.
 * @author Ed Hartnett
 */
int
rec_reattach_scales(NC_GRP_INFO_T *grp, int dimid, hid_t dimscaleid)
{
   NC_VAR_INFO_T *var;
   NC_GRP_INFO_T *child_grp;
   int d, i;
   int retval;

   assert(grp && grp->hdr.name && dimid >= 0 && dimscaleid >= 0);
   LOG((3, "%s: grp->hdr.name %s", __func__, grp->hdr.name));

   /* If there are any child groups, attach dimscale there, if needed. */
   for(i=0;i<ncindexsize(grp->children);i++) {
      child_grp = (NC_GRP_INFO_T*)ncindexith(grp->children,i);
      if(child_grp == NULL) continue;
      if ((retval = rec_reattach_scales(child_grp, dimid, dimscaleid)))
         return retval;
   }

   /* Find any vars that use this dimension id. */
   for(i=0;i<ncindexsize(grp->vars);i++)
   {
      var = (NC_VAR_INFO_T*)ncindexith(grp->vars,i);
      if(var == NULL) continue;
      for (d = 0; d < var->ndims; d++)
         if (var->dimids[d] == dimid && !var->dimscale)
         {
            LOG((2, "%s: attaching scale for dimid %d to var %s",
                 __func__, var->dimids[d], var->hdr.name));
            if (var->created)
            {
               if (H5DSattach_scale(var->hdf_datasetid, dimscaleid, d) < 0)
                  return NC_EHDFERR;
               var->dimscale_attached[d] = NC_TRUE;
            }
         }
   }
   return NC_NOERR;
}

/**
 * @internal This function is needed to handle one special case: what
 * if the user defines a dim, writes metadata, then goes back into
 * define mode and adds a coordinate var for the already existing
 * dim. In that case, I need to recreate the dim's dimension scale
 * dataset, and then I need to go to every var in the file which uses
 * that dimension, and attach the new dimension scale.
 *
 * @param grp Pointer to group info struct.
 * @param dimid Dimension ID.
 * @param dimscaleid HDF5 dimension scale ID.
 *
 * @returns NC_NOERR No error.
 * @returns NC_EHDFERR HDF5 returned an error.
 * @author Ed Hartnett
 */
int
rec_detach_scales(NC_GRP_INFO_T *grp, int dimid, hid_t dimscaleid)
{
   NC_VAR_INFO_T *var;
   NC_GRP_INFO_T *child_grp;
   int d, i;
   int retval;

   assert(grp && grp->hdr.name && dimid >= 0 && dimscaleid >= 0);
   LOG((3, "%s: grp->hdr.name %s", __func__, grp->hdr.name));

   /* If there are any child groups, detach dimscale there, if needed. */
   for(i=0;i<ncindexsize(grp->children);i++) {
      child_grp = (NC_GRP_INFO_T*)ncindexith(grp->children,i);
      if(child_grp == NULL) continue;
      if ((retval = rec_detach_scales(child_grp, dimid, dimscaleid)))
         return retval;
   }

   /* Find any vars that use this dimension id. */
   for(i=0;i<ncindexsize(grp->vars);i++) {
      var = (NC_VAR_INFO_T*)ncindexith(grp->vars,i);
      if(var == NULL) continue;
      for (d = 0; d < var->ndims; d++)
         if (var->dimids[d] == dimid && !var->dimscale)
         {
            LOG((2, "%s: detaching scale for dimid %d to var %s",
                 __func__, var->dimids[d], var->hdr.name));
            if (var->created)
               if (var->dimscale_attached && var->dimscale_attached[d])
               {
                  if (H5DSdetach_scale(var->hdf_datasetid, dimscaleid, d) < 0)
                     return NC_EHDFERR;
                  var->dimscale_attached[d] = NC_FALSE;
               }
         }
   }
   return NC_NOERR;
}

/**
 * @internal Open a HDF5 dataset and leave it open.
 *
 * @param grp Pointer to group info struct.
 * @param varid Variable ID.
 * @param dataset Pointer that gets the HDF5 dataset ID.
 *
 * @returns NC_NOERR No error.
 * @returns NC_EHDFERR HDF5 returned an error.
 * @author Ed Hartnett
 */
int
nc4_open_var_grp2(NC_GRP_INFO_T *grp, int varid, hid_t *dataset)
{
   NC_VAR_INFO_T *var;

   assert(grp && grp->format_grp_info && dataset);

   /* Find the requested varid. */
   var = (NC_VAR_INFO_T*)ncindexith(grp->vars,varid);
   if (!var) return NC_ENOTVAR;
   assert(var->hdr.id == varid);

   /* Open this dataset if necessary. */
   if (!var->hdf_datasetid)
   {
      NC_HDF5_GRP_INFO_T *hdf5_grp;
      hdf5_grp = (NC_HDF5_GRP_INFO_T *)grp->format_grp_info;

      if ((var->hdf_datasetid = H5Dopen2(hdf5_grp->hdf_grpid,
                                         var->hdr.name, H5P_DEFAULT)) < 0)
         return NC_ENOTVAR;
   }

   *dataset = var->hdf_datasetid;

   return NC_NOERR;
}

/**
 * @internal Get the default fill value for an atomic type. Memory for
 * fill_value must already be allocated, or you are DOOMED!
 *
 * @param type_info Pointer to type info struct.
 * @param fill_value Pointer that gets the default fill value.
 *
 * @returns NC_NOERR No error.
 * @returns NC_EINVAL Can't find atomic type.
 * @author Ed Hartnett
 */
int
nc4_get_default_fill_value(const NC_TYPE_INFO_T *type_info, void *fill_value)
{
   switch (type_info->hdr.id)
   {
   case NC_CHAR:
      *(char *)fill_value = NC_FILL_CHAR;
      break;

   case NC_STRING:
      *(char **)fill_value = strdup(NC_FILL_STRING);
      break;

   case NC_BYTE:
      *(signed char *)fill_value = NC_FILL_BYTE;
      break;

   case NC_SHORT:
      *(short *)fill_value = NC_FILL_SHORT;
      break;

   case NC_INT:
      *(int *)fill_value = NC_FILL_INT;
      break;

   case NC_UBYTE:
      *(unsigned char *)fill_value = NC_FILL_UBYTE;
      break;

   case NC_USHORT:
      *(unsigned short *)fill_value = NC_FILL_USHORT;
      break;

   case NC_UINT:
      *(unsigned int *)fill_value = NC_FILL_UINT;
      break;

   case NC_INT64:
      *(long long *)fill_value = NC_FILL_INT64;
      break;

   case NC_UINT64:
      *(unsigned long long *)fill_value = NC_FILL_UINT64;
      break;

   case NC_FLOAT:
      *(float *)fill_value = NC_FILL_FLOAT;
      break;

   case NC_DOUBLE:
      *(double *)fill_value = NC_FILL_DOUBLE;
      break;

   default:
      return NC_EINVAL;
   }

   return NC_NOERR;
}

/**
 * @internal What fill value should be used for a variable?
 *
 * @param h5 Pointer to HDF5 file info struct.
 * @param var Pointer to variable info struct.
 * @param fillp Pointer that gets pointer to fill value.
 *
 * @returns NC_NOERR No error.
 * @returns NC_ENOMEM Out of memory.
 * @author Ed Hartnett
 */
int
nc4_get_fill_value(NC_FILE_INFO_T *h5, NC_VAR_INFO_T *var, void **fillp)
{
   size_t size;
   int retval;

   /* Find out how much space we need for this type's fill value. */
   if (var->type_info->nc_type_class == NC_VLEN)
      size = sizeof(nc_vlen_t);
   else if (var->type_info->nc_type_class == NC_STRING)
      size = sizeof(char *);
   else
   {
      if ((retval = nc4_get_typelen_mem(h5, var->type_info->hdr.id, &size)))
         return retval;
   }
   assert(size);

   /* Allocate the space. */
   if (!((*fillp) = calloc(1, size)))
      return NC_ENOMEM;

   /* If the user has set a fill_value for this var, use, otherwise
    * find the default fill value. */
   if (var->fill_value)
   {
      LOG((4, "Found a fill value for var %s", var->hdr.name));
      if (var->type_info->nc_type_class == NC_VLEN)
      {
         nc_vlen_t *in_vlen = (nc_vlen_t *)(var->fill_value), *fv_vlen = (nc_vlen_t *)(*fillp);
	 size_t basetypesize = 0;

	 if((retval=nc4_get_typelen_mem(h5, var->type_info->u.v.base_nc_typeid, &basetypesize)))
	     return retval;

         fv_vlen->len = in_vlen->len;
         if (!(fv_vlen->p = malloc(basetypesize * in_vlen->len)))
         {
            free(*fillp);
            *fillp = NULL;
            return NC_ENOMEM;
         }
         memcpy(fv_vlen->p, in_vlen->p, in_vlen->len * basetypesize);
      }
      else if (var->type_info->nc_type_class == NC_STRING)
      {
         if (*(char **)var->fill_value)
            if (!(**(char ***)fillp = strdup(*(char **)var->fill_value)))
            {
               free(*fillp);
               *fillp = NULL;
               return NC_ENOMEM;
            }
      }
      else
         memcpy((*fillp), var->fill_value, size);
   }
   else
   {
      if (nc4_get_default_fill_value(var->type_info, *fillp))
      {
         /* Note: release memory, but don't return error on failure */
         free(*fillp);
         *fillp = NULL;
      }
   }

   return NC_NOERR;
}

/**
 * @internal Given a netcdf type, return appropriate HDF typeid.  (All
 * hdf_typeid's returned from this routine must be H5Tclosed by the
 * caller).
 *
 * @param h5 Pointer to HDF5 file info struct.
 * @param xtype NetCDF type ID.
 * @param hdf_typeid Pointer that gets the HDF5 type ID.
 * @param endianness Desired endianness in HDF5 type.
 *
 * @return NC_NOERR No error.
 * @return NC_ECHAR Conversions of NC_CHAR forbidden.
 * @return NC_EVARMETA HDF5 returning error creating datatype.
 * @return NC_EHDFERR HDF5 returning error.
 * @return NC_EBADTYE Type not found.
 * @author Ed Hartnett
 */
int
nc4_get_hdf_typeid(NC_FILE_INFO_T *h5, nc_type xtype,
                   hid_t *hdf_typeid, int endianness)
{
   NC_TYPE_INFO_T *type;
   hid_t typeid = 0;
   int retval = NC_NOERR;

   assert(hdf_typeid && h5);

   *hdf_typeid = -1;

   /* Determine an appropriate HDF5 datatype */
   if (xtype == NC_NAT)
      return NC_EBADTYPE;
   else if (xtype == NC_CHAR || xtype == NC_STRING)
   {
      /* NC_CHAR & NC_STRING types create a new HDF5 datatype */
      if (xtype == NC_CHAR)
      {
         if ((typeid = H5Tcopy(H5T_C_S1)) < 0)
            return NC_EHDFERR;
         if (H5Tset_strpad(typeid, H5T_STR_NULLTERM) < 0)
            BAIL(NC_EVARMETA);
         if(H5Tset_cset(typeid, H5T_CSET_ASCII) < 0)
            BAIL(NC_EVARMETA);

         /* Take ownership of the newly created HDF5 datatype */
         *hdf_typeid = typeid;
         typeid = 0;
      }
      else
      {
         if ((typeid = H5Tcopy(H5T_C_S1)) < 0)
            return NC_EHDFERR;
         if (H5Tset_size(typeid, H5T_VARIABLE) < 0)
            BAIL(NC_EVARMETA);
         if(H5Tset_cset(typeid, H5T_CSET_UTF8) < 0)
            BAIL(NC_EVARMETA);

         /* Take ownership of the newly created HDF5 datatype */
         *hdf_typeid = typeid;
         typeid = 0;
      }
   }
   else
   {
      /* All other types use an existing HDF5 datatype */
      switch (xtype)
      {
      case NC_BYTE: /* signed 1 byte integer */
         if (endianness == NC_ENDIAN_LITTLE)
            typeid = H5T_STD_I8LE;
         else if (endianness == NC_ENDIAN_BIG)
            typeid = H5T_STD_I8BE;
         else
            typeid = H5T_NATIVE_SCHAR;
         break;

      case NC_SHORT: /* signed 2 byte integer */
         if (endianness == NC_ENDIAN_LITTLE)
            typeid = H5T_STD_I16LE;
         else if (endianness == NC_ENDIAN_BIG)
            typeid = H5T_STD_I16BE;
         else
            typeid = H5T_NATIVE_SHORT;
         break;

      case NC_INT:
         if (endianness == NC_ENDIAN_LITTLE)
            typeid = H5T_STD_I32LE;
         else if (endianness == NC_ENDIAN_BIG)
            typeid = H5T_STD_I32BE;
         else
            typeid = H5T_NATIVE_INT;
         break;

      case NC_UBYTE:
         if (endianness == NC_ENDIAN_LITTLE)
            typeid = H5T_STD_U8LE;
         else if (endianness == NC_ENDIAN_BIG)
            typeid = H5T_STD_U8BE;
         else
            typeid = H5T_NATIVE_UCHAR;
         break;

      case NC_USHORT:
         if (endianness == NC_ENDIAN_LITTLE)
            typeid = H5T_STD_U16LE;
         else if (endianness == NC_ENDIAN_BIG)
            typeid = H5T_STD_U16BE;
         else
            typeid = H5T_NATIVE_USHORT;
         break;

      case NC_UINT:
         if (endianness == NC_ENDIAN_LITTLE)
            typeid = H5T_STD_U32LE;
         else if (endianness == NC_ENDIAN_BIG)
            typeid = H5T_STD_U32BE;
         else
            typeid = H5T_NATIVE_UINT;
         break;

      case NC_INT64:
         if (endianness == NC_ENDIAN_LITTLE)
            typeid = H5T_STD_I64LE;
         else if (endianness == NC_ENDIAN_BIG)
            typeid = H5T_STD_I64BE;
         else
            typeid = H5T_NATIVE_LLONG;
         break;

      case NC_UINT64:
         if (endianness == NC_ENDIAN_LITTLE)
            typeid = H5T_STD_U64LE;
         else if (endianness == NC_ENDIAN_BIG)
            typeid = H5T_STD_U64BE;
         else
            typeid = H5T_NATIVE_ULLONG;
         break;

      case NC_FLOAT:
         if (endianness == NC_ENDIAN_LITTLE)
            typeid = H5T_IEEE_F32LE;
         else if (endianness == NC_ENDIAN_BIG)
            typeid = H5T_IEEE_F32BE;
         else
            typeid = H5T_NATIVE_FLOAT;
         break;

      case NC_DOUBLE:
         if (endianness == NC_ENDIAN_LITTLE)
            typeid = H5T_IEEE_F64LE;
         else if (endianness == NC_ENDIAN_BIG)
            typeid = H5T_IEEE_F64BE;
         else
            typeid = H5T_NATIVE_DOUBLE;
         break;

      default:
         /* Maybe this is a user defined type? */
         if (nc4_find_type(h5, xtype, &type))
            return NC_EBADTYPE;
         if (!type)
            return NC_EBADTYPE;
         typeid = type->hdf_typeid;
         break;
      }
      assert(typeid);

      /* Copy the HDF5 datatype, so the function operates uniformly */
      if ((*hdf_typeid = H5Tcopy(typeid)) < 0)
         return NC_EHDFERR;
      typeid = 0;
   }
   assert(*hdf_typeid != -1);

exit:
   if (typeid > 0 && H5Tclose(typeid) < 0)
      BAIL2(NC_EHDFERR);
   return retval;
}

/**
 * @internal Write an attribute.
 *
 * @param grp Pointer to group info struct.
 * @param varid Variable ID or NC_GLOBAL.
 * @param att Pointer to attribute info struct.
 *
 * @returns ::NC_NOERR No error.
 * @returns ::NC_ENOTVAR Variable not found.
 * @returns ::NC_EPERM Read-only file.
 * @returns ::NC_EHDFERR HDF5 returned error.
 * @returns ::NC_EATTMETA HDF5 returned error with attribute calls.
 * @author Ed Hartnett
 */
static int
put_att_grpa(NC_GRP_INFO_T *grp, int varid, NC_ATT_INFO_T *att)
{
   NC_HDF5_GRP_INFO_T *hdf5_grp;
   hid_t datasetid = 0, locid;
   hid_t attid = 0, spaceid = 0, file_typeid = 0;
   hid_t existing_att_typeid = 0, existing_attid = 0, existing_spaceid = 0;
   hsize_t dims[1]; /* netcdf attributes always 1-D. */
   htri_t attr_exists;
   int reuse_att = 0; /* Will be true if we can re-use an existing att. */
   void *data;
   int phoney_data = 99;
   int retval = NC_NOERR;

   assert(att->hdr.name && grp && grp->format_grp_info);
   LOG((3, "%s: varid %d att->hdr.id %d att->hdr.name %s att->nc_typeid %d "
        "att->len %d", __func__, varid, att->hdr.id, att->hdr.name,
        att->nc_typeid, att->len));

   /* Get HDF5-specific group info. */
   hdf5_grp = (NC_HDF5_GRP_INFO_T *)grp->format_grp_info;

   /* If the file is read-only, return an error. */
   if (grp->nc4_info->no_write)
      BAIL(NC_EPERM);

   /* Get the hid to attach the attribute to, or read it from. */
   if (varid == NC_GLOBAL)
      locid = hdf5_grp->hdf_grpid;
   else
   {
      if ((retval = nc4_open_var_grp2(grp, varid, &datasetid)))
         BAIL(retval);
      locid = datasetid;
   }

   /* Get the length ready, and find the HDF type we'll be
    * writing. */
   dims[0] = att->len;
   if ((retval = nc4_get_hdf_typeid(grp->nc4_info, att->nc_typeid,
                                    &file_typeid, 0)))
      BAIL(retval);

   /* Even if the length is zero, HDF5 won't let me write with a
    * NULL pointer. So if the length of the att is zero, point to
    * some phoney data (which won't be written anyway.)*/
   if (!dims[0])
      data = &phoney_data;
   else if (att->data)
      data = att->data;
   else if (att->stdata)
      data = att->stdata;
   else
      data = att->vldata;

   /* NC_CHAR types require some extra work. The space ID is set to
    * scalar, and the type is told how long the string is. If it's
    * really zero length, set the size to 1. (The fact that it's
    * really zero will be marked by the NULL dataspace, but HDF5
    * doesn't allow me to set the size of the type to zero.)*/
   if (att->nc_typeid == NC_CHAR)
   {
      size_t string_size = dims[0];
      if (!string_size)
      {
         string_size = 1;
         if ((spaceid = H5Screate(H5S_NULL)) < 0)
            BAIL(NC_EATTMETA);
      }
      else
      {
         if ((spaceid = H5Screate(H5S_SCALAR)) < 0)
            BAIL(NC_EATTMETA);
      }
      if (H5Tset_size(file_typeid, string_size) < 0)
         BAIL(NC_EATTMETA);
      if (H5Tset_strpad(file_typeid, H5T_STR_NULLTERM) < 0)
         BAIL(NC_EATTMETA);
   }
   else
   {
      if (!att->len)
      {
         if ((spaceid = H5Screate(H5S_NULL)) < 0)
            BAIL(NC_EATTMETA);
      }
      else
      {
         if ((spaceid = H5Screate_simple(1, dims, NULL)) < 0)
            BAIL(NC_EATTMETA);
      }
   }

   /* Does the att exists already? */
   if ((attr_exists = H5Aexists(locid, att->hdr.name)) < 0)
      BAIL(NC_EHDFERR);
   if (attr_exists)
   {
      hssize_t npoints;

      /* Open the attribute. */
      if ((existing_attid = H5Aopen(locid, att->hdr.name, H5P_DEFAULT)) < 0)
         BAIL(NC_EATTMETA);

      /* Find the type of the existing attribute. */
      if ((existing_att_typeid = H5Aget_type(existing_attid)) < 0)
         BAIL(NC_EATTMETA);

      /* How big is the attribute? */
      if ((existing_spaceid = H5Aget_space(existing_attid)) < 0)
         BAIL(NC_EATTMETA);
      if ((npoints = H5Sget_simple_extent_npoints(existing_spaceid)) < 0)
         BAIL(NC_EATTMETA);

      /* Delete the attribute. */
      if (file_typeid != existing_att_typeid || npoints != att->len)
      {
         if (H5Adelete(locid, att->hdr.name) < 0)
            BAIL(NC_EHDFERR);
      }
      else
      {
         reuse_att++;
      }
   }

   /* Create the attribute. */
   if ((attid = H5Acreate(locid, att->hdr.name, file_typeid, spaceid,
                          H5P_DEFAULT)) < 0)
      BAIL(NC_EATTMETA);

   /* Write the values, (even if length is zero). */
   if (H5Awrite(attid, file_typeid, data) < 0)
      BAIL(NC_EATTMETA);

exit:
   if (file_typeid && H5Tclose(file_typeid))
      BAIL2(NC_EHDFERR);
   if (attid > 0 && H5Aclose(attid) < 0)
      BAIL2(NC_EHDFERR);
   if (existing_att_typeid && H5Tclose(existing_att_typeid))
      BAIL2(NC_EHDFERR);
   if (existing_attid > 0 && H5Aclose(existing_attid) < 0)
      BAIL2(NC_EHDFERR);
   if (spaceid > 0 && H5Sclose(spaceid) < 0)
      BAIL2(NC_EHDFERR);
   if (existing_spaceid > 0 && H5Sclose(existing_spaceid) < 0)
      BAIL2(NC_EHDFERR);
   return retval;
}

/**
 * @internal Write all the dirty atts in an attlist.
 *
 * @param attlist Pointer to the list if attributes.
 * @param varid Variable ID.
 * @param grp Pointer to group info struct.
 *
 * @returns NC_NOERR No error.
 * @returns NC_EHDFERR HDF5 returned an error.
 * @author Ed Hartnett
 */
static int
write_attlist(NCindex *attlist, int varid, NC_GRP_INFO_T *grp)
{
   NC_ATT_INFO_T *att;
   int retval;
   int i;

   for(i = 0; i < ncindexsize(attlist); i++)
   {
      att = (NC_ATT_INFO_T *)ncindexith(attlist, i);
      assert(att);
      if (att->dirty)
      {
         LOG((4, "%s: writing att %s to varid %d", __func__, att->hdr.name, varid));
         if ((retval = put_att_grpa(grp, varid, att)))
            return retval;
         att->dirty = NC_FALSE;
         att->created = NC_TRUE;
      }
   }
   return NC_NOERR;
}

/**
 * @internal HDF5 dimension scales cannot themselves have scales
 * attached. This is a problem for multidimensional coordinate
 * variables. So this function writes a special attribute for such a
 * variable, which has the ids of all the dimensions for that
 * coordinate variable.
 *
 * @param var Pointer to var info struct.
 *
 * @returns NC_NOERR No error.
 * @returns NC_EHDFERR HDF5 returned an error.
 * @author Ed Hartnett
 */
static int
write_coord_dimids(NC_VAR_INFO_T *var)
{
   hsize_t coords_len[1];
   hid_t c_spaceid = -1, c_attid = -1;
   int retval = NC_NOERR;

   /* Set up space for attribute. */
   coords_len[0] = var->ndims;
   if ((c_spaceid = H5Screate_simple(1, coords_len, coords_len)) < 0)
      BAIL(NC_EHDFERR);

   /* Create the attribute. */
   if ((c_attid = H5Acreate(var->hdf_datasetid, COORDINATES, H5T_NATIVE_INT,
                            c_spaceid, H5P_DEFAULT)) < 0)
      BAIL(NC_EHDFERR);

   /* Write our attribute. */
   if (H5Awrite(c_attid, H5T_NATIVE_INT, var->dimids) < 0)
      BAIL(NC_EHDFERR);

exit:
   if (c_spaceid >= 0 && H5Sclose(c_spaceid) < 0)
      BAIL2(NC_EHDFERR);
   if (c_attid >= 0 && H5Aclose(c_attid) < 0)
      BAIL2(NC_EHDFERR);
   return retval;
}

/**
 * @internal Write a special attribute for the netCDF-4 dimension ID.
 *
 * @param datasetid HDF5 datasset ID.
 * @param dimid NetCDF dimension ID.
 *
 * @returns NC_NOERR No error.
 * @returns NC_EHDFERR HDF5 returned an error.
 * @author Ed Hartnett
 */
static int
write_netcdf4_dimid(hid_t datasetid, int dimid)
{
   hid_t dimid_spaceid, dimid_attid;
   htri_t attr_exists;

   /* Create the space. */
   if ((dimid_spaceid = H5Screate(H5S_SCALAR)) < 0)
      return NC_EHDFERR;

   /* Does the attribute already exist? If so, don't try to create it. */
   if ((attr_exists = H5Aexists(datasetid, NC_DIMID_ATT_NAME)) < 0)
      return NC_EHDFERR;
   if (attr_exists)
      dimid_attid = H5Aopen_by_name(datasetid, ".", NC_DIMID_ATT_NAME,
                                    H5P_DEFAULT, H5P_DEFAULT);
   else
      /* Create the attribute if needed. */
      dimid_attid = H5Acreate(datasetid, NC_DIMID_ATT_NAME,
                              H5T_NATIVE_INT, dimid_spaceid, H5P_DEFAULT);
   if (dimid_attid  < 0)
      return NC_EHDFERR;


   /* Write it. */
   LOG((4, "%s: writing secret dimid %d", __func__, dimid));
   if (H5Awrite(dimid_attid, H5T_NATIVE_INT, &dimid) < 0)
      return NC_EHDFERR;

   /* Close stuff*/
   if (H5Sclose(dimid_spaceid) < 0)
      return NC_EHDFERR;
   if (H5Aclose(dimid_attid) < 0)
      return NC_EHDFERR;

   return NC_NOERR;
}

/**
 * @internal This function creates the HDF5 dataset for a variable.
 *
 * @param grp Pointer to group info struct.
 * @param var Pointer to variable info struct.
 * @param write_dimid True to write dimid.
 *
 * @return ::NC_NOERR
 * @returns NC_ECHAR Conversions of NC_CHAR forbidden.
 * @returns NC_EVARMETA HDF5 returning error creating datatype.
 * @returns NC_EHDFERR HDF5 returning error.
 * @returns NC_EBADTYE Type not found.
 * @author Ed Hartnett
 */
static int
var_create_dataset(NC_GRP_INFO_T *grp, NC_VAR_INFO_T *var, nc_bool_t write_dimid)
{
   NC_HDF5_GRP_INFO_T *hdf5_grp;
   hid_t plistid = 0, access_plistid = 0, typeid = 0, spaceid = 0;
   hsize_t chunksize[H5S_MAX_RANK], dimsize[H5S_MAX_RANK], maxdimsize[H5S_MAX_RANK];
   int d;
   void *fillp = NULL;
   NC_DIM_INFO_T *dim = NULL;
   char *name_to_use;
   int retval;

   assert(grp && grp->format_grp_info && var);

   LOG((3, "%s:: name %s", __func__, var->hdr.name));

   /* Get HDF5-specific group info. */
   hdf5_grp = (NC_HDF5_GRP_INFO_T *)grp->format_grp_info;

   /* Scalar or not, we need a creation property list. */
   if ((plistid = H5Pcreate(H5P_DATASET_CREATE)) < 0)
      BAIL(NC_EHDFERR);
   if ((access_plistid = H5Pcreate(H5P_DATASET_ACCESS)) < 0)
      BAIL(NC_EHDFERR);

   /* Turn off object tracking times in HDF5. */
   if (H5Pset_obj_track_times(plistid, 0) < 0)
      BAIL(NC_EHDFERR);

   /* Find the HDF5 type of the dataset. */
   if ((retval = nc4_get_hdf_typeid(grp->nc4_info, var->type_info->hdr.id, &typeid,
                                    var->type_info->endianness)))
      BAIL(retval);

   /* Figure out what fill value to set, if any. */
   if (var->no_fill)
   {
      /* Required to truly turn HDF5 fill values off */
      if (H5Pset_fill_time(plistid, H5D_FILL_TIME_NEVER) < 0)
         BAIL(NC_EHDFERR);
   }
   else
   {
      if ((retval = nc4_get_fill_value(grp->nc4_info, var, &fillp)))
         BAIL(retval);

      /* If there is a fill value, set it. */
      if (fillp)
      {
         if (var->type_info->nc_type_class == NC_STRING)
         {
            if (H5Pset_fill_value(plistid, typeid, fillp) < 0)
               BAIL(NC_EHDFERR);
         }
         else
         {
            /* The fill value set in HDF5 must always be presented as
             * a native type, even if the endianness for this dataset
             * is non-native. HDF5 will translate the fill value to
             * the target endiannesss. */
            hid_t fill_typeid = 0;

            if ((retval = nc4_get_hdf_typeid(grp->nc4_info, var->type_info->hdr.id, &fill_typeid,
                                             NC_ENDIAN_NATIVE)))
               BAIL(retval);
            if (H5Pset_fill_value(plistid, fill_typeid, fillp) < 0)
            {
               if (H5Tclose(fill_typeid) < 0)
                  BAIL(NC_EHDFERR);
               BAIL(NC_EHDFERR);
            }
            if (H5Tclose(fill_typeid) < 0)
               BAIL(NC_EHDFERR);
         }
      }
   }

   /* If the user wants to shuffle the data, set that up now. */
   if (var->shuffle) {
      if (H5Pset_shuffle(plistid) < 0)
         BAIL(NC_EHDFERR);
   }

   /* If the user wants to deflate the data, set that up now. */
   if (var->deflate) {
      if (H5Pset_deflate(plistid, var->deflate_level) < 0)
         BAIL(NC_EHDFERR);
   } else if(var->filterid) {
      /* Handle szip case here */
      if(var->filterid == H5Z_FILTER_SZIP) {
         int options_mask;
         int bits_per_pixel;
         if(var->nparams != 2)
            BAIL(NC_EFILTER);
         options_mask = (int)var->params[0];
         bits_per_pixel = (int)var->params[1];
         if(H5Pset_szip(plistid, options_mask, bits_per_pixel) < 0)
            BAIL(NC_EFILTER);
      } else {
         herr_t code = H5Pset_filter(plistid, var->filterid, H5Z_FLAG_MANDATORY, var->nparams, var->params);
         if(code < 0) {
            BAIL(NC_EFILTER);
         }
      }
   }

   /* If the user wants to fletcher error correcton, set that up now. */
   if (var->fletcher32)
      if (H5Pset_fletcher32(plistid) < 0)
         BAIL(NC_EHDFERR);

   /* If ndims non-zero, get info for all dimensions. We look up the
      dimids and get the len of each dimension. We need this to create
      the space for the dataset. In netCDF a dimension length of zero
      means an unlimited dimension. */
   if (var->ndims)
   {
      int unlimdim = 0;

      /* Check to see if any unlimited dimensions are used in this var. */
      for (d = 0; d < var->ndims; d++) {
         dim = var->dim[d];
         assert(dim && dim->hdr.id == var->dimids[d]);
         if (dim->unlimited)
            unlimdim++;
      }

      /* If there are no unlimited dims, and no filters, and the user
       * has not specified chunksizes, use contiguous variable for
       * better performance. */
      if (!var->shuffle && !var->deflate && !var->fletcher32 &&
          (var->chunksizes == NULL || !var->chunksizes[0]) && !unlimdim)
         var->contiguous = NC_TRUE;

      /* Gather current & maximum dimension sizes, along with chunk sizes */
      for (d = 0; d < var->ndims; d++)
      {
         dim = var->dim[d];
         assert(dim && dim->hdr.id == var->dimids[d]);
         dimsize[d] = dim->unlimited ? NC_HDF5_UNLIMITED_DIMSIZE : dim->len;
         maxdimsize[d] = dim->unlimited ? H5S_UNLIMITED : (hsize_t)dim->len;
         if (!var->contiguous) {
            if (var->chunksizes[d])
               chunksize[d] = var->chunksizes[d];
            else
            {
               size_t type_size;
               if (var->type_info->nc_type_class == NC_STRING)
                  type_size = sizeof(char *);
               else
                  type_size = var->type_info->size;

               /* Unlimited dim always gets chunksize of 1. */
               if (dim->unlimited)
                  chunksize[d] = 1;
               else
                  chunksize[d] = pow((double)DEFAULT_CHUNK_SIZE/type_size,
                                     1/(double)(var->ndims - unlimdim));

               /* If the chunksize is greater than the dim
                * length, make it the dim length. */
               if (!dim->unlimited && chunksize[d] > dim->len)
                  chunksize[d] = dim->len;

               /* Remember the computed chunksize */
               var->chunksizes[d] = chunksize[d];
            }
         }
      }

      if (var->contiguous)
      {
         if (H5Pset_layout(plistid, H5D_CONTIGUOUS) < 0)
            BAIL(NC_EHDFERR);
      }
      else
      {
         if (H5Pset_chunk(plistid, var->ndims, chunksize) < 0)
            BAIL(NC_EHDFERR);
      }

      /* Create the dataspace. */
      if ((spaceid = H5Screate_simple(var->ndims, dimsize, maxdimsize)) < 0)
         BAIL(NC_EHDFERR);
   }
   else
   {
      if ((spaceid = H5Screate(H5S_SCALAR)) < 0)
         BAIL(NC_EHDFERR);
   }

   /* Turn on creation order tracking. */
   if (H5Pset_attr_creation_order(plistid, H5P_CRT_ORDER_TRACKED|
                                  H5P_CRT_ORDER_INDEXED) < 0)
      BAIL(NC_EHDFERR);

   /* Set per-var chunk cache, for chunked datasets. */
   if (!var->contiguous && var->chunk_cache_size)
      if (H5Pset_chunk_cache(access_plistid, var->chunk_cache_nelems,
                             var->chunk_cache_size, var->chunk_cache_preemption) < 0)
         BAIL(NC_EHDFERR);

   /* At long last, create the dataset. */
   name_to_use = var->hdf5_name ? var->hdf5_name : var->hdr.name;
   LOG((4, "%s: about to H5Dcreate2 dataset %s of type 0x%x", __func__,
        name_to_use, typeid));
   if ((var->hdf_datasetid = H5Dcreate2(hdf5_grp->hdf_grpid, name_to_use, typeid,
                                        spaceid, H5P_DEFAULT, plistid, access_plistid)) < 0)
      BAIL(NC_EHDFERR);
   var->created = NC_TRUE;
   var->is_new_var = NC_FALSE;

   /* If this is a dimscale, mark it as such in the HDF5 file. Also
    * find the dimension info and store the dataset id of the dimscale
    * dataset. */
   if (var->dimscale)
   {
      if (H5DSset_scale(var->hdf_datasetid, var->hdr.name) < 0)
         BAIL(NC_EHDFERR);

      /* If this is a multidimensional coordinate variable, write a
       * coordinates attribute. */
      if (var->ndims > 1)
         if ((retval = write_coord_dimids(var)))
            BAIL(retval);

      /* If desired, write the netCDF dimid. */
      if (write_dimid)
         if ((retval = write_netcdf4_dimid(var->hdf_datasetid, var->dimids[0])))
            BAIL(retval);
   }

   /* Write attributes for this var. */
   if ((retval = write_attlist(var->att, var->hdr.id, grp)))
      BAIL(retval);
   var->attr_dirty = NC_FALSE;

exit:
   if (typeid > 0 && H5Tclose(typeid) < 0)
      BAIL2(NC_EHDFERR);
   if (plistid > 0 && H5Pclose(plistid) < 0)
      BAIL2(NC_EHDFERR);
   if (access_plistid > 0 && H5Pclose(access_plistid) < 0)
      BAIL2(NC_EHDFERR);
   if (spaceid > 0 && H5Sclose(spaceid) < 0)
      BAIL2(NC_EHDFERR);
   if (fillp)
   {
      if (var->type_info->nc_type_class == NC_VLEN)
         nc_free_vlen((nc_vlen_t *)fillp);
      else if (var->type_info->nc_type_class == NC_STRING && *(char **)fillp)
         free(*(char **)fillp);
      free(fillp);
   }

   return retval;
}

/**
 * @internal Adjust the chunk cache of a var for better
 * performance.
 *
 * @param grp Pointer to group info struct.
 * @param var Pointer to var info struct.
 *
 * @return NC_NOERR No error.
 * @author Ed Hartnett
 */
int
nc4_adjust_var_cache(NC_GRP_INFO_T *grp, NC_VAR_INFO_T *var)
{
   size_t chunk_size_bytes = 1;
   int d;
   int retval;

   /* Nothing to be done. */
   if (var->contiguous)
      return NC_NOERR;
#ifdef USE_PARALLEL4
   return NC_NOERR;
#endif

   /* How many bytes in the chunk? */
   for (d = 0; d < var->ndims; d++)
      chunk_size_bytes *= var->chunksizes[d];
   if (var->type_info->size)
      chunk_size_bytes *= var->type_info->size;
   else
      chunk_size_bytes *= sizeof(char *);

   /* If the chunk cache is too small, and the user has not changed
    * the default value of the chunk cache size, then increase the
    * size of the cache. */
   if (var->chunk_cache_size == CHUNK_CACHE_SIZE)
      if (chunk_size_bytes > var->chunk_cache_size)
      {
         var->chunk_cache_size = chunk_size_bytes * DEFAULT_CHUNKS_IN_CACHE;
         if (var->chunk_cache_size > MAX_DEFAULT_CACHE_SIZE)
            var->chunk_cache_size = MAX_DEFAULT_CACHE_SIZE;
         if ((retval = nc4_reopen_dataset(grp, var)))
            return retval;
      }

   return NC_NOERR;
}

/**
 * @internal Create a HDF5 defined type from a NC_TYPE_INFO_T struct,
 * and commit it to the file.
 *
 * @param grp Pointer to group info struct.
 * @param type Pointer to type info struct.
 *
 * @return NC_NOERR No error.
 * @return NC_EHDFERR HDF5 error.
 * @return NC_ECHAR Conversions of NC_CHAR forbidden.
 * @return NC_EVARMETA HDF5 returning error creating datatype.
 * @return NC_EHDFERR HDF5 returning error.
 * @return NC_EBADTYE Type not found.
 * @author Ed Hartnett
 */
static int
commit_type(NC_GRP_INFO_T *grp, NC_TYPE_INFO_T *type)
{
   NC_HDF5_GRP_INFO_T *hdf5_grp;
   int retval;

   assert(grp && grp->format_grp_info && type);

   /* Get HDF5-specific group info. */
   hdf5_grp = (NC_HDF5_GRP_INFO_T *)grp->format_grp_info;

   /* Did we already record this type? */
   if (type->committed)
      return NC_NOERR;

   /* Is this a compound type? */
   if (type->nc_type_class == NC_COMPOUND)
   {
      NC_FIELD_INFO_T *field;
      hid_t hdf_base_typeid, hdf_typeid;
      int i;

      if ((type->hdf_typeid = H5Tcreate(H5T_COMPOUND, type->size)) < 0)
         return NC_EHDFERR;
      LOG((4, "creating compound type %s hdf_typeid 0x%x", type->hdr.name,
           type->hdf_typeid));

      for(i=0;i<nclistlength(type->u.c.field);i++)
      {
         field = (NC_FIELD_INFO_T *)nclistget(type->u.c.field, i);
         assert(field);
         if ((retval = nc4_get_hdf_typeid(grp->nc4_info, field->nc_typeid,
                                          &hdf_base_typeid, type->endianness)))
            return retval;

         /* If this is an array, create a special array type. */
         if (field->ndims)
         {
            int d;
            hsize_t dims[NC_MAX_VAR_DIMS];

            for (d = 0; d < field->ndims; d++)
               dims[d] = field->dim_size[d];
            if ((hdf_typeid = H5Tarray_create(hdf_base_typeid, field->ndims,
                                              dims, NULL)) < 0)
            {
               if (H5Tclose(hdf_base_typeid) < 0)
                  return NC_EHDFERR;
               return NC_EHDFERR;
            }
            if (H5Tclose(hdf_base_typeid) < 0)
               return NC_EHDFERR;
         }
         else
            hdf_typeid = hdf_base_typeid;
         LOG((4, "inserting field %s offset %d hdf_typeid 0x%x", field->hdr.name,
              field->offset, hdf_typeid));
         if (H5Tinsert(type->hdf_typeid, field->hdr.name, field->offset,
                       hdf_typeid) < 0)
            return NC_EHDFERR;
         if (H5Tclose(hdf_typeid) < 0)
            return NC_EHDFERR;
      }
   }
   else if (type->nc_type_class == NC_VLEN)
   {
      /* Find the HDF typeid of the base type of this vlen. */
      if ((retval = nc4_get_hdf_typeid(grp->nc4_info, type->u.v.base_nc_typeid,
                                       &type->u.v.base_hdf_typeid, type->endianness)))
         return retval;

      /* Create a vlen type. */
      if ((type->hdf_typeid = H5Tvlen_create(type->u.v.base_hdf_typeid)) < 0)
         return NC_EHDFERR;
   }
   else if (type->nc_type_class == NC_OPAQUE)
   {
      /* Create the opaque type. */
      if ((type->hdf_typeid = H5Tcreate(H5T_OPAQUE, type->size)) < 0)
         return NC_EHDFERR;
   }
   else if (type->nc_type_class == NC_ENUM)
   {
      NC_ENUM_MEMBER_INFO_T *enum_m;
      int i;

      if (nclistlength(type->u.e.enum_member) == 0)
         return NC_EINVAL;

      /* Find the HDF typeid of the base type of this enum. */
      if ((retval = nc4_get_hdf_typeid(grp->nc4_info, type->u.e.base_nc_typeid,
                                       &type->u.e.base_hdf_typeid, type->endianness)))
         return retval;

      /* Create an enum type. */
      if ((type->hdf_typeid =  H5Tenum_create(type->u.e.base_hdf_typeid)) < 0)
         return NC_EHDFERR;

      /* Add all the members to the HDF5 type. */
      for(i=0;i<nclistlength(type->u.e.enum_member);i++) {
         enum_m = (NC_ENUM_MEMBER_INFO_T*)nclistget(type->u.e.enum_member,i);
         if (H5Tenum_insert(type->hdf_typeid, enum_m->name, enum_m->value) < 0)
            return NC_EHDFERR;
      }
   }
   else
   {
      LOG((0, "Unknown class: %d", type->nc_type_class));
      return NC_EBADTYPE;
   }

   /* Commit the type. */
   if (H5Tcommit(hdf5_grp->hdf_grpid, type->hdr.name, type->hdf_typeid) < 0)
      return NC_EHDFERR;
   type->committed = NC_TRUE;
   LOG((4, "just committed type %s, HDF typeid: 0x%x", type->hdr.name,
        type->hdf_typeid));

   /* Later we will always use the native typeid. In this case, it is
    * a copy of the same type pointed to by hdf_typeid, but it's
    * easier to maintain a copy. */
   if ((type->native_hdf_typeid = H5Tget_native_type(type->hdf_typeid,
                                                     H5T_DIR_DEFAULT)) < 0)
      return NC_EHDFERR;

   return NC_NOERR;
}

/**
 * @internal Write an attribute, with value 1, to indicate that strict
 * NC3 rules apply to this file.
 *
 * @param hdf_grpid HDF5 group ID.
 *
 * @returns NC_NOERR No error.
 * @returns NC_EHDFERR HDF5 returned an error.
 * @author Ed Hartnett
 */
static int
write_nc3_strict_att(hid_t hdf_grpid)
{
   hid_t attid = 0, spaceid = 0;
   int one = 1;
   int retval = NC_NOERR;
   htri_t attr_exists;

   /* If the attribute already exists, call that a success and return
    * NC_NOERR. */
   if ((attr_exists = H5Aexists(hdf_grpid, NC3_STRICT_ATT_NAME)) < 0)
      return NC_EHDFERR;
   if (attr_exists)
      return NC_NOERR;

   /* Create the attribute to mark this as a file that needs to obey
    * strict netcdf-3 rules. */
   if ((spaceid = H5Screate(H5S_SCALAR)) < 0)
      BAIL(NC_EFILEMETA);
   if ((attid = H5Acreate(hdf_grpid, NC3_STRICT_ATT_NAME,
                          H5T_NATIVE_INT, spaceid, H5P_DEFAULT)) < 0)
      BAIL(NC_EFILEMETA);
   if (H5Awrite(attid, H5T_NATIVE_INT, &one) < 0)
      BAIL(NC_EFILEMETA);

exit:
   if (spaceid > 0 && (H5Sclose(spaceid) < 0))
      BAIL2(NC_EFILEMETA);
   if (attid > 0 && (H5Aclose(attid) < 0))
      BAIL2(NC_EFILEMETA);
   return retval;
}

/**
 * @internal Create a HDF5 group that is not the root group. HDF5
 * properties are set in the group to ensure that objects and
 * attributes are kept in creation order, instead of alphebetical
 * order (the HDF5 default).
 *
 * @param grp Pointer to group info struct.
 *
 * @return NC_NOERR No error.
 * @return NC_EHDFERR HDF5 error.
 * @author Ed Hartnett
 */
static int
create_group(NC_GRP_INFO_T *grp)
{
   NC_HDF5_GRP_INFO_T *hdf5_grp, *parent_hdf5_grp;
   hid_t gcpl_id = -1;
   int retval = NC_NOERR;;

   assert(grp && grp->format_grp_info && grp->parent &&
          grp->parent->format_grp_info);

   /* Get HDF5 specific group info for group and parent. */
   hdf5_grp = (NC_HDF5_GRP_INFO_T *)grp->format_grp_info;
   parent_hdf5_grp = (NC_HDF5_GRP_INFO_T *)grp->parent->format_grp_info;
   assert(parent_hdf5_grp->hdf_grpid);

   /* Create group, with link_creation_order set in the group
    * creation property list. */
   if ((gcpl_id = H5Pcreate(H5P_GROUP_CREATE)) < 0)
      BAIL(NC_EHDFERR);

   /* Set track_times to be FALSE. */
   if (H5Pset_obj_track_times(gcpl_id, 0) < 0)
      BAIL(NC_EHDFERR);

   /* Tell HDF5 to keep track of objects in creation order. */
   if (H5Pset_link_creation_order(gcpl_id, H5P_CRT_ORDER_TRACKED|H5P_CRT_ORDER_INDEXED) < 0)
      BAIL(NC_EHDFERR);

   /* Tell HDF5 to keep track of attributes in creation order. */
   if (H5Pset_attr_creation_order(gcpl_id, H5P_CRT_ORDER_TRACKED|H5P_CRT_ORDER_INDEXED) < 0)
      BAIL(NC_EHDFERR);

   /* Create the group. */
   if ((hdf5_grp->hdf_grpid = H5Gcreate2(parent_hdf5_grp->hdf_grpid, grp->hdr.name,
                                         H5P_DEFAULT, gcpl_id, H5P_DEFAULT)) < 0)
      BAIL(NC_EHDFERR);

exit:
   if (gcpl_id > -1 && H5Pclose(gcpl_id) < 0)
      BAIL2(NC_EHDFERR);
   if (retval)
      if (hdf5_grp->hdf_grpid > 0 && H5Gclose(hdf5_grp->hdf_grpid) < 0)
         BAIL2(NC_EHDFERR);
   return retval;
}

/**
 * @internal After all the datasets of the file have been read, it's
 * time to sort the wheat from the chaff. Which of the datasets are
 * netCDF dimensions, and which are coordinate variables, and which
 * are non-coordinate variables.
 *
 * @param grp Pointer to group info struct.
 *
 * @return ::NC_NOERR No error.
 * @author Ed Hartnett
 */
static int
attach_dimscales(NC_GRP_INFO_T *grp)
{
   NC_VAR_INFO_T *var;
   NC_DIM_INFO_T *dim1;
   int d, i;
   int retval = NC_NOERR;

   /* Attach dimension scales. */
   for(i=0;i<ncindexsize(grp->vars);i++)
   {
      var = (NC_VAR_INFO_T*)ncindexith(grp->vars,i);
      if (!var) continue;
      /* Scales themselves do not attach. But I really wish they
       * would. */
      if (var->dimscale)
      {
         /* If this is a multidimensional coordinate variable, it will
          * have a special coords attribute (read earlier) with a list
          * of the dimensions for this variable. */
      }
      else /* not a dimscale... */
      {
         /* Find the scale for each dimension, if any, and attach it. */
         for (d = 0; d < var->ndims; d++)
         {
            /* Is there a dimscale for this dimension? */
            if (var->dimscale_attached)
            {
               if (!var->dimscale_attached[d])
               {
                  NC_HDF5_DIM_INFO_T *hdf5_dim1;
                  hid_t dim_datasetid;  /* Dataset ID for dimension */
                  dim1 = var->dim[d];
                  assert(dim1 && dim1->hdr.id == var->dimids[d] && dim1->format_dim_info);
                  hdf5_dim1 = (NC_HDF5_DIM_INFO_T *)dim1->format_dim_info;

                  LOG((2, "%s: attaching scale for dimid %d to var %s",
                       __func__, var->dimids[d], var->hdr.name));

                  /* Find dataset ID for dimension */
                  if (dim1->coord_var)
                     dim_datasetid = dim1->coord_var->hdf_datasetid;
                  else
                     dim_datasetid = hdf5_dim1->hdf_dimscaleid;
                  if(!(dim_datasetid > 0))
                     assert(dim_datasetid > 0);
                  if (H5DSattach_scale(var->hdf_datasetid, dim_datasetid, d) < 0)
                     BAIL(NC_EHDFERR);
                  var->dimscale_attached[d] = NC_TRUE;
               }

               /* If we didn't find a dimscale to attach, that's a problem! */
               if (!var->dimscale_attached[d])
               {
                  LOG((0, "no dimscale found!"));
                  return NC_EDIMSCALE;
               }
            }
         }
      }
   }

exit:
   return retval;
}

/**
 * @internal Does a variable exist?
 *
 * @param grpid HDF5 group ID.
 * @param name Name of variable.
 * @param exists Pointer that gets 1 of the variable exists, 0 otherwise.
 *
 * @return ::NC_NOERR No error.
 * @author Ed Hartnett
 */
static int
var_exists(hid_t grpid, char *name, nc_bool_t *exists)
{
   htri_t link_exists;

   /* Reset the boolean */
   *exists = NC_FALSE;

   /* Check if the object name exists in the group */
   if ((link_exists = H5Lexists(grpid, name, H5P_DEFAULT)) < 0)
      return NC_EHDFERR;
   if (link_exists)
   {
      H5G_stat_t statbuf;

      /* Get info about the object */
      if (H5Gget_objinfo(grpid, name, 1, &statbuf) < 0)
         return NC_EHDFERR;

      if (H5G_DATASET == statbuf.type)
         *exists = NC_TRUE;
   }

   return NC_NOERR;
}

/**
 * @internal Convert a coordinate variable HDF5 dataset into one that
 * is not a coordinate variable. This happens during renaming of vars
 * and dims. This function removes the HDF5 NAME and CLASS attributes
 * associated with dimension scales, and also the NC_DIMID_ATT_NAME
 * attribute which may be present, and, if it does, holds the dimid of
 * the coordinate variable.
 *
 * @param hdf_datasetid The HDF5 dataset ID of the coordinate variable dataset.
 *
 * @return ::NC_NOERR No error.
 * @return ::NC_EHDFERR HDF5 error.
 * @author Ed Hartnett
 */
int
remove_coord_atts(hid_t hdf_datasetid)
{
   htri_t attr_exists;

   /* If the variable dataset has an optional NC_DIMID_ATT_NAME
    * attribute, delete it. */
   if ((attr_exists = H5Aexists(hdf_datasetid, NC_DIMID_ATT_NAME)) < 0)
      return NC_EHDFERR;
   if (attr_exists)
   {
      if (H5Adelete(hdf_datasetid, NC_DIMID_ATT_NAME) < 0)
         return NC_EHDFERR;
   }

   /* (We could do a better job here and verify that the attributes are
    * really dimension scale 'CLASS' & 'NAME' attributes, but that would be
    * poking about in the HDF5 DimScale internal data) */
   if ((attr_exists = H5Aexists(hdf_datasetid, HDF5_DIMSCALE_CLASS_ATT_NAME)) < 0)
      return NC_EHDFERR;
   if (attr_exists)
   {
      if (H5Adelete(hdf_datasetid, HDF5_DIMSCALE_CLASS_ATT_NAME) < 0)
         return NC_EHDFERR;
   }
   if ((attr_exists = H5Aexists(hdf_datasetid, HDF5_DIMSCALE_NAME_ATT_NAME)) < 0)
      return NC_EHDFERR;
   if (attr_exists)
   {
      if (H5Adelete(hdf_datasetid, HDF5_DIMSCALE_NAME_ATT_NAME) < 0)
         return NC_EHDFERR;
   }
   return NC_NOERR;
}

/**
 * @internal This function writes a variable. The principle difficulty
 * comes from the possibility that this is a coordinate variable, and
 * was already written to the file as a dimension-only dimscale. If
 * this occurs, then it must be deleted and recreated.
 *
 * @param var Pointer to variable info struct.
 * @param grp Pointer to group info struct.
 * @param write_dimid
 *
 * @returns NC_NOERR No error.
 * @returns NC_EHDFERR HDF5 returned an error.
 * @author Ed Hartnett, Quincey Koziol
 */
static int
write_var(NC_VAR_INFO_T *var, NC_GRP_INFO_T *grp, nc_bool_t write_dimid)
{
   NC_HDF5_GRP_INFO_T *hdf5_grp;
   nc_bool_t replace_existing_var = NC_FALSE;
   int retval;

   assert(var && grp && grp->format_grp_info);

   LOG((4, "%s: writing var %s", __func__, var->hdr.name));

   /* Get HDF5-specific group info. */
   hdf5_grp = (NC_HDF5_GRP_INFO_T *)grp->format_grp_info;

   /* If the variable has already been created & the fill value changed,
    * indicate that the existing variable should be replaced. */
   if (var->created && var->fill_val_changed)
   {
      replace_existing_var = NC_TRUE;
      var->fill_val_changed = NC_FALSE;
      /* If the variable is going to be replaced,
         we need to flag any other attributes associated
         with the variable as 'dirty', or else
         *only* the fill value attribute will be copied over
         and the rest will be lost.  See:

         * https://github.com/Unidata/netcdf-c/issues/239 */

      flag_atts_dirty(var->att);
   }

   /* Is this a coordinate var that has already been created in
    * the HDF5 file as a dimscale dataset? Check for dims with the
    * same name in this group. If there is one, check to see if
    * this object exists in the HDF group. */
   if (var->became_coord_var)
   {
      NC_DIM_INFO_T *d1;
      int i;

      for (i = 0; i < ncindexsize(grp->dim); i++)
      {
         d1 = (NC_DIM_INFO_T*)ncindexith(grp->dim,i);
         assert(d1);
         if (!strcmp(d1->hdr.name, var->hdr.name))
         {
            nc_bool_t exists;

            if ((retval = var_exists(hdf5_grp->hdf_grpid, var->hdr.name, &exists)))
               return retval;
            if (exists)
            {
               /* Indicate that the variable already exists, and should be replaced */
               replace_existing_var = NC_TRUE;
               flag_atts_dirty(var->att);
               break;
            }
         }
      }
   }

   /* Check dims if the variable will be replaced, so that the
    * dimensions will be de-attached and re-attached correctly. (Note:
    * There's a temptation to merge this loop over the dimensions with
    * the prior loop over dimensions, but that blurs the line over the
    * purpose of them, so they are currently separate.  If performance
    * becomes an issue here, it would be possible to merge them. -QAK)
    */
   if (replace_existing_var)
   {
      int i;

      for (i = 0; i < ncindexsize(grp->dim); i++)
      {
         NC_DIM_INFO_T *d1;
         NC_HDF5_DIM_INFO_T *hdf5_d1;

         /* Get info about the dim, including HDF5-specific info. */
         d1 = (NC_DIM_INFO_T *)ncindexith(grp->dim, i);
         assert(d1 && d1->format_dim_info && d1->hdr.name);
         hdf5_d1 = (NC_HDF5_DIM_INFO_T *)d1->format_dim_info;

         if (!strcmp(d1->hdr.name, var->hdr.name))
         {
            nc_bool_t exists;

            if ((retval = var_exists(hdf5_grp->hdf_grpid, var->hdr.name,
                                     &exists)))
               return retval;
            if (exists)
            {
               hid_t dim_datasetid;  /* Dataset ID for dimension */

               /* Find dataset ID for dimension */
               if (d1->coord_var)
                  dim_datasetid = d1->coord_var->hdf_datasetid;
               else
                  dim_datasetid = hdf5_d1->hdf_dimscaleid;
               assert(dim_datasetid > 0);

               /* If we're replacing an existing dimscale dataset, go to
                * every var in the file and detach this dimension scale,
                * because we have to delete it. */
               if ((retval = rec_detach_scales(grp->nc4_info->root_grp,
                                               var->dimids[0], dim_datasetid)))
                  return retval;
               break;
            }
         }
      }
   }

   /* If this is not a dimension scale, do this stuff. */
   if (var->was_coord_var && var->dimscale_attached)
   {
      /* If the variable already exists in the file, Remove any dimension scale
       * attributes from it, if they exist. */
      /* (The HDF5 Dimension Scale API should really have an API routine
       * for making a dataset not a scale. -QAK) */
      if (var->created)
      {
         if ((retval = remove_coord_atts(var->hdf_datasetid)))
            BAIL(retval);
      }

      if (var->dimscale_attached)
      {
         int d;

         /* If this is a regular var, detach all its dim scales. */
         for (d = 0; d < var->ndims; d++)
         {
            if (var->dimscale_attached[d])
            {
               hid_t dim_datasetid;  /* Dataset ID for dimension */
               NC_DIM_INFO_T *dim1 = var->dim[d];
               NC_HDF5_DIM_INFO_T *hdf5_dim1;
               assert(dim1 && dim1->hdr.id == var->dimids[d] && dim1->format_dim_info);
               hdf5_dim1 = (NC_HDF5_DIM_INFO_T *)dim1->format_dim_info;

               /* Find dataset ID for dimension */
               if (dim1->coord_var)
                  dim_datasetid = dim1->coord_var->hdf_datasetid;
               else
                  dim_datasetid = hdf5_dim1->hdf_dimscaleid;
               assert(dim_datasetid > 0);

               if (H5DSdetach_scale(var->hdf_datasetid, dim_datasetid, d) < 0)
                  BAIL(NC_EHDFERR);
               var->dimscale_attached[d] = NC_FALSE;
            }
         }
      }
   }

   /* Delete the HDF5 dataset that is to be replaced. */
   if (replace_existing_var)
   {
      /* Free the HDF5 dataset id. */
      if (var->hdf_datasetid && H5Dclose(var->hdf_datasetid) < 0)
         BAIL(NC_EHDFERR);
      var->hdf_datasetid = 0;

      /* Now delete the variable. */
      if (H5Gunlink(hdf5_grp->hdf_grpid, var->hdr.name) < 0)
         return NC_EDIMMETA;
   }

   /* Create the dataset. */
   if (var->is_new_var || replace_existing_var)
   {
      if ((retval = var_create_dataset(grp, var, write_dimid)))
         return retval;
   }
   else
   {
      if (write_dimid && var->ndims)
         if ((retval = write_netcdf4_dimid(var->hdf_datasetid, var->dimids[0])))
            BAIL(retval);
   }

   if (replace_existing_var)
   {
      /* If this is a dimension scale, reattach the scale everywhere it
       * is used. (Recall that netCDF dimscales are always 1-D). */
      if(var->dimscale)
      {
         if ((retval = rec_reattach_scales(grp->nc4_info->root_grp,
                                           var->dimids[0], var->hdf_datasetid)))
            return retval;
      }
      /* If it's not a dimension scale, clear the dimscale attached flags,
       * so the dimensions are re-attached. */
      else
      {
         if (var->dimscale_attached)
            memset(var->dimscale_attached, 0, sizeof(nc_bool_t) * var->ndims);
      }
   }

   /* Clear coord. var state transition flags */
   var->was_coord_var = NC_FALSE;
   var->became_coord_var = NC_FALSE;

   /* Now check the attributes for this var. */
   if (var->attr_dirty)
   {
      /* Write attributes for this var. */
      if ((retval = write_attlist(var->att, var->hdr.id, grp)))
         BAIL(retval);
      var->attr_dirty = NC_FALSE;
   }

   return NC_NOERR;
exit:
   return retval;
}

/**
 * @internal Write a dimension.
 *
 * @param dim Pointer to dim info struct.
 * @param grp Pointer to group info struct.
 * @param write_dimid
 *
 * @returns ::NC_NOERR No error.
 * @returns ::NC_EPERM Read-only file.
 * @returns ::NC_EHDFERR HDF5 returned error.
 * @author Ed Hartnett
 */
static int
write_dim(NC_DIM_INFO_T *dim, NC_GRP_INFO_T *grp, nc_bool_t write_dimid)
{
   NC_HDF5_GRP_INFO_T *hdf5_grp;
   NC_HDF5_DIM_INFO_T *hdf5_dim;
   int retval;

   assert(dim && dim->format_dim_info && grp && grp->format_grp_info);

   /* Get HDF5-specific dim and group info. */
   hdf5_dim = (NC_HDF5_DIM_INFO_T *)dim->format_dim_info;
   hdf5_grp = (NC_HDF5_GRP_INFO_T *)grp->format_grp_info;

   /* If there's no dimscale dataset for this dim, create one,
    * and mark that it should be hidden from netCDF as a
    * variable. (That is, it should appear as a dimension
    * without an associated variable.) */
   if (!hdf5_dim->hdf_dimscaleid)
   {
      hid_t spaceid, create_propid;
      hsize_t dims[1], max_dims[1], chunk_dims[1] = {1};
      char dimscale_wo_var[NC_MAX_NAME];

      LOG((4, "%s: creating dim %s", __func__, dim->hdr.name));

      /* Sanity check */
      assert(NULL == dim->coord_var);

      /* Create a property list. If this dimension scale is
       * unlimited (i.e. it's an unlimited dimension), then set
       * up chunking, with a chunksize of 1. */
      if ((create_propid = H5Pcreate(H5P_DATASET_CREATE)) < 0)
         BAIL(NC_EHDFERR);

      /* RJ: this suppose to be FALSE that is defined in H5 private.h as 0 */
      if (H5Pset_obj_track_times(create_propid,0)<0)
         BAIL(NC_EHDFERR);

      dims[0] = dim->len;
      max_dims[0] = dim->len;
      if (dim->unlimited)
      {
         max_dims[0] = H5S_UNLIMITED;
         if (H5Pset_chunk(create_propid, 1, chunk_dims) < 0)
            BAIL(NC_EHDFERR);
      }

      /* Set up space. */
      if ((spaceid = H5Screate_simple(1, dims, max_dims)) < 0)
         BAIL(NC_EHDFERR);

      if (H5Pset_attr_creation_order(create_propid, H5P_CRT_ORDER_TRACKED|
                                     H5P_CRT_ORDER_INDEXED) < 0)
         BAIL(NC_EHDFERR);

      /* Create the dataset that will be the dimension scale. */
      LOG((4, "%s: about to H5Dcreate1 a dimscale dataset %s", __func__, dim->hdr.name));
      if ((hdf5_dim->hdf_dimscaleid = H5Dcreate1(hdf5_grp->hdf_grpid, dim->hdr.name, H5T_IEEE_F32BE,
                                                 spaceid, create_propid)) < 0)
         BAIL(NC_EHDFERR);

      /* Close the spaceid and create_propid. */
      if (H5Sclose(spaceid) < 0)
         BAIL(NC_EHDFERR);
      if (H5Pclose(create_propid) < 0)
         BAIL(NC_EHDFERR);

      /* Indicate that this is a scale. Also indicate that not
       * be shown to the user as a variable. It is hidden. It is
       * a DIM WITHOUT A VARIABLE! */
      sprintf(dimscale_wo_var, "%s%10d", DIM_WITHOUT_VARIABLE, (int)dim->len);
      if (H5DSset_scale(hdf5_dim->hdf_dimscaleid, dimscale_wo_var) < 0)
         BAIL(NC_EHDFERR);
   }

   /* Did we extend an unlimited dimension? */
   if (dim->extended)
   {
      NC_VAR_INFO_T *v1 = NULL;

      assert(dim->unlimited);
      /* If this is a dimension without a variable, then update
       * the secret length information at the end of the NAME
       * attribute. */
      v1 = (NC_VAR_INFO_T*)ncindexlookup(grp->vars,dim->hdr.name);
      if (v1)
      {
         hsize_t *new_size = NULL;
         int d1;

         /* Extend the dimension scale dataset to reflect the new
          * length of the dimension. */
         if (!(new_size = malloc(v1->ndims * sizeof(hsize_t))))
            BAIL(NC_ENOMEM);
         for (d1 = 0; d1 < v1->ndims; d1++)
         {
            assert(v1->dim[d1] && v1->dim[d1]->hdr.id == v1->dimids[d1]);
            new_size[d1] = v1->dim[d1]->len;
         }
         if (H5Dset_extent(v1->hdf_datasetid, new_size) < 0) {
            free(new_size);
            BAIL(NC_EHDFERR);
         }
         free(new_size);
      }
   }

   /* If desired, write the secret dimid. This will be used instead of
    * the dimid that the dimension would otherwise receive based on
    * creation order. This can be necessary when dims and their
    * coordinate variables were created in different order. */
   if (write_dimid && hdf5_dim->hdf_dimscaleid)
      if ((retval = write_netcdf4_dimid(hdf5_dim->hdf_dimscaleid, dim->hdr.id)))
         BAIL(retval);

   return NC_NOERR;
exit:

   return retval;
}

/**
 * @internal Recursively determine if there is a mismatch between
 * order of coordinate creation and associated dimensions in this
 * group or any subgroups, to find out if we have to handle that
 * situation.  Also check if there are any multidimensional coordinate
 * variables defined, which require the same treatment to fix a
 * potential bug when such variables occur in subgroups.
 *
 * @param grp Pointer to group info struct.
 * @param bad_coord_orderp Pointer that gets 1 if there is a bad
 * coordinate order.
 *
 * @returns NC_NOERR No error.
 * @returns NC_EHDFERR HDF5 returned an error.
 * @author Ed Hartnett
 */
int
nc4_detect_preserve_dimids(NC_GRP_INFO_T *grp, nc_bool_t *bad_coord_orderp)
{
   NC_VAR_INFO_T *var;
   NC_GRP_INFO_T *child_grp;
   int last_dimid = -1;
   int retval;
   int i;

   /* Iterate over variables in this group */
   for (i=0; i < ncindexsize(grp->vars); i++)
   {
      var = (NC_VAR_INFO_T*)ncindexith(grp->vars,i);
      if (var == NULL) continue;
      /* Only matters for dimension scale variables, with non-scalar dimensionality */
      if (var->dimscale && var->ndims)
      {
         /* If the user writes coord vars in a different order then he
          * defined their dimensions, then, when the file is reopened, the
          * order of the dimids will change to match the order of the coord
          * vars. Detect if this is about to happen. */
         if (var->dimids[0] < last_dimid)
         {
            LOG((5, "%s: %s is out of order coord var", __func__, var->hdr.name));
            *bad_coord_orderp = NC_TRUE;
            return NC_NOERR;
         }
         last_dimid = var->dimids[0];

         /* If there are multidimensional coordinate variables defined, then
          * it's also necessary to preserve dimension IDs when the file is
          * reopened ... */
         if (var->ndims > 1)
         {
            LOG((5, "%s: %s is multidimensional coord var", __func__, var->hdr.name));
            *bad_coord_orderp = NC_TRUE;
            return NC_NOERR;
         }

         /* Did the user define a dimension, end define mode, reenter define
          * mode, and then define a coordinate variable for that dimension?
          * If so, dimensions will be out of order. */
         if (var->is_new_var || var->became_coord_var)
         {
            LOG((5, "%s: coord var defined after enddef/redef", __func__));
            *bad_coord_orderp = NC_TRUE;
            return NC_NOERR;
         }
      }
   }

   /* If there are any child groups, check them also for this condition. */
   for (i = 0; i < ncindexsize(grp->children); i++)
   {
      if (!(child_grp = (NC_GRP_INFO_T *)ncindexith(grp->children, i)))
         continue;
      if ((retval = nc4_detect_preserve_dimids(child_grp, bad_coord_orderp)))
         return retval;
   }
   return NC_NOERR;
}

/**
 * @internal Recursively write all the metadata in a group. Groups and
 * types have all already been written. Propagate bad cooordinate
 * order to subgroups, if detected.
 *
 * @param grp Pointer to group info struct.
 * @param bad_coord_order 1 if there is a bad coordinate order.
 *
 * @returns NC_NOERR No error.
 * @returns NC_EHDFERR HDF5 returned an error.
 * @author Ed Hartnett
 */
int
nc4_rec_write_metadata(NC_GRP_INFO_T *grp, nc_bool_t bad_coord_order)
{
   NC_DIM_INFO_T *dim = NULL;
   NC_VAR_INFO_T *var = NULL;
   NC_GRP_INFO_T *child_grp = NULL;
   int coord_varid = -1;
   int var_index = 0;
   int dim_index = 0;
   int retval;
   int i;

   assert(grp && grp->hdr.name &&
          ((NC_HDF5_GRP_INFO_T *)(grp->format_grp_info))->hdf_grpid);
   LOG((3, "%s: grp->hdr.name %s, bad_coord_order %d", __func__, grp->hdr.name,
        bad_coord_order));

   /* Write global attributes for this group. */
   if ((retval = write_attlist(grp->att, NC_GLOBAL, grp)))
      return retval;

   /* Set the pointers to the beginning of the list of dims & vars in this
    * group. */
   dim_index = 0;
   var_index = 0;
   /* prime the loop */
   dim = (NC_DIM_INFO_T*)ncindexith(grp->dim,dim_index);
   var = (NC_VAR_INFO_T*)ncindexith(grp->vars,var_index);
   /* Because of HDF5 ordering the dims and vars have to be stored in
    * this way to ensure that the dims and coordinate vars come out in
    * the correct order. */
   while (dim || var)
   {
      nc_bool_t found_coord, wrote_coord;

      /* Write non-coord dims in order, stopping at the first one that
       * has an associated coord var. */
      for (found_coord = NC_FALSE; dim && !found_coord; )
      {
         if (!dim->coord_var)
         {
            if ((retval = write_dim(dim, grp, bad_coord_order)))
               return retval;
         }
         else
         {
            coord_varid = dim->coord_var->hdr.id;
            found_coord = NC_TRUE;
         }
         dim = (NC_DIM_INFO_T*)ncindexith(grp->dim,++dim_index);
      }

      /* Write each var. When we get to the coord var we are waiting
       * for (if any), then we break after writing it. */
      for (wrote_coord = NC_FALSE; var && !wrote_coord; )
      {
         if ((retval = write_var(var, grp, bad_coord_order)))
            return retval;
         if (found_coord && var->hdr.id == coord_varid)
            wrote_coord = NC_TRUE;
         var = (NC_VAR_INFO_T*)ncindexith(grp->vars,++var_index);
      }
   } /* end while */

   if ((retval = attach_dimscales(grp)))
      return retval;

   /* If there are any child groups, write their metadata. */
   for(i=0;i<ncindexsize(grp->children);i++) {
      if((child_grp = (NC_GRP_INFO_T*)ncindexith(grp->children,i)) == NULL) continue;
      if ((retval = nc4_rec_write_metadata(child_grp, bad_coord_order)))
         return retval;
   }
   return NC_NOERR;
}

/**
 * @internal Recursively write all groups and types.
 *
 * @param grp Pointer to group info struct.
 *
 * @returns NC_NOERR No error.
 * @returns NC_EHDFERR HDF5 returned an error.
 * @author Ed Hartnett
 */
int
nc4_rec_write_groups_types(NC_GRP_INFO_T *grp)
{
   NC_GRP_INFO_T *child_grp;
   NC_HDF5_GRP_INFO_T *hdf5_grp;
   NC_TYPE_INFO_T *type;
   int retval;
   int i;

   assert(grp && grp->hdr.name && grp->format_grp_info);
   LOG((3, "%s: grp->hdr.name %s", __func__, grp->hdr.name));

   /* Get HDF5-specific group info. */
   hdf5_grp = (NC_HDF5_GRP_INFO_T *)grp->format_grp_info;

   /* Create the group in the HDF5 file if it doesn't exist. */
   if (!hdf5_grp->hdf_grpid)
      if ((retval = create_group(grp)))
         return retval;

   /* If this is the root group of a file with strict NC3 rules, write
    * an attribute. But don't leave the attribute open. */
   if (!grp->parent && (grp->nc4_info->cmode & NC_CLASSIC_MODEL))
      if ((retval = write_nc3_strict_att(hdf5_grp->hdf_grpid)))
         return retval;

   /* If there are any user-defined types, write them now. */
   for(i=0;i<ncindexsize(grp->type);i++) {
      type = (NC_TYPE_INFO_T *)ncindexith(grp->type, i);
      assert(type);
      if ((retval = commit_type(grp, type)))
         return retval;
   }

   /* If there are any child groups, write their groups and types. */
   for(i=0;i<ncindexsize(grp->children);i++) {
      if((child_grp = (NC_GRP_INFO_T*)ncindexith(grp->children,i)) == NULL) continue;
      if ((retval = nc4_rec_write_groups_types(child_grp)))
         return retval;
   }
   return NC_NOERR;
}

/**
 * @internal Copy data from one buffer to another, performing
 * appropriate data conversion.
 *
 * This function will copy data from one buffer to another, in
 * accordance with the types. Range errors will be noted, and the fill
 * value used (or the default fill value if none is supplied) for
 * values that overflow the type.
 *
 * @note I should be able to take this out when HDF5 does the right thing
 * with data type conversion. Ed Hartnett, 11/15/3
 *
 * @param src Pointer to source of data.
 * @param dest Pointer that gets data.
 * @param src_type Type ID of source data.
 * @param dest_type Type ID of destination data.
 * @param len Number of elements of data to copy.
 * @param range_error Pointer that gets 1 if there was a range error.
 * @param fill_value The fill value.
 * @param strict_nc3 Non-zero if strict model in effect.
 *
 * @returns NC_NOERR No error.
 * @returns NC_EBADTYPE Type not found.
 * @author Ed Hartnett, Dennis Heimbigner
 */
int
nc4_convert_type(const void *src, void *dest, const nc_type src_type,
                 const nc_type dest_type, const size_t len, int *range_error,
                 const void *fill_value, int strict_nc3)
{
   char *cp, *cp1;
   float *fp, *fp1;
   double *dp, *dp1;
   int *ip, *ip1;
   short *sp, *sp1;
   signed char *bp, *bp1;
   unsigned char *ubp, *ubp1;
   unsigned short *usp, *usp1;
   unsigned int *uip, *uip1;
   long long *lip, *lip1;
   unsigned long long *ulip, *ulip1;
   size_t count = 0;

   *range_error = 0;
   LOG((3, "%s: len %d src_type %d dest_type %d", __func__, len, src_type,
        dest_type));

   /* OK, this is ugly. If you can think of anything better, I'm open
      to suggestions!

      Note that we don't use a default fill value for type
      NC_BYTE. This is because Lord Voldemort cast a nofilleramous spell
      at Harry Potter, but it bounced off his scar and hit the netcdf-4
      code.
   */
   switch (src_type)
   {
   case NC_CHAR:
      switch (dest_type)
      {
      case NC_CHAR:
         for (cp = (char *)src, cp1 = dest; count < len; count++)
            *cp1++ = *cp++;
         break;
      default:
         LOG((0, "%s: Unknown destination type.", __func__));
      }
      break;

   case NC_BYTE:
      switch (dest_type)
      {
      case NC_BYTE:
         for (bp = (signed char *)src, bp1 = dest; count < len; count++)
            *bp1++ = *bp++;
         break;
      case NC_UBYTE:
         for (bp = (signed char *)src, ubp = dest; count < len; count++)
         {
            if (*bp < 0)
               (*range_error)++;
            *ubp++ = *bp++;
         }
         break;
      case NC_SHORT:
         for (bp = (signed char *)src, sp = dest; count < len; count++)
            *sp++ = *bp++;
         break;
      case NC_USHORT:
         for (bp = (signed char *)src, usp = dest; count < len; count++)
         {
            if (*bp < 0)
               (*range_error)++;
            *usp++ = *bp++;
         }
         break;
      case NC_INT:
         for (bp = (signed char *)src, ip = dest; count < len; count++)
            *ip++ = *bp++;
         break;
      case NC_UINT:
         for (bp = (signed char *)src, uip = dest; count < len; count++)
         {
            if (*bp < 0)
               (*range_error)++;
            *uip++ = *bp++;
         }
         break;
      case NC_INT64:
         for (bp = (signed char *)src, lip = dest; count < len; count++)
            *lip++ = *bp++;
         break;
      case NC_UINT64:
         for (bp = (signed char *)src, ulip = dest; count < len; count++)
         {
            if (*bp < 0)
               (*range_error)++;
            *ulip++ = *bp++;
         }
         break;
      case NC_FLOAT:
         for (bp = (signed char *)src, fp = dest; count < len; count++)
            *fp++ = *bp++;
         break;
      case NC_DOUBLE:
         for (bp = (signed char *)src, dp = dest; count < len; count++)
            *dp++ = *bp++;
         break;
      default:
         LOG((0, "%s: unexpected dest type. src_type %d, dest_type %d",
              __func__, src_type, dest_type));
         return NC_EBADTYPE;
      }
      break;

   case NC_UBYTE:
      switch (dest_type)
      {
      case NC_BYTE:
         for (ubp = (unsigned char *)src, bp = dest; count < len; count++)
         {
            if (!strict_nc3 && *ubp > X_SCHAR_MAX)
               (*range_error)++;
            *bp++ = *ubp++;
         }
         break;
      case NC_SHORT:
         for (ubp = (unsigned char *)src, sp = dest; count < len; count++)
            *sp++ = *ubp++;
         break;
      case NC_UBYTE:
         for (ubp = (unsigned char *)src, ubp1 = dest; count < len; count++)
            *ubp1++ = *ubp++;
         break;
      case NC_USHORT:
         for (ubp = (unsigned char *)src, usp = dest; count < len; count++)
            *usp++ = *ubp++;
         break;
      case NC_INT:
         for (ubp = (unsigned char *)src, ip = dest; count < len; count++)
            *ip++ = *ubp++;
         break;
      case NC_UINT:
         for (ubp = (unsigned char *)src, uip = dest; count < len; count++)
            *uip++ = *ubp++;
         break;
      case NC_INT64:
         for (ubp = (unsigned char *)src, lip = dest; count < len; count++)
            *lip++ = *ubp++;
         break;
      case NC_UINT64:
         for (ubp = (unsigned char *)src, ulip = dest; count < len; count++)
            *ulip++ = *ubp++;
         break;
      case NC_FLOAT:
         for (ubp = (unsigned char *)src, fp = dest; count < len; count++)
            *fp++ = *ubp++;
         break;
      case NC_DOUBLE:
         for (ubp = (unsigned char *)src, dp = dest; count < len; count++)
            *dp++ = *ubp++;
         break;
      default:
         LOG((0, "%s: unexpected dest type. src_type %d, dest_type %d",
              __func__, src_type, dest_type));
         return NC_EBADTYPE;
      }
      break;

   case NC_SHORT:
      switch (dest_type)
      {
      case NC_UBYTE:
         for (sp = (short *)src, ubp = dest; count < len; count++)
         {
            if (*sp > X_UCHAR_MAX || *sp < 0)
               (*range_error)++;
            *ubp++ = *sp++;
         }
         break;
      case NC_BYTE:
         for (sp = (short *)src, bp = dest; count < len; count++)
         {
            if (*sp > X_SCHAR_MAX || *sp < X_SCHAR_MIN)
               (*range_error)++;
            *bp++ = *sp++;
         }
         break;
      case NC_SHORT:
         for (sp = (short *)src, sp1 = dest; count < len; count++)
            *sp1++ = *sp++;
         break;
      case NC_USHORT:
         for (sp = (short *)src, usp = dest; count < len; count++)
         {
            if (*sp < 0)
               (*range_error)++;
            *usp++ = *sp++;
         }
         break;
      case NC_INT:
         for (sp = (short *)src, ip = dest; count < len; count++)
            *ip++ = *sp++;
         break;
      case NC_UINT:
         for (sp = (short *)src, uip = dest; count < len; count++)
         {
            if (*sp < 0)
               (*range_error)++;
            *uip++ = *sp++;
         }
         break;
      case NC_INT64:
         for (sp = (short *)src, lip = dest; count < len; count++)
            *lip++ = *sp++;
         break;
      case NC_UINT64:
         for (sp = (short *)src, ulip = dest; count < len; count++)
         {
            if (*sp < 0)
               (*range_error)++;
            *ulip++ = *sp++;
         }
         break;
      case NC_FLOAT:
         for (sp = (short *)src, fp = dest; count < len; count++)
            *fp++ = *sp++;
         break;
      case NC_DOUBLE:
         for (sp = (short *)src, dp = dest; count < len; count++)
            *dp++ = *sp++;
         break;
      default:
         LOG((0, "%s: unexpected dest type. src_type %d, dest_type %d",
              __func__, src_type, dest_type));
         return NC_EBADTYPE;
      }
      break;

   case NC_USHORT:
      switch (dest_type)
      {
      case NC_UBYTE:
         for (usp = (unsigned short *)src, ubp = dest; count < len; count++)
         {
            if (*usp > X_UCHAR_MAX)
               (*range_error)++;
            *ubp++ = *usp++;
         }
         break;
      case NC_BYTE:
         for (usp = (unsigned short *)src, bp = dest; count < len; count++)
         {
            if (*usp > X_SCHAR_MAX)
               (*range_error)++;
            *bp++ = *usp++;
         }
         break;
      case NC_SHORT:
         for (usp = (unsigned short *)src, sp = dest; count < len; count++)
         {
            if (*usp > X_SHORT_MAX)
               (*range_error)++;
            *sp++ = *usp++;
         }
         break;
      case NC_USHORT:
         for (usp = (unsigned short *)src, usp1 = dest; count < len; count++)
            *usp1++ = *usp++;
         break;
      case NC_INT:
         for (usp = (unsigned short *)src, ip = dest; count < len; count++)
            *ip++ = *usp++;
         break;
      case NC_UINT:
         for (usp = (unsigned short *)src, uip = dest; count < len; count++)
            *uip++ = *usp++;
         break;
      case NC_INT64:
         for (usp = (unsigned short *)src, lip = dest; count < len; count++)
            *lip++ = *usp++;
         break;
      case NC_UINT64:
         for (usp = (unsigned short *)src, ulip = dest; count < len; count++)
            *ulip++ = *usp++;
         break;
      case NC_FLOAT:
         for (usp = (unsigned short *)src, fp = dest; count < len; count++)
            *fp++ = *usp++;
         break;
      case NC_DOUBLE:
         for (usp = (unsigned short *)src, dp = dest; count < len; count++)
            *dp++ = *usp++;
         break;
      default:
         LOG((0, "%s: unexpected dest type. src_type %d, dest_type %d",
              __func__, src_type, dest_type));
         return NC_EBADTYPE;
      }
      break;

   case NC_INT:
      switch (dest_type)
      {
      case NC_UBYTE:
         for (ip = (int *)src, ubp = dest; count < len; count++)
         {
            if (*ip > X_UCHAR_MAX || *ip < 0)
               (*range_error)++;
            *ubp++ = *ip++;
         }
         break;
      case NC_BYTE:
         for (ip = (int *)src, bp = dest; count < len; count++)
         {
            if (*ip > X_SCHAR_MAX || *ip < X_SCHAR_MIN)
               (*range_error)++;
            *bp++ = *ip++;
         }
         break;
      case NC_SHORT:
         for (ip = (int *)src, sp = dest; count < len; count++)
         {
            if (*ip > X_SHORT_MAX || *ip < X_SHORT_MIN)
               (*range_error)++;
            *sp++ = *ip++;
         }
         break;
      case NC_USHORT:
         for (ip = (int *)src, usp = dest; count < len; count++)
         {
            if (*ip > X_USHORT_MAX || *ip < 0)
               (*range_error)++;
            *usp++ = *ip++;
         }
         break;
      case NC_INT: /* src is int */
         for (ip = (int *)src, ip1 = dest; count < len; count++)
         {
            if (*ip > X_INT_MAX || *ip < X_INT_MIN)
               (*range_error)++;
            *ip1++ = *ip++;
         }
         break;
      case NC_UINT:
         for (ip = (int *)src, uip = dest; count < len; count++)
         {
            if (*ip > X_UINT_MAX || *ip < 0)
               (*range_error)++;
            *uip++ = *ip++;
         }
         break;
      case NC_INT64:
         for (ip = (int *)src, lip = dest; count < len; count++)
            *lip++ = *ip++;
         break;
      case NC_UINT64:
         for (ip = (int *)src, ulip = dest; count < len; count++)
         {
            if (*ip < 0)
               (*range_error)++;
            *ulip++ = *ip++;
         }
         break;
      case NC_FLOAT:
         for (ip = (int *)src, fp = dest; count < len; count++)
            *fp++ = *ip++;
         break;
      case NC_DOUBLE:
         for (ip = (int *)src, dp = dest; count < len; count++)
            *dp++ = *ip++;
         break;
      default:
         LOG((0, "%s: unexpected dest type. src_type %d, dest_type %d",
              __func__, src_type, dest_type));
         return NC_EBADTYPE;
      }
      break;

   case NC_UINT:
      switch (dest_type)
      {
      case NC_UBYTE:
         for (uip = (unsigned int *)src, ubp = dest; count < len; count++)
         {
            if (*uip > X_UCHAR_MAX)
               (*range_error)++;
            *ubp++ = *uip++;
         }
         break;
      case NC_BYTE:
         for (uip = (unsigned int *)src, bp = dest; count < len; count++)
         {
            if (*uip > X_SCHAR_MAX)
               (*range_error)++;
            *bp++ = *uip++;
         }
         break;
      case NC_SHORT:
         for (uip = (unsigned int *)src, sp = dest; count < len; count++)
         {
            if (*uip > X_SHORT_MAX)
               (*range_error)++;
            *sp++ = *uip++;
         }
         break;
      case NC_USHORT:
         for (uip = (unsigned int *)src, usp = dest; count < len; count++)
         {
            if (*uip > X_USHORT_MAX)
               (*range_error)++;
            *usp++ = *uip++;
         }
         break;
      case NC_INT:
         for (uip = (unsigned int *)src, ip = dest; count < len; count++)
         {
            if (*uip > X_INT_MAX)
               (*range_error)++;
            *ip++ = *uip++;
         }
         break;
      case NC_UINT:
         for (uip = (unsigned int *)src, uip1 = dest; count < len; count++)
         {
            if (*uip > X_UINT_MAX)
               (*range_error)++;
            *uip1++ = *uip++;
         }
         break;
      case NC_INT64:
         for (uip = (unsigned int *)src, lip = dest; count < len; count++)
            *lip++ = *uip++;
         break;
      case NC_UINT64:
         for (uip = (unsigned int *)src, ulip = dest; count < len; count++)
            *ulip++ = *uip++;
         break;
      case NC_FLOAT:
         for (uip = (unsigned int *)src, fp = dest; count < len; count++)
            *fp++ = *uip++;
         break;
      case NC_DOUBLE:
         for (uip = (unsigned int *)src, dp = dest; count < len; count++)
            *dp++ = *uip++;
         break;
      default:
         LOG((0, "%s: unexpected dest type. src_type %d, dest_type %d",
              __func__, src_type, dest_type));
         return NC_EBADTYPE;
      }
      break;

   case NC_INT64:
      switch (dest_type)
      {
      case NC_UBYTE:
         for (lip = (long long *)src, ubp = dest; count < len; count++)
         {
            if (*lip > X_UCHAR_MAX || *lip < 0)
               (*range_error)++;
            *ubp++ = *lip++;
         }
         break;
      case NC_BYTE:
         for (lip = (long long *)src, bp = dest; count < len; count++)
         {
            if (*lip > X_SCHAR_MAX || *lip < X_SCHAR_MIN)
               (*range_error)++;
            *bp++ = *lip++;
         }
         break;
      case NC_SHORT:
         for (lip = (long long *)src, sp = dest; count < len; count++)
         {
            if (*lip > X_SHORT_MAX || *lip < X_SHORT_MIN)
               (*range_error)++;
            *sp++ = *lip++;
         }
         break;
      case NC_USHORT:
         for (lip = (long long *)src, usp = dest; count < len; count++)
         {
            if (*lip > X_USHORT_MAX || *lip < 0)
               (*range_error)++;
            *usp++ = *lip++;
         }
         break;
      case NC_UINT:
         for (lip = (long long *)src, uip = dest; count < len; count++)
         {
            if (*lip > X_UINT_MAX || *lip < 0)
               (*range_error)++;
            *uip++ = *lip++;
         }
         break;
      case NC_INT:
         for (lip = (long long *)src, ip = dest; count < len; count++)
         {
            if (*lip > X_INT_MAX || *lip < X_INT_MIN)
               (*range_error)++;
            *ip++ = *lip++;
         }
         break;
      case NC_INT64:
         for (lip = (long long *)src, lip1 = dest; count < len; count++)
            *lip1++ = *lip++;
         break;
      case NC_UINT64:
         for (lip = (long long *)src, ulip = dest; count < len; count++)
         {
            if (*lip < 0)
               (*range_error)++;
            *ulip++ = *lip++;
         }
         break;
      case NC_FLOAT:
         for (lip = (long long *)src, fp = dest; count < len; count++)
            *fp++ = *lip++;
         break;
      case NC_DOUBLE:
         for (lip = (long long *)src, dp = dest; count < len; count++)
            *dp++ = *lip++;
         break;
      default:
         LOG((0, "%s: unexpected dest type. src_type %d, dest_type %d",
              __func__, src_type, dest_type));
         return NC_EBADTYPE;
      }
      break;

   case NC_UINT64:
      switch (dest_type)
      {
      case NC_UBYTE:
         for (ulip = (unsigned long long *)src, ubp = dest; count < len; count++)
         {
            if (*ulip > X_UCHAR_MAX)
               (*range_error)++;
            *ubp++ = *ulip++;
         }
         break;
      case NC_BYTE:
         for (ulip = (unsigned long long *)src, bp = dest; count < len; count++)
         {
            if (*ulip > X_SCHAR_MAX)
               (*range_error)++;
            *bp++ = *ulip++;
         }
         break;
      case NC_SHORT:
         for (ulip = (unsigned long long *)src, sp = dest; count < len; count++)
         {
            if (*ulip > X_SHORT_MAX)
               (*range_error)++;
            *sp++ = *ulip++;
         }
         break;
      case NC_USHORT:
         for (ulip = (unsigned long long *)src, usp = dest; count < len; count++)
         {
            if (*ulip > X_USHORT_MAX)
               (*range_error)++;
            *usp++ = *ulip++;
         }
         break;
      case NC_UINT:
         for (ulip = (unsigned long long *)src, uip = dest; count < len; count++)
         {
            if (*ulip > X_UINT_MAX)
               (*range_error)++;
            *uip++ = *ulip++;
         }
         break;
      case NC_INT:
         for (ulip = (unsigned long long *)src, ip = dest; count < len; count++)
         {
            if (*ulip > X_INT_MAX)
               (*range_error)++;
            *ip++ = *ulip++;
         }
         break;
      case NC_INT64:
         for (ulip = (unsigned long long *)src, lip = dest; count < len; count++)
         {
            if (*ulip > X_INT64_MAX)
               (*range_error)++;
            *lip++ = *ulip++;
         }
         break;
      case NC_UINT64:
         for (ulip = (unsigned long long *)src, ulip1 = dest; count < len; count++)
            *ulip1++ = *ulip++;
         break;
      case NC_FLOAT:
         for (ulip = (unsigned long long *)src, fp = dest; count < len; count++)
            *fp++ = *ulip++;
         break;
      case NC_DOUBLE:
         for (ulip = (unsigned long long *)src, dp = dest; count < len; count++)
            *dp++ = *ulip++;
         break;
      default:
         LOG((0, "%s: unexpected dest type. src_type %d, dest_type %d",
              __func__, src_type, dest_type));
         return NC_EBADTYPE;
      }
      break;

   case NC_FLOAT:
      switch (dest_type)
      {
      case NC_UBYTE:
         for (fp = (float *)src, ubp = dest; count < len; count++)
         {
            if (*fp > X_UCHAR_MAX || *fp < 0)
               (*range_error)++;
            *ubp++ = *fp++;
         }
         break;
      case NC_BYTE:
         for (fp = (float *)src, bp = dest; count < len; count++)
         {
            if (*fp > (double)X_SCHAR_MAX || *fp < (double)X_SCHAR_MIN)
               (*range_error)++;
            *bp++ = *fp++;
         }
         break;
      case NC_SHORT:
         for (fp = (float *)src, sp = dest; count < len; count++)
         {
            if (*fp > (double)X_SHORT_MAX || *fp < (double)X_SHORT_MIN)
               (*range_error)++;
            *sp++ = *fp++;
         }
         break;
      case NC_USHORT:
         for (fp = (float *)src, usp = dest; count < len; count++)
         {
            if (*fp > X_USHORT_MAX || *fp < 0)
               (*range_error)++;
            *usp++ = *fp++;
         }
         break;
      case NC_UINT:
         for (fp = (float *)src, uip = dest; count < len; count++)
         {
            if (*fp > X_UINT_MAX || *fp < 0)
               (*range_error)++;
            *uip++ = *fp++;
         }
         break;
      case NC_INT:
         for (fp = (float *)src, ip = dest; count < len; count++)
         {
            if (*fp > (double)X_INT_MAX || *fp < (double)X_INT_MIN)
               (*range_error)++;
            *ip++ = *fp++;
         }
         break;
      case NC_INT64:
         for (fp = (float *)src, lip = dest; count < len; count++)
         {
            if (*fp > X_INT64_MAX || *fp <X_INT64_MIN)
               (*range_error)++;
            *lip++ = *fp++;
         }
         break;
      case NC_UINT64:
         for (fp = (float *)src, lip = dest; count < len; count++)
         {
            if (*fp > X_UINT64_MAX || *fp < 0)
               (*range_error)++;
            *lip++ = *fp++;
         }
         break;
      case NC_FLOAT:
         for (fp = (float *)src, fp1 = dest; count < len; count++)
         {
            /*                if (*fp > X_FLOAT_MAX || *fp < X_FLOAT_MIN)
                              (*range_error)++;*/
            *fp1++ = *fp++;
         }
         break;
      case NC_DOUBLE:
         for (fp = (float *)src, dp = dest; count < len; count++)
            *dp++ = *fp++;
         break;
      default:
         LOG((0, "%s: unexpected dest type. src_type %d, dest_type %d",
              __func__, src_type, dest_type));
         return NC_EBADTYPE;
      }
      break;

   case NC_DOUBLE:
      switch (dest_type)
      {
      case NC_UBYTE:
         for (dp = (double *)src, ubp = dest; count < len; count++)
         {
            if (*dp > X_UCHAR_MAX || *dp < 0)
               (*range_error)++;
            *ubp++ = *dp++;
         }
         break;
      case NC_BYTE:
         for (dp = (double *)src, bp = dest; count < len; count++)
         {
            if (*dp > X_SCHAR_MAX || *dp < X_SCHAR_MIN)
               (*range_error)++;
            *bp++ = *dp++;
         }
         break;
      case NC_SHORT:
         for (dp = (double *)src, sp = dest; count < len; count++)
         {
            if (*dp > X_SHORT_MAX || *dp < X_SHORT_MIN)
               (*range_error)++;
            *sp++ = *dp++;
         }
         break;
      case NC_USHORT:
         for (dp = (double *)src, usp = dest; count < len; count++)
         {
            if (*dp > X_USHORT_MAX || *dp < 0)
               (*range_error)++;
            *usp++ = *dp++;
         }
         break;
      case NC_UINT:
         for (dp = (double *)src, uip = dest; count < len; count++)
         {
            if (*dp > X_UINT_MAX || *dp < 0)
               (*range_error)++;
            *uip++ = *dp++;
         }
         break;
      case NC_INT:
         for (dp = (double *)src, ip = dest; count < len; count++)
         {
            if (*dp > X_INT_MAX || *dp < X_INT_MIN)
               (*range_error)++;
            *ip++ = *dp++;
         }
         break;
      case NC_INT64:
         for (dp = (double *)src, lip = dest; count < len; count++)
         {
            if (*dp > X_INT64_MAX || *dp < X_INT64_MIN)
               (*range_error)++;
            *lip++ = *dp++;
         }
         break;
      case NC_UINT64:
         for (dp = (double *)src, lip = dest; count < len; count++)
         {
            if (*dp > X_UINT64_MAX || *dp < 0)
               (*range_error)++;
            *lip++ = *dp++;
         }
         break;
      case NC_FLOAT:
         for (dp = (double *)src, fp = dest; count < len; count++)
         {
            if (*dp > X_FLOAT_MAX || *dp < X_FLOAT_MIN)
               (*range_error)++;
            *fp++ = *dp++;
         }
         break;
      case NC_DOUBLE:
         for (dp = (double *)src, dp1 = dest; count < len; count++)
         {
            /* if (*dp > X_DOUBLE_MAX || *dp < X_DOUBLE_MIN) */
            /*    (*range_error)++; */
            *dp1++ = *dp++;
         }
         break;
      default:
         LOG((0, "%s: unexpected dest type. src_type %d, dest_type %d",
              __func__, src_type, dest_type));
         return NC_EBADTYPE;
      }
      break;

   default:
      LOG((0, "%s: unexpected src type. src_type %d, dest_type %d",
           __func__, src_type, dest_type));
      return NC_EBADTYPE;
   }
   return NC_NOERR;
}

/**
 * @internal In our first pass through the data, we may have
 * encountered variables before encountering their dimscales, so go
 * through the vars in this file and make sure we've got a dimid for
 * each.
 *
 * @param grp Pointer to group info struct.
 *
 * @returns NC_NOERR No error.
 * @returns NC_EHDFERR HDF5 returned an error.
 * @returns NC_ENOMEM Out of memory.
 * @author Ed Hartnett
 */
int
nc4_rec_match_dimscales(NC_GRP_INFO_T *grp)
{
   NC_GRP_INFO_T *g;
   NC_VAR_INFO_T *var;
   NC_DIM_INFO_T *dim;
   int retval = NC_NOERR;
   int i;

   assert(grp && grp->hdr.name);
   LOG((4, "%s: grp->hdr.name %s", __func__, grp->hdr.name));

   /* Perform var dimscale match for child groups. */
   for(i=0;i<ncindexsize(grp->children);i++) {
      if((g = (NC_GRP_INFO_T*)ncindexith(grp->children,i)) == NULL) continue;
      if ((retval = nc4_rec_match_dimscales(g)))
         return retval;
   }
   /* Check all the vars in this group. If they have dimscale info,
    * try and find a dimension for them. */
   for(i=0;i<ncindexsize(grp->vars);i++)
   {
      int ndims;
      int d;
      if((var = (NC_VAR_INFO_T*)ncindexith(grp->vars,i)) == NULL) continue;
      /* Check all vars and see if dim[i] != NULL if dimids[i] valid. */
      /* This loop is very odd. Under normal circumstances, var->dimid[d] is zero
         (from the initial calloc) which is a legitimate dimid. The code does not
         distinquish this case from the dimscale case where the id might actually
         be defined.
         The original nc4_find_dim searched up the group tree looking for the given
         dimid in one of the dim lists associated with each ancestor group.
         I changed nc4_fnd_dim to use the dimid directly using h5->alldims.
         However, here that is incorrect because it will find the dimid 0 always
         (if any dimensions were defined). Except that when dimscale dimids have
         been defined, one or more of the values in var->dimids will have a
         legitimate value.
         The solution I choose is to modify nc4_var_list_add to initialize dimids to
         illegal values (-1). This is another example of the problems with dimscales.
      */
      ndims = var->ndims;
      for (d = 0; d < ndims; d++)
      {
         if (var->dim[d] == NULL) {
            nc4_find_dim(grp, var->dimids[d], &var->dim[d], NULL);
         }
         /*       assert(var->dim[d] && var->dim[d]->hdr.id == var->dimids[d]); */
      }

      /* Skip dimension scale variables */
      if (!var->dimscale)
      {
         int d;
         int j;

         /* Are there dimscales for this variable? */
         if (var->dimscale_hdf5_objids)
         {
            for (d = 0; d < var->ndims; d++)
            {
               nc_bool_t finished = NC_FALSE;
               LOG((5, "%s: var %s has dimscale info...", __func__, var->hdr.name));

               /* Check this and parent groups. */
               for (g = grp; g && !finished; g = g->parent)
               {
                  /* Check all dims in this group. */
                  for (j = 0; j < ncindexsize(g->dim); j++)
                  {
                     /* Get the HDF5 specific dim info. */
                     NC_HDF5_DIM_INFO_T *hdf5_dim;
                     dim = (NC_DIM_INFO_T *)ncindexith(g->dim, j);
                     assert(dim && dim->format_dim_info);
                     hdf5_dim = (NC_HDF5_DIM_INFO_T *)dim->format_dim_info;

                     /* Check for exact match of fileno/objid arrays
                      * to find identical objects in HDF5 file. */
                     if (var->dimscale_hdf5_objids[d].fileno[0] == hdf5_dim->hdf5_objid.fileno[0] &&
                         var->dimscale_hdf5_objids[d].objno[0] == hdf5_dim->hdf5_objid.objno[0] &&
                         var->dimscale_hdf5_objids[d].fileno[1] == hdf5_dim->hdf5_objid.fileno[1] &&
                         var->dimscale_hdf5_objids[d].objno[1] == hdf5_dim->hdf5_objid.objno[1])
                     {
                        LOG((4, "%s: for dimension %d, found dim %s", __func__,
                             d, dim->hdr.name));
                        var->dimids[d] = dim->hdr.id;
                        var->dim[d] = dim;
                        finished = NC_TRUE;
                        break;
                     }
                  } /* next dim */
               } /* next grp */
               LOG((5, "%s: dimid for this dimscale is %d", __func__,
                    var->type_info->hdr.id));
            } /* next var->dim */
         }
         /* No dimscales for this var! Invent phony dimensions. */
         else
         {
            hid_t spaceid = 0;
            hsize_t *h5dimlen = NULL, *h5dimlenmax = NULL;
            int dataset_ndims;

            /* Find the space information for this dimension. */
            if ((spaceid = H5Dget_space(var->hdf_datasetid)) < 0)
               return NC_EHDFERR;

            /* Get the len of each dim in the space. */
            if (var->ndims)
            {
               if (!(h5dimlen = malloc(var->ndims * sizeof(hsize_t))))
                  return NC_ENOMEM;
               if (!(h5dimlenmax = malloc(var->ndims * sizeof(hsize_t))))
               {
                  free(h5dimlen);
                  return NC_ENOMEM;
               }
               if ((dataset_ndims = H5Sget_simple_extent_dims(spaceid, h5dimlen,
                                                              h5dimlenmax)) < 0) {
                  free(h5dimlenmax);
                  free(h5dimlen);
                  return NC_EHDFERR;
               }
               if (dataset_ndims != var->ndims) {
                  free(h5dimlenmax);
                  free(h5dimlen);
                  return NC_EHDFERR;
               }
            }
            else
            {
               /* Make sure it's scalar. */
               if (H5Sget_simple_extent_type(spaceid) != H5S_SCALAR)
                  return NC_EHDFERR;
            }

            /* Release the space object. */
            if (H5Sclose(spaceid) < 0) {
               free(h5dimlen);
               free(h5dimlenmax);
               return NC_EHDFERR;
            }

            /* Create a phony dimension for each dimension in the
             * dataset, unless there already is one the correct
             * size. */
            for (d = 0; d < var->ndims; d++)
            {
               int k;
               int match;
               /* Is there already a phony dimension of the correct size? */
               for(match=-1,k=0;k<ncindexsize(grp->dim);k++) {
                  if((dim = (NC_DIM_INFO_T*)ncindexith(grp->dim,k)) == NULL) continue;
                  if ((dim->len == h5dimlen[d]) &&
                      ((h5dimlenmax[d] == H5S_UNLIMITED && dim->unlimited) ||
                       (h5dimlenmax[d] != H5S_UNLIMITED && !dim->unlimited)))
                  {match = k; break;}
               }

               /* Didn't find a phony dim? Then create one. */
               if (match < 0)
               {
                  char phony_dim_name[NC_MAX_NAME + 1];
                  sprintf(phony_dim_name, "phony_dim_%d", grp->nc4_info->next_dimid);
                  LOG((3, "%s: creating phony dim for var %s", __func__, var->hdr.name));
                  if ((retval = nc4_dim_list_add(grp, phony_dim_name, h5dimlen[d], -1, &dim)))
                  {
                     free(h5dimlenmax);
                     free(h5dimlen);
                     return retval;
                  }
                  /* Create struct for HDF5-specific dim info. */
                  if (!(dim->format_dim_info = calloc(1, sizeof(NC_HDF5_DIM_INFO_T))))
                     return NC_ENOMEM;
                  if (h5dimlenmax[d] == H5S_UNLIMITED)
                     dim->unlimited = NC_TRUE;
               }

               /* The variable must remember the dimid. */
               var->dimids[d] = dim->hdr.id;
               var->dim[d] = dim;
            } /* next dim */

              /* Free the memory we malloced. */
            free(h5dimlen);
            free(h5dimlenmax);
         }
      }
   }

   return retval;
}

/**
 * @internal Get the length, in bytes, of one element of a type in
 * memory.
 *
 * @param h5 Pointer to HDF5 file info struct.
 * @param xtype NetCDF type ID.
 * @param len Pointer that gets length in bytes.
 *
 * @returns NC_NOERR No error.
 * @returns NC_EBADTYPE Type not found
 * @author Ed Hartnett
 */
int
nc4_get_typelen_mem(NC_FILE_INFO_T *h5, nc_type xtype, size_t *len)
{
   NC_TYPE_INFO_T *type;
   int retval;

   LOG((4, "%s xtype: %d", __func__, xtype));
   assert(len);

   /* If this is an atomic type, the answer is easy. */
   switch (xtype)
   {
   case NC_BYTE:
   case NC_CHAR:
   case NC_UBYTE:
      *len = sizeof(char);
      return NC_NOERR;
   case NC_SHORT:
   case NC_USHORT:
      *len = sizeof(short);
      return NC_NOERR;
   case NC_INT:
   case NC_UINT:
      *len = sizeof(int);
      return NC_NOERR;
   case NC_FLOAT:
      *len = sizeof(float);
      return NC_NOERR;
   case NC_DOUBLE:
      *len = sizeof(double);
      return NC_NOERR;
   case NC_INT64:
   case NC_UINT64:
      *len = sizeof(long long);
      return NC_NOERR;
   case NC_STRING:
      *len = sizeof(char *);
      return NC_NOERR;
   }

   /* See if var is compound type. */
   if ((retval = nc4_find_type(h5, xtype, &type)))
      return retval;

   if (!type)
      return NC_EBADTYPE;

   *len = type->size;

   LOG((5, "type->size: %d", type->size));

   return NC_NOERR;
}

/**
 * @internal Get the class of a type
 *
 * @param h5 Pointer to the HDF5 file info struct.
 * @param xtype NetCDF type ID.
 * @param type_class Pointer that gets class of type, NC_INT,
 * NC_FLOAT, NC_CHAR, or NC_STRING, NC_ENUM, NC_VLEN, NC_COMPOUND, or
 * NC_OPAQUE.
 *
 * @return ::NC_NOERR No error.
 * @author Ed Hartnett
 */
int
nc4_get_typeclass(const NC_FILE_INFO_T *h5, nc_type xtype, int *type_class)
{
   int retval = NC_NOERR;

   LOG((4, "%s xtype: %d", __func__, xtype));
   assert(type_class);

   /* If this is an atomic type, the answer is easy. */
   if (xtype <= NC_STRING)
   {
      switch (xtype)
      {
      case NC_BYTE:
      case NC_UBYTE:
      case NC_SHORT:
      case NC_USHORT:
      case NC_INT:
      case NC_UINT:
      case NC_INT64:
      case NC_UINT64:
         /* NC_INT is class used for all integral types */
         *type_class = NC_INT;
         break;

      case NC_FLOAT:
      case NC_DOUBLE:
         /* NC_FLOAT is class used for all floating-point types */
         *type_class = NC_FLOAT;
         break;

      case NC_CHAR:
         *type_class = NC_CHAR;
         break;

      case NC_STRING:
         *type_class = NC_STRING;
         break;

      default:
         BAIL(NC_EBADTYPE);
      }
   }
   else
   {
      NC_TYPE_INFO_T *type;

      /* See if it's a used-defined type */
      if ((retval = nc4_find_type(h5, xtype, &type)))
         BAIL(retval);
      if (!type)
         BAIL(NC_EBADTYPE);

      *type_class = type->nc_type_class;
   }

exit:
   return retval;
}

/**
 * @internal Report information about an open HDF5 object. This is
 * called on any still-open objects when a HDF5 file close is
 * attempted.
 *
 * @param uselog If true, send output to LOG not stderr.
 * @param id HDF5 ID of open object.
 * @param type Type of HDF5 object, file, dataset, etc.
 *
 * @author Dennis Heimbigner
 */
void
reportobject(int uselog, hid_t id, unsigned int type)
{
   char name[NC_HDF5_MAX_NAME];
   ssize_t len;
   const char* typename = NULL;
   long long printid = (long long)id;

   len = H5Iget_name(id, name, NC_HDF5_MAX_NAME);
   if(len < 0) return;
   name[len] = '\0';

   switch (type) {
   case H5F_OBJ_FILE: typename = "File"; break;
   case H5F_OBJ_DATASET: typename = "Dataset"; break;
   case H5F_OBJ_GROUP: typename = "Group"; break;
   case H5F_OBJ_DATATYPE: typename = "Datatype"; break;
   case H5F_OBJ_ATTR:
      typename = "Attribute";
      len = H5Aget_name(id, NC_HDF5_MAX_NAME, name);
      if(len < 0) len = 0;
      name[len] = '\0';
      break;
   default: typename = "<unknown>"; break;
   }
#ifdef LOGGING
   if(uselog) {
      LOG((0,"Type = %s(%lld) name='%s'",typename,printid,name));
   } else
#endif
   {
      fprintf(stderr,"Type = %s(%lld) name='%s'",typename,printid,name);
   }

}

/**
 * @internal
 *
 * @param uselog
 * @param fid HDF5 ID.
 * @param ntypes Number of types.
 * @param otypes Pointer that gets number of open types.
 *
 * @author Dennis Heimbigner
 */
static void
reportopenobjectsT(int uselog, hid_t fid, int ntypes, unsigned int* otypes)
{
   int t,i;
   ssize_t ocount;
   size_t maxobjs = -1;
   hid_t* idlist = NULL;

   /* Always report somehow */
#ifdef LOGGING
   if(uselog)
      LOG((0,"\nReport: open objects on %lld",(long long)fid));
   else
#endif
      fprintf(stdout,"\nReport: open objects on %lld\n",(long long)fid);
   maxobjs = H5Fget_obj_count(fid,H5F_OBJ_ALL);
   if(idlist != NULL) free(idlist);
   idlist = (hid_t*)malloc(sizeof(hid_t)*maxobjs);
   for(t=0;t<ntypes;t++) {
      unsigned int ot = otypes[t];
      ocount = H5Fget_obj_ids(fid,ot,maxobjs,idlist);
      for(i=0;i<ocount;i++) {
         hid_t o = idlist[i];
         reportobject(uselog,o,ot);
      }
   }
   if(idlist != NULL) free(idlist);
}

/**
 * @internal Report open objects.
 *
 * @param uselog
 * @param fid HDF5 file ID.
 *
 * @author Dennit Heimbigner
 */
void
reportopenobjects(int uselog, hid_t fid)
{
   unsigned int OTYPES[5] = {H5F_OBJ_FILE, H5F_OBJ_DATASET, H5F_OBJ_GROUP,
                             H5F_OBJ_DATATYPE, H5F_OBJ_ATTR};

   reportopenobjectsT(uselog, fid ,5, OTYPES);
}

/**
 * @internal Report open objects given a pointer to NC_FILE_INFO_T object
 *
 * @param h5 file object
 *
 * @author Dennis Heimbigner
 */
void
showopenobjects5(NC_FILE_INFO_T* h5)
{
   NC_HDF5_FILE_INFO_T *hdf5_info;

   assert(h5 && h5->format_file_info);
   hdf5_info = (NC_HDF5_FILE_INFO_T *)h5->format_file_info;

   fprintf(stderr,"===== begin showopenobjects =====\n");
   reportopenobjects(0,hdf5_info->hdfid);
   fprintf(stderr,"===== end showopenobjects =====\n");
   fflush(stderr);
}

/**
 * @internal Report open objects given an ncid
 * Defined to support user or gdb level call.
 *
 * @param ncid file id
 *
 * @author Dennis Heimbigner
 */
void
showopenobjects(int ncid)
{
   NC_FILE_INFO_T* h5 = NULL;

   /* Find our metadata for this file. */
   if (nc4_find_nc_grp_h5(ncid, NULL, NULL, &h5) != NC_NOERR)
      fprintf(stderr,"failed\n");
   else
      showopenobjects5(h5);
   fflush(stderr);
}

/**
 * @internal Get HDF5 library version.
 *
 * @param major Pointer that gets major version number.
 * @param minor Pointer that gets minor version number.
 * @param release Pointer that gets release version number.
 *
 * @returns NC_NOERR No error.
 * @returns NC_EHDFERR HDF5 returned error.
 * @author Dennis Heimbigner
 */
int
NC4_hdf5get_libversion(unsigned* major,unsigned* minor,unsigned* release)
{
   if(H5get_libversion(major,minor,release) < 0)
      return NC_EHDFERR;
   return NC_NOERR;
}

/**
 * @internal Get HDF5 superblock version.
 *
 * @param h5 Pointer to HDF5 file info struct.
 * @param idp Pointer that gets superblock version.
 *
 * @returns NC_NOERR No error.
 * @returns NC_EHDFERR HDF5 returned error.
 * @author Dennis Heimbigner
 */
int
NC4_hdf5get_superblock(struct NC_FILE_INFO* h5, int* idp)
{
   NC_HDF5_FILE_INFO_T *hdf5_info;
   int stat = NC_NOERR;
   unsigned super;
   hid_t plist = -1;

   assert(h5 && h5->format_file_info);
   hdf5_info = (NC_HDF5_FILE_INFO_T *)h5->format_file_info;

   if((plist = H5Fget_create_plist(hdf5_info->hdfid)) < 0)
   {stat = NC_EHDFERR; goto done;}
   if(H5Pget_version(plist, &super, NULL, NULL, NULL) < 0)
   {stat = NC_EHDFERR; goto done;}
   if(idp) *idp = (int)super;
done:
   if(plist >= 0) H5Pclose(plist);
   return stat;
}

static int NC4_get_strict_att(NC_FILE_INFO_T*);
static int NC4_walk(hid_t, int*);

/**
 * @internal Determine whether file is netCDF-4.
 *
 * We define a file as being from netcdf-4 if any of the following
 * are true:
 * 1. NCPROPS attribute exists in root group
 * 2. NC3_STRICT_ATT_NAME exists in root group
 * 3. any of NC_ATT_REFERENCE_LIST, NC_ATT_CLASS,
 * NC_ATT_DIMENSION_LIST, NC_ATT_NAME,
 * NC_ATT_COORDINATES, NC_DIMID_ATT_NAME
 * exist anywhere in the file; note that this
 * requires walking the file.

 * @note WARNINGS:
 *   1. False negatives are possible for a small subset of netcdf-4
 *   created files.
 *   2. Deliberate falsification in the file can be used to cause
 *   a false positive.
 *
 * @param h5 Pointer to HDF5 file info struct.
 *
 * @returns NC_NOERR No error.
 * @author Dennis Heimbigner.
 */
int
NC4_isnetcdf4(struct NC_FILE_INFO* h5)
{
   int stat;
   int isnc4 = 0;
   int count;

   /* Look for NC3_STRICT_ATT_NAME */
   isnc4 = NC4_get_strict_att(h5);
   if(isnc4 > 0)
      goto done;
   /* attribute did not exist */
   /* => last resort: walk the HDF5 file looking for markers */
   count = 0;
   stat = NC4_walk(((NC_HDF5_GRP_INFO_T *)(h5->root_grp->format_grp_info))->hdf_grpid,
                   &count);
   if(stat != NC_NOERR)
      isnc4 = 0;
   else /* Threshold is at least two matches */
      isnc4 = (count >= 2);

done:
   return isnc4;
}

/**
 * @internal Get the NC3 strict attribute.
 *
 * @param h5 Pointer to HDF5 file info struct.
 *
 * @returns NC_NOERR No error.
 * @author Dennis Heimbigner.
 */
static int
NC4_get_strict_att(NC_FILE_INFO_T *h5)
{
   hid_t grpid = -1;
   hid_t attid = -1;

   /* Get root group ID. */
   grpid = ((NC_HDF5_GRP_INFO_T *)(h5->root_grp->format_grp_info))->hdf_grpid;

   /* Try to extract the NC3_STRICT_ATT_NAME attribute */
   attid = H5Aopen_name(grpid, NC3_STRICT_ATT_NAME);
   H5Aclose(attid);
   return attid;
}

/**
 * @internal Walk group struct.
 *
 * @param gid HDF5 ID of starting group.
 * @param countp Pointer that gets count.
 *
 * @returns NC_NOERR No error.
 * @author Dennis Heimbigner
 */
static int
NC4_walk(hid_t gid, int* countp)
{
   int ncstat = NC_NOERR;
   int i,j,na;
   ssize_t len;
   hsize_t nobj;
   herr_t err;
   int otype;
   hid_t grpid, dsid;
   char name[NC_HDF5_MAX_NAME];

   /* walk group members of interest */
   err = H5Gget_num_objs(gid, &nobj);
   if(err < 0) return err;

   for(i = 0; i < nobj; i++) {
      /* Get name & kind of object in the group */
      len = H5Gget_objname_by_idx(gid,(hsize_t)i,name,(size_t)NC_HDF5_MAX_NAME);
      if(len < 0) return len;

      otype =  H5Gget_objtype_by_idx(gid,(size_t)i);
      switch(otype) {
      case H5G_GROUP:
         grpid = H5Gopen(gid,name);
         NC4_walk(grpid,countp);
         H5Gclose(grpid);
         break;
      case H5G_DATASET: /* variables */
         /* Check for phony_dim */
         if(strcmp(name,"phony_dim")==0)
            *countp = *countp + 1;
         dsid = H5Dopen(gid,name);
         na = H5Aget_num_attrs(dsid);
         for(j = 0; j < na; j++) {
            hid_t aid =  H5Aopen_idx(dsid,(unsigned int)    j);
            if(aid >= 0) {
               const NC_reservedatt* ra;
               ssize_t len = H5Aget_name(aid, NC_HDF5_MAX_NAME, name);
               if(len < 0) return len;
               /* Is this a netcdf-4 marker attribute */
               /* Is this a netcdf-4 marker attribute */
               ra = NC_findreserved(name);
               if(ra != NULL)
                  *countp = *countp + 1;
            }
            H5Aclose(aid);
         }
         H5Dclose(dsid);
         break;
      default:/* ignore */
         break;
      }
   }
   return ncstat;
}


/* Copyright 2003-2018, University Corporation for Atmospheric
 * Research. See COPYRIGHT file for copying and redistribution
 * conditions. */
/**
 * @file
 *
 * @internal This file is part of netcdf-4, a netCDF-like interface
 * for HDF5, or a HDF5 backend for netCDF, depending on your point of
 * view.
 *
 * This file handles the nc4 attribute functions.
 *
 * Remember that with atts, type conversion can take place when
 * writing them, and when reading them.
 *
 * @author Ed Hartnett
 */

#include "nc.h"
#include "nc4internal.h"
#include "hdf5internal.h"
#include "nc4dispatch.h"
#include "ncdispatch.h"

/**
 * @internal Get special informatation about the attribute.
 *
 * @param h5 Pointer to HDF5 file info struct.
 * @param name Name of attribute.
 * @param filetypep Pointer that gets type of the attribute data in
 * file.
 * @param mem_type Type of attribute data in memory.
 * @param lenp Pointer that gets length of attribute array.
 * @param attnump Pointer that gets the attribute number.
 * @param data Attribute data.
 *
 * @return ::NC_NOERR No error.
 * @return ::NC_EBADID Bad ncid.
 * @return ::NC_ERANGE Data conversion out of range.
 * @author Dennis Heimbigner
 */
static int
nc4_get_att_special(NC_FILE_INFO_T* h5, const char* name,
                    nc_type* filetypep, nc_type mem_type, size_t* lenp,
                    int* attnump, void* data)
{
   /* Fail if asking for att id */
   if(attnump)
      return NC_EATTMETA;

   if(strcmp(name,NCPROPS)==0) {
      char* propdata = NULL;
      int stat = NC_NOERR;
      int len;
      if(h5->provenance->propattr.version == 0)
         return NC_ENOTATT;
      if(mem_type == NC_NAT) mem_type = NC_CHAR;
      if(mem_type != NC_CHAR)
         return NC_ECHAR;
      if(filetypep) *filetypep = NC_CHAR;
      stat = NC4_buildpropinfo(&h5->provenance->propattr, &propdata);
      if(stat != NC_NOERR) return stat;
      len = strlen(propdata);
      if(lenp) *lenp = len;
      if(data) strncpy((char*)data,propdata,len+1);
      free(propdata);
   } else if(strcmp(name,ISNETCDF4ATT)==0
             || strcmp(name,SUPERBLOCKATT)==0) {
      unsigned long long iv = 0;
      if(filetypep) *filetypep = NC_INT;
      if(lenp) *lenp = 1;
      if(strcmp(name,SUPERBLOCKATT)==0)
         iv = (unsigned long long)h5->provenance->superblockversion;
      else /* strcmp(name,ISNETCDF4ATT)==0 */
         iv = NC4_isnetcdf4(h5);
      if(mem_type == NC_NAT) mem_type = NC_INT;
      if(data)
         switch (mem_type) {
         case NC_BYTE: *((char*)data) = (char)iv; break;
         case NC_SHORT: *((short*)data) = (short)iv; break;
         case NC_INT: *((int*)data) = (int)iv; break;
         case NC_UBYTE: *((unsigned char*)data) = (unsigned char)iv; break;
         case NC_USHORT: *((unsigned short*)data) = (unsigned short)iv; break;
         case NC_UINT: *((unsigned int*)data) = (unsigned int)iv; break;
         case NC_INT64: *((long long*)data) = (long long)iv; break;
         case NC_UINT64: *((unsigned long long*)data) = (unsigned long long)iv; break;
         default:
            return NC_ERANGE;
         }
   }
   return NC_NOERR;
}

/**
 * @internal Get or put attribute metadata from our linked list of
 * file info. Always locate the attribute by name, never by attnum.
 * The mem_type is ignored if data=NULL.
 *
 * @param ncid File and group ID.
 * @param varid Variable ID.
 * @param name Name of attribute.
 * @param xtype Pointer that gets (file) type of attribute.
 * @param mem_type The type of attribute data in memory.
 * @param lenp Pointer that gets length of attribute array.
 * @param attnum Pointer that gets the index number of this attribute.
 * @param data Pointer that gets attribute data.
 *
 * @return ::NC_NOERR No error.
 * @return ::NC_EBADID Bad ncid.
 * @author Ed Hartnett
 */
static int
nc4_get_att(int ncid, int varid, const char *name, nc_type *xtype,
            nc_type mem_type, size_t *lenp, int *attnum, void *data)
{
   NC *nc;
   NC_GRP_INFO_T *grp;
   NC_FILE_INFO_T *h5;
   NC_ATT_INFO_T *att = NULL;
   NC_VAR_INFO_T *var;
   int my_attnum = -1;
   int need_to_convert = 0;
   int range_error = NC_NOERR;
   void *bufr = NULL;
   size_t type_size;
   char norm_name[NC_MAX_NAME + 1];
   int i;
   int retval;

   if (attnum)
      my_attnum = *attnum;

   LOG((3, "%s: ncid 0x%x varid %d name %s attnum %d mem_type %d",
        __func__, ncid, varid, name, my_attnum, mem_type));

   /* Find info for this file, group, and h5 info. */
   if ((retval = nc4_find_nc_grp_h5(ncid, &nc, &grp, &h5)))
      return retval;

   /* Check varid */
   if (varid != NC_GLOBAL)
   {
      if (!(var = (NC_VAR_INFO_T*)ncindexith(grp->vars,varid)))
         return NC_ENOTVAR;
      assert(var->hdr.id == varid);
   }

   if (name == NULL)
      BAIL(NC_EBADNAME);

   /* Normalize name. */
   if ((retval = nc4_normalize_name(name, norm_name)))
      BAIL(retval);

   /* Read the atts for this group/var, if they have not been read. */
   if (varid == NC_GLOBAL)
   {
      if (grp->atts_not_read)
         if ((retval = nc4_read_atts(grp, NULL)))
            return retval;
   }
   else
   {
      if (var->atts_not_read)
         if ((retval = nc4_read_atts(grp, var)))
            return retval;
   }

   /* If this is one of the reserved atts, use nc_get_att_special. */
   if (nc->ext_ncid == ncid && varid == NC_GLOBAL) {
      const NC_reservedatt* ra = NC_findreserved(norm_name);
      if(ra != NULL && (ra->flags & NAMEONLYFLAG))
	return nc4_get_att_special(h5, norm_name, xtype, mem_type, lenp, attnum, data);
   }

   /* Find the attribute, if it exists. */
   if ((retval = nc4_find_grp_att(grp, varid, norm_name, my_attnum, &att)))
      return retval;

   /* If mem_type is NC_NAT, it means we want to use the attribute's
    * file type as the mem type as well. */
   if (mem_type == NC_NAT)
      mem_type = att->nc_typeid;

   /* If the attribute is NC_CHAR, and the mem_type isn't, or vice
    * versa, that's a freakish attempt to convert text to
    * numbers. Some pervert out there is trying to pull a fast one!
    * Send him an NC_ECHAR error. */
   if (data && att->len)
      if ((att->nc_typeid == NC_CHAR && mem_type != NC_CHAR) ||
           (att->nc_typeid != NC_CHAR && mem_type == NC_CHAR))
         BAIL(NC_ECHAR); /* take that, you freak! */

   /* Copy the info. */
   if (lenp)
      *lenp = att->len;
   if (xtype)
      *xtype = att->nc_typeid;
   if (attnum) {
      *attnum = att->hdr.id;
   }

   /* Zero len attributes are easy to read! */
   if (!att->len)
      BAIL(NC_NOERR);

   /* Later on, we will need to know the size of this type. */
   if ((retval = nc4_get_typelen_mem(h5, mem_type, &type_size)))
      BAIL(retval);

   /* We may have to convert data. Treat NC_CHAR the same as
    * NC_UBYTE. If the mem_type is NAT, don't try any conversion - use
    * the attribute's type. */
   if (data && att->len && mem_type != att->nc_typeid &&
       mem_type != NC_NAT &&
       !(mem_type == NC_CHAR &&
         (att->nc_typeid == NC_UBYTE || att->nc_typeid == NC_BYTE)))
   {
      if (!(bufr = malloc((size_t)(att->len * type_size))))
         BAIL(NC_ENOMEM);
      need_to_convert++;
      if ((retval = nc4_convert_type(att->data, bufr, att->nc_typeid,
                                     mem_type, (size_t)att->len, &range_error,
                                     NULL, (h5->cmode & NC_CLASSIC_MODEL))))
         BAIL(retval);

      /* For strict netcdf-3 rules, ignore erange errors between UBYTE
       * and BYTE types. */
      if ((h5->cmode & NC_CLASSIC_MODEL) &&
          (att->nc_typeid == NC_UBYTE || att->nc_typeid == NC_BYTE) &&
          (mem_type == NC_UBYTE || mem_type == NC_BYTE) &&
          range_error)
         range_error = 0;
   }
   else
   {
      bufr = att->data;
   }

   /* If the caller wants data, copy it for him. If he hasn't
      allocated enough memory for it, he will burn in segmentation
      fault hell, writhing with the agony of undiscovered memory
      bugs! */
   if (data)
   {
      if (att->vldata)
      {
         size_t base_typelen;
         hvl_t *vldest = data;
         NC_TYPE_INFO_T *type;

         /* Get the type object for the attribute's type */
         if ((retval = nc4_find_type(h5, att->nc_typeid, &type)))
            BAIL(retval);

         /* Retrieve the size of the base type */
         if ((retval = nc4_get_typelen_mem(h5, type->u.v.base_nc_typeid, &base_typelen)))
            BAIL(retval);

         for (i = 0; i < att->len; i++)
         {
            vldest[i].len = att->vldata[i].len;
            if (!(vldest[i].p = malloc(vldest[i].len * base_typelen)))
               BAIL(NC_ENOMEM);
            memcpy(vldest[i].p, att->vldata[i].p, vldest[i].len * base_typelen);
         }
      }
      else if (att->stdata)
      {
         for (i = 0; i < att->len; i++)
         {
            /* Check for NULL pointer for string (valid in HDF5) */
            if(att->stdata[i])
            {
               if (!(((char **)data)[i] = strdup(att->stdata[i])))
                  BAIL(NC_ENOMEM);
            }
            else
               ((char **)data)[i] = att->stdata[i];
         }
      }
      else
      {
         memcpy(data, bufr, (size_t)(att->len * type_size));
      }
   }

exit:
   if (need_to_convert)
      free(bufr);
   if (range_error)
      retval = NC_ERANGE;
   return retval;
}

/**
 * @internal Learn about an att. All the nc4 nc_inq_ functions just
 * call nc4_get_att to get the metadata on an attribute.
 *
 * @param ncid File and group ID.
 * @param varid Variable ID.
 * @param name Name of attribute.
 * @param xtypep Pointer that gets type of attribute.
 * @param lenp Pointer that gets length of attribute data array.
 *
 * @return ::NC_NOERR No error.
 * @return ::NC_EBADID Bad ncid.
 * @author Ed Hartnett
 */
int
NC4_inq_att(int ncid, int varid, const char *name, nc_type *xtypep,
            size_t *lenp)
{
   LOG((2, "%s: ncid 0x%x varid %d name %s", __func__, ncid, varid, name));
   return nc4_get_att(ncid, varid, name, xtypep, NC_NAT, lenp, NULL, NULL);
}

/**
 * @internal Learn an attnum, given a name.
 *
 * @param ncid File and group ID.
 * @param varid Variable ID.
 * @param name Name of attribute.
 * @param attnump Pointer that gets the attribute index number.
 *
 * @return ::NC_NOERR No error.
 * @author Ed Hartnett
 */
int
NC4_inq_attid(int ncid, int varid, const char *name, int *attnump)
{
   LOG((2, "%s: ncid 0x%x varid %d name %s", __func__, ncid, varid, name));
   return nc4_get_att(ncid, varid, name, NULL, NC_NAT, NULL, attnump, NULL);
}


/**
 * @internal Given an attnum, find the att's name.
 *
 * @param ncid File and group ID.
 * @param varid Variable ID.
 * @param attnum The index number of the attribute.
 * @param name Pointer that gets name of attrribute.
 *
 * @return ::NC_NOERR No error.
 * @return ::NC_EBADID Bad ncid.
 * @author Ed Hartnett
 */
int
NC4_inq_attname(int ncid, int varid, int attnum, char *name)
{
   NC_ATT_INFO_T *att;
   int retval;

   LOG((2, "nc_inq_attname: ncid 0x%x varid %d attnum %d", ncid, varid,
        attnum));

   /* Find the attribute metadata. */
   if ((retval = nc4_find_nc_att(ncid, varid, NULL, attnum, &att)))
      return retval;

   /* Get the name. */
   if (name)
      strcpy(name, att->hdr.name);

   return NC_NOERR;
}

/**
 * @internal Get an attribute.
 *
 * @param ncid File and group ID.
 * @param varid Variable ID.
 * @param name Name of attribute.
 * @param value Pointer that gets attribute data.
 * @param memtype The type the data should be converted to as it is read.
 *
 * @return ::NC_NOERR No error.
 * @return ::NC_EBADID Bad ncid.
 * @author Ed Hartnett
 */
int
NC4_get_att(int ncid, int varid, const char *name, void *value, nc_type memtype)
{
   return nc4_get_att(ncid, varid, name, NULL, memtype, NULL, NULL, value);
}

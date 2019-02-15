/* Copyright 2018, UCAR/Unidata See netcdf/COPYRIGHT file for copying
 * and redistribution conditions.*/
/**
 * @file @internal The HDF4 file functions.
 *
 * @author Ed Hartnett
 */

#include "config.h"
#include "nc4internal.h"
#include "hdf4dispatch.h"
#include <mfhdf.h>

#define NUM_TYPES 12 /**< Number of netCDF atomic types. */

extern int nc4_vararray_add(NC_GRP_INFO_T *grp, NC_VAR_INFO_T *var);

/** @internal These flags may not be set for open mode. */
static const int
ILLEGAL_OPEN_FLAGS = (NC_MMAP|NC_64BIT_OFFSET|NC_DISKLESS|NC_WRITE);

/** @internal NetCDF atomic type names. */
static const char*
nc_type_name_g[NUM_TYPES] = {"char", "byte", "short", "int", "float", "double",
                             "ubyte", "ushort", "uint", "int64", "uint64",
                             "string"};

/** @internal NetCDF atomic type sizes. */
static const size_t
nc_type_size_g[NUM_TYPES] = {sizeof(char), sizeof(char), sizeof(short),
                             sizeof(int), sizeof(float), sizeof(double),
                             sizeof(unsigned char), sizeof(unsigned short),
                             sizeof(unsigned int), sizeof(long long),
                             sizeof(unsigned long long), sizeof(char *)};

/**
 * @internal Recursively delete the data for a group (and everything
 * it contains) in our internal metadata store.
 *
 * @param grp Pointer to group info struct.
 *
 * @return ::NC_NOERR No error.
 * @author Ed Hartnett
 */
static int
hdf4_rec_grp_del(NC_GRP_INFO_T *grp)
{
   NC_VAR_INFO_T *var;
   int i;

   assert(grp);
   LOG((3, "%s: grp->name %s", __func__, grp->hdr.name));

   /* Delete all vars. */
   for (i = 0; i < ncindexsize(grp->vars); i++)
   {
      NC_VAR_HDF4_INFO_T *hdf4_var;
      
      var = (NC_VAR_INFO_T *)ncindexith(grp->vars, i);
      if (var == NULL) continue;
      LOG((4, "%s: deleting var %s", __func__, var->hdr.name));

      /* Get the HDF4 specific var metadata. */
      hdf4_var = (NC_VAR_HDF4_INFO_T *)var->format_var_info;

      /* Close HDF4 dataset associated with this var, unless it's a
       * scale. */
      if (hdf4_var->sdsid && SDendaccess(hdf4_var->sdsid) < 0)
         return NC_EHDFERR;
      free(hdf4_var);
   }

   return NC_NOERR;
}

/**
 * @internal Given an HDF4 type, set a pointer to netcdf type.
 *
 * See http://www.hdfgroup.org/training/HDFtraining/UsersGuide/Fundmtls.fm3.html
 * for more information re: HDF4 types.
 *
 * @param h5 Pointer to HDF5 file info struct.
 * @param hdf4_typeid Type ID for hdf4 datatype.
 * @param xtype Pointer that gets netcdf type. Ignored if NULL.
 * @param endniannessp Pointer that gets endianness. Ignored if NULL.
 * @param type_sizep Pointer that gets type size. Ignored if NULL.
 * @param type_name Pointer that gets the type name. Ignored if NULL.
 *
 * @return ::NC_NOERR No error.
 * @author Ed Hartnett
 */
static int
hdf4_type_info(NC_FILE_INFO_T *h5, int32 hdf4_typeid, nc_type* xtypep,
               int *endiannessp, size_t *type_sizep, char *type_name)
{
   int t = 0;
   int endianness = NC_ENDIAN_BIG;
   nc_type xtype;

   assert(h5);

   switch(hdf4_typeid)
   {
   case DFNT_CHAR:
      xtype = NC_CHAR;
      t = 0;
      break;
   case DFNT_UCHAR:
   case DFNT_UINT8:
      xtype = NC_UBYTE;
      t = 6;
      break;
   case DFNT_LUINT8:
      xtype = NC_UBYTE;
      t = 6;
      endianness = NC_ENDIAN_LITTLE;
      break;
   case DFNT_INT8:
      xtype = NC_BYTE;
      t = 1;
      break;
   case DFNT_LINT8:
      xtype = NC_BYTE;
      t = 1;
      endianness = NC_ENDIAN_LITTLE;
      break;
   case DFNT_INT16:
      xtype = NC_SHORT;
      t = 2;
      break;
   case DFNT_LINT16:
      xtype = NC_SHORT;
      t = 2;
      endianness = NC_ENDIAN_LITTLE;
      break;
   case DFNT_UINT16:
      xtype = NC_USHORT;
      t = 7;
      break;
   case DFNT_LUINT16:
      xtype = NC_USHORT;
      t = 7;
      endianness = NC_ENDIAN_LITTLE;
      break;
   case DFNT_INT32:
      xtype = NC_INT;
      t = 3;
      break;
   case DFNT_LINT32:
      xtype = NC_INT;
      t = 3;
      endianness = NC_ENDIAN_LITTLE;
      break;
   case DFNT_UINT32:
      xtype = NC_UINT;
      t = 8;
      break;
   case DFNT_LUINT32:
      xtype = NC_UINT;
      t = 8;
      endianness = NC_ENDIAN_LITTLE;
      break;
   case DFNT_FLOAT32:
      xtype = NC_FLOAT;
      t = 4;
      break;
   case DFNT_LFLOAT32:
      xtype = NC_FLOAT;
      t = 4;
      endianness = NC_ENDIAN_LITTLE;
      break;
   case DFNT_FLOAT64:
      xtype = NC_DOUBLE;
      t = 5;
      break;
   case DFNT_LFLOAT64:
      xtype = NC_DOUBLE;
      t = 5;
      endianness = NC_ENDIAN_LITTLE;
      break;
   default:
      return NC_EBADTYPID;
   }

   /* Return results to caller. */
   if (xtypep)
      *xtypep = xtype;
   if (endiannessp)
      *endiannessp = endianness;
   if (type_sizep)
      *type_sizep = nc_type_size_g[t];
   if (type_name)
      strncpy(type_name, nc_type_name_g[t], NC_MAX_NAME);

   return NC_NOERR;
}

/**
 * @internal Set the type of a netCDF-4 variable.
 *
 * @param xtype A netcdf type.
 * @param endianness The endianness of the data.
 * @param type_size The size in bytes of one element of this type.
 * @param type_name A name for the type.
 * @param typep Pointer to a pointer that gets the TYPE_INFO_T struct.
 *
 * @return ::NC_NOERR No error.
 * @author Ed Hartnett
 */
static int
nc4_set_var_type(nc_type xtype, int endianness, size_t type_size, char *type_name,
                 NC_TYPE_INFO_T **typep)
{
   NC_TYPE_INFO_T *type;

   /* Check inputs. */
   assert(typep);

   /* Allocate space for the type info struct. */
   if (!(type = calloc(1, sizeof(NC_TYPE_INFO_T))))
      return NC_ENOMEM;
   if (!(type->hdr.name = strdup(type_name)))
   {
      free(type);
      return NC_ENOMEM;
   }
   type->hdr.sort = NCTYP;

   /* Determine the type class. */
   if (xtype == NC_FLOAT)
      type->nc_type_class = NC_FLOAT;
   else if (xtype == NC_DOUBLE)
      type->nc_type_class = NC_DOUBLE;
   else if (xtype == NC_CHAR)
      type->nc_type_class = NC_STRING;
   else
      type->nc_type_class = NC_INT;

   /* Set other type info values. */
   type->endianness = endianness;
   type->size = type_size;
   type->hdr.id = (size_t)xtype;
   type->hdr.hashkey = NC_hashmapkey(type->hdr.name, strlen(type->hdr.name));

   /* Return to caller. */
   *typep = type;

   return NC_NOERR;
}

/**
 * @internal Read an attribute from a HDF4 file.
 *
 * @param h5 Pointer to the file metadata struct.
 * @param var Pointer to variable metadata struct or NULL for global
 * attributes.
 * @param a Index of attribute to read.
 *
 * @return ::NC_NOERR No error.
 * @return ::NC_EHDFERR HDF4 error.
 * @return ::NC_EATTMETA Error reading HDF4 attribute.
 * @return ::NC_ENOMEM Out of memory.
 * @author Ed Hartnett
 */
static int
hdf4_read_att(NC_FILE_INFO_T *h5, NC_VAR_INFO_T *var, int a)
{
   NC_HDF4_FILE_INFO_T *hdf4_file;
   NC_ATT_INFO_T *att;   
   NCindex *att_list;   
   int32 att_data_type, att_count;
   size_t att_type_size;
   char name[NC_MAX_HDF4_NAME+1];
   int sd_id;
   int retval;
   
   LOG((3, "%s: a %d var %s", __func__, a, var ? var->hdr.name : "global"));

   /* Check inputs. */
   assert(h5 && h5->format_file_info);

   /* Get the HDF4 file info. */
   hdf4_file = h5->format_file_info;

   /* Decide what att list to use, global or from a var. */
   if (var)
   {
      NC_VAR_HDF4_INFO_T *hdf4_var;
      assert(var->format_var_info);
      att_list = var->att;
      hdf4_var = var->format_var_info;
      sd_id = hdf4_var->sdsid;
   } else {
      att_list = h5->root_grp->att;
      sd_id = hdf4_file->sdid;
   }

   /* Learn about this attribute. */
   if (SDattrinfo(sd_id, a, name, &att_data_type, &att_count))
      return NC_EATTMETA;
   
   /* Get information about the attribute type. */
   nc_type xtype;
   if ((retval = hdf4_type_info(h5, att_data_type, &xtype, NULL,
                                &att_type_size, NULL)))
      return retval;

   /* Add to the end of the list of atts for this var. */
   if ((retval = nc4_att_list_add(att_list, name, &att)))
      return retval;
   att->nc_typeid = xtype;
   att->created = NC_TRUE;
   att->len = att_count;
   
   /* Allocate memory to hold the data. */
   if (att->len)
      if (!(att->data = malloc(att_type_size * att->len)))
         return NC_ENOMEM;
   
   /* Read the data. */
   if (SDreadattr(sd_id, a, att->data))
      return NC_EHDFERR;

   return NC_NOERR;
}

/**
 * @internal Read a HDF4 dimension. As new dimensions are found, add
 * them to the metadata list of dimensions.
 *
 * @param h5 Pointer to the file metadata struct.
 * @param var Pointer to variable metadata struct or NULL for global
 * attributes.
 * @param rec_dim_len Actual length of first dim for this SD.
 * @param d Dimension index for this SD.
 * 
 * @return ::NC_NOERR No error.
 * @return ::NC_EHDFERR HDF4 error.
 * @return ::NC_EDIMMETA Error reading HDF4 dimension info.
 * @return ::NC_ENOMEM Out of memory.
 * @return ::NC_EMAXNAME Name too long.
 * @author Ed Hartnett
 */
static int
hdf4_read_dim(NC_FILE_INFO_T *h5, NC_VAR_INFO_T *var, int rec_dim_len, int d)
{
   NC_VAR_HDF4_INFO_T *hdf4_var;
   NC_DIM_INFO_T *dim = NULL;
   int32 dimid, dim_len, dim_data_type, dim_num_attrs;
   char dim_name[NC_MAX_NAME + 1];
   int i;
   int retval;
   
   assert(h5 && h5->format_file_info && var && var->format_var_info);
   hdf4_var = var->format_var_info;         

   /* Get information about the dimension. */
   if ((dimid = SDgetdimid(hdf4_var->sdsid, d)) == FAIL)
      return NC_EDIMMETA;
   if (SDdiminfo(dimid, dim_name, &dim_len, &dim_data_type, &dim_num_attrs))
      return NC_EDIMMETA;
   if (strlen(dim_name) > NC_MAX_HDF4_NAME)
      return NC_EMAXNAME;
   
   /* Do we already have this dimension? HDF4 explicitly uses
    * the name to tell. */
   for (i = 0; i < ncindexsize(h5->root_grp->dim); i++)
   {
      dim = (NC_DIM_INFO_T*)ncindexith(h5->root_grp->dim, i);
      if (!strcmp(dim->hdr.name, dim_name))
         break;
      dim = NULL;
   }
   
   /* If we didn't find this dimension, add one. */
   if (!dim)
   {
      LOG((4, "adding dim %s for dataset %s", dim_name, var->hdr.name));
      if ((retval = nc4_dim_list_add(h5->root_grp, dim_name,
                                     (dim_len ? dim_len : rec_dim_len), -1, &dim)))
         return retval;
   }
   
   /* Tell the variable the id of this dimension. */
   var->dimids[d] = dim->hdr.id;
   var->dim[d] = dim;

   return NC_NOERR;
}

/**
 * @internal Create a new variable and insert int relevant lists
 *
 * @param grp the containing group
 * @param name the name for the new variable
 * @param ndims the rank of the new variable
 * @param format_var_info Pointer to format-specific var info struct.
 * @param var Pointer in which to return a pointer to the new var.
 *
 * @return ::NC_NOERR No error.
 * @return ::NC_ENOMEM Out of memory.
 * @author Ed Hartnett
 */
static int
nc4_var_list_add_full(NC_GRP_INFO_T* grp, const char* name, int ndims, nc_type xtype,
                      int endianness, size_t type_size, char *type_name, void *fill_value,
                      int contiguous, size_t *chunksizes, void *format_var_info,
                      NC_VAR_INFO_T **var)
{
   int d;
   int retval;

   /* Add the VAR_INFO_T struct to our list of vars. */
   if ((retval = nc4_var_list_add(grp, name, ndims, var)))
      return retval;
   (*var)->created = NC_TRUE;
   (*var)->written_to = NC_TRUE;
   (*var)->format_var_info = format_var_info;

   /* Fill special type_info struct for variable type information. */
   if ((retval = nc4_set_var_type(xtype, endianness, type_size, type_name,
                                  &(*var)->type_info)))
      return retval;

   (*var)->type_info->rc++;
   
   /* Handle fill value, if provided. */
   if (fill_value)
   {
      if (!((*var)->fill_value = malloc(type_size)))
         return NC_ENOMEM;
      memcpy((*var)->fill_value, fill_value, type_size);
   }

   /* Var contiguous or chunked? */
   (*var)->contiguous = contiguous;
   
   /* Were chunksizes provided? */
   if (chunksizes)
   {
      if (!((*var)->chunksizes = malloc(ndims * sizeof(size_t))))
         return NC_ENOMEM;
      for (d = 0; d < ndims; d++) 
         (*var)->chunksizes[d] = chunksizes[d];
   }

   return NC_NOERR;
}

/**
 * @internal Read a HDF4 variable, including its associated dimensions
 * and attributes.
 *
 * @param h5 Pointer to the file metadata struct.
 * @param v Index of variable to read.
 * 
 * @return ::NC_NOERR No error.
 * @return ::NC_EHDFERR HDF4 error.
 * @return ::NC_EDIMMETA Error reading HDF4 dimension info.
 * @return ::NC_EVARMETA Error reading HDF4 dataset or att.
 * @return ::NC_EATTMETA Error reading HDF4 attribute.
 * @return ::NC_ENOMEM Out of memory.
 * @return ::NC_EMAXNAME Name too long.
 * @author Ed Hartnett
 */
static int
hdf4_read_var(NC_FILE_INFO_T *h5, int v)
{
   NC_HDF4_FILE_INFO_T *hdf4_file;
   NC_VAR_INFO_T *var;
   NC_VAR_HDF4_INFO_T *hdf4_var;
   HDF_CHUNK_DEF chunkdefs;
   int32 data_type, num_atts;
   int32 dimsize[NC_MAX_HDF4_DIMS];
   int32 rec_dim_len;
   int32 rank;
   int32 sdsid;
   int contiguous;
   int d, a;
   int flag;
   char name[NC_MAX_HDF4_NAME+1];
   int xtype;
   char type_name[NC_MAX_NAME + 1];
   int endianness;
   size_t type_size;
   void *fill_value;
   size_t *chunksizes = NULL;
   int retval;

   /* Check inputs. */
   assert(h5 && h5->format_file_info);

   /* Get HDF4 file metadata. */
   hdf4_file = h5->format_file_info;

   /* Open this dataset in HDF4 file. */
   if ((sdsid = SDselect(hdf4_file->sdid, v)) == FAIL)
      return  NC_EVARMETA;

   /* Learn about this dataset. */
   if (SDgetinfo(sdsid, name, &rank, dimsize, &data_type, &num_atts))
      return NC_EVARMETA;
   rec_dim_len = dimsize[0];

   /* Get chunking info from HDF4 file. */
   if (SDgetchunkinfo(sdsid, &chunkdefs, &flag))
      return NC_EVARMETA;

   /* Learn about the HDF4 type. */
   if ((retval = hdf4_type_info(h5, data_type, &xtype, &endianness, &type_size,
                                type_name)))
      return retval;

   /* Get the fill value. */
   if (!(fill_value = malloc(type_size)))
      return NC_ENOMEM;
   if (SDgetfillvalue(sdsid, fill_value))
   {
      /* Whoops! No fill value! */
      free(fill_value);
      fill_value = NULL;
   }

   /* Is variable chunked or contiguous? */
   if (flag == HDF_NONE)
      contiguous = NC_TRUE;
   else if (flag & HDF_CHUNK)
   {
      contiguous = NC_FALSE;
      if (!(chunksizes = malloc(rank * sizeof(size_t))))
         return NC_ENOMEM;
      for (d = 0; d < rank; d++) 
         chunksizes[d] = chunkdefs.chunk_lengths[d];
   }

   /* Malloc a struct to hold HDF4-specific variable
    * information. */
   if (!(hdf4_var = malloc(sizeof(NC_VAR_HDF4_INFO_T))))
      return NC_ENOMEM;

   /* Remember these values. */
   hdf4_var->hdf4_data_type = data_type;
   hdf4_var->sdsid = sdsid;   

   /* Add a variable to metadata structures. */
   LOG((3, "adding var for HDF4 dataset %s, rank %d netCDF type %d", name,
        rank, xtype));
   retval = nc4_var_list_add_full(h5->root_grp, name, (int)rank,
                                  xtype, endianness, type_size, type_name,                                       
                                  fill_value, contiguous, chunksizes, hdf4_var,
                                  &var);

   /* Free resources. */
   if (chunksizes)
      free(chunksizes);
   if (fill_value)
      free(fill_value);

   /* Did the add fail? */
   if (retval)
   {
      free(hdf4_var);
      return retval;
   }

   /* Find the variable's dimensions. */
   for (d = 0; d < var->ndims; d++)
      if ((retval = hdf4_read_dim(h5, var, rec_dim_len, d)))
         return retval;

   /* Read the variable's attributes. */
   for (a = 0; a < num_atts; a++)
      if ((retval = hdf4_read_att(h5, var, a)))
         return retval;

   return NC_NOERR;
}

/**
 * @internal Open a HDF4 file.
 *
 * @param path The file name of the file.
 * @param mode The open mode flag.
 * @param basepe Ignored by this function.
 * @param chunksizehintp Ignored by this function.
 * @param parameters pointer to struct holding extra data (e.g. for
 * parallel I/O) layer. Ignored if NULL.
 * @param dispatch Pointer to the dispatch table for this file.
 * @param nc_file Pointer to an instance of NC.
 *
 * @return ::NC_NOERR No error.
 * @return ::NC_EINVAL Invalid input.
 * @author Ed Hartnett
 */
int
NC_HDF4_open(const char *path, int mode, int basepe, size_t *chunksizehintp,
             void *parameters, NC_Dispatch *dispatch, NC *nc_file)
{
   NC_FILE_INFO_T *h5;
   NC_HDF4_FILE_INFO_T *hdf4_file;
   int32 num_datasets, num_gatts;
   int32 sdid;
   int v, a;
   NC_FILE_INFO_T* nc4_info = NULL;
   int retval = NC_NOERR;

   /* Check inputs. */
   assert(nc_file && path);

   LOG((1, "%s: path %s mode %d params %x", __func__, path, mode, parameters));

   /* Check the mode for validity */
   if (mode & ILLEGAL_OPEN_FLAGS)
      return NC_EINVAL;

   /* Open the file. */
   nc_file->int_ncid = nc_file->ext_ncid;

   /* Open the file and initialize SD interface. */
   if ((sdid = SDstart(path, DFACC_READ)) == FAIL)
      return NC_EHDFERR;

   /* Learn how many datasets and global atts we have. */
   if (SDfileinfo(sdid, &num_datasets, &num_gatts))
      return NC_EHDFERR;

   /* Add necessary structs to hold netcdf-4 file data. */
   if ((retval = nc4_nc4f_list_add(nc_file, path, mode)))
      return retval;
   nc4_info = NC4_DATA(nc_file);
   assert(nc4_info && nc4_info->root_grp);
   h5 = nc4_info;
   h5->no_write = NC_TRUE;

   /* Allocate data to hold HDF4 specific file data. */
   if (!(hdf4_file = malloc(sizeof(NC_HDF4_FILE_INFO_T))))
      return NC_ENOMEM;
   h5->format_file_info = hdf4_file;
   hdf4_file->sdid = sdid;   

   /* Read the global atts. */
   for (a = 0; a < num_gatts; a++)
      if ((retval = hdf4_read_att(h5, NULL, a)))
         break;

   /* Read each dataset. */
   if (!retval)
      for (v = 0; v < num_datasets; v++)
         if ((retval = hdf4_read_var(h5, v)))
            break;

   /* If there is an error, free resources. */
   if (retval)
      free(hdf4_file);
   
#ifdef LOGGING
   /* This will print out the names, types, lens, etc of the vars and
      atts in the file, if the logging level is 2 or greater. */
   log_metadata_nc(h5);
#endif

   return retval;
}

/**
 * @internal Close the HDF4 file.
 *
 * @param ncid File ID.
 *
 * @return ::NC_NOERR No error.
 * @return ::NC_EBADID Bad ncid.
 * @author Ed Hartnett
 */
int
NC_HDF4_close(int ncid)
{
   NC_GRP_INFO_T *grp;
   NC *nc;
   NC_FILE_INFO_T *h5;
   NC_HDF4_FILE_INFO_T *hdf4_file;
   int retval;

   LOG((1, "%s: ncid 0x%x", __func__, ncid));

   /* Find our metadata for this file. */
   if ((retval = nc4_find_nc_grp_h5(ncid, &nc, &grp, &h5)))
      return retval;
   assert(nc && h5 && grp && !grp->parent);

   /* Clean up HDF4 specific allocations. */
   if ((retval = hdf4_rec_grp_del(h5->root_grp)))
      return retval;

   /* Delete all the list contents for vars, dims, and atts, in each
    * group. */
   if ((retval = nc4_rec_grp_del(h5->root_grp)))
      return retval;

   /* Close hdf4 file and free HDF4 file info. */
   hdf4_file = (NC_HDF4_FILE_INFO_T *)h5->format_file_info;
   if (SDend(hdf4_file->sdid))
      return NC_EHDFERR;
   free(hdf4_file);

   /* Misc. Cleanup */
   nclistfree(h5->alldims);
   nclistfree(h5->allgroups);
   nclistfree(h5->alltypes);

   /* Free the nc4_info struct; above code should have reclaimed
      everything else */
   free(h5);

   return NC_NOERR;
}




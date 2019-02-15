/* Copyright 2018, UCAR/Unidata See netcdf/COPYRIGHT file for copying
 * and redistribution conditions.*/
/**
 * @file
 * @internal Dispatch code for HDF4. HDF4 access is read-only.
 *
 * @author Ed Hartnett
 */

#include "config.h"
#include <stdlib.h>
#include "hdf4dispatch.h"
#include "nc4dispatch.h"

/* This is the dispatch object that holds pointers to all the
 * functions that make up the HDF4 dispatch interface. */
static NC_Dispatch HDF4_dispatcher = {

NC_FORMATX_NC_HDF4,

NC_RO_create,
NC_HDF4_open,

NC_RO_redef,
NC_RO__enddef,
NC_RO_sync,
NC_HDF4_close,
NC_HDF4_close,
NC_RO_set_fill,
NC_NOTNC3_inq_base_pe,
NC_NOTNC3_set_base_pe,
NC_HDF4_inq_format,
NC_HDF4_inq_format_extended,

NC4_inq,
NC4_inq_type,

NC_RO_def_dim,
NC4_inq_dimid,
NC4_inq_dim,
NC4_inq_unlimdim,
NC_RO_rename_dim,

NC4_inq_att,
NC4_inq_attid,
NC4_inq_attname,
NC_RO_rename_att,
NC_RO_del_att,
NC4_get_att,
NC_RO_put_att,

NC_RO_def_var,
NC4_inq_varid,
NC_RO_rename_var,
NC_HDF4_get_vara,
NC_RO_put_vara,
NCDEFAULT_get_vars,
NCDEFAULT_put_vars,
NCDEFAULT_get_varm,
NCDEFAULT_put_varm,

NC4_inq_var_all,

NC_NOTNC4_var_par_access,
NC_RO_def_var_fill,

NC4_show_metadata,
NC4_inq_unlimdims,

NC4_inq_ncid,
NC4_inq_grps,
NC4_inq_grpname,
NC4_inq_grpname_full,
NC4_inq_grp_parent,
NC4_inq_grp_full_ncid,
NC4_inq_varids,
NC4_inq_dimids,
NC4_inq_typeids,
NC4_inq_type_equal,
NC_NOTNC4_def_grp,
NC_NOTNC4_rename_grp,
NC4_inq_user_type,
NC4_inq_typeid,

NC_NOTNC4_def_compound,
NC_NOTNC4_insert_compound,
NC_NOTNC4_insert_array_compound,
NC_NOTNC4_inq_compound_field,
NC_NOTNC4_inq_compound_fieldindex,
NC_NOTNC4_def_vlen,
NC_NOTNC4_put_vlen_element,
NC_NOTNC4_get_vlen_element,
NC_NOTNC4_def_enum,
NC_NOTNC4_insert_enum,
NC_NOTNC4_inq_enum_member,
NC_NOTNC4_inq_enum_ident,
NC_NOTNC4_def_opaque,
NC_NOTNC4_def_var_deflate,
NC_NOTNC4_def_var_fletcher32,
NC_NOTNC4_def_var_chunking,
NC_NOTNC4_def_var_endian,
NC_NOTNC4_def_var_filter,
NC_NOTNC4_set_var_chunk_cache,
NC_NOTNC4_get_var_chunk_cache
};

NC_Dispatch* HDF4_dispatch_table = NULL;

/**
 * @internal Initialize HDF4 dispatch layer.
 *
 * @return ::NC_NOERR No error.
 * @author Ed Hartnett
 */
int
NC_HDF4_initialize(void)
{
    HDF4_dispatch_table = &HDF4_dispatcher;
    return NC_NOERR;
}

/**
 * @internal Finalize HDF4 dispatch layer.
 *
 * @return ::NC_NOERR No error.
 * @author Ed Hartnett
 */
int
NC_HDF4_finalize(void)
{
    return NC_NOERR;
}

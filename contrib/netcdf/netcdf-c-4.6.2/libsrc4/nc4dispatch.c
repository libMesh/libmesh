/*********************************************************************
 *   Copyright 1993, UCAR/Unidata
 *   See netcdf/COPYRIGHT file for copying and redistribution conditions.
 *   $Header: /upc/share/CVS/netcdf-3/libsrc4/nc4dispatch.c,v 1.5 2010/05/27 02:19:37 dmh Exp $
 *********************************************************************/

#include "config.h"
#include <stdlib.h>
#include "nc4internal.h"
#include "nc4dispatch.h"
#include "nc.h"

/* If user-defined formats are in use, we need to declare their
 * dispatch tables. */
#ifdef USE_UDF0
extern NC_Dispatch UDF0_DISPATCH;
#endif /* USE_UDF0 */
#ifdef USE_UDF1
extern NC_Dispatch UDF1_DISPATCH;
#endif /* USE_UDF1 */

#ifdef USE_NETCDF4
/* Pointers to dispatch tables for user-defined formats. */
extern NC_Dispatch *UDF0_dispatch_table;
extern NC_Dispatch *UDF1_dispatch_table;
#endif /* USE_NETCDF4 */

static NC_Dispatch NC4_dispatcher = {

NC_FORMATX_NC4,

NC4_create,
NC4_open,

NC4_redef,
NC4__enddef,
NC4_sync,
NC4_abort,
NC4_close,
NC4_set_fill,
NC_NOTNC3_inq_base_pe,
NC_NOTNC3_set_base_pe,
NC4_inq_format,
NC4_inq_format_extended,

NC4_inq,
NC4_inq_type,

NC4_def_dim,
NC4_inq_dimid,
NC4_inq_dim,
NC4_inq_unlimdim,
NC4_rename_dim,

NC4_inq_att,
NC4_inq_attid,
NC4_inq_attname,
NC4_rename_att,
NC4_del_att,
NC4_get_att,
NC4_put_att,

NC4_def_var,
NC4_inq_varid,
NC4_rename_var,
NC4_get_vara,
NC4_put_vara,
NC4_get_vars,
NC4_put_vars,
NCDEFAULT_get_varm,
NCDEFAULT_put_varm,

NC4_inq_var_all,

NC4_var_par_access,
NC4_def_var_fill,

#ifdef USE_NETCDF4
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
NC4_def_grp,
NC4_rename_grp,
NC4_inq_user_type,
NC4_inq_typeid,

NC4_def_compound,
NC4_insert_compound,
NC4_insert_array_compound,
NC4_inq_compound_field,
NC4_inq_compound_fieldindex,
NC4_def_vlen,
NC4_put_vlen_element,
NC4_get_vlen_element,
NC4_def_enum,
NC4_insert_enum,
NC4_inq_enum_member,
NC4_inq_enum_ident,
NC4_def_opaque,
NC4_def_var_deflate,
NC4_def_var_fletcher32,
NC4_def_var_chunking,
NC4_def_var_endian,
NC4_def_var_filter,
NC4_set_var_chunk_cache,
NC4_get_var_chunk_cache,
#endif

};

NC_Dispatch* NC4_dispatch_table = NULL; /* moved here from ddispatch.c */

/**
 * @internal Initialize netCDF-4. If user-defined format(s) have been
 * specified in configure, load their dispatch table(s).
 *
 * @return ::NC_NOERR No error.
 * @author Dennis Heimbigner
 */
int
NC4_initialize(void)
{
   int ret = NC_NOERR;
   
   NC4_dispatch_table = &NC4_dispatcher;

   /* This needs some kind of conditional on if libhdf5 is enabled */
   if(!nc4_hdf5_initialized)
      nc4_hdf5_initialize();

#ifdef USE_UDF0
   /* If user-defined format 0 was specified during configure, set up
    * it's dispatch table. */
   if ((ret = nc_def_user_format(NC_UDF0, UDF0_DISPATCH_FUNC, NULL)))
      return ret;
#endif /* USE_UDF0 */
    
#ifdef USE_UDF1
   /* If user-defined format 0 was specified during configure, set up
    * it's dispatch table. */
   if ((ret = nc_def_user_format(NC_UDF1F, &UDF1_DISPATCH_FUNC, NULL)))
      return ret;
#endif /* USE_UDF0 */
    
#ifdef LOGGING
   if(getenv(NCLOGLEVELENV) != NULL) {
   char* slevel = getenv(NCLOGLEVELENV);
   long level = atol(slevel);
   if(level >= 0)
      nc_set_log_level((int)level);
}
#endif
   return ret;
}

/**
 * @internal Finalize netCDF-4.
 *
 * @return ::NC_NOERR No error.
 * @author Dennis Heimbigner
 */
int
NC4_finalize(void)
{
    nc4_hdf5_finalize();
    return NC_NOERR;
}

/*********************************************************************
 *   Copyright 1993, UCAR/Unidata
 *   See netcdf/COPYRIGHT file for copying and redistribution conditions.
 *********************************************************************/

/* WARNING: Order of mpi.h, nc.h, and pnetcdf.h is important */
#include "config.h"
#include <stdlib.h>
#include <mpi.h>
#include "nc.h"
#include "ncdispatch.h"
/* Must follow netcdf.h */
#include <pnetcdf.h>

typedef struct NCP_INFO
{
   /* pnetcdf_file will be true if the file is created/opened with the
    * parallel-netcdf library. pnetcdf_access_mode keeps track of
    * whether independpent or collective mode is
    * desired. pnetcdf_ndims keeps track of how many dims each var
    * has, which I need to know to convert start, count, and stride
    * arrays from size_t to MPI_Offset. (I can't use an inq function
    * to find out the number of dims, because these are collective in
    * pnetcdf.) */
   int pnetcdf_access_mode;
} NCP_INFO;

/* Define accessors for the dispatchdata */
#define NCP_DATA(nc) ((NCP_INFO*)(nc)->dispatchdata)
#define NCP_DATA_SET(nc,data) ((nc)->dispatchdata = (void*)(data))

/* Cannot have NC_MPIPOSIX flag, ignore NC_MPIIO as PnetCDF use MPIIO */
static const int LEGAL_CREATE_FLAGS = (NC_NOCLOBBER | NC_64BIT_OFFSET | NC_CLASSIC_MODEL | NC_SHARE | NC_LOCK | NC_64BIT_DATA | NC_MPIIO);

static const int LEGAL_OPEN_FLAGS = (NC_WRITE | NC_NOCLOBBER | NC_SHARE | NC_LOCK | NC_CLASSIC_MODEL | NC_64BIT_OFFSET | NC_64BIT_DATA | NC_MPIIO);


/**************************************************/

static int
NCP_create(const char *path, int cmode,
	  size_t initialsz, int basepe, size_t *chunksizehintp,
	  int use_parallel, void* mpidata,
	  struct NC_Dispatch* table, NC* nc)
{
    int res, default_format;
    NCP_INFO* nc5;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Info info = MPI_INFO_NULL;

    /* Check the cmode for only valid flags*/
    if(cmode & ~LEGAL_CREATE_FLAGS)
	return NC_EINVAL;

    /* Cannot have both MPIO flags */
    if((cmode & (NC_MPIIO|NC_MPIPOSIX)) == (NC_MPIIO|NC_MPIPOSIX))
	return NC_EINVAL;

    /* Cannot have both NC_64BIT_OFFSET & NC_64BIT_DATA */
    if((cmode & (NC_64BIT_OFFSET|NC_64BIT_DATA)) == (NC_64BIT_OFFSET|NC_64BIT_DATA))
 	return NC_EINVAL;

    default_format = nc_get_default_format();
    /* if (default_format == NC_FORMAT_CLASSIC) then we respect the format set in cmode */
    if (default_format == NC_FORMAT_64BIT_OFFSET) {
        if (! (cmode & NC_64BIT_OFFSET)) /* check if cmode has NC_64BIT_OFFSET already */
            cmode |= NC_64BIT_OFFSET;
    }
    else if (default_format == NC_FORMAT_CDF5) {
        if (! (cmode & NC_64BIT_DATA)) /* check if cmode has NC_64BIT_DATA already */
            cmode |= NC_64BIT_DATA;
    }

    /* No MPI environment initialized */
    if (mpidata == NULL) return NC_ENOPAR;

    comm = ((NC_MPI_INFO *)mpidata)->comm;
    info = ((NC_MPI_INFO *)mpidata)->info;

    /* Create our specific NCP_INFO instance */
    nc5 = (NCP_INFO*)calloc(1,sizeof(NCP_INFO));
    if(nc5 == NULL) return NC_ENOMEM;

    /* Link nc5 and nc */
    NCP_DATA_SET(nc,nc5);

    /* Fix up the cmode by keeping only essential flags;
       these are the flags that are the same in netcf.h and pnetcdf.h
    */
    /* It turns out that pnetcdf.h defines a flag called
       NC_64BIT_DATA (not to be confused with NC_64BIT_OFFSET).
       This flag is essential to getting ncmpi_create to create
       a proper pnetcdf format file.
       We have set the value of NC_64BIT_DATA to be the same as in pnetcdf.h
       (as of pnetcdf version 1.6.0) to avoid conflicts.
       In any case, this flag must be set.
    */
    /* PnetCDF recognizes the flags below for create and ignores NC_LOCK and  NC_SHARE */
    cmode &= (NC_WRITE | NC_NOCLOBBER | NC_SHARE | NC_64BIT_OFFSET | NC_64BIT_DATA);
    res = ncmpi_create(comm, path, cmode, info, &(nc->int_ncid));

    if(res && nc5 != NULL) free(nc5); /* reclaim allocated space */
    return res;
}

static int
NCP_open(const char *path, int cmode,
	    int basepe, size_t *chunksizehintp,
	    int use_parallel, void* mpidata,
	    struct NC_Dispatch* table, NC* nc)
{
    int res;
    NCP_INFO* nc5;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Info info = MPI_INFO_NULL;

    /* Check the cmode for only valid flags*/
    if(cmode & ~LEGAL_OPEN_FLAGS)
	return NC_EINVAL;

    /* Cannot have both MPIO flags */
    if((cmode & (NC_MPIIO|NC_MPIPOSIX)) == (NC_MPIIO|NC_MPIPOSIX))
	return NC_EINVAL;

    /* Appears that this comment is wrong; allow 64 bit offset*/
    /* Cannot have 64 bit offset flag */
    /* if(cmode & (NC_64BIT_OFFSET)) return NC_EINVAL; */

    if(mpidata != NULL) {
        comm = ((NC_MPI_INFO *)mpidata)->comm;
        info = ((NC_MPI_INFO *)mpidata)->info;
    } else {
	comm = MPI_COMM_WORLD;
	info = MPI_INFO_NULL;
    }

   /* PnetCDF recognizes the flags NC_WRITE and NC_NOCLOBBER for file open
     * and ignores NC_LOCK, NC_SHARE, NC_64BIT_OFFSET, and NC_64BIT_DATA.
     * Ignoring the NC_64BIT_OFFSET and NC_64BIT_DATA flags is because the
     * file is already in one of the CDF-formats, and setting these 2 flags
     * will not change the format of that file.
     */

    cmode &= (NC_WRITE | NC_NOCLOBBER);

    /* Create our specific NCP_INFO instance */
    nc5 = (NCP_INFO*)calloc(1,sizeof(NCP_INFO));
    if(nc5 == NULL) return NC_ENOMEM;

    /* Link nc5 and nc */
    NCP_DATA_SET(nc,nc5);

    res = ncmpi_open(comm, path, cmode, info, &(nc->int_ncid));

    /* Default to independent access, like netCDF-4/HDF5 files. */
    if(!res) {
	res = ncmpi_begin_indep_data(nc->int_ncid);
	nc5->pnetcdf_access_mode = NC_INDEPENDENT;
    }

    return res;
}

static int
NCP_redef(int ncid)
{
    NC* nc;
    int status = NC_check_id(ncid, &nc);
    if(status != NC_NOERR) return status;
    return ncmpi_redef(nc->int_ncid);
}

static int
NCP__enddef(int ncid, size_t h_minfree, size_t v_align, size_t v_minfree, size_t r_align)
{
    int status;
    NC* nc;
    NCP_INFO* nc5;
    MPI_Offset mpi_h_minfree = h_minfree;
    MPI_Offset mpi_v_align   = v_align;
    MPI_Offset mpi_v_minfree = v_minfree;
    MPI_Offset mpi_r_align   = r_align;

    status = NC_check_id(ncid, &nc);
    if(status != NC_NOERR)
	return status;

    nc5 = NCP_DATA(nc);
    assert(nc5);

    /* causes implicitly defined warning; may be because of old installed pnetcdf? */
#if 0
    /* In PnetCDF ncmpi__enddef() is only implemented in v1.5.0 and later */
    status = ncmpi__enddef(nc->int_ncid, mpi_h_minfree, mpi_v_align,
                           mpi_v_minfree, mpi_r_align);
#else
    status = ncmpi_enddef(nc->int_ncid);
#endif

    if(!status) {
	if (nc5->pnetcdf_access_mode == NC_INDEPENDENT)
	    status = ncmpi_begin_indep_data(nc->int_ncid);
    }
    return status;
}

static int
NCP_sync(int ncid)
{
    NC* nc;
    int status = NC_check_id(ncid, &nc);
    if(status != NC_NOERR) return status;
    return ncmpi_sync(nc->int_ncid);
}

static int
NCP_abort(int ncid)
{
    NC* nc;
    NCP_INFO* nc5;
    int status = NC_check_id(ncid, &nc);
    if(status != NC_NOERR) goto done;

    status = ncmpi_abort(nc->int_ncid);

done:
    nc5 = NCP_DATA(nc);
    if(nc5 != NULL) free(nc5); /* reclaim allocated space */
    return status;
}


static int
NCP_close(int ncid)
{
    NC* nc;
    NCP_INFO* nc5;
    int status = NC_check_id(ncid, &nc);
    if(status != NC_NOERR) goto done;

    status = ncmpi_close(nc->int_ncid);

done:
    nc5 = NCP_DATA(nc);
    if(nc5 != NULL) free(nc5); /* reclaim allocated space */
    return status;
}

static int
NCP_set_fill(int ncid, int fillmode, int *old_mode_ptr)
{
    NC* nc;
    int status = NC_check_id(ncid, &nc);
    if(status != NC_NOERR) return status;
    return ncmpi_set_fill(nc->int_ncid,fillmode,old_mode_ptr);
}

static int
NCP_inq_base_pe(int ncid, int* pep)
{
    if(pep) *pep = 0;
    return NC_NOERR;
}

static int
NCP_set_base_pe(int ncid, int pe)
{
    return NC_NOERR;
}

static int
NCP_inq_format(int ncid, int* formatp)
{
    NC* nc;
    int status = NC_check_id(ncid, &nc);
    if(status != NC_NOERR) return status;
    status = ncmpi_inq_format(nc->int_ncid,formatp);
    return status;
}

static int
NCP_inq_format_extended(int ncid, int* formatp, int *modep)
{
    NC* nc;
    int status = NC_check_id(ncid, &nc);
    if(status != NC_NOERR) return status;
    if(modep) *modep = nc->mode;
    /* Note that we do not use NC_FORMAT_CDF5 because PNETCDF has a dispatch table */
    if(formatp) *formatp = NC_FORMATX_PNETCDF;
    return NC_NOERR;
}

static int
NCP_inq(int ncid,
	int *ndimsp,
	int *nvarsp,
	int *nattsp,
	int *unlimp)
{
    NC* nc;
    int status = NC_check_id(ncid, &nc);
    if(status != NC_NOERR) return status;
    return ncmpi_inq(nc->int_ncid,ndimsp,nvarsp,nattsp,unlimp);
}


static int
NCP_inq_type(int ncid, nc_type typeid, char* name, size_t* size)
{
    /* Assert mode & NC_FORMAT_CDF5 */
    if (typeid < NC_BYTE || typeid >= NC_STRING)
           return NC_EBADTYPE;
    if(name)
	strcpy(name, NC_atomictypename(typeid));
    if(size)
        *size = NC_atomictypelen(typeid);
   return NC_NOERR;
}

static int
NCP_def_dim(int ncid, const char* name, size_t len, int* idp)
{
    int status;
    NCP_INFO* nc5;
    NC* nc;

    status = NC_check_id(ncid, &nc);
    if(status != NC_NOERR)
	return status;

    nc5 = NCP_DATA(nc);
    assert(nc5);

    return ncmpi_def_dim(nc->int_ncid, name, len, idp);
}

static int
NCP_inq_dimid(int ncid, const char *name, int *idp)
{
    NC* nc;
    int status = NC_check_id(ncid, &nc);
    if(status != NC_NOERR) return status;
    return ncmpi_inq_dimid(nc->int_ncid,name,idp);
}

static int
NCP_inq_dim(int ncid, int dimid, char *name, size_t* lenp)
{
    int status;
    NC* nc;
    MPI_Offset mpilen;
    status = NC_check_id(ncid, &nc);
    if(status != NC_NOERR) return status;
    status = ncmpi_inq_dim(nc->int_ncid,dimid,name,&mpilen);
    if(lenp) *lenp = mpilen;
    return status;
}

static int
NCP_inq_unlimdim(int ncid,  int *unlimdimidp)
{
    NC* nc;
    int status = NC_check_id(ncid, &nc);
    if(status != NC_NOERR) return status;
    return ncmpi_inq_unlimdim(nc->int_ncid,unlimdimidp);
}

static int
NCP_rename_dim(int ncid, int dimid, const char* newname)
{
    NC* nc;
    int status = NC_check_id(ncid, &nc);
    if(status != NC_NOERR) return status;
    return ncmpi_rename_dim(nc->int_ncid,dimid,newname);
}

static int
NCP_inq_att(int ncid, int varid, const char* name, nc_type* xtypep, size_t* lenp)
{
    NC* nc;
    MPI_Offset mpilen;
    int status = NC_check_id(ncid, &nc);
    if(status != NC_NOERR) return status;
    status = ncmpi_inq_att(nc->int_ncid,varid,name,xtypep,&mpilen);
    if(status != NC_NOERR) return status;
    if(lenp) *lenp = mpilen;
    return status;
}

static int
NCP_inq_attid(int ncid, int varid, const char *name, int *idp)
{
    NC* nc;
    int status = NC_check_id(ncid, &nc);
    if(status != NC_NOERR) return status;
    return ncmpi_inq_attid(nc->int_ncid,varid,name,idp);
}

static int
NCP_inq_attname(int ncid, int varid, int attnum, char *name)
{
    NC* nc;
    int status = NC_check_id(ncid, &nc);
    if(status != NC_NOERR) return status;
    return ncmpi_inq_attname(nc->int_ncid,varid,attnum,name);

}

static int
NCP_rename_att(int ncid, int varid, const char *name,
		const char *newname)
{
    NC* nc;
    int status = NC_check_id(ncid, &nc);
    if(status != NC_NOERR) return status;
    return ncmpi_rename_att(nc->int_ncid,varid,name,newname);
}

static int
NCP_del_att(int ncid, int varid, const char *name)
{
    NC* nc;
    int status = NC_check_id(ncid, &nc);
    if(status != NC_NOERR) return status;
    return ncmpi_del_att(nc->int_ncid,varid,name);
}

int
NCP_get_att(
	int ncid,
	int varid,
	const char *name,
	void *ip,
	nc_type memtype)
{
    NC* nc;
    int status;
    nc_type xtype;

    status = NC_check_id(ncid, &nc);
    if(status != NC_NOERR) return status;

    status = NCP_inq_att(ncid,varid,name,&xtype,NULL);

    if(memtype == NC_NAT) memtype = xtype;

    switch (memtype) {
    case NC_CHAR:
        return ncmpi_get_att_text(nc->int_ncid, varid, name, (char*)ip);
    case NC_BYTE:
        return ncmpi_get_att_schar(nc->int_ncid, varid, name, (signed char*)ip);
    case NC_SHORT:
        return ncmpi_get_att_short(nc->int_ncid, varid, name, (short*)ip);
    case NC_INT:
        return ncmpi_get_att_int(nc->int_ncid, varid, name, (int*)ip);
    case NC_FLOAT:
        return ncmpi_get_att_float(nc->int_ncid, varid, name, (float*)ip);
    case NC_DOUBLE:
        return ncmpi_get_att_double(nc->int_ncid, varid, name, (double*)ip);
    case NC_UBYTE:
         return ncmpi_get_att_uchar(nc->int_ncid, varid, name, (unsigned char*)ip);
    case NC_USHORT:
        return ncmpi_get_att_ushort(nc->int_ncid, varid, name, (unsigned short*)ip);
    case NC_UINT:
        return ncmpi_get_att_uint(nc->int_ncid, varid, name, (unsigned int*)ip);
    case NC_INT64:
        return ncmpi_get_att_longlong(nc->int_ncid, varid, name, (long long*)ip);
    case NC_UINT64:
        return ncmpi_get_att_ulonglong(nc->int_ncid, varid, name, (unsigned long long*)ip);
    default:
	break;
    }
    return NC_EBADTYPE;
}

int
NCP_put_att(
	int ncid,
	int varid,
	const char *name,
	nc_type xtype,
	size_t len,
	const void *ip,
	nc_type memtype)
{
    NC* nc;
    int status;
    MPI_Offset mpilen;

    if (!name || (strlen(name) > NC_MAX_NAME))
	return NC_EBADNAME;

    /* The length needs to be positive (cast needed for braindead
       systems with signed size_t). */
    if(((unsigned long) len) > X_INT_MAX)
	return NC_EINVAL;

    status = NC_check_id(ncid, &nc);
    if(status != NC_NOERR) return status;

    mpilen = len;

    switch (memtype) {
    case NC_CHAR:
        return ncmpi_put_att_text(nc->int_ncid, varid, name, mpilen, (char*)ip);
    case NC_BYTE:
        return ncmpi_put_att_schar(nc->int_ncid, varid, name, xtype, mpilen, (signed char*)ip);
    case NC_SHORT:
        return ncmpi_put_att_short(nc->int_ncid, varid, name, xtype, mpilen, (short*)ip);
    case NC_INT:
        return ncmpi_put_att_int(nc->int_ncid, varid, name, xtype, mpilen, (int*)ip);
    case NC_FLOAT:
        return ncmpi_put_att_float(nc->int_ncid, varid, name, xtype, mpilen, (float*)ip);
    case NC_DOUBLE:
        return ncmpi_put_att_double(nc->int_ncid, varid, name, xtype, mpilen, (double*)ip);
    case NC_UBYTE:
         return ncmpi_put_att_uchar(nc->int_ncid, varid, name, xtype, mpilen, (unsigned char*)ip);
    case NC_USHORT:
        return ncmpi_put_att_ushort(nc->int_ncid, varid, name, xtype, mpilen, (unsigned short*)ip);
    case NC_UINT:
        return ncmpi_put_att_uint(nc->int_ncid, varid, name, xtype, mpilen, (unsigned int*)ip);
    case NC_INT64:
        return ncmpi_put_att_longlong(nc->int_ncid, varid, name, xtype, mpilen, (long long*)ip);
    case NC_UINT64:
        return ncmpi_put_att_ulonglong(nc->int_ncid, varid, name, xtype, mpilen, (unsigned long long*)ip);
    default:
	break;
    }
    return NC_EBADTYPE;
}

static int
NCP_def_var(int ncid, const char *name, nc_type xtype,
	    int ndims, const int *dimidsp, int *varidp)
{
    NC* nc;
    NCP_INFO* nc5;
    int status;

    status = NC_check_id(ncid, &nc);
    if(status != NC_NOERR) return status;
    nc5 = NCP_DATA(nc);
    assert(nc5);

    status = ncmpi_def_var(nc->int_ncid,name,xtype,ndims,dimidsp,varidp);
    return status;
}

static int
NCP_inq_varid(int ncid, const char *name, int *varidp)
{
    NC* nc;
    int status = NC_check_id(ncid, &nc);
    if(status != NC_NOERR) return status;
    return ncmpi_inq_varid(nc->int_ncid,name,varidp);
}

static int
NCP_rename_var(int ncid, int varid, const char *name)
{
    NC* nc;
    int status = NC_check_id(ncid, &nc);
    if(status != NC_NOERR) return status;
    return ncmpi_rename_var(nc->int_ncid,varid,name);
}

static int
NCP_get_vara(int ncid,
		int varid,
		const size_t* startp,
		const size_t* countp,
		void* ip,
		nc_type memtype)
{
    NC* nc;
    NCP_INFO* nc5;
    int status;
    MPI_Offset mpi_start[NC_MAX_VAR_DIMS], mpi_count[NC_MAX_VAR_DIMS];
    int d;
    int rank = 0;

    status = NC_check_id(ncid, &nc);
    if(status != NC_NOERR) return status;

    nc5 = NCP_DATA(nc);
    assert(nc5);

    /* get variable's rank */
    status= ncmpi_inq_varndims(nc->int_ncid, varid, &rank);
    if(status) return status;

    /* We must convert the start and count arrays to MPI_Offset type. */
    for (d = 0; d < rank; d++) {
	 mpi_start[d] = startp[d];
	 mpi_count[d] = countp[d];
    }

    if (memtype == NC_NAT) {
        status = ncmpi_inq_vartype(nc->int_ncid, varid, &memtype);
        if (status) return status;
    }

    if(nc5->pnetcdf_access_mode == NC_INDEPENDENT) {
	switch(memtype) {
	case NC_BYTE:
	    status=ncmpi_get_vara_schar(nc->int_ncid, varid, mpi_start, mpi_count, ip); break;
	case NC_CHAR:
	    status=ncmpi_get_vara_text(nc->int_ncid, varid, mpi_start, mpi_count, ip); break;
	case NC_SHORT:
	    status=ncmpi_get_vara_short(nc->int_ncid, varid, mpi_start, mpi_count, ip); break;
	case NC_INT:
	    status=ncmpi_get_vara_int(nc->int_ncid, varid, mpi_start, mpi_count, ip); break;
	case NC_FLOAT:
	    status=ncmpi_get_vara_float(nc->int_ncid, varid, mpi_start, mpi_count, ip); break;
	case NC_DOUBLE:
	    status=ncmpi_get_vara_double(nc->int_ncid, varid, mpi_start, mpi_count, ip); break;
	case NC_UBYTE:
 	    status=ncmpi_get_vara_uchar(nc->int_ncid, varid, mpi_start, mpi_count, ip); break;
	case NC_USHORT:
	    status=ncmpi_get_vara_ushort(nc->int_ncid, varid, mpi_start, mpi_count, ip); break;
	case NC_UINT:
	    status=ncmpi_get_vara_uint(nc->int_ncid, varid, mpi_start, mpi_count, ip); break;
	case NC_INT64:
	    status=ncmpi_get_vara_longlong(nc->int_ncid, varid, mpi_start, mpi_count, ip); break;
	case NC_UINT64:
	    status=ncmpi_get_vara_ulonglong(nc->int_ncid, varid, mpi_start, mpi_count, ip); break;
	default:
	    status = NC_EBADTYPE;
	}
      } else {
	switch(memtype) {
	case NC_BYTE:
	    status=ncmpi_get_vara_schar_all(nc->int_ncid, varid, mpi_start, mpi_count, ip); break;
	case NC_CHAR:
	    status=ncmpi_get_vara_text_all(nc->int_ncid, varid, mpi_start, mpi_count, ip); break;
	case NC_SHORT:
	    status=ncmpi_get_vara_short_all(nc->int_ncid, varid, mpi_start, mpi_count, ip); break;
	case NC_INT:
	    status=ncmpi_get_vara_int_all(nc->int_ncid, varid, mpi_start, mpi_count, ip); break;
	case NC_FLOAT:
	    status=ncmpi_get_vara_float_all(nc->int_ncid, varid, mpi_start, mpi_count, ip); break;
	case NC_DOUBLE:
	    status=ncmpi_get_vara_double_all(nc->int_ncid, varid, mpi_start, mpi_count, ip); break;
	case NC_UBYTE:
	    status=ncmpi_get_vara_uchar_all(nc->int_ncid, varid, mpi_start, mpi_count, ip); break;
	case NC_USHORT:
	    status=ncmpi_get_vara_ushort_all(nc->int_ncid, varid, mpi_start, mpi_count, ip); break;
	case NC_UINT:
	    status=ncmpi_get_vara_uint_all(nc->int_ncid, varid, mpi_start, mpi_count, ip); break;
 	case NC_INT64:
	    status=ncmpi_get_vara_longlong_all(nc->int_ncid, varid, mpi_start, mpi_count, ip); break;
	case NC_UINT64:
	    status=ncmpi_get_vara_ulonglong_all(nc->int_ncid, varid, mpi_start, mpi_count, ip); break;
	default:
	    status = NC_EBADTYPE;
	}
    }
    return status;
}

static int
NCP_put_vara(int ncid,
	int varid,
	const size_t* startp,
	const size_t* countp,
	const void*ip,
	nc_type memtype)
{
    NC* nc;
    NCP_INFO* nc5;
    int status;
    MPI_Offset mpi_start[NC_MAX_VAR_DIMS], mpi_count[NC_MAX_VAR_DIMS];
    int d;
    int rank;

    status = NC_check_id(ncid, &nc);
    if(status != NC_NOERR) return status;

    nc5 = NCP_DATA(nc);
    assert(nc5);

    /* get variable's rank */
    status = ncmpi_inq_varndims(nc->int_ncid, varid, &rank);
    if(status) return status;

    /* We must convert the start and count arrays to MPI_Offset type. */
    for (d = 0; d < rank; d++) {
	 mpi_start[d] = startp[d];
	 mpi_count[d] = countp[d];
    }

    if (memtype == NC_NAT) {
        status = ncmpi_inq_vartype(nc->int_ncid, varid, &memtype);
        if (status) return status;
    }

    if(nc5->pnetcdf_access_mode == NC_INDEPENDENT) {
	switch(memtype) {
	case NC_BYTE:
	    status = ncmpi_put_vara_schar(nc->int_ncid, varid, mpi_start, mpi_count, ip); break;
	case NC_CHAR:
	    status = ncmpi_put_vara_text(nc->int_ncid, varid, mpi_start, mpi_count, ip); break;
	case NC_SHORT:
	    status = ncmpi_put_vara_short(nc->int_ncid, varid, mpi_start, mpi_count, ip); break;
	case NC_INT:
	    status = ncmpi_put_vara_int(nc->int_ncid, varid, mpi_start, mpi_count, ip); break;
	case NC_FLOAT:
	    status = ncmpi_put_vara_float(nc->int_ncid, varid, mpi_start, mpi_count, ip); break;
	case NC_DOUBLE:
	    status = ncmpi_put_vara_double(nc->int_ncid, varid, mpi_start, mpi_count, ip); break;
	case NC_UBYTE:
	    status = ncmpi_put_vara_uchar(nc->int_ncid, varid, mpi_start, mpi_count, ip); break;
	case NC_USHORT:
	    status = ncmpi_put_vara_ushort(nc->int_ncid, varid, mpi_start, mpi_count, ip); break;
	case NC_UINT:
	    status = ncmpi_put_vara_uint(nc->int_ncid, varid, mpi_start, mpi_count, ip); break;
	case NC_INT64:
      status = ncmpi_put_vara_longlong(nc->int_ncid, varid, mpi_start, mpi_count, ip); break;
	case NC_UINT64:
	    status = ncmpi_put_vara_ulonglong(nc->int_ncid, varid, mpi_start, mpi_count, ip); break;
	default:
	    status = NC_EBADTYPE;
	}
      } else {
	switch(memtype) {
	case NC_BYTE:
	    status = ncmpi_put_vara_schar_all(nc->int_ncid, varid, mpi_start, mpi_count, ip); break;
	case NC_CHAR:
	    status = ncmpi_put_vara_text_all(nc->int_ncid, varid, mpi_start, mpi_count, ip); break;
	case NC_SHORT:
	    status = ncmpi_put_vara_short_all(nc->int_ncid, varid, mpi_start, mpi_count, ip); break;
	case NC_INT:
	    status = ncmpi_put_vara_int_all(nc->int_ncid, varid, mpi_start, mpi_count, ip); break;
	case NC_FLOAT:
	    status = ncmpi_put_vara_float_all(nc->int_ncid, varid, mpi_start, mpi_count, ip); break;
	case NC_DOUBLE:
	    status = ncmpi_put_vara_double_all(nc->int_ncid, varid, mpi_start, mpi_count, ip); break;
	case NC_UBYTE:
	    status = ncmpi_put_vara_uchar_all(nc->int_ncid, varid, mpi_start, mpi_count, ip); break;
	case NC_USHORT:
	    status = ncmpi_put_vara_ushort_all(nc->int_ncid, varid, mpi_start, mpi_count, ip); break;
	case NC_UINT:
	    status = ncmpi_put_vara_uint_all(nc->int_ncid, varid, mpi_start, mpi_count, ip); break;
 	case NC_INT64:
	    status = ncmpi_put_vara_longlong_all(nc->int_ncid, varid, mpi_start, mpi_count, ip); break;
	case NC_UINT64:
	    status = ncmpi_put_vara_ulonglong_all(nc->int_ncid, varid, mpi_start, mpi_count, ip); break;
	default:
	    status = NC_EBADTYPE;
	}
    }
    return status;
}

static int
NCP_get_vars(int ncid,
		int varid,
		const size_t* startp,
		const size_t* countp,
		const ptrdiff_t* stridep,
		void* ip,
		nc_type memtype)
{
    NC* nc;
    NCP_INFO* nc5;
    int status;
    MPI_Offset mpi_start[NC_MAX_VAR_DIMS], mpi_count[NC_MAX_VAR_DIMS], mpi_stride[NC_MAX_VAR_DIMS];
    int d;
    int rank = 0;

    status = NC_check_id(ncid, &nc);
    if(status != NC_NOERR) return status;

    nc5 = NCP_DATA(nc);
    assert(nc5);

    /* get variable's rank */
    status= ncmpi_inq_varndims(nc->int_ncid, varid, &rank);
    if(status) return status;

    /* We must convert the start, count, and stride arrays to MPI_Offset type. */
    for (d = 0; d < rank; d++) {
	 mpi_start[d] = startp[d];
	 mpi_count[d] = countp[d];
	 mpi_stride[d] = stridep[d];
    }

    if (memtype == NC_NAT) {
        status = ncmpi_inq_vartype(nc->int_ncid, varid, &memtype);
        if (status) return status;
    }

    if(nc5->pnetcdf_access_mode == NC_INDEPENDENT) {
	switch(memtype) {
	case NC_BYTE:
	    status=ncmpi_get_vars_schar(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, ip); break;
	case NC_CHAR:
	    status=ncmpi_get_vars_text(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, ip); break;
	case NC_SHORT:
	    status=ncmpi_get_vars_short(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, ip); break;
	case NC_INT:
	    status=ncmpi_get_vars_int(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, ip); break;
	case NC_FLOAT:
	    status=ncmpi_get_vars_float(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, ip); break;
	case NC_DOUBLE:
	    status=ncmpi_get_vars_double(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, ip); break;
	case NC_UBYTE:
	    status=ncmpi_get_vars_uchar(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, ip); break;
	case NC_USHORT:
	    status=ncmpi_get_vars_ushort(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, ip); break;
	case NC_UINT:
	    status=ncmpi_get_vars_uint(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, ip); break;
	case NC_INT64:
	    status=ncmpi_get_vars_longlong(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, ip); break;
	case NC_UINT64:
	    status=ncmpi_get_vars_ulonglong(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, ip); break;
	default:
	    status = NC_EBADTYPE;
	}
      } else {
	switch(memtype) {
	case NC_BYTE:
	    status=ncmpi_get_vars_schar_all(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, ip); break;
	case NC_CHAR:
	    status=ncmpi_get_vars_text_all(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, ip); break;
	case NC_SHORT:
	    status=ncmpi_get_vars_short_all(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, ip); break;
	case NC_INT:
	    status=ncmpi_get_vars_int_all(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, ip); break;
	case NC_FLOAT:
	    status=ncmpi_get_vars_float_all(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, ip); break;
	case NC_DOUBLE:
	    status=ncmpi_get_vars_double_all(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, ip); break;
	case NC_UBYTE:
	    status=ncmpi_get_vars_uchar_all(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, ip); break;
	case NC_USHORT:
	    status=ncmpi_get_vars_ushort_all(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, ip); break;
	case NC_UINT:
	    status=ncmpi_get_vars_uint_all(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, ip); break;
	case NC_INT64:
	    status=ncmpi_get_vars_longlong_all(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, ip); break;
	case NC_UINT64:
	    status=ncmpi_get_vars_ulonglong_all(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, ip); break;
	default:
	    status = NC_EBADTYPE;
	}
    }
    return status;
}

static int
NCP_put_vars(int ncid,
	int varid,
	const size_t* startp,
	const size_t* countp,
	const ptrdiff_t* stridep,
	const void*ip,
	nc_type memtype)
{
    NC* nc;
    NCP_INFO* nc5;
    int status;
    MPI_Offset mpi_start[NC_MAX_VAR_DIMS], mpi_count[NC_MAX_VAR_DIMS], mpi_stride[NC_MAX_VAR_DIMS];
    int d;
    int rank;

    status = NC_check_id(ncid, &nc);
    if(status != NC_NOERR) return status;

    nc5 = NCP_DATA(nc);
    assert(nc5);

    /* get variable's rank */
    status = ncmpi_inq_varndims(nc->int_ncid, varid, &rank);
    if(status) return status;

    /* We must convert the start, count, and stride arrays to MPI_Offset type. */
    for (d = 0; d < rank; d++) {
	 mpi_start[d] = startp[d];
	 mpi_count[d] = countp[d];
	 mpi_stride[d] = stridep[d];
    }

    if (memtype == NC_NAT) {
        status = ncmpi_inq_vartype(nc->int_ncid, varid, &memtype);
        if (status) return status;
    }

    if(nc5->pnetcdf_access_mode == NC_INDEPENDENT) {
	switch(memtype) {
	case NC_BYTE:
	    status = ncmpi_put_vars_schar(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, ip); break;
	case NC_CHAR:
	    status = ncmpi_put_vars_text(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, ip); break;
	case NC_SHORT:
	    status = ncmpi_put_vars_short(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, ip); break;
	case NC_INT:
	    status = ncmpi_put_vars_int(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, ip); break;
	case NC_FLOAT:
	    status = ncmpi_put_vars_float(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, ip); break;
	case NC_DOUBLE:
	    status = ncmpi_put_vars_double(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, ip); break;
	case NC_UBYTE:
	    status = ncmpi_put_vars_uchar(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, ip); break;
	case NC_USHORT:
	    status = ncmpi_put_vars_ushort(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, ip); break;
	case NC_UINT:
	    status = ncmpi_put_vars_uint(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, ip); break;
	case NC_INT64:
	    status = ncmpi_put_vars_longlong(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, ip); break;
	case NC_UINT64:
	    status = ncmpi_put_vars_ulonglong(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, ip); break;
	default:
	    status = NC_EBADTYPE;
	}
      } else {
	switch(memtype) {
	case NC_BYTE:
	    status = ncmpi_put_vars_schar_all(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, ip); break;
	case NC_CHAR:
	    status = ncmpi_put_vars_text_all(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, ip); break;
	case NC_SHORT:
	    status = ncmpi_put_vars_short_all(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, ip); break;
	case NC_INT:
	    status = ncmpi_put_vars_int_all(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, ip); break;
	case NC_FLOAT:
	    status = ncmpi_put_vars_float_all(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, ip); break;
	case NC_DOUBLE:
	    status = ncmpi_put_vars_double_all(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, ip); break;
	case NC_UBYTE:
	    status = ncmpi_put_vars_uchar_all(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, ip); break;
	case NC_USHORT:
	    status = ncmpi_put_vars_ushort_all(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, ip); break;
	case NC_UINT:
	    status = ncmpi_put_vars_uint_all(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, ip); break;
	case NC_INT64:
	    status = ncmpi_put_vars_longlong_all(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, ip); break;
	case NC_UINT64:
	    status = ncmpi_put_vars_ulonglong_all(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, ip); break;
	default:
	    status = NC_EBADTYPE;
	}
    }
    return status;
}

static int
NCP_get_varm(int ncid,
		int varid,
		const size_t* startp,
		const size_t* countp,
		const ptrdiff_t* stridep,
		const ptrdiff_t* imapp,
		void* ip,
		nc_type memtype)
{
    NC* nc;
    NCP_INFO* nc5;
    int status;
    MPI_Offset mpi_start[NC_MAX_VAR_DIMS], mpi_count[NC_MAX_VAR_DIMS], mpi_stride[NC_MAX_VAR_DIMS], mpi_imap[NC_MAX_VAR_DIMS];
    int d;
    int rank = 0;

    status = NC_check_id(ncid, &nc);
    if(status != NC_NOERR) return status;

    nc5 = NCP_DATA(nc);
    assert(nc5);

    /* get variable's rank */
    status= ncmpi_inq_varndims(nc->int_ncid, varid, &rank);
    if(status) return status;

    /* We must convert the start, count, stride, and imap arrays to MPI_Offset type. */
    for (d = 0; d < rank; d++) {
	 mpi_start[d] = startp[d];
	 mpi_count[d] = countp[d];
	 mpi_stride[d] = stridep[d];
	 mpi_imap[d] = imapp[d];
    }

    if (memtype == NC_NAT) {
        status = ncmpi_inq_vartype(nc->int_ncid, varid, &memtype);
        if (status) return status;
    }

    if(nc5->pnetcdf_access_mode == NC_INDEPENDENT) {
	switch(memtype) {
	case NC_BYTE:
	    status=ncmpi_get_varm_schar(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, mpi_imap, ip); break;
	case NC_CHAR:
	    status=ncmpi_get_varm_text(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, mpi_imap, ip); break;
	case NC_SHORT:
	    status=ncmpi_get_varm_short(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, mpi_imap, ip); break;
	case NC_INT:
	    status=ncmpi_get_varm_int(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, mpi_imap, ip); break;
	case NC_FLOAT:
	    status=ncmpi_get_varm_float(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, mpi_imap, ip); break;
	case NC_DOUBLE:
	    status=ncmpi_get_varm_double(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, mpi_imap, ip); break;
	case NC_UBYTE:
	    status=ncmpi_get_varm_uchar(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, mpi_imap, ip); break;
	case NC_USHORT:
	    status=ncmpi_get_varm_ushort(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, mpi_imap, ip); break;
	case NC_UINT:
	    status=ncmpi_get_varm_uint(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, mpi_imap, ip); break;
	case NC_INT64:
	    status=ncmpi_get_varm_longlong(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, mpi_imap, ip); break;
	case NC_UINT64:
	    status=ncmpi_get_varm_ulonglong(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, mpi_imap, ip); break;
	default:
	    status = NC_EBADTYPE;
	}
      } else {
	switch(memtype) {
	case NC_BYTE:
	    status=ncmpi_get_varm_schar_all(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, mpi_imap, ip); break;
	case NC_CHAR:
	    status=ncmpi_get_varm_text_all(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, mpi_imap, ip); break;
	case NC_SHORT:
	    status=ncmpi_get_varm_short_all(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, mpi_imap, ip); break;
	case NC_INT:
	    status=ncmpi_get_varm_int_all(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, mpi_imap, ip); break;
	case NC_FLOAT:
	    status=ncmpi_get_varm_float_all(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, mpi_imap, ip); break;
	case NC_DOUBLE:
	    status=ncmpi_get_varm_double_all(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, mpi_imap, ip); break;
	case NC_UBYTE:
	    status=ncmpi_get_varm_uchar_all(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, mpi_imap, ip); break;
	case NC_USHORT:
	    status=ncmpi_get_varm_ushort_all(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, mpi_imap, ip); break;
	case NC_UINT:
	    status=ncmpi_get_varm_uint_all(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, mpi_imap, ip); break;
	case NC_INT64:
	    status=ncmpi_get_varm_longlong_all(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, mpi_imap, ip); break;
	case NC_UINT64:
	    status=ncmpi_get_varm_ulonglong_all(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, mpi_imap, ip); break;
	default:
	    status = NC_EBADTYPE;
	}
    }
    return status;
}

static int
NCP_put_varm(int ncid,
	int varid,
	const size_t* startp,
	const size_t* countp,
	const ptrdiff_t* stridep,
	const ptrdiff_t* imapp,
	const void*ip,
	nc_type memtype)
{
    NC* nc;
    NCP_INFO* nc5;
    int status;
    MPI_Offset mpi_start[NC_MAX_VAR_DIMS], mpi_count[NC_MAX_VAR_DIMS], mpi_stride[NC_MAX_VAR_DIMS], mpi_imap[NC_MAX_VAR_DIMS];
    int d;
    int rank;

    status = NC_check_id(ncid, &nc);
    if(status != NC_NOERR) return status;

    nc5 = NCP_DATA(nc);
    assert(nc5);

    /* get variable's rank */
    status = ncmpi_inq_varndims(nc->int_ncid, varid, &rank);
    if(status) return status;

    /* We must convert the start, count, stride, and imap arrays to MPI_Offset type. */
    for (d = 0; d < rank; d++) {
	 mpi_start[d] = startp[d];
	 mpi_count[d] = countp[d];
	 mpi_stride[d] = stridep[d];
	 mpi_imap[d] = imapp[d];
    }

    if (memtype == NC_NAT) {
        status = ncmpi_inq_vartype(nc->int_ncid, varid, &memtype);
        if (status) return status;
    }

    if(nc5->pnetcdf_access_mode == NC_INDEPENDENT) {
	switch(memtype) {
	case NC_BYTE:
	    status = ncmpi_put_varm_schar(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, mpi_imap, ip); break;
	case NC_CHAR:
	    status = ncmpi_put_varm_text(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, mpi_imap, ip); break;
	case NC_SHORT:
	    status = ncmpi_put_varm_short(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, mpi_imap, ip); break;
	case NC_INT:
	    status = ncmpi_put_varm_int(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, mpi_imap, ip); break;
	case NC_FLOAT:
	    status = ncmpi_put_varm_float(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, mpi_imap, ip); break;
	case NC_DOUBLE:
	    status = ncmpi_put_varm_double(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, mpi_imap, ip); break;
	case NC_UBYTE:
	    status = ncmpi_put_varm_uchar(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, mpi_imap, ip); break;
	case NC_USHORT:
	    status = ncmpi_put_varm_ushort(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, mpi_imap, ip); break;
	case NC_UINT:
	    status = ncmpi_put_varm_uint(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, mpi_imap, ip); break;
	case NC_INT64:
	    status = ncmpi_put_varm_longlong(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, mpi_imap, ip); break;
	case NC_UINT64:
	    status = ncmpi_put_varm_ulonglong(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, mpi_imap, ip); break;
	default:
	    status = NC_EBADTYPE;
	}
      } else {
	switch(memtype) {
	case NC_BYTE:
	    status = ncmpi_put_varm_schar_all(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, mpi_imap, ip); break;
	case NC_CHAR:
	    status = ncmpi_put_varm_text_all(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, mpi_imap, ip); break;
	case NC_SHORT:
	    status = ncmpi_put_varm_short_all(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, mpi_imap, ip); break;
	case NC_INT:
	    status = ncmpi_put_varm_int_all(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, mpi_imap, ip); break;
	case NC_FLOAT:
	    status = ncmpi_put_varm_float_all(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, mpi_imap, ip); break;
	case NC_DOUBLE:
	    status = ncmpi_put_varm_double_all(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, mpi_imap, ip); break;
	case NC_UBYTE:
	    status = ncmpi_put_varm_uchar_all(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, mpi_imap, ip); break;
	case NC_USHORT:
	    status = ncmpi_put_varm_ushort_all(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, mpi_imap, ip); break;
	case NC_UINT:
	    status = ncmpi_put_varm_uint_all(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, mpi_imap, ip); break;
	case NC_INT64:
	    status = ncmpi_put_varm_longlong_all(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, mpi_imap, ip); break;
	case NC_UINT64:
	    status = ncmpi_put_varm_ulonglong_all(nc->int_ncid, varid, mpi_start, mpi_count, mpi_stride, mpi_imap, ip); break;
	default:
	    status = NC_EBADTYPE;
	}
    }
    return status;
}

static int
NCP_inq_var_all(int ncid, int varid, char *name, nc_type *xtypep,
               int *ndimsp, int *dimidsp, int *nattsp,
               int *shufflep, int *deflatep, int *deflate_levelp,
               int *fletcher32p, int *contiguousp, size_t *chunksizesp,
               int *no_fill, void *fill_valuep, int *endiannessp,
	       int *options_maskp, int *pixels_per_blockp)
{
    int status;
    NC* nc;

    status = NC_check_id(ncid, &nc);
    if(status != NC_NOERR) return status;

    status = ncmpi_inq_var(nc->int_ncid, varid, name, xtypep, ndimsp, dimidsp, nattsp);
    if(status) return status;
    if(shufflep) *shufflep = 0;
    if(deflatep) *deflatep = 0;
    if(fletcher32p) *fletcher32p = 0;
    if(contiguousp) *contiguousp = NC_CONTIGUOUS;
    if(no_fill) *no_fill = 1;
    if(endiannessp) return NC_ENOTNC4;
    if(options_maskp) return NC_ENOTNC4;
    return NC_NOERR;
}

static int
NCP_var_par_access(int ncid, int varid, int par_access)
{
    NC *nc;
    NCP_INFO* nc5;
    int status;

    if (par_access != NC_INDEPENDENT && par_access != NC_COLLECTIVE)
	return NC_EINVAL;

    status = NC_check_id(ncid, &nc);
    if(status != NC_NOERR) return status;

    nc5 = NCP_DATA(nc);
    assert(nc5);

    if(par_access == nc5->pnetcdf_access_mode)
	return NC_NOERR;
    nc5->pnetcdf_access_mode = par_access;
    if (par_access == NC_INDEPENDENT)
	return ncmpi_begin_indep_data(nc->int_ncid);
    else
	return ncmpi_end_indep_data(nc->int_ncid);
}

#ifdef USE_NETCDF4

static int
NCP_show_metadata(int ncid)
{
    return NC_NOERR;
}

static int
NCP_inq_unlimdims(int ncid, int *ndimsp, int *unlimdimidsp)
{
    int retval;
    int unlimid;

    if((retval = NCP_inq_unlimdim(ncid, &unlimid)))
        return retval;
    if(unlimid != -1) {
        if(ndimsp) *ndimsp = 1;
        if(unlimdimidsp)
            unlimdimidsp[0] = unlimid;
    } else
        if(ndimsp) *ndimsp = 0;
    return NC_NOERR;
}

static int
NCP_inq_type_equal(int ncid1, nc_type typeid1, int ncid2, nc_type typeid2, int* equalp)
{
    /* Check input. */
    if(equalp == NULL) return NC_NOERR;

    if (typeid1 <= NC_NAT || typeid2 <= NC_NAT)
       return NC_EINVAL;

    *equalp = 0; /* assume */

    /* If one is atomic, and the other user-defined, the types are not equal */
    if ((typeid1 <= NC_STRING && typeid2 > NC_STRING) ||
        (typeid2 <= NC_STRING && typeid1 > NC_STRING)) {
        if (equalp) *equalp = 0;
        return NC_NOERR;
    }

    /* If both are atomic types, the answer is easy. */
    if (typeid1 <= ATOMICTYPEMAX5) {
        if (equalp) {
            if (typeid1 == typeid2)
                *equalp = 1;
            else
                *equalp = 0;
        }
        return NC_NOERR;
    }
    return NC_NOERR;
}

static int
NCP_def_grp(int parent_ncid, const char *name, int *new_ncid)
{
    return NC_ENOTNC4;
}

static int
NCP_rename_grp(int ncid, const char *name)
{
    return NC_ENOTNC4;
}

static int
NCP_inq_ncid(int ncid, const char *name, int *grp_ncid)
{
    if(grp_ncid) *grp_ncid = ncid;
    return NC_NOERR;
}

static int
NCP_inq_grps(int ncid, int *numgrps, int *ncids)
{
    if(numgrps)
       *numgrps = 0;
    return NC_NOERR;
}

static int
NCP_inq_grpname(int ncid, char *name)
{
    if(name)
        strcpy(name, "/");
    return NC_NOERR;
}

static int
NCP_inq_grpname_full(int ncid, size_t *lenp, char *full_name)
{
    if(full_name)
        strcpy(full_name, "/");
    if(lenp) *lenp = 1;
    return NC_NOERR;
}

static int
NCP_inq_grp_parent(int ncid, int *parent_ncid)
{
    return NC_ENOGRP;
}

static int
NCP_inq_grp_full_ncid(int ncid, const char *full_name, int *grp_ncid)
{
    return NC_ENOGRP;
}

static int
NCP_inq_varids(int ncid, int *nvarsp, int *varids)
{
    int retval,v,nvars;
    /* This is, effectively, a netcdf-3 file, there is only one group, the root
        group, and its vars have ids 0 thru nvars - 1. */
    if((retval = NCP_inq(ncid, NULL, &nvars, NULL, NULL)))
        return retval;
    if(nvarsp) *nvarsp = nvars;
    if(varids)
        for (v = 0; v < nvars; v++)
            varids[v] = v;
    return NC_NOERR;
}

static int
NCP_inq_dimids(int ncid, int *ndimsp, int *dimids, int include_parents)
{
    int retval,d,ndims;
    /* If this is like a netcdf-3 file, then the dimids are going to be 0
       thru ndims-1, so just provide them. */
    if((retval = NCP_inq(ncid, &ndims,  NULL, NULL, NULL)))
        return retval;
    if(ndimsp) *ndimsp = ndims;
    if(dimids)
        for (d = 0; d < ndims; d++)
            dimids[d] = d;
    return NC_NOERR;
}

static int
NCP_inq_typeid(int ncid, const char *name, nc_type *typeidp)
{
    int i;
    for (i = 0; i <= ATOMICTYPEMAX5; i++)
        if(!strcmp(name, NC_atomictypename(i))) {
            if(typeidp) *typeidp = i;
                return NC_NOERR;
        }
    return NC_ENOTNC4;
}

static int
NCP_inq_typeids(int ncid, int *ntypes, int *typeids)
{
    if(ntypes) *ntypes = 0;
    return NC_NOERR;
}

static int
NCP_inq_user_type(int ncid, nc_type typeid, char *name, size_t *size,
		 nc_type *base_nc_typep, size_t *nfieldsp, int *classp)
{
    return NC_ENOTNC4;
}

static int
NCP_def_compound(int ncid, size_t size, const char *name, nc_type *typeidp)
{
    return NC_ENOTNC4;
}

static int
NCP_insert_compound(int ncid, nc_type typeid, const char *name, size_t offset,
                    nc_type field_typeid)
{
    return NC_ENOTNC4;
}

static int
NCP_insert_array_compound(int ncid, nc_type typeid, const char *name,
			 size_t offset, nc_type field_typeid,
			 int ndims, const int *dim_sizes)
{
    return NC_ENOTNC4;
}


static int
NCP_inq_compound_field(int ncid, nc_type typeid, int fieldid, char *name,
		      size_t *offsetp, nc_type *field_typeidp, int *ndimsp,
		      int *dim_sizesp)
{
    return NC_ENOTNC4;
}

static int
NCP_inq_compound_fieldindex(int ncid, nc_type typeid, const char *name, int *fieldidp)
{
    return NC_ENOTNC4;
}

static int
NCP_def_opaque(int ncid, size_t datum_size, const char *name, nc_type* xtypep)
{
    return NC_ENOTNC4;
}

static int
NCP_def_vlen(int ncid, const char *name, nc_type base_typeid, nc_type* xtypep)
{
    return NC_ENOTNC4;
}

static int
NCP_def_enum(int ncid, nc_type base_typeid, const char *name,
	    nc_type *typeidp)
{
    return NC_ENOTNC4;
}

static int
NCP_inq_enum_ident(int ncid, nc_type xtype, long long value, char *identifier)
{
    return NC_ENOTNC4;
}

static int
NCP_inq_enum_member(int ncid, nc_type typeid, int idx, char *identifier,
		   void *value)
{
    return NC_ENOTNC4;
}

static int
NCP_insert_enum(int ncid, nc_type typeid, const char *identifier,
	       const void *value)
{
    return NC_ENOTNC4;
}

static int
NCP_put_vlen_element(int ncid, int typeid, void *vlen_element,
		    size_t len, const void *data)
{
    return NC_ENOTNC4;
}

static int
NCP_get_vlen_element(int ncid, int typeid, const void *vlen_element,
		    size_t *len, void *data)
{
    return NC_ENOTNC4;
}

static int
NCP_set_var_chunk_cache(int ncid, int varid, size_t size, size_t nelems, float preemption)
{
    return NC_ENOTNC4;
}

static int
NCP_get_var_chunk_cache(int ncid, int varid, size_t *sizep, size_t *nelemsp, float *preemptionp)
{
    return NC_ENOTNC4;
}

static int
NCP_def_var_deflate(int ncid, int varid, int shuffle, int deflate,
		   int deflate_level)
{
    return NC_ENOTNC4;
}

static int
NCP_def_var_fletcher32(int ncid, int varid, int fletcher32)
{
    return NC_ENOTNC4;
}

static int
NCP_def_var_chunking(int ncid, int varid, int contiguous, const size_t *chunksizesp)
{
    return NC_ENOTNC4;
}

static int
NCP_def_var_fill(int ncid, int varid, int no_fill, const void *fill_value)
{
    return NC_ENOTNC4;
}

static int
NCP_def_var_endian(int ncid, int varid, int endianness)
{
    return NC_ENOTNC4;
}

#endif /*USE_NETCDF4*/

/**************************************************/
/* Pnetcdf Dispatch table */

NC_Dispatch NCP_dispatcher = {

NC_FORMATX_PNETCDF,

NCP_create,
NCP_open,

NCP_redef,
NCP__enddef,
NCP_sync,
NCP_abort,
NCP_close,
NCP_set_fill,
NCP_inq_base_pe,
NCP_set_base_pe,
NCP_inq_format,
NCP_inq_format_extended,

NCP_inq,
NCP_inq_type,

NCP_def_dim,
NCP_inq_dimid,
NCP_inq_dim,
NCP_inq_unlimdim,
NCP_rename_dim,

NCP_inq_att,
NCP_inq_attid,
NCP_inq_attname,
NCP_rename_att,
NCP_del_att,
NCP_get_att,
NCP_put_att,

NCP_def_var,
NCP_inq_varid,
NCP_rename_var,
NCP_get_vara,
NCP_put_vara,
NCP_get_vars,
NCP_put_vars,
NCP_get_varm,
NCP_put_varm,

NCP_inq_var_all,

NCP_var_par_access,

#ifdef USE_NETCDF4
NCP_show_metadata,
NCP_inq_unlimdims,

NCP_inq_ncid,
NCP_inq_grps,
NCP_inq_grpname,
NCP_inq_grpname_full,
NCP_inq_grp_parent,
NCP_inq_grp_full_ncid,
NCP_inq_varids,
NCP_inq_dimids,
NCP_inq_typeids,
NCP_inq_type_equal,
NCP_def_grp,
NCP_rename_grp,
NCP_inq_user_type,
NCP_inq_typeid,

NCP_def_compound,
NCP_insert_compound,
NCP_insert_array_compound,
NCP_inq_compound_field,
NCP_inq_compound_fieldindex,
NCP_def_vlen,
NCP_put_vlen_element,
NCP_get_vlen_element,
NCP_def_enum,
NCP_insert_enum,
NCP_inq_enum_member,
NCP_inq_enum_ident,
NCP_def_opaque,
NCP_def_var_deflate,
NCP_def_var_fletcher32,
NCP_def_var_chunking,
NCP_def_var_fill,
NCP_def_var_endian,
NCP_set_var_chunk_cache,
NCP_get_var_chunk_cache,
#endif /*USE_NETCDF4*/

};

NC_Dispatch* NCP_dispatch_table = NULL; /* moved here from ddispatch.c */

int
NCP_initialize(void)
{
    NCP_dispatch_table = &NCP_dispatcher;
    return NC_NOERR;
}

int
NCP_finalize(void)
{
    return NC_NOERR;
}

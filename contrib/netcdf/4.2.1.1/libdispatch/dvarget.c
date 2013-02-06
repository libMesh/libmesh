/*! \file
Functions for getting data from variables.

Copyright 2011 University Corporation for Atmospheric
Research/Unidata. See \ref COPYRIGHT file for more info. */

#include "ncdispatch.h"


/** \internal
\ingroup variables 

 */
int
NC_get_vara(int ncid, int varid,
	    const size_t *start, const size_t *edges,
            void *value, nc_type memtype)
{
   NC* ncp;
   int stat = NC_check_id(ncid, &ncp);
   if(stat != NC_NOERR) return stat;
#ifdef USE_NETCDF4
   if(memtype >= NC_FIRSTUSERTYPEID) memtype = NC_NAT;
#endif
   if(edges == NULL) {
      size_t shape[NC_MAX_VAR_DIMS];
      int ndims;
      stat = nc_inq_varndims(ncid, varid, &ndims); 
      if(stat != NC_NOERR) return stat;
      stat = NC_getshape(ncid,varid,ndims,shape);
      if(stat != NC_NOERR) return stat;
      return ncp->dispatch->get_vara(ncid,varid,start,shape,value,memtype);
   } else
      return ncp->dispatch->get_vara(ncid,varid,start,edges,value,memtype);
}

/** \ingroup variables 
\internal
 */
static int
NC_get_var(int ncid, int varid, void *value, nc_type memtype)
{
   int ndims;
   size_t shape[NC_MAX_VAR_DIMS];
   int stat = nc_inq_varndims(ncid,varid, &ndims);
   if(stat) return stat;
   stat = NC_getshape(ncid,varid, ndims, shape);
   if(stat) return stat;
   return NC_get_vara(ncid, varid, NC_coord_zero, shape, value, memtype);
}

/** \internal
\ingroup variables 
 Most dispatch tables will use the default procedures
*/
int
NCDEFAULT_get_vars(int ncid, int varid, const size_t * start,
	    const size_t * edges, const ptrdiff_t * stride,
	    void *value, nc_type memtype)
{
   NC* ncp;
   int stat = NC_check_id(ncid, &ncp);

   if(stat != NC_NOERR) return stat;
   return ncp->dispatch->get_varm(ncid,varid,start,edges,stride,NULL,value,memtype);
}

/** \internal
\ingroup variables 
 */
static int
NC_get_var1(int ncid, int varid, const size_t *coord, void* value, 
	    nc_type memtype)
{
   return NC_get_vara(ncid, varid, coord, NC_coord_one, value, memtype);
}

/** \internal
\ingroup variables 
 */
int
NCDEFAULT_get_varm(int ncid, int varid, const size_t *start,
	    const size_t *edges, const ptrdiff_t *stride,
	    const ptrdiff_t *imapp, void *value0, nc_type memtype)
{
   int status = NC_NOERR;
   nc_type vartype = NC_NAT;
   int varndims,maxidim;
   NC* ncp;
   size_t memtypelen;
   ptrdiff_t cvtmap[NC_MAX_VAR_DIMS];
   char* value = (char*)value0;

   status = NC_check_id (ncid, &ncp);
   if(status != NC_NOERR) return status;

/*
  if(NC_indef(ncp)) return NC_EINDEFINE;
*/

   status = nc_inq_vartype(ncid, varid, &vartype); 
   if(status != NC_NOERR) return status;
   /* Check that this is an atomic type */
   if(vartype >= NC_MAX_ATOMIC_TYPE)
	return NC_EMAPTYPE;

   status = nc_inq_varndims(ncid, varid, &varndims); 
   if(status != NC_NOERR) return status;

   if(memtype == NC_NAT) {
      if(imapp != NULL && varndims != 0) {
	 /*
	  * convert map units from bytes to units of sizeof(type)
	  */
	 size_t ii;
	 const ptrdiff_t szof = (ptrdiff_t) nctypelen(vartype);
	 for(ii = 0; ii < varndims; ii++) {
	    if(imapp[ii] % szof != 0) {
	       /*free(cvtmap);*/
	       return NC_EINVAL;
	    }
	    cvtmap[ii] = imapp[ii] / szof;
	 }
	 imapp = cvtmap;
      }
      memtype = vartype;
   }

   if(memtype == NC_CHAR && vartype != NC_CHAR)
      return NC_ECHAR;
   else if(memtype != NC_CHAR && vartype == NC_CHAR)  
      return NC_ECHAR;

   memtypelen = nctypelen(memtype);

   maxidim = (int) varndims - 1;

   if (maxidim < 0)
   {
      /*
       * The variable is a scalar; consequently,
       * there s only one thing to get and only one place to put it.
       * (Why was I called?)
       */
      size_t edge1[1] = {1};
      return NC_get_vara(ncid, varid, start, edge1, value, memtype);
   }

   /*
    * else
    * The variable is an array.
    */
   {
      int idim;
      size_t *mystart = NULL;
      size_t *myedges;
      size_t *iocount;    /* count vector */
      size_t *stop;   /* stop indexes */
      size_t *length; /* edge lengths in bytes */
      ptrdiff_t *mystride;
      ptrdiff_t *mymap;
      size_t varshape[NC_MAX_VAR_DIMS];
      int isrecvar;
      size_t numrecs;

      /* Compute some dimension related values */
      isrecvar = NC_is_recvar(ncid,varid,&numrecs);
      NC_getshape(ncid,varid,varndims,varshape);	

      /*
       * Verify stride argument; also see if stride is all ones
       */
      if(stride != NULL) {
	 int stride1 = 1;
	 for (idim = 0; idim <= maxidim; ++idim)
	 {
            if (stride[idim] == 0
		/* cast needed for braindead systems with signed size_t */
                || ((unsigned long) stride[idim] >= X_INT_MAX))
            {
	       return NC_ESTRIDE;
            }
	    if(stride[idim] != 1) stride1 = 0;
	 }
         /* If stride1 is true, and there is no imap 
            then call get_vara directly.
         */
         if(stride1 && imapp == NULL) {
	     return NC_get_vara(ncid, varid, start, edges, value, memtype);
	 }
      }

      /* assert(sizeof(ptrdiff_t) >= sizeof(size_t)); */
      /* Allocate space for mystart,mystride,mymap etc.all at once */
      mystart = (size_t *)calloc(varndims * 7, sizeof(ptrdiff_t));
      if(mystart == NULL) return NC_ENOMEM;
      myedges = mystart + varndims;
      iocount = myedges + varndims;
      stop = iocount + varndims;
      length = stop + varndims;
      mystride = (ptrdiff_t *)(length + varndims);
      mymap = mystride + varndims;

      /*
       * Initialize I/O parameters.
       */
      for (idim = maxidim; idim >= 0; --idim)
      {
	 mystart[idim] = start != NULL
	    ? start[idim]
	    : 0;

	 if (edges != NULL && edges[idim] == 0)
	 {
	    status = NC_NOERR;    /* read/write no data */
	    goto done;
	 }

#ifdef COMPLEX
	 myedges[idim] = edges != NULL
	    ? edges[idim]
	    : idim == 0 && isrecvar
	    ? numrecs - mystart[idim]
	    : varshape[idim] - mystart[idim];
#else
	 if(edges != NULL)
	    myedges[idim] = edges[idim];
	 else if (idim == 0 && isrecvar)
	    myedges[idim] = numrecs - mystart[idim];
	 else
	    myedges[idim] = varshape[idim] - mystart[idim];
#endif

	 mystride[idim] = stride != NULL
	    ? stride[idim]
	    : 1;

	 /* Remember: imapp is byte oriented, not index oriented */
#ifdef COMPLEX
	 mymap[idim] = (imapp != NULL
			? imapp[idim]
			: (idim == maxidim ? 1
			   : mymap[idim + 1] * (ptrdiff_t) myedges[idim + 1]));
#else
	 if(imapp != NULL)
	    mymap[idim] = imapp[idim];
	 else if (idim == maxidim)
	    mymap[idim] = 1;
	 else
	    mymap[idim] = 
	       mymap[idim + 1] * (ptrdiff_t) myedges[idim + 1];
#endif
	 iocount[idim] = 1;
	 length[idim] = mymap[idim] * myedges[idim];
	 stop[idim] = mystart[idim] + myedges[idim] * mystride[idim];
      }

      /*
       * Check start, edges
       */
      for (idim = maxidim; idim >= 0; --idim)
      {
	 size_t dimlen = 
	    idim == 0 && isrecvar
	    ? numrecs
	    : varshape[idim];
	 if (mystart[idim] >= dimlen)
	 {
	    status = NC_EINVALCOORDS;
	    goto done;
	 }

	 if (mystart[idim] + myedges[idim] > dimlen)
	 {
	    status = NC_EEDGE;
	    goto done;
	 }

      }


      /* Lower body */
      /*
       * As an optimization, adjust I/O parameters when the fastest 
       * dimension has unity stride both externally and internally.
       * In this case, the user could have called a simpler routine
       * (i.e. ncvar$1()
       */
      if (mystride[maxidim] == 1
	  && mymap[maxidim] == 1)
      {
	 iocount[maxidim] = myedges[maxidim];
	 mystride[maxidim] = (ptrdiff_t) myedges[maxidim];
	 mymap[maxidim] = (ptrdiff_t) length[maxidim];
      }

      /* 
       * Perform I/O.  Exit when done.
       */
      for (;;)
      {
	 /* TODO: */
	 int lstatus = NC_get_vara(ncid, varid, mystart, iocount,
				   value, memtype);
	 if (lstatus != NC_NOERR) {
	    if(status == NC_NOERR || lstatus != NC_ERANGE)
	       status = lstatus;
	 }
	 /*
	  * The following code permutes through the variable s
	  * external start-index space and it s internal address
	  * space.  At the UPC, this algorithm is commonly
	  * called "odometer code".
	  */
	 idim = maxidim;
        carry:
	 value += (mymap[idim] * memtypelen);
	 mystart[idim] += mystride[idim];
	 if (mystart[idim] == stop[idim])
	 {
	    mystart[idim] = start[idim];
	    value -= (length[idim] * memtypelen);
	    if (--idim < 0)
	       break; /* normal return */
	    goto carry;
	 }
      } /* I/O loop */
     done:
      free(mystart);
   } /* variable is array */
   return status;
}

/** \ingroup variables 
\internal
Called by externally visible nc_get_vars_xxx routines */
static int
NC_get_vars(int ncid, int varid, const size_t *start, 
	    const size_t *edges, const ptrdiff_t *stride, void *value,
	    nc_type memtype)
{
   NC* ncp;
   int stat = NC_check_id(ncid, &ncp);

   if(stat != NC_NOERR) return stat;
#ifdef USE_NETCDF4
   if(memtype >= NC_FIRSTUSERTYPEID) memtype = NC_NAT;
#endif
   return ncp->dispatch->get_vars(ncid,varid,start,edges,stride,value,memtype);
}

/** \ingroup variables 
\internal
Called by externally visible nc_get_varm_xxx routines 
 */
static int
NC_get_varm(int ncid, int varid, const size_t *start, 
	    const size_t *edges, const ptrdiff_t *stride, const ptrdiff_t* map,
	    void *value, nc_type memtype)
{
   NC* ncp;
   int stat = NC_check_id(ncid, &ncp);

   if(stat != NC_NOERR) return stat;
#ifdef USE_NETCDF4
   if(memtype >= NC_FIRSTUSERTYPEID) memtype = NC_NAT;
#endif
   return ncp->dispatch->get_varm(ncid,varid,start,edges,stride,map,value,memtype);
}

/** \name Reading Data from Variables

Functions to read data from variables. */
/*! \{ */ /* All these functions are part of this named group... */

/** \ingroup variables
Read an array of values from a variable. 

The array to be read is specified by giving a corner and a vector of
edge lengths to \ref specify_hyperslab. 

The data values are read into consecutive locations with the last
dimension varying fastest. The netCDF dataset must be in data mode
(for netCDF-4/HDF5 files, the switch to data mode will happen
automatically, unless the classic model is used).

The nc_get_vara() function will read a variable of any type,
including user defined type. For this function, the type of the data
in memory must match the type of the variable - no data conversion is
done.

Other nc_get_vara_ functions will convert data to the desired output
type as needed.

\param ncid NetCDF or group ID, from a previous call to nc_open(),
nc_create(), nc_def_grp(), or associated inquiry functions such as 
nc_inq_ncid().

\param varid Variable ID

\param startp Start vector with one element for each dimension to \ref
specify_hyperslab.

\param countp Count vector with one element for each dimension to \ref
specify_hyperslab.

\param ip Pointer where the data will be copied. Memory must be
allocated by the user before this function is called.

\returns ::NC_NOERR No error.
\returns ::NC_ENOTVAR Variable not found.
\returns ::NC_EINVALCOORDS Index exceeds dimension bound.
\returns ::NC_EEDGE Start+count exceeds dimension bound.
\returns ::NC_ERANGE One or more of the values are out of range.
\returns ::NC_EINDEFINE Operation not allowed in define mode.
\returns ::NC_EBADID Bad ncid.

\section Example

Here is an example using nc_get_vara_double() to read all the values of
the variable named rh from an existing netCDF dataset named
foo.nc. For simplicity in this example, we assume that we know that rh
is dimensioned with time, lat, and lon, and that there are three time
values, five lat values, and ten lon values.

\code
     #include <netcdf.h>
        ...
     #define TIMES 3
     #define LATS 5
     #define LONS 10
     int  status;                      
     int ncid;                         
     int rh_id;                        
     static size_t start[] = {0, 0, 0};
     static size_t count[] = {TIMES, LATS, LONS};
     double rh_vals[TIMES*LATS*LONS]; 
        ...
     status = nc_open("foo.nc", NC_NOWRITE, &ncid);
     if (status != NC_NOERR) handle_error(status);
        ...
     status = nc_inq_varid (ncid, "rh", &rh_id);
     if (status != NC_NOERR) handle_error(status);
        ...
     status = nc_get_vara_double(ncid, rh_id, start, count, rh_vals);
     if (status != NC_NOERR) handle_error(status);
\endcode
 */
/**@{*/
int
nc_get_vara(int ncid, int varid, const size_t *startp, 
	    const size_t *countp, void *ip)
{
   NC* ncp = NULL;
   nc_type xtype = NC_NAT;
   int stat = NC_check_id(ncid, &ncp);
   if(stat != NC_NOERR) return stat;
   stat = nc_inq_vartype(ncid, varid, &xtype);
   if(stat != NC_NOERR) return stat;
   return NC_get_vara(ncid, varid, startp, countp, ip, xtype);
}

int
nc_get_vara_text(int ncid, int varid, const size_t *startp, 
		 const size_t *countp, char *ip)
{
   NC* ncp;
   int stat = NC_check_id(ncid, &ncp);
   if(stat != NC_NOERR) return stat;
   return NC_get_vara(ncid, varid, startp, countp, 
		      (void *)ip, NC_CHAR);
}

int
nc_get_vara_schar(int ncid, int varid, const size_t *startp, 
		  const size_t *countp, signed char *ip)
{
   NC* ncp;
   int stat = NC_check_id(ncid, &ncp);
   if(stat != NC_NOERR) return stat;
   return NC_get_vara(ncid, varid, startp, countp, 
		      (void *)ip, NC_BYTE);
}

int
nc_get_vara_uchar(int ncid, int varid, const size_t *startp, 
		  const size_t *countp, unsigned char *ip)
{
   NC* ncp;
   int stat = NC_check_id(ncid, &ncp);
   if(stat != NC_NOERR) return stat;
   return NC_get_vara(ncid, varid, startp, countp, 
		      (void *)ip, T_uchar);
}

int
nc_get_vara_short(int ncid, int varid, const size_t *startp, 
		  const size_t *countp, short *ip)
{
   NC* ncp;
   int stat = NC_check_id(ncid, &ncp);
   if(stat != NC_NOERR) return stat;
   return NC_get_vara(ncid, varid, startp, countp, 
		      (void *)ip, NC_SHORT);
}

int
nc_get_vara_int(int ncid, int varid,
		const size_t *startp, const size_t *countp, int *ip)
{
   NC* ncp;
   int stat = NC_check_id(ncid, &ncp);
   if(stat != NC_NOERR) return stat;
   return NC_get_vara(ncid,varid,startp,countp, (void *)ip,NC_INT);
}

int
nc_get_vara_long(int ncid, int varid,
		 const size_t *startp, const size_t *countp, long *ip)
{
   NC* ncp;
   int stat = NC_check_id(ncid, &ncp);
   if(stat != NC_NOERR) return stat;
   return NC_get_vara(ncid,varid,startp,countp, (void *)ip,T_long);
}

int
nc_get_vara_float(int ncid, int varid,
		  const size_t *startp, const size_t *countp, float *ip)
{
   NC* ncp;
   int stat = NC_check_id(ncid, &ncp);
   if(stat != NC_NOERR) return stat;
   return NC_get_vara(ncid,varid,startp,countp, (void *)ip,T_float);
}


int
nc_get_vara_double(int ncid, int varid, const size_t *startp, 
		   const size_t *countp, double *ip)
{
   NC* ncp;
   int stat = NC_check_id(ncid, &ncp);
   if(stat != NC_NOERR) return stat;
   return NC_get_vara(ncid,varid,startp,countp, (void *)ip,T_double);
}

int
nc_get_vara_ubyte(int ncid, int varid,
		  const size_t *startp, const size_t *countp, unsigned char *ip)
{
   NC* ncp;
   int stat = NC_check_id(ncid, &ncp);
   if(stat != NC_NOERR) return stat;
   return NC_get_vara(ncid,varid,startp,countp, (void *)ip,T_ubyte);
}

int
nc_get_vara_ushort(int ncid, int varid,
		   const size_t *startp, const size_t *countp, unsigned short *ip)
{
   NC* ncp;
   int stat = NC_check_id(ncid, &ncp);
   if(stat != NC_NOERR) return stat;
   return NC_get_vara(ncid,varid,startp,countp, (void *)ip,T_ushort);
}

int
nc_get_vara_uint(int ncid, int varid,
		 const size_t *startp, const size_t *countp, unsigned int *ip)
{
   NC* ncp;
   int stat = NC_check_id(ncid, &ncp);
   if(stat != NC_NOERR) return stat;
   return NC_get_vara(ncid,varid,startp,countp, (void *)ip,T_uint);
}

int
nc_get_vara_longlong(int ncid, int varid,
		     const size_t *startp, const size_t *countp, long long *ip)
{
   NC* ncp;
   int stat = NC_check_id(ncid, &ncp);
   if(stat != NC_NOERR) return stat;
   return NC_get_vara(ncid,varid,startp,countp, (void *)ip,T_longlong);
}

int
nc_get_vara_ulonglong(int ncid, int varid,
		      const size_t *startp, const size_t *countp, unsigned long long *ip)
{
   NC* ncp;
   int stat = NC_check_id(ncid, &ncp);
   if(stat != NC_NOERR) return stat;
   return NC_get_vara(ncid,varid,startp,countp, (void *)ip,NC_UINT64);
}

#ifdef USE_NETCDF4
int
nc_get_vara_string(int ncid, int varid,
		   const size_t *startp, const size_t *countp, char* *ip)
{
   NC* ncp;
   int stat = NC_check_id(ncid, &ncp);
   if(stat != NC_NOERR) return stat;
   return NC_get_vara(ncid,varid,startp,countp, (void *)ip,NC_STRING);
}

#endif /*USE_NETCDF4*/
/**@}*/

/** \ingroup variables
Read a single datum from a variable. 

Inputs are the netCDF ID, the variable ID, a multidimensional index
that specifies which value to get, and the address of a location into
which the data value will be read. The value is converted from the
external data type of the variable, if necessary.

The nc_get_var1() function will read a variable of any type, including
user defined type. For this function, the type of the data in memory
must match the type of the variable - no data conversion is done.

Other nc_get_var1_ functions will convert data to the desired output
type as needed.

\param ncid NetCDF or group ID, from a previous call to nc_open(),
nc_create(), nc_def_grp(), or associated inquiry functions such as 
nc_inq_ncid().

\param varid Variable ID

\param indexp Index vector with one element for each dimension.

\param ip Pointer where the data will be copied. Memory must be
allocated by the user before this function is called.

\returns ::NC_NOERR No error.
\returns ::NC_ENOTVAR Variable not found.
\returns ::NC_EINVALCOORDS Index exceeds dimension bound.
\returns ::NC_ERANGE One or more of the values are out of range.
\returns ::NC_EINDEFINE Operation not allowed in define mode.
\returns ::NC_EBADID Bad ncid.
*/
/** \{ */
int
nc_get_var1(int ncid, int varid, const size_t *indexp, void *ip)
{
   return NC_get_var1(ncid, varid, indexp, ip, NC_NAT);
}

int
nc_get_var1_text(int ncid, int varid, const size_t *indexp, char *ip)
{
   NC* ncp;
   int stat = NC_check_id(ncid, &ncp);
   if(stat != NC_NOERR) return stat;
   return NC_get_var1(ncid, varid, indexp, (void *)ip, NC_CHAR);
}

int
nc_get_var1_schar(int ncid, int varid, const size_t *indexp, signed char *ip)
{
   NC* ncp;
   int stat = NC_check_id(ncid, &ncp);
   if(stat != NC_NOERR) return stat;
   return NC_get_var1(ncid, varid, indexp, (void *)ip, NC_BYTE);
}

int
nc_get_var1_uchar(int ncid, int varid, const size_t *indexp, unsigned char *ip)
{
   NC* ncp;
   int stat = NC_check_id(ncid, &ncp);
   if(stat != NC_NOERR) return stat;
   return NC_get_var1(ncid, varid, indexp, (void *)ip, NC_UBYTE);
}

int
nc_get_var1_short(int ncid, int varid, const size_t *indexp, short *ip)
{
   NC* ncp;
   int stat = NC_check_id(ncid, &ncp);
   if(stat != NC_NOERR) return stat;
   return NC_get_var1(ncid, varid, indexp, (void *)ip, NC_SHORT);
}

int
nc_get_var1_int(int ncid, int varid, const size_t *indexp, int *ip)
{
   NC* ncp;
   int stat = NC_check_id(ncid, &ncp);
   if(stat != NC_NOERR) return stat;
   return NC_get_var1(ncid, varid, indexp, (void *)ip, NC_INT);
}

int
nc_get_var1_long(int ncid, int varid, const size_t *indexp, 
		 long *ip)
{
   NC* ncp;
   int stat = NC_check_id(ncid, &ncp);
   if(stat != NC_NOERR) return stat;
   return NC_get_var1(ncid, varid, indexp, (void *)ip, longtype);
}

int
nc_get_var1_float(int ncid, int varid, const size_t *indexp, 
		  float *ip)
{
   NC* ncp;
   int stat = NC_check_id(ncid, &ncp);
   if(stat != NC_NOERR) return stat;
   return NC_get_var1(ncid, varid, indexp, (void *)ip, NC_FLOAT);
}

int
nc_get_var1_double(int ncid, int varid, const size_t *indexp, 
		   double *ip)
{
   NC* ncp;
   int stat = NC_check_id(ncid, &ncp);
   if(stat != NC_NOERR) return stat;
   return NC_get_var1(ncid, varid, indexp, (void *)ip, NC_DOUBLE);
}

int
nc_get_var1_ubyte(int ncid, int varid, const size_t *indexp, 
		  unsigned char *ip)
{
   NC* ncp;
   int stat = NC_check_id(ncid, &ncp);
   if(stat != NC_NOERR) return stat;
   return NC_get_var1(ncid, varid, indexp, (void *)ip, NC_UBYTE);
}

int
nc_get_var1_ushort(int ncid, int varid, const size_t *indexp, 
		   unsigned short *ip)
{
   NC* ncp;
   int stat = NC_check_id(ncid, &ncp);
   if(stat != NC_NOERR) return stat;
   return NC_get_var1(ncid, varid, indexp, (void *)ip, NC_USHORT);
}

int
nc_get_var1_uint(int ncid, int varid, const size_t *indexp, 
		 unsigned int *ip)
{
   NC* ncp;
   int stat = NC_check_id(ncid, &ncp);
   if(stat != NC_NOERR) return stat;
   return NC_get_var1(ncid, varid, indexp, (void *)ip, NC_INT);
}

int
nc_get_var1_longlong(int ncid, int varid, const size_t *indexp, 
		     long long *ip)
{
   NC* ncp;
   int stat = NC_check_id(ncid, &ncp);
   if(stat != NC_NOERR) return stat;
   return NC_get_var1(ncid, varid, indexp, (void *)ip, NC_INT64);
}

int
nc_get_var1_ulonglong(int ncid, int varid, const size_t *indexp, 
		      unsigned long long *ip)
{
   NC* ncp;
   int stat = NC_check_id(ncid, &ncp);
   if(stat != NC_NOERR) return stat;
   return NC_get_var1(ncid, varid, indexp, (void *)ip, NC_UINT64);
}

#ifdef USE_NETCDF4
int
nc_get_var1_string(int ncid, int varid, const size_t *indexp, char* *ip)
{
   NC* ncp;
   int stat = NC_check_id(ncid, &ncp);
   if(stat != NC_NOERR) return stat;
   return NC_get_var1(ncid, varid, indexp, (void *)ip, NC_STRING);
}
#endif /*USE_NETCDF4*/
/** \} */

/** \ingroup variables
Read an entire variable in one call. 

This function will read all the values from a netCDF variable of an
open netCDF dataset.

This is the simplest interface to use for reading the value of a
scalar variable or when all the values of a multidimensional variable
can be read at once. The values are read into consecutive locations
with the last dimension varying fastest. The netCDF dataset must be in
data mode.

Take care when using this function with record variables (variables
that use the ::NC_UNLIMITED dimension). If you try to read all the
values of a record variable into an array but there are more records
in the file than you assume, more data will be read than you expect,
which may cause a segmentation violation. To avoid such problems, it
is better to use the nc_get_vara interfaces for variables that use the
::NC_UNLIMITED dimension.

The functions for types ubyte, ushort, uint, longlong, ulonglong, and
string are only available for netCDF-4/HDF5 files.

The nc_get_var() function will read a variable of any type, including
user defined type. For this function, the type of the data in memory
must match the type of the variable - no data conversion is done.

\param ncid NetCDF or group ID, from a previous call to nc_open(),
nc_create(), nc_def_grp(), or associated inquiry functions such as 
nc_inq_ncid().

\param varid Variable ID

\param ip Pointer where the data will be copied. Memory must be
allocated by the user before this function is called.

\returns ::NC_NOERR No error.
\returns ::NC_ENOTVAR Variable not found.
\returns ::NC_ERANGE One or more of the values are out of range.
\returns ::NC_EINDEFINE Operation not allowed in define mode.
\returns ::NC_EBADID Bad ncid.
*/
/** \{ */
int
nc_get_var(int ncid, int varid, void *ip)
{
   return NC_get_var(ncid, varid, ip, NC_NAT);
}

int
nc_get_var_text(int ncid, int varid, char *ip)
{
   NC *ncp;
   int stat = NC_check_id(ncid, &ncp);
   if(stat != NC_NOERR) return stat;
   return NC_get_var(ncid, varid, (void *)ip, NC_CHAR);
}

int
nc_get_var_schar(int ncid, int varid, signed char *ip)
{
   NC *ncp;
   int stat = NC_check_id(ncid, &ncp);
   if(stat != NC_NOERR) return stat;
   return NC_get_var(ncid, varid, (void *)ip, NC_BYTE);
}

int
nc_get_var_uchar(int ncid, int varid, unsigned char *ip)
{
   NC *ncp;
   int stat = NC_check_id(ncid, &ncp);
   if(stat != NC_NOERR) return stat;
   return NC_get_var(ncid,varid, (void *)ip, NC_UBYTE);
}

int
nc_get_var_short(int ncid, int varid, short *ip)
{
   NC* ncp;
   int stat = NC_check_id(ncid, &ncp);
   if(stat != NC_NOERR) return stat;
   return NC_get_var(ncid, varid, (void *)ip, NC_SHORT);
}

int
nc_get_var_int(int ncid, int varid, int *ip)
{
   NC* ncp;
   int stat = NC_check_id(ncid, &ncp);
   if(stat != NC_NOERR) return stat;
   return NC_get_var(ncid,varid, (void *)ip, NC_INT);
}

int
nc_get_var_long(int ncid, int varid, long *ip)
{
   NC* ncp;
   int stat = NC_check_id(ncid, &ncp);
   if(stat != NC_NOERR) return stat;
   return NC_get_var(ncid,varid, (void *)ip, longtype);
}

int
nc_get_var_float(int ncid, int varid, float *ip)
{
   NC* ncp;
   int stat = NC_check_id(ncid, &ncp);
   if(stat != NC_NOERR) return stat;
   return NC_get_var(ncid,varid, (void *)ip, NC_FLOAT);
}

int
nc_get_var_double(int ncid, int varid, double *ip)
{
   NC* ncp;
   int stat = NC_check_id(ncid, &ncp);
   if(stat != NC_NOERR) return stat;
   return NC_get_var(ncid,varid, (void *)ip, NC_DOUBLE);
}

int
nc_get_var_ubyte(int ncid, int varid, unsigned char *ip)
{
   NC* ncp;
   int stat = NC_check_id(ncid, &ncp);
   if(stat != NC_NOERR) return stat;
   return NC_get_var(ncid,varid, (void *)ip, NC_UBYTE);
}

int
nc_get_var_ushort(int ncid, int varid, unsigned short *ip)
{
   NC* ncp;
   int stat = NC_check_id(ncid, &ncp);
   if(stat != NC_NOERR) return stat;
   return NC_get_var(ncid,varid, (void *)ip, NC_USHORT);
}

int
nc_get_var_uint(int ncid, int varid, unsigned int *ip)
{
   NC* ncp;
   int stat = NC_check_id(ncid, &ncp);
   if(stat != NC_NOERR) return stat;
   return NC_get_var(ncid,varid, (void *)ip, NC_UINT);
}

int
nc_get_var_longlong(int ncid, int varid, long long *ip)
{
   NC* ncp;
   int stat = NC_check_id(ncid, &ncp);
   if(stat != NC_NOERR) return stat;
   return NC_get_var(ncid,varid, (void *)ip, NC_INT64);
}

int
nc_get_var_ulonglong(int ncid, int varid, unsigned long long *ip)
{
   NC* ncp;
   int stat = NC_check_id(ncid, &ncp);
   if(stat != NC_NOERR) return stat;
   return NC_get_var(ncid,varid, (void *)ip,NC_UINT64);
}

#ifdef USE_NETCDF4
int
nc_get_var_string(int ncid, int varid, char* *ip)
{
   NC* ncp;
   int stat = NC_check_id(ncid, &ncp);
   if(stat != NC_NOERR) return stat;
   return NC_get_var(ncid,varid, (void *)ip,NC_STRING);
}
#endif /*USE_NETCDF4*/
/** \} */

/** \ingroup variables
Read a strided array from a variable. 

This function reads a subsampled (strided) array section of values
from a netCDF variable of an open netCDF dataset. The subsampled array
section is specified by giving a corner, a vector of edge lengths, and
a stride vector. The values are read with the last dimension of the
netCDF variable varying fastest. The netCDF dataset must be in data
mode.

The nc_get_vars() function will read a variable of any type, including
user defined type. For this function, the type of the data in memory
must match the type of the variable - no data conversion is done.

\param ncid NetCDF or group ID, from a previous call to nc_open(),
nc_create(), nc_def_grp(), or associated inquiry functions such as 
nc_inq_ncid().

\param varid Variable ID

\param startp Start vector with one element for each dimension to \ref
specify_hyperslab.

\param countp Count vector with one element for each dimension to \ref
specify_hyperslab.

\param stridep Stride vector with one element for each dimension to
\ref specify_hyperslab.

\param ip Pointer where the data will be copied. Memory must be
allocated by the user before this function is called.

\returns ::NC_NOERR No error.
\returns ::NC_ENOTVAR Variable not found.
\returns ::NC_EINVALCOORDS Index exceeds dimension bound.
\returns ::NC_ERANGE One or more of the values are out of range.
\returns ::NC_EINDEFINE Operation not allowed in define mode.
\returns ::NC_EBADID Bad ncid.
*/
/** \{ */
int
nc_get_vars (int ncid, int varid, const size_t * startp,
	     const size_t * countp, const ptrdiff_t * stridep,
	     void *ip)
{
   NC* ncp;
   int stat = NC_NOERR;

   if ((stat = NC_check_id(ncid, &ncp)))
       return stat;
   return ncp->dispatch->get_vars(ncid, varid, startp, countp, stridep,
		      ip, NC_NAT);
}

int
nc_get_vars_text(int ncid, int varid, const size_t *startp, 
		 const size_t *countp, const ptrdiff_t * stridep,
		 char *ip)
{
   NC* ncp;
   int stat = NC_check_id(ncid, &ncp);
   if(stat != NC_NOERR) return stat;
   return NC_get_vars(ncid,varid,startp, countp, stridep,
		      (void *)ip, NC_CHAR);
}

int
nc_get_vars_schar(int ncid, int varid, const size_t *startp, 
		  const size_t *countp, const ptrdiff_t * stridep,
		  signed char *ip)
{
   NC* ncp;
   int stat = NC_check_id(ncid, &ncp);
   if(stat != NC_NOERR) return stat;
   return NC_get_vars(ncid,varid,startp, countp, stridep,
		      (void *)ip, NC_BYTE);
}

int
nc_get_vars_uchar(int ncid, int varid, const size_t *startp, 
		  const size_t *countp, const ptrdiff_t * stridep,
		  unsigned char *ip)
{
   NC* ncp;
   int stat = NC_check_id(ncid, &ncp);
   if(stat != NC_NOERR) return stat;
   return NC_get_vars(ncid,varid,startp, countp, stridep,
		      (void *)ip, T_uchar);
}

int
nc_get_vars_short(int ncid, int varid, const size_t *startp, 
		  const size_t *countp, const ptrdiff_t *stridep,
		  short *ip)
{
   NC* ncp;
   int stat = NC_check_id(ncid, &ncp);
   if(stat != NC_NOERR) return stat;
   return NC_get_vars(ncid,varid,startp, countp, stridep, 
		      (void *)ip, NC_SHORT);
}

int
nc_get_vars_int(int ncid, int varid, const size_t *startp, 
		const size_t *countp, const ptrdiff_t * stridep,
		int *ip)
{
   NC* ncp;
   int stat = NC_check_id(ncid, &ncp);
   if(stat != NC_NOERR) return stat;
   return NC_get_vars(ncid,varid,startp, countp, stridep,
		      (void *)ip, NC_INT);
}

int
nc_get_vars_long(int ncid, int varid, const size_t *startp, 
		 const size_t *countp, const ptrdiff_t * stridep,
		 long *ip)
{
   NC* ncp;
   int stat = NC_check_id(ncid, &ncp);
   if(stat != NC_NOERR) return stat;
   return NC_get_vars(ncid,varid,startp, countp, stridep,
		      (void *)ip, T_long);
}

int
nc_get_vars_float(int ncid, int varid, const size_t *startp, 
		  const size_t *countp, const ptrdiff_t * stridep,
		  float *ip)
{
   NC* ncp;
   int stat = NC_check_id(ncid, &ncp);
   if(stat != NC_NOERR) return stat;
   return NC_get_vars(ncid,varid,startp, countp, stridep,
		      (void *)ip, T_float);
}

int
nc_get_vars_double(int ncid, int varid, const size_t *startp, 
		   const size_t *countp, const ptrdiff_t * stridep,
		   double *ip)
{
   NC* ncp;
   int stat = NC_check_id(ncid, &ncp);
   if(stat != NC_NOERR) return stat;
   return NC_get_vars(ncid,varid,startp, countp, stridep,
		      (void *)ip, T_double);
}

int
nc_get_vars_ubyte(int ncid, int varid, const size_t *startp, 
		  const size_t *countp, const ptrdiff_t * stridep,
		  unsigned char *ip)
{
   NC* ncp;
   int stat = NC_check_id(ncid, &ncp);
   if(stat != NC_NOERR) return stat;
   return NC_get_vars(ncid,varid, startp, countp, stridep,
		      (void *)ip, T_ubyte);
}

int
nc_get_vars_ushort(int ncid, int varid, const size_t *startp, 
		   const size_t *countp, const ptrdiff_t * stridep,
		   unsigned short *ip)
{
   NC* ncp;
   int stat = NC_check_id(ncid, &ncp);
   if(stat != NC_NOERR) return stat;
   return NC_get_vars(ncid,varid,startp,countp, stridep,
		      (void *)ip, T_ushort);
}

int
nc_get_vars_uint(int ncid, int varid, const size_t *startp, 
		 const size_t *countp, const ptrdiff_t * stridep,
		 unsigned int *ip)
{
   NC* ncp;
   int stat = NC_check_id(ncid, &ncp);
   if(stat != NC_NOERR) return stat;
   return NC_get_vars(ncid,varid,startp, countp, stridep,
		      (void *)ip, T_uint);
}

int
nc_get_vars_longlong(int ncid, int varid, const size_t *startp, 
		     const size_t *countp, const ptrdiff_t * stridep,
		     long long *ip)
{
   NC* ncp;
   int stat = NC_check_id(ncid, &ncp);
   if(stat != NC_NOERR) return stat;
   return NC_get_vars(ncid, varid, startp, countp, stridep,
		      (void *)ip, T_longlong);
}

int
nc_get_vars_ulonglong(int ncid, int varid, const size_t *startp, 
		      const size_t *countp, const ptrdiff_t * stridep,
		      unsigned long long *ip)
{
   NC* ncp;
   int stat = NC_check_id(ncid, &ncp);
   if(stat != NC_NOERR) return stat;
   return NC_get_vars(ncid, varid, startp, countp, stridep,
		      (void *)ip, NC_UINT64);
}

#ifdef USE_NETCDF4
int
nc_get_vars_string(int ncid, int varid,
		   const size_t *startp, const size_t *countp,
		   const ptrdiff_t * stridep,
		   char* *ip)
{
   NC* ncp;
   int stat = NC_check_id(ncid, &ncp);
   if(stat != NC_NOERR) return stat;
   return NC_get_vars(ncid, varid, startp, countp, stridep, 
		      (void *)ip, NC_STRING);
}
#endif /*USE_NETCDF4*/
/** \} */

/** \ingroup variables
Read a mapped array from a variable. 

The nc_get_varm_ type family of functions reads a mapped array section
of values from a netCDF variable of an open netCDF dataset. The mapped
array section is specified by giving a corner, a vector of edge
lengths, a stride vector, and an index mapping vector. The index
mapping vector is a vector of integers that specifies the mapping
between the dimensions of a netCDF variable and the in-memory
structure of the internal data array. No assumptions are made about
the ordering or length of the dimensions of the data array. The netCDF
dataset must be in data mode.

The functions for types ubyte, ushort, uint, longlong, ulonglong, and
string are only available for netCDF-4/HDF5 files.

The nc_get_varm() function will read a variable of any type, including
user defined type. For this function, the type of the data in memory
must match the type of the variable - no data conversion is done.

\param ncid NetCDF or group ID, from a previous call to nc_open(),
nc_create(), nc_def_grp(), or associated inquiry functions such as 
nc_inq_ncid().

\param varid Variable ID

\param startp Start vector with one element for each dimension to \ref
specify_hyperslab.

\param countp Count vector with one element for each dimension to \ref
specify_hyperslab.

\param stridep Stride vector with one element for each dimension to
\ref specify_hyperslab.

\param imapp Mapping vector with one element for each dimension to
\ref specify_hyperslab.

\param ip Pointer where the data will be copied. Memory must be
allocated by the user before this function is called.

\returns ::NC_NOERR No error.
\returns ::NC_ENOTVAR Variable not found.
\returns ::NC_EINVALCOORDS Index exceeds dimension bound.
\returns ::NC_ERANGE One or more of the values are out of range.
\returns ::NC_EINDEFINE Operation not allowed in define mode.
\returns ::NC_EBADID Bad ncid.
*/
/** \{ */
int
nc_get_varm(int ncid, int varid, const size_t * startp,
	    const size_t * countp, const ptrdiff_t * stridep,
	    const ptrdiff_t * imapp, void *ip)
{
   NC* ncp;
   int stat = NC_NOERR;

   if ((stat = NC_check_id(ncid, &ncp)))
       return stat;
   return ncp->dispatch->get_varm(ncid, varid, startp, countp, 
				  stridep, imapp, ip, NC_NAT);
}

int
nc_get_varm_schar(int ncid, int varid,
		  const size_t *startp, const size_t *countp,
		  const ptrdiff_t *stridep, 
		  const ptrdiff_t *imapp, signed char *ip)
{
   NC *ncp;
   int stat = NC_check_id(ncid, &ncp);
   if(stat != NC_NOERR) return stat;
   return NC_get_varm(ncid, varid, startp, countp,
		      stridep, imapp, (void *)ip, NC_BYTE);
}

int
nc_get_varm_uchar(int ncid, int varid,
		  const size_t *startp, const size_t *countp,
		  const ptrdiff_t *stridep, const ptrdiff_t *imapp,
		  unsigned char *ip)
{
   NC *ncp;
   int stat = NC_check_id(ncid, &ncp);
   if(stat != NC_NOERR) return stat;
   return NC_get_varm(ncid,varid,startp,countp,stridep,imapp, (void *)ip,T_uchar);
}

int
nc_get_varm_short(int ncid, int varid, const size_t *startp, 
		  const size_t *countp, const ptrdiff_t *stridep, 
		  const ptrdiff_t *imapp, short *ip)
{
   NC *ncp;
   int stat = NC_check_id(ncid, &ncp);
   if(stat != NC_NOERR) return stat;
   return NC_get_varm(ncid,varid,startp,countp,stridep,imapp, (void *)ip,NC_SHORT);
}

int
nc_get_varm_int(int ncid, int varid,
		const size_t *startp, const size_t *countp,
		const ptrdiff_t *stridep, const ptrdiff_t *imapp,
		int *ip)
{
   NC *ncp;
   int stat = NC_check_id(ncid, &ncp);
   if(stat != NC_NOERR) return stat;
   return NC_get_varm(ncid,varid,startp,countp,stridep,imapp, (void *)ip,NC_INT);
}

int
nc_get_varm_long(int ncid, int varid,
		 const size_t *startp, const size_t *countp,
		 const ptrdiff_t *stridep, const ptrdiff_t *imapp,
		 long *ip)
{
   NC *ncp;
   int stat = NC_check_id(ncid, &ncp);
   if(stat != NC_NOERR) return stat;
   return NC_get_varm(ncid,varid,startp,countp,stridep,imapp, (void *)ip,T_long);
}

int
nc_get_varm_float(int ncid, int varid,
		  const size_t *startp, const size_t *countp,
		  const ptrdiff_t *stridep, const ptrdiff_t *imapp,
		  float *ip)
{
   NC *ncp;
   int stat = NC_check_id(ncid, &ncp);
   if(stat != NC_NOERR) return stat;
   return NC_get_varm(ncid,varid,startp,countp,stridep,imapp, (void *)ip,T_float);
}

int
nc_get_varm_double(int ncid, int varid,
		   const size_t *startp, const size_t *countp,
		   const ptrdiff_t *stridep, const ptrdiff_t *imapp,
		   double *ip)
{
   NC *ncp;
   int stat = NC_check_id(ncid, &ncp);
   if(stat != NC_NOERR) return stat;
   return NC_get_varm(ncid,varid,startp,countp,stridep,imapp, (void *)ip,T_double);
}

int
nc_get_varm_ubyte(int ncid, int varid,
		  const size_t *startp, const size_t *countp,
		  const ptrdiff_t *stridep, const ptrdiff_t *imapp,
		  unsigned char *ip)
{
   NC *ncp;
   int stat = NC_check_id(ncid, &ncp);
   if(stat != NC_NOERR) return stat;
   return NC_get_varm(ncid,varid,startp,countp,stridep,
		      imapp, (void *)ip, T_ubyte);
}

int
nc_get_varm_ushort(int ncid, int varid,
		   const size_t *startp, const size_t *countp,
		   const ptrdiff_t *stridep, const ptrdiff_t *imapp,
		   unsigned short *ip)
{
   NC *ncp;
   int stat = NC_check_id(ncid, &ncp);
   if(stat != NC_NOERR) return stat;
   return NC_get_varm(ncid, varid, startp, countp, stridep,
		      imapp, (void *)ip, T_ushort);
}

int
nc_get_varm_uint(int ncid, int varid,
		 const size_t *startp, const size_t *countp,
		 const ptrdiff_t *stridep, const ptrdiff_t *imapp,
		 unsigned int *ip)
{
   NC *ncp;
   int stat = NC_check_id(ncid, &ncp);
   if(stat != NC_NOERR) return stat;
   return NC_get_varm(ncid, varid, startp, countp,
		      stridep, imapp, (void *)ip, T_uint);
}

int
nc_get_varm_longlong(int ncid, int varid, const size_t *startp, 
		     const size_t *countp, const ptrdiff_t *stridep, 
		     const ptrdiff_t *imapp, long long *ip)
{
   NC *ncp;
   int stat = NC_check_id(ncid, &ncp);
   if(stat != NC_NOERR) return stat;
   return NC_get_varm(ncid, varid, startp, countp, stridep, imapp,
		      (void *)ip, T_longlong);
}

int
nc_get_varm_ulonglong(int ncid, int varid,
		      const size_t *startp, const size_t *countp,
		      const ptrdiff_t *stridep, const ptrdiff_t *imapp,
		      unsigned long long *ip)
{
   NC *ncp;
   int stat = NC_check_id(ncid, &ncp);
   if(stat != NC_NOERR) return stat;
   return NC_get_varm(ncid, varid, startp, countp, stridep, imapp,
		      (void *)ip, NC_UINT64);
}

int
nc_get_varm_text(int ncid, int varid, const size_t *startp, 
		 const size_t *countp, const ptrdiff_t *stridep, 
		 const ptrdiff_t *imapp, char *ip)
{
   NC *ncp;
   int stat = NC_check_id(ncid, &ncp);
   if(stat != NC_NOERR) return stat;
   return NC_get_varm(ncid, varid, startp, countp, stridep, imapp,
		      (void *)ip, NC_CHAR);
}

#ifdef USE_NETCDF4
int
nc_get_varm_string(int ncid, int varid, const size_t *startp, 
		   const size_t *countp, const ptrdiff_t *stridep, 
		   const ptrdiff_t *imapp, char **ip)
{
   NC *ncp;
   int stat = NC_check_id(ncid, &ncp);
   if(stat != NC_NOERR) return stat;
   return NC_get_varm(ncid, varid, startp, countp, stridep, imapp,
		      (void *)ip, NC_STRING);
}
/** \} */
#endif /*USE_NETCDF4*/

/*! \} */ /* End of named group... */


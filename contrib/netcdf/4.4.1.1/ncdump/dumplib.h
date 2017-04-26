/*********************************************************************
 *   Copyright 1993, University Corporation for Atmospheric Research
 *   See netcdf/COPYRIGHT file for copying and redistribution conditions.
 *   $Header: /upc/share/CVS/netcdf-3/ncdump/dumplib.h,v 1.28 2009/08/13 21:06:13 russ Exp $
 *********************************************************************/
#ifndef _DUMPLIB_H_
#define _DUMPLIB_H_

#include "config.h"

extern char *progname;		/* for error messages */

#ifndef NO_NETCDF_2
#define NO_NETCDF_2		/* assert we aren't using any netcdf-2 stuff */
#endif

#ifndef EXIT_FAILURE
#ifndef vms
#define EXIT_SUCCESS 0
#define EXIT_FAILURE 1
#else
#define EXIT_SUCCESS 1
#define EXIT_FAILURE 0
#endif
#endif

/* more than enough characters needed to represent any single numeric
 * primitive value (e.g. int, float, double, long long, ...) using
 * printf functions */
#define PRIM_LEN 100

#define FLT_DIGITS 7		/* default sig. digits for float data */
#define DBL_DIGITS 15		/* default sig. digits for double data */

extern int float_precision_specified; /* -p option specified float precision */
extern int double_precision_specified; /* -p option specified double precision */
extern char float_var_fmt[];
extern char double_var_fmt[];
extern char float_att_fmt[];
extern char float_attx_fmt[];
extern char double_att_fmt[];

/* Display for netCDF-4 and HDF5 string values representing NULL
 * pointers rather than empty strings.  HDF5 distinguishes these two
 * kinds of string values so NULL pointers can be used as fill values
 * for lists of strings that might include empty strings. */
#define NIL_STRING "NIL"

#ifdef __cplusplus
extern "C" {
#endif

#ifndef HAVE_STRLCAT
/* Append src to dst of size siz */
extern size_t strlcat(char *dst, const char *src, size_t siz);
#endif /* ! HAVE_STRLCAT */

/* In case different formats specified with -d option, set them here. */
extern void	set_formats ( int flt_digs, int dbl_digs );

/* Determine print format to use for each value for this variable. */
const char *	get_fmt ( int ncid, int varid, nc_type typeid );

/* Add type info to type list */
extern void	typeadd ( nctype_t *typep );

/* From type id, get full type info (in current typlelist context)  */
extern nctype_t *get_typeinfo ( int typeid );

/* From type id, get type name */
extern void get_type_name(int ncid, nc_type type, char *name);

/* set tostring member function */
extern void set_tostring_func ( ncvar_t *varp);

/* helper function for output of opaque attributes */
extern int ncopaque_val_as_hex ( size_t size, char *buf, const void *valp );

/* Initialize list of types, only need primitive types for netCDF-3 */
extern void init_types ( int ncid );

/* Deallocate type list  */
/* extern void     xfree_typeinfo ( ); */

/* Test if variable is a coordinate variable */
extern int      iscoordvar ( int ncid, int varid );

/* Test if user-defined type */
extern int  is_user_defined_type ( nc_type type );

/* Initialize global constants used in slightly fuzzy float comparisons */
extern void init_epsilons ( void );

/* Initialize string buffer */
safebuf_t *sbuf_new();

/* Copy string s2 to buffer in sbuf, growing if necessary */
void sbuf_cpy(safebuf_t *sbuf, const char *s2);

/* Concatenate string s2 to end of buffer in sbuf, growing if necessary */
void sbuf_cat(safebuf_t *sbuf, const char *s2);

/* Concatenate sbuf s2 to end of sbuf s1, growing if necessary */
void sbuf_catb(safebuf_t *s1, const safebuf_t *s2);

/* Return length of the string in sbuf */
size_t sbuf_len(const safebuf_t *sbuf);

/* Return the C string inside an sbuf */
char *sbuf_str(const safebuf_t *sbuf);

/* Free string buffer */
void sbuf_free(safebuf_t *sbuf);

/* Print simple type name symbol or as absolute path, if necessary.
 * Escape any special chars. */
void print_type_name(int locid, int typeid);

int nctime_val_tostring(const ncvar_t *varp, safebuf_t *sfbf, const void *valp);

/* Return true if dimid is an unlimited dimension */
extern bool_t is_unlim_dim(int ncid, int dimid);

#ifdef __cplusplus
}
#endif

#endif	/*_DUMPLIB_H_ */

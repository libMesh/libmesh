/*********************************************************************
 *   Copyright 1996, UCAR/Unidata
 *   See netcdf/COPYRIGHT file for copying and redistribution conditions.
 *   $Id: tests.h 2785 2014-10-26 05:21:20Z wkliao $
 *********************************************************************/

#include <config.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <limits.h>
#include <float.h>
#define NO_NETCDF_2 1
#include "netcdf.h"
#ifdef USE_PARALLEL
#include "netcdf_par.h"
#endif
#include "error.h"

#if defined(_CRAY) && !defined(_CRAYIEEE) && !defined(__crayx1)
#define CRAYFLOAT 1 /* CRAY Floating point */
#elif defined(_SX) && defined(_FLOAT2)	/* NEC SUPER-UX in CRAY mode */
#define CRAYFLOAT 1 /* CRAY Floating point */
#endif

    /* Limits of external types (based on those in ncx.h) */

#ifdef __CHAR_UNSIGNED__
#define X_CHAR_MIN	SCHAR_MIN
#define X_CHAR_MAX	SCHAR_MAX
#else
#define X_CHAR_MIN	CHAR_MIN
#define X_CHAR_MAX	CHAR_MAX
#endif //__unsigned_char__

#define X_BYTE_MIN	(-128)
#define X_BYTE_MAX	127
#define X_SHORT_MIN	(-32768)
#define X_SHORT_MAX	32767
#define X_INT_MIN	(-2147483647-1)
#define X_INT_MAX	2147483647
#if defined(FLT_MAX_EXP) && FLT_MAX_EXP < 128
/* FLT_MAX < X_FLOAT_MAX */
#define X_FLOAT_MAX	FLT_MAX
#else
#ifdef WIN32 /* Windows, of course, has to be a *little* different. */
#define X_FLOAT_MAX	3.402823466e+38f
#else
#define X_FLOAT_MAX	3.40282347e+38f
#endif /* WIN32 */
#endif
#define X_FLOAT_MIN	(-X_FLOAT_MAX)
#if CRAYFLOAT
/* ldexp(1. - ldexp(.5 , -46), 1024) */
#define X_DOUBLE_MAX    1.79769313486230e+308
#else
/* scalb(1. - scalb(.5 , -52), 1024) */
#define X_DOUBLE_MAX	1.7976931348623157e+308
#endif
#define X_DOUBLE_MIN	(-X_DOUBLE_MAX)

#define X_SCHAR_MAX     SCHAR_MAX
#define X_SCHAR_MIN     SCHAR_MIN
#define X_UCHAR_MAX     UCHAR_MAX
#define X_UCHAR_MIN     0
#define X_UBYTE_MAX     X_UCHAR_MAX
#define X_UBYTE_MIN     X_UCHAR_MIN
#define X_USHORT_MAX    USHRT_MAX
#define X_USHORT_MIN    0
#define X_UINT_MAX      UINT_MAX
#define X_UINT_MIN      0

#ifndef LLONG_MAX
#define LLONG_MAX  0x7fffffffffffffffLL
#endif
#ifndef LLONG_MIN
#define LLONG_MIN (-0x7fffffffffffffffLL-1)
#endif
#ifndef ULLONG_MAX
#define ULLONG_MAX  0xffffffffffffffffULL
#endif

#ifndef X_INT64_MAX
#define X_INT64_MAX    LLONG_MAX
#endif
#ifndef X_INT64_MIN
#define X_INT64_MIN    LLONG_MIN
#endif
#ifndef X_UINT64_MAX
#define X_UINT64_MAX  ULLONG_MAX
#endif
#ifndef X_UINT64_MIN
#define X_UINT64_MIN  ULLONG_MIN
#endif


#if _SX /* NEC SUPER UX */
#if _INT64
#undef  INT_MAX /* workaround cpp bug */
#define INT_MAX  X_INT_MAX
#undef  INT_MIN /* workaround cpp bug */
#define INT_MIN  X_INT_MIN
#undef  LONG_MAX /* workaround cpp bug */
#define LONG_MAX  X_INT_MAX
#undef  LONG_MIN /* workaround cpp bug */
#define LONG_MIN  X_INT_MIN
#elif _LONG64
#undef  LONG_MAX /* workaround cpp bug */
#define LONG_MAX  4294967295L
#undef  LONG_MIN /* workaround cpp bug */
#define LONG_MIN -4294967295L
#endif
#endif /* _SX */


#ifndef MAX
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#endif /* MAX */

#ifndef MIN
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#endif /* MIN */

#ifndef ABS
#define ABS(x)  ((x) < 0 ? -(x) : (x))
#endif /* ABS */


    /* Parameters of test data */

#define NTYPES 11
#define NDIMS 5
#define NVARS 166
#define NRECS 2
#define NGATTS NTYPES
#define RECDIM 0
#define MAX_RANK 3
#define MAX_NELS 64
#define MAX_DIM_LEN 4
#define MAX_NATTS 3

extern int numGatts;  /* number of global attributes */
extern int numVars;   /* number of variables */
extern int numTypes;  /* number of netCDF data types to test */


    /* Limits of internal types */

#define text_min SCHAR_MIN
#define uchar_min 0
#define schar_min SCHAR_MIN
#define short_min SHRT_MIN
#define int_min INT_MIN
#define long_min LONG_MIN
#define float_min (-FLT_MAX)
#define double_min (-DBL_MAX)
#define ushort_min 0
#define uint_min 0
#define ulong_min 0
#define int64_min LLONG_MIN
#define longlong_min int64_min
#define uint64_min 0
#define ulonglong_min uint64_min

#define text_max SCHAR_MAX
#define uchar_max UCHAR_MAX
#define schar_max SCHAR_MAX
#define short_max SHRT_MAX
#define int_max INT_MAX
#define long_max LONG_MAX
#define float_max FLT_MAX
#define double_max DBL_MAX
#define ushort_max USHRT_MAX
#define uint_max UINT_MAX
#define ulong_max ULONG_MAX
#define int64_max LLONG_MAX
#define longlong_max int64_max
#define uint64_max ULLONG_MAX
#define ulonglong_max uint64_max



    /* Examples of invalid argument values */

#define BAD_ID -1               /* invalid netCDF ID */
#define BAD_DIMID -1            /* invalid dim ID */
#define BAD_VARID -2            /* invalid var ID */
#define BAD_ATTNUM -1           /* invalid att number */
#define BAD_TYPE (nc_type) 0    /* invalid data type */
#define BAD_FILLMODE -1         /* invalid fill mode */
#define BAD_NAME "a/b"		/* invalid name */
#define BAD_DEFAULT_FORMAT 12   /* invalid default format */

#define LEN_OF(array) ((sizeof array) / (sizeof array[0]))

#ifdef __cplusplus
extern "C" {
#endif


    /* Non-standard internal types */

#ifdef HAVE_SYS_TYPES_H
# include <sys/types.h>
#endif

typedef signed char text;
typedef signed char schar;

#ifndef HAVE_USHORT
typedef unsigned short int  ushort;
#endif

#ifndef HAVE_UINT
typedef unsigned       int  uint;
#endif

#ifndef HAVE_INT64
typedef          long long  int64;
#endif

#ifndef HAVE_UINT64
typedef unsigned long long  uint64;
#endif


    /* Global variables - filenames */

extern char testfile[];		/* netCDF read-only test data */
extern char scratch[];		/* netCDF test file for writing */

    /* Global variables - command-line arguments */

extern int  read_only;		/* if 1, don't try to change files */
extern int  verbose;		/* if 1, print details of tests */
extern int  nfails;		/* number of failures in specific test */
extern int  max_nmpt;		/* max. number of messages per test */

    /* Global variables - test data */

extern char dim_name[NDIMS][3];
extern size_t dim_len[NDIMS];
extern char var_name[NVARS][2+MAX_RANK];
extern nc_type var_type[NVARS];
extern size_t var_rank[NVARS];
extern int var_dimid[NVARS][MAX_RANK];
extern size_t var_shape[NVARS][MAX_RANK];
extern size_t var_nels[NVARS];
extern size_t var_natts[NVARS];
extern char att_name[NVARS][MAX_NATTS][2];
extern char gatt_name[NGATTS][3];
extern nc_type att_type[NVARS][NGATTS];
extern nc_type gatt_type[NGATTS];
extern size_t att_len[NVARS][MAX_NATTS];
extern size_t gatt_len[NGATTS];


    /* Macros for accessing attribute test data */
    /* varid is -1 for NC_GLOBAL so can do global atts in same loop */

#define VARID(varid)      (varid < 0 ? NC_GLOBAL : varid)
#define NATTS(varid)      (varid < 0 ? numGatts : var_natts[varid])
#define ATT_NAME(varid,j) (varid < 0 ? gatt_name[j] : att_name[varid][j])
#define ATT_TYPE(varid,j) (varid < 0 ? gatt_type[j] : att_type[varid][j])
#define ATT_LEN(varid,j)  (varid < 0 ? gatt_len[j] : att_len[varid][j])

extern const char *s_nc_type(nc_type);

extern void test_nc_strerror(void);
extern void test_nc_open(void);
extern void test_nc_close(void);

extern void test_nc_inq(void);
extern void test_nc_inq_natts(void);
extern void test_nc_inq_ndims(void);
extern void test_nc_inq_nvars(void);
extern void test_nc_inq_unlimdim(void);

extern void test_nc_inq_dimid(void);
extern void test_nc_inq_dim(void);
extern void test_nc_inq_dimlen(void);
extern void test_nc_inq_dimname(void);

extern void test_nc_inq_varid(void);
extern void test_nc_inq_vardimid(void);
extern void test_nc_inq_varname(void);
extern void test_nc_inq_varnatts(void);
extern void test_nc_inq_varndims(void);
extern void test_nc_inq_vartype(void);
extern void test_nc_inq_var(void);

extern void test_nc_get_var_double(void);
extern void test_nc_get_var_float(void);
extern void test_nc_get_var_int(void);
extern void test_nc_get_var_long(void);
extern void test_nc_get_var_schar(void);
extern void test_nc_get_var_short(void);
extern void test_nc_get_var_text(void);
extern void test_nc_get_var_uchar(void);
extern void test_nc_get_var(void);
extern void test_nc_get_var_ushort(void);
extern void test_nc_get_var_uint(void);
extern void test_nc_get_var_longlong(void);
extern void test_nc_get_var_ulonglong(void);

extern void test_nc_get_var1_double(void);
extern void test_nc_get_var1_float(void);
extern void test_nc_get_var1_int(void);
extern void test_nc_get_var1_long(void);
extern void test_nc_get_var1_schar(void);
extern void test_nc_get_var1_short(void);
extern void test_nc_get_var1_text(void);
extern void test_nc_get_var1_uchar(void);
extern void test_nc_get_var1(void);
extern void test_nc_get_var1_ushort(void);
extern void test_nc_get_var1_uint(void);
extern void test_nc_get_var1_longlong(void);
extern void test_nc_get_var1_ulonglong(void);

extern void test_nc_get_vara_double(void);
extern void test_nc_get_vara_float(void);
extern void test_nc_get_vara_int(void);
extern void test_nc_get_vara_long(void);
extern void test_nc_get_vara_schar(void);
extern void test_nc_get_vara_short(void);
extern void test_nc_get_vara_text(void);
extern void test_nc_get_vara_uchar(void);
extern void test_nc_get_vara(void);
extern void test_nc_get_vara_ushort(void);
extern void test_nc_get_vara_uint(void);
extern void test_nc_get_vara_longlong(void);
extern void test_nc_get_vara_ulonglong(void);

extern void test_nc_get_vars(void);
extern void test_nc_get_vars_double(void);
extern void test_nc_get_vars_float(void);
extern void test_nc_get_vars_int(void);
extern void test_nc_get_vars_long(void);
extern void test_nc_get_vars_schar(void);
extern void test_nc_get_vars_short(void);
extern void test_nc_get_vars_text(void);
extern void test_nc_get_vars_uchar(void);
extern void test_nc_get_vars(void);
extern void test_nc_get_vars_ushort(void);
extern void test_nc_get_vars_uint(void);
extern void test_nc_get_vars_longlong(void);
extern void test_nc_get_vars_ulonglong(void);

extern void test_nc_get_varm(void);
extern void test_nc_get_varm_double(void);
extern void test_nc_get_varm_float(void);
extern void test_nc_get_varm_int(void);
extern void test_nc_get_varm_long(void);
extern void test_nc_get_varm_schar(void);
extern void test_nc_get_varm_short(void);
extern void test_nc_get_varm_text(void);
extern void test_nc_get_varm_uchar(void);
extern void test_nc_get_varm(void);
extern void test_nc_get_varm_ushort(void);
extern void test_nc_get_varm_uint(void);
extern void test_nc_get_varm_longlong(void);
extern void test_nc_get_varm_ulonglong(void);

extern void test_nc_get_att(void);
extern void test_nc_get_att_double(void);
extern void test_nc_get_att_float(void);
extern void test_nc_get_att_int(void);
extern void test_nc_get_att_long(void);
extern void test_nc_get_att_schar(void);
extern void test_nc_get_att_short(void);
extern void test_nc_get_att_text(void);
extern void test_nc_get_att_uchar(void);
extern void test_nc_get_att_ushort(void);
extern void test_nc_get_att_uint(void);
extern void test_nc_get_att_longlong(void);
extern void test_nc_get_att_ulonglong(void);

extern void test_nc_put_att(void);
extern void test_nc_put_var_double(void);
extern void test_nc_put_var_float(void);
extern void test_nc_put_var_int(void);
extern void test_nc_put_var_long(void);
extern void test_nc_put_var_schar(void);
extern void test_nc_put_var_short(void);
extern void test_nc_put_var_text(void);
extern void test_nc_put_var_uchar(void);
extern void test_nc_put_var(void);
extern void test_nc_put_var_ushort(void);
extern void test_nc_put_var_uint(void);
extern void test_nc_put_var_longlong(void);
extern void test_nc_put_var_ulonglong(void);

extern void test_nc_put_var1_double(void);
extern void test_nc_put_var1_float(void);
extern void test_nc_put_var1_int(void);
extern void test_nc_put_var1_long(void);
extern void test_nc_put_var1_schar(void);
extern void test_nc_put_var1_short(void);
extern void test_nc_put_var1_text(void);
extern void test_nc_put_var1_uchar(void);
extern void test_nc_put_var1(void);
extern void test_nc_put_var1_ushort(void);
extern void test_nc_put_var1_uint(void);
extern void test_nc_put_var1_longlong(void);
extern void test_nc_put_var1_ulonglong(void);

extern void test_nc_put_vara_double(void);
extern void test_nc_put_vara_float(void);
extern void test_nc_put_vara_int(void);
extern void test_nc_put_vara_long(void);
extern void test_nc_put_vara_schar(void);
extern void test_nc_put_vara_short(void);
extern void test_nc_put_vara_text(void);
extern void test_nc_put_vara_uchar(void);
extern void test_nc_put_vara(void);
extern void test_nc_put_vara_ushort(void);
extern void test_nc_put_vara_uint(void);
extern void test_nc_put_vara_longlong(void);
extern void test_nc_put_vara_ulonglong(void);

extern void test_nc_put_vars_double(void);
extern void test_nc_put_vars_float(void);
extern void test_nc_put_vars_int(void);
extern void test_nc_put_vars_long(void);
extern void test_nc_put_vars_schar(void);
extern void test_nc_put_vars_short(void);
extern void test_nc_put_vars_text(void);
extern void test_nc_put_vars_uchar(void);
extern void test_nc_put_vars(void);
extern void test_nc_put_vars_ushort(void);
extern void test_nc_put_vars_uint(void);
extern void test_nc_put_vars_longlong(void);
extern void test_nc_put_vars_ulonglong(void);

extern void test_nc_put_varm_double(void);
extern void test_nc_put_varm_float(void);
extern void test_nc_put_varm_int(void);
extern void test_nc_put_varm_long(void);
extern void test_nc_put_varm_schar(void);
extern void test_nc_put_varm_short(void);
extern void test_nc_put_varm_text(void);
extern void test_nc_put_varm_uchar(void);
extern void test_nc_put_varm(void);
extern void test_nc_put_varm_ushort(void);
extern void test_nc_put_varm_uint(void);
extern void test_nc_put_varm_longlong(void);
extern void test_nc_put_varm_ulonglong(void);

extern void test_nc_put_att_double(void);
extern void test_nc_put_att_float(void);
extern void test_nc_put_att_int(void);
extern void test_nc_put_att_long(void);
extern void test_nc_put_att_schar(void);
extern void test_nc_put_att_short(void);
extern void test_nc_put_att_text(void);
extern void test_nc_put_att_uchar(void);
extern void test_nc_put_att_ushort(void);
extern void test_nc_put_att_uint(void);
extern void test_nc_put_att_longlong(void);
extern void test_nc_put_att_ulonglong(void);

extern void test_nc_create(void);
extern void test_nc_redef(void);
extern void test_nc_enddef(void);
extern void test_nc_sync(void);
extern void test_nc_abort(void);
extern void test_nc_def_dim(void);
extern void test_nc_rename_dim(void);
extern void test_nc_def_var(void);
extern void test_nc_rename_var(void);
extern void test_nc_copy_att(void);

extern void test_nc_inq_att(void);
extern void test_nc_inq_attname(void);
extern void test_nc_inq_attid(void);
extern void test_nc_inq_attlen(void);
extern void test_nc_inq_atttype(void);

extern void test_nc_rename_att(void);
extern void test_nc_del_att(void);
extern void test_nc_set_fill(void);
extern void test_nc_set_default_format(void);

void print_nok(int nok);

int inRange(const double value, const nc_type datatype);

/*
 * internal types
 */
typedef enum {
	NCT_UNSPECIFIED = 0,
	NCT_UCHAR =	1,	/* unsigned char */
	NCT_TEXT =	16,	/* char */
#define NCT_CHAR NCT_TEXT
	NCT_SCHAR =	17,	/* signed char */
	NCT_SHORT =	18,	/* short */
	NCT_INT =	20,	/* int */
	NCT_LONG =	22,	/* long */
	NCT_FLOAT =	36,	/* float */
	NCT_DOUBLE =	40,	/* double */
        NCT_USHORT =    41,
        NCT_UINT =      42,
        NCT_INT64 =     43,
#define NCT_LONGLONG NCT_INT64
        NCT_UINT64 =    44
#define NCT_ULONGLONG NCT_UINT64
} nct_itype;

int inRange3(const double value, const nc_type datatype, const nct_itype itype);

int equal(const double x, const double y, nc_type extType, nct_itype itype);

int int_vec_eq(const int *v1, const int *v2, const int n);

int roll( int n );

int
toMixedBase(
    size_t number,        /* number to be converted to mixed base */
    size_t length,
    const size_t base[],        /* dimensioned [length], base[0] ignored */
    size_t result[]);      /* dimensioned [length] */

size_t
fromMixedBase(
    size_t length,
    size_t number[],      /* dimensioned [length] */
    size_t base[]);        /* dimensioned [length], base[0] ignored */

int nc2dbl ( const nc_type datatype, const void *p, double *result);

int dbl2nc ( const double d, const nc_type datatype, void *p);

double hash( const nc_type type, const int rank, const size_t *index );

double hash4(
    const nc_type type,
    const int rank,
    const size_t *index,
    const nct_itype itype);

void init_gvars(void);

void def_dims(int ncid);

void def_vars(int ncid);

void put_atts(int ncid);

void put_vars(int ncid);

void write_file(char *filename);

void check_dims(int  ncid);

void check_vars(int  ncid);

void check_atts(int  ncid);

void check_file(char *filename);

int file_create(const char *filename, int cmode, int *ncid);

int file_open(const char *filename, int omode, int *ncid);

#ifdef __cplusplus
}
#endif

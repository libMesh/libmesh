/*********************************************************************
 *   Copyright 1993, University Corporation for Atmospheric Research
 *   See netcdf/README file for copying and redistribution conditions.
 *   $Header: /upc/share/CVS/netcdf-3/ncdump/dumplib.c,v 1.85 2010/05/05 22:15:39 dmh Exp $
 *********************************************************************/

/*
 * We potentially include <stdarg.h> before <stdio.h> in order to obtain a
 * definition for va_list from the GNU C compiler.
 */

#include <config.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#ifndef NO_FLOAT_H
#include <float.h>		/* for FLT_EPSILON, DBL_EPSILON */
#endif /* NO_FLOAT_H */
#include <math.h>
#include <netcdf.h>
#include "utils.h"
#include "nccomps.h"
#include "dumplib.h"
#include "ncdump.h"
#include "isnan.h"
#include "nctime0.h"

static float float_eps;
static double double_eps;

extern fspec_t formatting_specs; /* set from command-line options */

static float
float_epsilon(void)
{
    float float_eps;
#ifndef NO_FLOAT_H
    float_eps = FLT_EPSILON;
#else /* NO_FLOAT_H */
    {
	float etop, ebot, eps;
	float one = 1.0;
	float two = 2.0;
	etop = 1.0;
	ebot = 0.0;
	eps = ebot + (etop - ebot)/two;
	while (eps != ebot && eps != etop) {
	    float epsp1;

	    epsp1 = one + eps;
	    if (epsp1 > one)
		etop = eps;
	    else
		ebot = eps;
	    eps = ebot + (etop - ebot)/two;
	}
	float_eps = two * etop;
    }
#endif /* NO_FLOAT_H */
    return float_eps;
}


static double
double_epsilon(void)
{
    double double_eps;
#ifndef NO_FLOAT_H
    double_eps = DBL_EPSILON;
#else /* NO_FLOAT_H */
    {
	double etop, ebot, eps;
	double one = 1.0;
	double two = 2.0;
	etop = 1.0;
	ebot = 0.0;
	eps = ebot + (etop - ebot)/two;
	while (eps != ebot && eps != etop) {
	    double epsp1;

	    epsp1 = one + eps;
	    if (epsp1 > one)
		etop = eps;
	    else
		ebot = eps;
	    eps = ebot + (etop - ebot)/two;
	}
	double_eps = two * etop;
    }
#endif /* NO_FLOAT_H */
    return double_eps;
}


void
init_epsilons(void)
{
    float_eps = float_epsilon();
    double_eps = double_epsilon();
}


static char* has_c_format_att(int ncid, int varid);

int float_precision_specified = 0; /* -p option specified float precision */
int double_precision_specified = 0; /* -p option specified double precision */
char float_var_fmt[] = "%.NNg";
char double_var_fmt[] = "%.NNg";
char float_att_fmt[] = "%#.NNgf";
char float_attx_fmt[] = "%#.NNg";
char double_att_fmt[] = "%#.NNg";

#ifndef HAVE_STRLCAT
/*	$OpenBSD: strlcat.c,v 1.12 2005/03/30 20:13:52 otto Exp $	*/

/*
 * Copyright (c) 1998 Todd C. Miller <Todd.Miller@courtesan.com>
 *
 * Permission to use, copy, modify, and distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 * WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
 * ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
 * ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
 * OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 */

/*
 * Appends src to string dst of size siz (unlike strncat, siz is the
 * full size of dst, not space left).  At most siz-1 characters
 * will be copied.  Always NUL terminates (unless siz <= strlen(dst)).
 * Returns strlen(src) + MIN(siz, strlen(initial dst)).
 * If retval >= siz, truncation occurred.
 */
size_t
strlcat(char *dst, const char *src, size_t siz)
{
	char *d = dst;
	const char *s = src;
	size_t n = siz;
	size_t dlen;

	/* Find the end of dst and adjust bytes left but don't go past end */
	while (n-- != 0 && *d != '\0')
		d++;
	dlen = d - dst;
	n = siz - dlen;

	if (n == 0)
		return(dlen + strlen(s));
	while (*s != '\0') {
		if (n != 1) {
			*d++ = *s;
			n--;
		}
		s++;
	}
	*d = '\0';

	return(dlen + (s - src));	/* count does not include NUL */
}
#endif /* ! HAVE_STRLCAT */


/* magic number stored in a safebuf and checked, hoping it will be
 * changed if buffer was overwritten inadvertently */
#define SAFEBUF_CERT 2147114711

/* expression for where SAFEBUF_CERT is stored within safebuf (at end
 * of buffer, after data) */
#define SAFEBUF_EXPR(sbuf) (*(int *)((sbuf)->buf + (sbuf)->len))

/* expression to be checked whenever a safebuf is used */
#define SAFEBUF_CHECK(sbuf) (SAFEBUF_EXPR(sbuf) == SAFEBUF_CERT)

/* somewhat arbitrary initial size of safebufs, grow as needed */
#define SAFEBUF_INIT_LEN 128

/* initialize safe buffer */
safebuf_t *
sbuf_new() {
    size_t len = SAFEBUF_INIT_LEN;
    safebuf_t *sb;
    sb = (safebuf_t *) emalloc(sizeof(safebuf_t));
    sb->buf = (char *)emalloc(len + sizeof(int));
    sb->len = len;
    /* write a "stamp" in last 4 bytes of buffer for id and to check for overflow */
    SAFEBUF_EXPR(sb) = SAFEBUF_CERT;
    sb->buf[0] = 0;
    sb->cl = strlen(sb->buf);
    assert(SAFEBUF_CHECK(sb));
    return sb;
}


/* grow buffer to at least len bytes, copying previous contents if
 * necessary */
void
sbuf_grow(safebuf_t *sb, size_t len) {
    size_t m = sb->len;
    void *tmp;
    assert(SAFEBUF_CHECK(sb));
    if (len <= m)
	return;

    /* Make sure we at least double size of buffer to get what's
     * needed.  If we just used realloc(), no guarantee that length
     * would be expanded by a multiple, which we want. */
    while(len > m) {
	m *= 2;
    }
    tmp = emalloc(m + sizeof(int));
    memcpy(tmp, sb->buf, sb->len);
    sb->len = m;
    free(sb->buf);
    sb->buf = tmp;
    SAFEBUF_EXPR(sb) = SAFEBUF_CERT;
    assert(SAFEBUF_CHECK(sb));
}

/* Copy string s2 to safe buffer, growing if necessary */
void
sbuf_cpy(safebuf_t *sb, const char *s2) {
    size_t s2len;
    assert(SAFEBUF_CHECK(sb));
    s2len = strlen(s2);
    sbuf_grow(sb, 1 + s2len);
    strncpy(sb->buf, s2, sb->len);
    sb->cl = s2len;
    assert(SAFEBUF_CHECK(sb));
}

/* Concatenate string s2 to end of string in safe buffer, growing if necessary */
void
sbuf_cat(safebuf_t *sb, const char *s2) {
    size_t s2len;
    size_t res;
    assert(SAFEBUF_CHECK(sb));
    s2len = strlen(s2);
    sbuf_grow(sb, 1 + sb->cl + s2len);
    res = strlcat(sb->buf + sb->cl, s2, sb->len);
    assert( res < sb->len );
    sb->cl += s2len;
    assert(SAFEBUF_CHECK(sb));
}

/* Concatenate string in safebuf s2 to end of string in safebuf s1,
 * growing if necessary */
void
sbuf_catb(safebuf_t *s1, const safebuf_t *s2) {
    size_t s2len;
    size_t res;
    assert(SAFEBUF_CHECK(s1));
    assert(SAFEBUF_CHECK(s2));
    s2len = sbuf_len(s2);
    sbuf_grow(s1, 1 + s1->cl + s2len);
    res = strlcat(s1->buf + s1->cl, s2->buf, s1->len);
    assert( res < s1->len );
    s1->cl += s2len;
    assert(SAFEBUF_CHECK(s1));
}

/* Return length of string in sbuf */
size_t
sbuf_len(const safebuf_t *sb) {
    assert(SAFEBUF_CHECK(sb));
    return sb->cl;
}

/* Return C string in an sbuf */
char *
sbuf_str(const safebuf_t *sb) {
    assert(SAFEBUF_CHECK(sb));
    return sb->buf;
}

/* free safe buffer */
void
sbuf_free(safebuf_t *sb) {
    assert(SAFEBUF_CHECK(sb));
    free(sb->buf);
    free(sb);
}


/* In case different formats specified with -d option, set them here. */
void
set_formats(int float_digits, int double_digits)
{
    int res;
    res = snprintf(float_var_fmt, strlen(float_var_fmt) + 1, "%%.%dg",
		   float_digits) + 1;
    assert(res <= sizeof(float_var_fmt));
    res = snprintf(double_var_fmt, strlen(double_var_fmt) + 1, "%%.%dg",
		   double_digits) + 1;
    assert(res <= sizeof(double_var_fmt));
    res = snprintf(float_att_fmt, strlen(float_att_fmt) + 1, "%%#.%dgf",
		   float_digits) + 1;
    assert(res <= sizeof(float_att_fmt));
    res = snprintf(float_attx_fmt, strlen(float_attx_fmt) + 1, "%%#.%dg",
		   float_digits) + 1;
    assert(res <= sizeof(float_attx_fmt));
    res = snprintf(double_att_fmt, strlen(double_att_fmt) + 1, "%%#.%dg",
		   double_digits) + 1;
    assert(res <= sizeof(double_att_fmt));
}


static char *
has_c_format_att(
    int ncid,			/* netcdf id */
    int varid			/* variable id */
    )
{
    nc_type cfmt_type;
    size_t cfmt_len;
#define C_FMT_NAME	"C_format" /* name of C format attribute */
#define	MAX_CFMT_LEN	100	/* max length of C format attribute */
    static char cfmt[MAX_CFMT_LEN];

    /* we expect nc_inq_att to fail if there is no "C_format" attribute */
    int nc_stat = nc_inq_att(ncid, varid, "C_format", &cfmt_type, &cfmt_len);

    switch(nc_stat) {
    case NC_NOERR:
	if (cfmt_type == NC_CHAR && cfmt_len != 0 && cfmt_len < MAX_CFMT_LEN) {
	    int nc_stat = nc_get_att_text(ncid, varid, "C_format", cfmt);
	    if(nc_stat != NC_NOERR) {
		fprintf(stderr, "Getting 'C_format' attribute %s\n",
			nc_strerror(nc_stat));
		(void) fflush(stderr);
	    }
	    cfmt[cfmt_len] = '\0';
	    return &cfmt[0];
	}
	break;
    case NC_ENOTATT:
	break;
    default:
	fprintf(stderr, "Inquiring about 'C_format' attribute %s\n",
		nc_strerror(nc_stat));
	(void) fflush(stderr);
	break;
    }
    return 0;
}


/* Return default format to use for a primitive type */
const char *
get_default_fmt(nc_type typeid) {
    /* Otherwise return sensible default. */
    switch (typeid) {
    case NC_BYTE:
	return "%d";
    case NC_CHAR:
	return "%s";
    case NC_SHORT:
	return "%d";
    case NC_INT:
	return "%d";
    case NC_FLOAT:
	return float_var_fmt;
    case NC_DOUBLE:
	return double_var_fmt;
    case NC_UBYTE:
	return "%u";
    case NC_USHORT:
	return "%u";
    case NC_UINT:
	return "%u";
    case NC_INT64:
	return "%lld";
    case NC_UINT64:
	return "%llu";
    case NC_STRING:
	return "\"%s\"";
    default:
	break;
    }
    return "";		/* user-defined types don't use fmt member */
}

/*
 * Determine print format to use for each primitive value for this
 * variable.  Use value of attribute C_format if it exists, otherwise
 * a sensible default.
 */
const char *
get_fmt(
     int ncid,
     int varid,
     nc_type typeid
     )
{
    char *c_format_att;

    /* float or double precision specified with -p option overrides any
       C_format attribute value, so check for that first. */
    if (float_precision_specified && typeid == NC_FLOAT)
	return float_var_fmt;
    if (double_precision_specified && typeid == NC_DOUBLE)
	return double_var_fmt;
    /* If C_format attribute exists, return it */
    c_format_att = has_c_format_att(ncid, varid);
    if (c_format_att)
      return c_format_att;
    return get_default_fmt(typeid);
}

/* Return primitive type name */
static const char *
prim_type_name(nc_type type)
{
    switch (type) {
      case NC_BYTE:
	return "byte";
      case NC_CHAR:
	return "char";
      case NC_SHORT:
	return "short";
      case NC_INT:
	return "int";
      case NC_FLOAT:
	return "float";
      case NC_DOUBLE:
	return "double";
      case NC_UBYTE:
	return "ubyte";
      case NC_USHORT:
	return "ushort";
      case NC_UINT:
	return "uint";
      case NC_INT64:
	return "int64";
      case NC_UINT64:
	return "uint64";
      case NC_STRING:
	return "string";
      default:
	error("prim_type_name: bad type %d", type);
	return "bogus";
    }
}

static int max_type = 0;
static int max_atomic_type = 0;
static nctype_t **nctypes = 0;	/* holds all types in a netCDF dataset */


#ifdef USE_NETCDF4
/* return number of user-defined types in a group and all its subgroups */
static int
count_udtypes(int ncid) {
    int ntypes = 0;
    int numgrps;
    int *ncids;
    int i;
    int format;

    NC_CHECK( nc_inq_format(ncid, &format) );

    if (format == NC_FORMAT_NETCDF4) {
	/* Get number of types in this group */
	NC_CHECK( nc_inq_typeids(ncid, &ntypes, NULL) ) ;
	NC_CHECK( nc_inq_grps(ncid, &numgrps, NULL) ) ;
	ncids = (int *) emalloc(sizeof(int) * (numgrps + 1));
	NC_CHECK( nc_inq_grps(ncid, NULL, ncids) ) ;
	/* Add number of types in each subgroup, if any */
	for (i=0; i < numgrps; i++) {
	    ntypes += count_udtypes(ncids[i]);
	}
	free(ncids);
    }
    return ntypes;
}
#endif /*USE_NETCDF4*/

/* This routine really is intended to return the max atomic typeid */
static int
max_typeid(int ncid) {
    int maxtypes = NC_NAT;
    int maxatomictypes = NC_NAT;
    int format = 0;
    int err = NC_NOERR;

    /* get the file type */
    err = nc_inq_format(ncid,&format);
    if(err) {
	fprintf(stderr,"%s: Cannot get file format.\n",nc_strerror(err));
	return 0;
    }
    switch (format) {
    case NC_FORMAT_CLASSIC:
    case NC_FORMAT_NETCDF4_CLASSIC:
    case NC_FORMAT_64BIT_OFFSET:
        maxatomictypes = (maxtypes = NC_DOUBLE); /*ignore NC_NAT?*/
	break;
    case NC_FORMAT_64BIT_DATA:
	maxatomictypes = (maxtypes = NC_UINT64);
	break;
    case NC_FORMAT_NETCDF4:
#ifdef USE_NETCDF4
	{
        int nuser = 0;
        maxatomictypes = (maxtypes = NC_STRING); /* extra netCDF-4 primitive types */
        maxtypes += 4;		/* user-defined classes */
        nuser = count_udtypes(ncid);
        if(nuser > 0)
	    maxtypes = NC_FIRSTUSERTYPEID + (nuser - 1);
	} break;
#else
	/* fallthru */
#endif
    default:
	fprintf(stderr,"Unexpected file format: %d\n",format);
	return 0;
    }
    max_type = maxtypes;
    max_atomic_type = maxatomictypes;
    return maxtypes;
}

void typeadd(nctype_t *typep) {
    nctypes[typep->tid] = typep;
}

/* From type id, get full type info */
nctype_t *
get_typeinfo ( int typeid ) {
    if(typeid < 0 || typeid > max_type)
	error("ncdump: %d is an invalid type id", typeid);
    return nctypes[typeid];
}

/* void */
/* xfree_typeinfo(int ncid) { */
/*     int i; */
/*     for (i = 0; i < number_of_types; i++) { */
/* 	nctype_t *tinfop = nctypes[i]; */
/* 	if (tinfop) { */
/* 	    if(tinfop->name) */
/* 		free(tinfop->name); */
/* 	    if(tinfop->grps) */
/* 		free(tinfop->grps); */
/* 	    free(tinfop); */
/* 	} */
/*     } */
/* } */


bool_t
ncbyte_val_equals(const nctype_t *this,
		  const void *v1p, const void *v2p) {
    return ( *(signed char* )v1p == *(signed char* )v2p);
}

bool_t
ncchar_val_equals(const nctype_t *this,
		  const void *v1p, const void *v2p) {
    return ( *(char* )v1p == *(char* )v2p);
}

bool_t
ncshort_val_equals(const nctype_t *this,
		   const void *v1p, const void *v2p) {
    return ( *(short* )v1p == *(short* )v2p);
}

bool_t
ncint_val_equals(const nctype_t *this,
		 const void *v1p, const void *v2p) {
    return ( *(int* )v1p == *(int* )v2p);
}

#define absval(x)  ( (x) < 0 ? -(x) : (x) )

/*
 * Return ( *(float* )v1p == *(float* )v2p);
 * except use floating epsilon to compare very close vals as equal
 * and handle IEEE NaNs and infinities.
 */
bool_t
ncfloat_val_equals(const nctype_t *this,
		   const void *v1p, const void *v2p) {
    float v1 = *(float* )v1p;
    float v2 = *(float* )v2p;
    if((v1 > 0.0f) != (v2 > 0.0f))	/* avoid overflow */
	return false;
    if(isfinite(v1) && isfinite(v2))
	return (absval(v1 - v2) <= absval(float_eps * v2)) ;
    if(isnan(v1) && isnan(v2))
	return true;
    if(isinf(v1) && isinf(v2))
	return true;
    return false;
}

/*
 * Return ( *(double* )v1p == *(double* )v2p);
 * except use floating epsilon to compare very close vals as equal
 * and handle IEEE NaNs and infinities.
 */
bool_t
ncdouble_val_equals(const nctype_t *this,
		    const void *v1p, const void *v2p) {
    double v1 = *(double* )v1p;
    double v2 = *(double* )v2p;
    if((v1 > 0.0) != (v2 > 0.0))	/* avoid overflow */
	return false;
    if(isfinite(v1) && isfinite(v2))
	return (absval(v1 - v2) <= absval(double_eps * v2)) ;
    if(isnan(v1) && isnan(v2))
	return true;
    if(isinf(v1) && isinf(v2))
	return true;
    return false;
}

bool_t
ncubyte_val_equals(const nctype_t *this,
		   const void *v1p, const void *v2p) {
    return ( *(unsigned char* )v1p == *(unsigned char* )v2p);
}

bool_t
ncushort_val_equals(const nctype_t *this,
		    const void *v1p, const void *v2p) {
    return ( *(unsigned short* )v1p == *(unsigned short* )v2p);
}

bool_t
ncuint_val_equals(const nctype_t *this,
		  const void *v1p, const void *v2p) {
    return ( *(unsigned int* )v1p == *(unsigned int* )v2p);
}

bool_t
ncint64_val_equals(const nctype_t *this,
		   const void *v1p, const void *v2p) {
    return ( *(long long* )v1p == *(long long* )v2p);
}

bool_t
ncuint64_val_equals(const nctype_t *this,
		    const void *v1p, const void *v2p) {
    return ( *(unsigned long long* )v1p == *(unsigned long long* )v2p);
}

bool_t
ncstring_val_equals(const nctype_t *this,
		    const void *v1p, const void *v2p) {
    if (NULL == *((char **)v1p) && NULL == *((char **)v2p))
        return(1);
    else if (NULL != *((char **)v1p) && NULL == *((char **)v2p))
        return(0);
    else if (NULL == *((char **)v1p) && NULL != *((char **)v2p))
        return(0);
    return (strcmp(*((char **)v1p), *((char **)v2p)) == 0);
}

#ifdef USE_NETCDF4
bool_t
ncopaque_val_equals(const nctype_t *this,
		    const void *v1p, const void *v2p) {
    size_t nbytes = this->size;
    const char *c1p = (const char *) v1p;
    const char *c2p = (const char *) v2p;
    int i;
    for (i=0; i < nbytes; i++) {
	if (*c1p++ != *c2p++)
	    return false;
    }
    return true;
}

bool_t
ncvlen_val_equals(const nctype_t *this,
		  const void *v1p, const void *v2p) {
    size_t v1len = ((nc_vlen_t *)v1p)->len;
    size_t v2len = ((nc_vlen_t *)v2p)->len;
    if (v1len != v2len)
	return false;
    {
	size_t base_size = this->size;
	nc_type base_type = this->base_tid;
	nctype_t *base_info = get_typeinfo(base_type);
	val_equals_func base_val_equals = base_info->val_equals;
	const char *v1dat = ((nc_vlen_t *)v1p)->p;
	const char *v2dat = ((nc_vlen_t *)v2p)->p;
	size_t i;
	for(i = 0; i < v1len; i++) {
	    if (base_val_equals(base_info, (const void *)v1dat,
				(const void *)v2dat) != true)
		return false;
	    v1dat += base_size;
	    v2dat += base_size;
	}
    }
    return true;
}

/* Determine if two compound values are equal, by testing eqaulity of
 * each member field. */
bool_t
nccomp_val_equals(const nctype_t *this,
		  const void *v1p, const void *v2p) {
    int nfields = this->nfields;
    int fidx;			/* field id */

    for (fidx = 0; fidx < nfields; fidx++) {
	size_t offset = this->offsets[fidx];
	nc_type fid = this->fids[fidx];		/* field type id */
	nctype_t *finfo = get_typeinfo(fid);
	if(finfo->ranks == 0 || finfo->ranks[fidx] == 0) {
	    if(! finfo->val_equals(finfo,
				   (char *)v1p + offset, (char *)v2p + offset))
		return false;
	} else {		/* this field is an array */
	    int i;		/* array element counter when rank > 0 */
	    void *v1elem = (char *)v1p + offset;
	    void *v2elem = (char *)v2p + offset;
	    for(i = 0; i < finfo->nvals[fidx]; i++) {
		if(! finfo->val_equals(finfo, v1elem, v2elem))
		    return false;
		v1elem = (char *)v1elem + finfo->size;
		v2elem = (char *)v1elem + finfo->size;
	    }
	}
    }
    return true;
}
#endif /* USE_NETCDF4 */

int
ncbyte_typ_tostring(const nctype_t *typ, safebuf_t *sfbf, const void *valp) {
    char sout[PRIM_LEN];
    int res;
    res = snprintf(sout, PRIM_LEN, typ->fmt, *(signed char *)valp);
    assert(res < PRIM_LEN);
    sbuf_cpy(sfbf, sout);
    return sbuf_len(sfbf);
}

int
ncchar_typ_tostring(const nctype_t *typ, safebuf_t *sfbf, const void *valp) {
    char sout[PRIM_LEN];
    int res;
    res = snprintf(sout, PRIM_LEN, typ->fmt, *(char *)valp);
    assert(res < PRIM_LEN);
    sbuf_cpy(sfbf, sout);
    return sbuf_len(sfbf);
}

int
ncshort_typ_tostring(const nctype_t *typ, safebuf_t *sfbf, const void *valp) {
    char sout[PRIM_LEN];
    int res;
    res = snprintf(sout, PRIM_LEN, typ->fmt, *(short *)valp);
    assert(res < PRIM_LEN);
    sbuf_cpy(sfbf, sout);
    return sbuf_len(sfbf);
}

int
ncint_typ_tostring(const nctype_t *typ, safebuf_t *sfbf, const void *valp) {
    char sout[PRIM_LEN];
    int res;
    res = snprintf(sout, PRIM_LEN, typ->fmt, *(int *)valp);
    assert(res < PRIM_LEN);
    sbuf_cpy(sfbf, sout);
    return sbuf_len(sfbf);
}

/* CDL canonical representations of some special floating point values */
#define NCDL_NANF "NaNf"
#define NCDL_NAN "NaN"
#define NCDL_INFF "Infinityf"
#define NCDL_INF "Infinity"

/* Convert a float NaN or Infinity to an allocated string large enough
 * to hold it (at least PRIM_LEN chars) */
static void
float_special_tostring(float vv, char *sout) {
    if(isnan(vv)) {
	snprintf(sout, PRIM_LEN, "%s", NCDL_NANF);
    } else if(isinf(vv)) {
	if(vv < 0.0) {
	    snprintf(sout, PRIM_LEN, "-%s", NCDL_INFF);
	} else {
	    snprintf(sout, PRIM_LEN, "%s", NCDL_INFF);
	}
    } else
	assert(false);		/* vv was finite */
}

/* Convert a double NaN or Infinity to an allocated string large enough
 * to hold it (at least PRIM_LEN chars) */
static void
double_special_tostring(double vv, char *sout) {
    if(isnan(vv)) {
	    snprintf(sout, PRIM_LEN, "%s", NCDL_NAN);
    } else if(isinf(vv)) {
	if(vv < 0.0) {
	    snprintf(sout, PRIM_LEN, "-%s", NCDL_INF);
	} else {
	    snprintf(sout, PRIM_LEN, "%s", NCDL_INF);
	}
    } else
	assert(false);		/* vv was finite */
}

int
ncfloat_typ_tostring(const nctype_t *typ, safebuf_t *sfbf, const void *valp) {
    char sout[PRIM_LEN];
    float vv = *(float *)valp;
    if(isfinite(vv)) {
	int res;
	res = snprintf(sout, PRIM_LEN, typ->fmt, vv);
	assert(res < PRIM_LEN);
    } else {
	float_special_tostring(vv, sout);
    }
    sbuf_cpy(sfbf, sout);
    return sbuf_len(sfbf);
}

int
ncdouble_typ_tostring(const nctype_t *typ, safebuf_t *sfbf, const void *valp) {
    char sout[PRIM_LEN];
    double vv = *(double *)valp;
    if(isfinite(vv)) {
	int res;
	res = snprintf(sout, PRIM_LEN, typ->fmt, vv);
	assert(res < PRIM_LEN);
    } else {
	double_special_tostring(vv, sout);
    }
    sbuf_cpy(sfbf, sout);
    return sbuf_len(sfbf);
}

int
ncubyte_typ_tostring(const nctype_t *typ, safebuf_t *sfbf, const void *valp) {
    char sout[PRIM_LEN];
    int res;
    res = snprintf(sout, PRIM_LEN, typ->fmt, *(unsigned char *)valp);
    assert(res < PRIM_LEN);
    sbuf_cpy(sfbf, sout);
    return sbuf_len(sfbf);
}

int
ncushort_typ_tostring(const nctype_t *typ, safebuf_t *sfbf, const void *valp) {
    char sout[PRIM_LEN];
    int res;
    res = snprintf(sout, PRIM_LEN, typ->fmt, *(unsigned short *)valp);
    assert(res < PRIM_LEN);
    sbuf_cpy(sfbf, sout);
    return sbuf_len(sfbf);
}

int
ncuint_typ_tostring(const nctype_t *typ, safebuf_t *sfbf, const void *valp) {
    char sout[PRIM_LEN];
    int res;
    res = snprintf(sout, PRIM_LEN, typ->fmt, *(unsigned int *)valp);
    assert(res < PRIM_LEN);
    sbuf_cpy(sfbf, sout);
    return sbuf_len(sfbf);
}

int
ncint64_typ_tostring(const nctype_t *typ, safebuf_t *sfbf, const void *valp) {
    char sout[PRIM_LEN];
    int res;
    res = snprintf(sout, PRIM_LEN, typ->fmt, *(long long *)valp);
    assert(res < PRIM_LEN);
    sbuf_cpy(sfbf, sout);
    return sbuf_len(sfbf);
}

int
ncuint64_typ_tostring(const nctype_t *typ, safebuf_t *sfbf, const void *valp) {
    char sout[PRIM_LEN];
    int res;
    res = snprintf(sout, PRIM_LEN, typ->fmt, *(unsigned long long *)valp);
    assert(res < PRIM_LEN);
    sbuf_cpy(sfbf, sout);
    return sbuf_len(sfbf);
}

int ncstring_typ_tostring(const nctype_t *typ, safebuf_t *sfbf, const void *valp) {
    const char *cp;

    cp = ((char **)valp)[0];
    if(cp) {
        size_t slen;
        char *sout;
        char *sp;
        unsigned char uc;

        slen = 3 + 5 * strlen(cp); /* need "'s around string, and extra space to escape control characters */
        sout = emalloc(slen);
        sp = sout;
        *sp++ = '"' ;
        while(*cp) {
            switch (uc = *cp++ & 0377) {
            case '\b':
                *sp++ = '\\';
                *sp++ = 'b' ;
                break;
            case '\f':
                *sp++ = '\\';
                *sp++ = 'f';
                break;
            case '\n':
                *sp++ = '\\';
                *sp++ = 'n';
                break;
            case '\r':
                *sp++ = '\\';
                *sp++ = 'r';
                break;
            case '\t':
                *sp++ = '\\';
                *sp++ = 't';
                break;
            case '\v':
                *sp++ = '\\';
                *sp++ = 'n';
                break;
            case '\\':
                *sp++ = '\\';
                *sp++ = '\\';
                break;
            case '\'':
                *sp++ = '\\';
                *sp++ = '\'';
                break;
            case '\"':
                *sp++ = '\\';
                *sp++ = '\"';
                break;
            default:
                if (iscntrl(uc)) {
                    snprintf(sp,3,"\\%03o",uc);
                    sp += 4;
                }
                else
                    *sp++ = uc;
                break;
            }
        }
        *sp++ = '"' ;
        *sp = '\0' ;
        sbuf_cpy(sfbf, sout);
        free(sout);
    }
    else {
        sbuf_cpy(sfbf, "NIL");
    }
    return sbuf_len(sfbf);
}

#ifdef USE_NETCDF4
int
ncenum_typ_tostring(const nctype_t *typ, safebuf_t *sfbf, const void *valp) {
    char symbol[NC_MAX_NAME + 1];
    long long val;

    switch (typ->base_tid) {
    case NC_BYTE:
	val = *(signed char *)valp;
	break;
    case NC_UBYTE:
	val = *(unsigned char *)valp;
	break;
    case NC_SHORT:
	val = *(short *)valp;
	break;
    case NC_USHORT:
	val = *(unsigned short *)valp;
	break;
    case NC_INT:
	val = *(int *)valp;
	break;
    case NC_UINT:
	val = *(unsigned int *)valp;
	break;
    case NC_INT64:
	val = *(long long *)valp;
	break;
    case NC_UINT64:
	val = *(long long *)valp;
	break;
    default:
	error("bad base type for enum");
	break;
    }
    NC_CHECK( nc_inq_enum_ident(typ->ncid, typ->tid, val, symbol));
    sbuf_cpy(sfbf, symbol);
    return sbuf_len(sfbf);
}

/* Given an opaque type size and opaque value, convert to a string,
 * represented as hexadecimal characters, returning number of chars in
 * output string */
int
ncopaque_val_as_hex(size_t size, char *sout, const void *valp) {
    const unsigned char *cp = valp;
    char *sp = sout;
    int i;
    char *prefix = "0X";
    int prelen = strlen(prefix);

    snprintf(sp, prelen + 1, "%s", prefix);
    sp += prelen;
    for(i = 0; i < size; i++) {
	int res;
	res = snprintf(sp, prelen + 1, "%.2X", *cp++);
	assert (res == 2);
	sp += 2;
    }
    *sp = '\0';
    return 2*size + prelen;
}

/* Convert an opaque value to a string, represented as hexadecimal
 * characters */
int
ncopaque_typ_tostring(const nctype_t *typ, safebuf_t *sfbf,
		      const void *valp) {
    char* sout = (char *) emalloc(2 * typ->size + strlen("0X") + 1);
    (void) ncopaque_val_as_hex(typ->size, sout, valp);
    sbuf_cpy(sfbf, sout);
    free(sout);
    return sbuf_len(sfbf);
}

/* Convert a vlen value to a string, by using tostring function for base type */
int
ncvlen_typ_tostring(const nctype_t *tinfo, safebuf_t *sfbf, const void *valp) {
    nc_type base_type = tinfo->base_tid;
    nctype_t *base_info = get_typeinfo(base_type);
    size_t base_size = base_info->size;
    size_t len = ((nc_vlen_t *)valp)->len;
    typ_tostring_func base_typ_tostring = base_info->typ_tostring;
    size_t i;
    const char *vp;		/* instead of void* so can increment to next */
    safebuf_t* sout2 = sbuf_new();

    sbuf_cpy(sfbf, "{");
    /* put each val in sout2, then append sout2 to sfbf */
    vp = ((nc_vlen_t *)valp)->p;
    for(i = 0; i < len; i++) {
	(void) base_typ_tostring(base_info, sout2, vp);
	sbuf_catb(sfbf, sout2);
	if(i < len - 1) {
	    sbuf_cat(sfbf, ", ");
	}
	vp += base_size;
    }
    sbuf_cat(sfbf, "}");
    sbuf_free(sout2);
    return sbuf_len(sfbf);
}

/*
 * Print a number of char values as a text string.
 */
static int
chars_tostring(
    safebuf_t *sbuf,		/* for output */
    size_t len,			/* number of characters */
    const char *vals		/* pointer to block of values */
    )
{
    long iel;
    const char *sp;
    char *sout = (char *)emalloc(4*len + 5); /* max len of string */
    char *cp = sout;
    *cp++ = '"';

    /* adjust len so trailing nulls don't get printed */
    sp = vals + len;
    while (len != 0 && *--sp == '\0')
	len--;
    for (iel = 0; iel < len; iel++) {
	unsigned char uc;
	switch (uc = *vals++ & 0377) {
	case '\b':
	case '\f':
	case '\n':
	case '\r':
	case '\t':
	case '\v':
	case '\\':
	case '\'':
	case '\"':
	    *cp++ = '\\';
	    *cp++ = *(char *)&uc; /* just copy, even if char is signed */
	    break;
	default:
	    if (isprint(uc))
		*cp++ = *(char *)&uc; /* just copy, even if char is signed */
	    else {
		sprintf(cp,"\\%.3o",uc);
		cp += 4;
	    }
	    break;
	}
    }
    *cp++ = '"';
    *cp = '\0';
    sbuf_cpy(sbuf, sout);
    free(sout);
    return sbuf_len(sbuf);
}


/* Convert a compound value to a string, by using tostring function for
   each member field */
int
nccomp_typ_tostring(const nctype_t *tinfo, safebuf_t *sfbf, const void *valp) {
    int nfields = tinfo->nfields;
    int fidx;			/* field id */
    safebuf_t* sout2 = sbuf_new();

    sbuf_cpy(sfbf, "{");
    /* put each val in sout2, then append sout2 to sfbf if enough room */
    for (fidx = 0; fidx < nfields; fidx++) {
	size_t offset = tinfo->offsets[fidx];
	nc_type fid = tinfo->fids[fidx];		/* field type id */
	nctype_t *finfo = get_typeinfo(fid);

	if(tinfo->ranks[fidx] == 0) {
	    if(finfo->tid == NC_CHAR) { /* aggregate char rows into strings */
		chars_tostring(sout2, 1, ((char *)valp + offset));
	    } else {
		finfo->typ_tostring(finfo, sout2, ((char *)valp + offset));
	    }
	} else {		/* this field is an array */
	    int i;		/* array element counter when rank > 0 */
	    void *vp = (char *)valp + offset;
	    safebuf_t *sout3 = sbuf_new();
	    sbuf_cpy(sout2, "{");
	    if(finfo->tid == NC_CHAR) { /* aggregate char rows into strings */
		int rank = tinfo->ranks[fidx];
		size_t nstrings;
		size_t slen;
		int j;
		slen = tinfo->sides[fidx][rank-1];
		nstrings = 1;	/* product of all but last array dimension */
		for(j=0; j < rank-1; j++) {
		    nstrings *= tinfo->sides[fidx][j];
		}
		for(i=0; i < nstrings; i++) { /* loop on product of all but
						 last index of array */
		    chars_tostring(sout3, slen, (char *)vp);
		    vp = (char *)vp + slen;
		    if(i < nstrings - 1) {
			sbuf_cat(sout3, ", ");
		    }
		    sbuf_catb(sout2, sout3);
		}
	    } else {
		for(i = 0; i < tinfo->nvals[fidx]; i++) {
		    (void) finfo->typ_tostring(finfo, sout3, vp);
		    vp = (char *)vp + finfo->size;
		    if(i < tinfo->nvals[fidx] - 1) {
			sbuf_cat(sout3, ", ");
		    }
		    sbuf_catb(sout2, sout3);
		}
	    }
	    sbuf_cat(sout2, "}");
	    sbuf_free(sout3);
	}
	sbuf_catb(sfbf, sout2);
	if(fidx < nfields - 1) {
	    sbuf_cat(sfbf, ", ");
	}
    }
    sbuf_cat(sfbf, "}");
    sbuf_free(sout2);
    return sbuf_len(sfbf);
}
#endif	/* USE_NETCDF4 */

int
ncbyte_val_tostring(const ncvar_t *varp, safebuf_t *sfbf, const void *valp) {
    char sout[PRIM_LEN];
    int res;
    res = snprintf(sout, PRIM_LEN, varp->fmt, *(signed char *)valp);
    assert(res < PRIM_LEN);
    sbuf_cpy(sfbf, sout);
    return sbuf_len(sfbf);
}

int
ncchar_val_tostring(const ncvar_t *varp, safebuf_t *sfbf, const void *valp) {
    char sout[PRIM_LEN];
    int res;
    res = snprintf(sout, PRIM_LEN, varp->fmt, *(char *)valp);
    assert(res < PRIM_LEN);
    sbuf_cpy(sfbf, sout);
    return sbuf_len(sfbf);
}

int
ncshort_val_tostring(const ncvar_t *varp, safebuf_t *sfbf, const void *valp) {
    char sout[PRIM_LEN];
    int res;
    res = snprintf(sout, PRIM_LEN, varp->fmt, *(short *)valp);
    assert(res < PRIM_LEN);
    sbuf_cpy(sfbf, sout);
    return sbuf_len(sfbf);
}

int
ncint_val_tostring(const ncvar_t *varp, safebuf_t *sfbf, const void *valp) {
    char sout[PRIM_LEN];
    int res;
    res = snprintf(sout, PRIM_LEN, varp->fmt, *(int *)valp);
    assert(res < PRIM_LEN);
    sbuf_cpy(sfbf, sout);
    return sbuf_len(sfbf);
}

int
ncfloat_val_tostring(const ncvar_t *varp, safebuf_t *sfbf, const void *valp) {
    char sout[PRIM_LEN];
    float vv = *(float *)valp;
    if(isfinite(vv)) {
	int res;
	res = snprintf(sout, PRIM_LEN, varp->fmt, vv);
	assert(res < PRIM_LEN);
    } else {
	float_special_tostring(vv, sout);
    }
    sbuf_cpy(sfbf, sout);
    return sbuf_len(sfbf);
}

int
ncdouble_val_tostring(const ncvar_t *varp, safebuf_t *sfbf, const void *valp) {
    char sout[PRIM_LEN];
    double vv = *(double *)valp;
    if(isfinite(vv)) {
	int res;
	res = snprintf(sout, PRIM_LEN, varp->fmt, vv);
	assert(res < PRIM_LEN);
    } else {
	double_special_tostring(vv, sout);
    }
    sbuf_cpy(sfbf, sout);
    return sbuf_len(sfbf);
}

/* Convert value of any numeric type to a double.  Beware, this may
 * lose precision for values of type NC_INT64 or NC_UINT64 */
static
double to_double(const ncvar_t *varp, const void *valp) {
    double dd;
    switch (varp->type) {
    case NC_BYTE:
	dd = *(signed char *)valp;
	break;
    case NC_SHORT:
	dd = *(short *)valp;
	break;
    case NC_INT:
	dd = *(int *)valp;
	break;
    case NC_FLOAT:
	dd = *(float *)valp;
	break;
    case NC_DOUBLE:
	dd = *(double *)valp;
	break;
    case NC_UBYTE:
	dd = *(unsigned char *)valp;
	break;
    case NC_USHORT:
	dd = *(unsigned short *)valp;
	break;
    case NC_UINT:
	dd = *(unsigned int *)valp;
	break;
    case NC_INT64:
	dd = *(long long *)valp;
	break;
    case NC_UINT64:
	dd = *(unsigned long long *)valp;
	break;
    default:
	error("to_double: type not numeric primitive");
    }
    return dd;
}

int
nctime_val_tostring(const ncvar_t *varp, safebuf_t *sfbf, const void *valp) {
    char sout[PRIM_LEN];
    double vv = to_double(varp, valp);
    int separator = formatting_specs.iso_separator ? 'T' : ' ';
    if(isfinite(vv)) {
	int res;
	sout[0]='"';
	cdRel2Iso(varp->timeinfo->calendar, varp->timeinfo->units, separator, vv, &sout[1]);
	res = strlen(sout);
	sout[res++] = '"';
	sout[res] = '\0';
	assert(res < PRIM_LEN);
    } else {
	double_special_tostring(vv, sout);
    }
    sbuf_cpy(sfbf, sout);
    return sbuf_len(sfbf);
}

int
ncubyte_val_tostring(const ncvar_t *varp, safebuf_t *sfbf, const void *valp) {
    char sout[PRIM_LEN];
    int res;
    res = snprintf(sout, PRIM_LEN, varp->fmt, *(unsigned char *)valp);
    assert(res < PRIM_LEN);
    sbuf_cpy(sfbf, sout);
    return sbuf_len(sfbf);
}

int
ncushort_val_tostring(const ncvar_t *varp, safebuf_t *sfbf, const void *valp) {
    char sout[PRIM_LEN];
    int res;
    res = snprintf(sout, PRIM_LEN, varp->fmt, *(unsigned short *)valp);
    assert(res < PRIM_LEN);
    sbuf_cpy(sfbf, sout);
    return sbuf_len(sfbf);
}

int
ncuint_val_tostring(const ncvar_t *varp, safebuf_t *sfbf, const void *valp) {
    char sout[PRIM_LEN];
    int res;
    res = snprintf(sout, PRIM_LEN, varp->fmt, *(unsigned int *)valp);
    assert(res < PRIM_LEN);
    sbuf_cpy(sfbf, sout);
    return sbuf_len(sfbf);
}

int
ncint64_val_tostring(const ncvar_t *varp, safebuf_t *sfbf, const void *valp) {
    char sout[PRIM_LEN];
    int res;
    res = snprintf(sout, PRIM_LEN, varp->fmt, *(long long *)valp);
    assert(res < PRIM_LEN);
    sbuf_cpy(sfbf, sout);
    return sbuf_len(sfbf);
}

int
ncuint64_val_tostring(const ncvar_t *varp, safebuf_t *sfbf, const void *valp) {
    char sout[PRIM_LEN];
    int res;
    res = snprintf(sout, PRIM_LEN, varp->fmt, *(unsigned long long *)valp);
    assert(res < PRIM_LEN);
    sbuf_cpy(sfbf, sout);
    return sbuf_len(sfbf);
}

int
ncstring_val_tostring(const ncvar_t *varp, safebuf_t *sfbf, const void *valp) {
    return ncstring_typ_tostring(varp->tinfo, sfbf, valp);
}

#ifdef USE_NETCDF4
int
ncenum_val_tostring(const ncvar_t *varp, safebuf_t *sfbf, const void *valp) {
    return ncenum_typ_tostring(varp->tinfo, sfbf, valp);
}

/* Convert an opaque value to a string, represented as hexadecimal
 * characters */
int
ncopaque_val_tostring(const ncvar_t *varp, safebuf_t *sfbf, const void *valp) {
    return ncopaque_typ_tostring(varp->tinfo, sfbf, valp);
}

/* Convert a vlen value to a string, by using tostring function for base type */
int
ncvlen_val_tostring(const ncvar_t *varp, safebuf_t *sfbf, const void *valp) {
    return ncvlen_typ_tostring(varp->tinfo, sfbf, valp);
}

int
nccomp_val_tostring(const ncvar_t *varp, safebuf_t *sfbf, const void *valp) {
    return nccomp_typ_tostring(varp->tinfo, sfbf, valp);
}
#endif /*USE_NETCDF4*/

static val_equals_func eq_funcs[] = {
	ncbyte_val_equals,
	ncchar_val_equals,
	ncshort_val_equals,
	ncint_val_equals,
	ncfloat_val_equals,
	ncdouble_val_equals,
	ncubyte_val_equals,
	ncushort_val_equals,
	ncuint_val_equals,
	ncint64_val_equals,
	ncuint64_val_equals,
	ncstring_val_equals
    };

static typ_tostring_func ts_funcs[] = {
	ncbyte_typ_tostring,
	ncchar_typ_tostring,
	ncshort_typ_tostring,
	ncint_typ_tostring,
	ncfloat_typ_tostring,
	ncdouble_typ_tostring,
	ncubyte_typ_tostring,
	ncushort_typ_tostring,
	ncuint_typ_tostring,
	ncint64_typ_tostring,
	ncuint64_typ_tostring,
	ncstring_typ_tostring
    };


/* Set function pointer of function to convert a value to a string for
 * the variable pointed to by varp. */
void
set_tostring_func(ncvar_t *varp) {
    val_tostring_func tostring_funcs[] = {
	ncbyte_val_tostring,
	ncchar_val_tostring,
	ncshort_val_tostring,
	ncint_val_tostring,
	ncfloat_val_tostring,
	ncdouble_val_tostring,
	ncubyte_val_tostring,
	ncushort_val_tostring,
	ncuint_val_tostring,
	ncint64_val_tostring,
	ncuint64_val_tostring,
	ncstring_val_tostring
    };
    if(varp->has_timeval && formatting_specs.string_times) {
	varp->val_tostring = (val_tostring_func) nctime_val_tostring;
	return;
    }
    if( !is_user_defined_type(varp->type) ) {
	varp->val_tostring = tostring_funcs[varp->type - 1];
	return;
    }
#ifdef USE_NETCDF4
    switch(varp->tinfo->class) {
    case NC_VLEN:
	varp->val_tostring = (val_tostring_func) ncvlen_val_tostring;
	break;
    case NC_OPAQUE:
	varp->val_tostring = (val_tostring_func) ncopaque_val_tostring;
	break;
    case NC_ENUM:
	varp->val_tostring = (val_tostring_func) ncenum_val_tostring;
	break;
    case NC_COMPOUND:
	varp->val_tostring = (val_tostring_func) nccomp_val_tostring;
	break;
    default:
	error("unrecognized class of user defined type: %d",
	      varp->tinfo->class);
    }
#endif /* USE_NETCDF4 */
    return;
}


/* Initialize typelist with primitive types.  For netCDF-3 only need primitive
   types. */
static void
init_prim_types(int ncid) {
    nctype_t *tp;
    int i;
    int types[] =
{
	NC_BYTE,
	NC_CHAR,
	NC_SHORT,
	NC_INT,
	NC_FLOAT,
	NC_DOUBLE,
	NC_UBYTE,
	NC_USHORT,
	NC_UINT,
	NC_INT64,
	NC_UINT64,
	NC_STRING
    };
    size_t sizes[] = {
	sizeof(char),
	sizeof(char),
	sizeof(short),
	sizeof(int),
	sizeof(float),
	sizeof(double),
	sizeof(unsigned char),
	sizeof(unsigned short),
	sizeof(unsigned int),
	sizeof(long long),
	sizeof(unsigned long long),
	sizeof(char **)
    };

#if 0
    for(i=0; i < sizeof(types)/sizeof(int); i++) {
#else
    for(i=0; i < max_atomic_type; i++) {
#endif
	tp = (nctype_t *)emalloc(sizeof(nctype_t));
	tp->ncid = ncid;
	tp->tid = types[i];
	tp->name = strdup(prim_type_name(tp->tid));
	tp->grps = 0;
	tp->class = 0;		/* primitive type */
	tp->size = sizes[i];
	tp->base_tid = NC_NAT;	/* not used for primitive types */
	tp->nfields = 0;	/* not used for primitive types */
	tp->fmt = get_default_fmt(types[i]);
	tp->fids = 0;		/* not used for primitive types */
	tp->offsets = 0;	/* not used for primitive types */
	tp->ranks = 0;		/* not used for primitive types */
	tp->sides = 0;		/* not used for primitive types */
	tp->nvals = 0;		/* not used for primitive types */
	tp->val_equals = (val_equals_func) eq_funcs[i];
	tp->typ_tostring = (typ_tostring_func) ts_funcs[i];
	typeadd(tp);
    }
}

/* Initialize typelist.
 *
 * This must be done over all groups in netCDF-4, because
 * variables in one group may be declared using types in a
 * different group.  For netCDF-3, this is just the info about
 * primitive types.
 */
void
init_types(int ncid) {
#ifdef USE_NETCDF4
    int ntypes;
#endif
    if (max_type == 0) {	/* if called for first time */
	int maxtype = max_typeid(ncid);
	int i;
	nctypes = (nctype_t **) emalloc((maxtype + 2) * sizeof(nctype_t *));
	for(i=0; i < maxtype+1; i++)
	    nctypes[i] = NULL;	/* so can later skip over unused type slots */
	init_prim_types(ncid);
    }

#ifdef USE_NETCDF4
   /* Are there any user defined types in this group? */
   NC_CHECK( nc_inq_typeids(ncid, &ntypes, NULL) );
   if (ntypes)
   {
      int t;
      int *typeids = emalloc((ntypes + 1) * sizeof(int));
      NC_CHECK( nc_inq_typeids(ncid, NULL, typeids) );
      for (t = 0; t < ntypes; t++) {
	  nctype_t *tinfo;	/* details about the type */
	  char type_name[NC_MAX_NAME + 1];
	  size_t group_name_len;
	  char* group_name;
	  int fidx;		/* for compound type, field index */

	  tinfo = (nctype_t *) emalloc(sizeof(nctype_t));

	  NC_CHECK( nc_inq_user_type(ncid, typeids[t], type_name, &tinfo->size,
		                     &tinfo->base_tid, &tinfo->nfields,
				     &tinfo->class) );
	  tinfo->tid = typeids[t];
	  tinfo->ncid = ncid;
	  tinfo->name = strdup(type_name);
	  tinfo->grps = 0;
	  if(tinfo->class == NC_VLEN) {
	      tinfo->size = sizeof(nc_vlen_t); /* not size of base type */
	  }
	  NC_CHECK( nc_inq_grpname_full(ncid, &group_name_len, NULL) );
	  group_name = (char *) emalloc(group_name_len + 1);
	  NC_CHECK( nc_inq_grpname_full(ncid, &group_name_len, group_name) );

	  tinfo->grps = strdup(group_name);
	  free(group_name);
	  switch(tinfo->class) {
	  case NC_ENUM:
	      tinfo->val_equals = eq_funcs[tinfo->base_tid-1];
	      tinfo->typ_tostring = (typ_tostring_func) ncenum_typ_tostring;
	      break;
	  case NC_COMPOUND:
	      tinfo->val_equals = (val_equals_func) nccomp_val_equals;
	      tinfo->typ_tostring = (typ_tostring_func) nccomp_typ_tostring;
	      tinfo->fids = (nc_type *) emalloc((tinfo->nfields + 1)
						  * sizeof(nc_type));
	      tinfo->offsets = (size_t *) emalloc((tinfo->nfields + 1)
						  * sizeof(size_t));
	      tinfo->ranks = (int *) emalloc((tinfo->nfields + 1)
					     * sizeof(int));
	      tinfo->sides = (int **) emalloc((tinfo->nfields + 1)
						 * sizeof(int *));
	      tinfo->nvals = (int *) emalloc((tinfo->nfields + 1)
					     * sizeof(int));
	      for (fidx = 0; fidx < tinfo->nfields; fidx++) {
		  size_t offset;
		  nc_type ftype;
		  int rank;
		  int *sides;
		  int i;
		  sides = NULL;
		  NC_CHECK( nc_inq_compound_field(ncid, tinfo->tid, fidx, NULL,
						  &offset, &ftype, &rank,
						  sides) );
		  if(rank > 0) sides = (int *) emalloc(rank * sizeof(int));
		  NC_CHECK( nc_inq_compound_field(ncid, tinfo->tid, fidx, NULL,
						  NULL, NULL, NULL, sides) );
		  tinfo->fids[fidx] = ftype;
		  tinfo->offsets[fidx] = offset;
		  tinfo->ranks[fidx] = rank;
		  if (rank > 0)
		      tinfo->sides[fidx] = (int *) emalloc(rank * sizeof(int));
		  tinfo->nvals[fidx] = 1;
		  for(i = 0; i < rank; i++) {
		      tinfo->sides[fidx][i] = sides[i];
		      tinfo->nvals[fidx] *= sides[i];
		  }
		  if (rank > 0)
		      free(sides);
	      }
	      break;
	  case NC_VLEN:
	      tinfo->val_equals = (val_equals_func) ncvlen_val_equals;
	      tinfo->typ_tostring = (typ_tostring_func) ncvlen_typ_tostring;
	      break;
	  case NC_OPAQUE:
	      tinfo->val_equals = (val_equals_func) ncopaque_val_equals;
	      tinfo->typ_tostring = (typ_tostring_func) ncopaque_typ_tostring;
	      break;
	  default:
	      error("bad class: %d", tinfo->class);
	      break;
	  }

	  typeadd(tinfo);
      }
      free(typeids);
   }
   /* For netCDF-4, check to see if this group has any subgroups and call
    * recursively on each of them. */
   {
      int g, numgrps, *ncids;

      /* See how many groups there are. */
      NC_CHECK( nc_inq_grps(ncid, &numgrps, NULL) );
      if (numgrps > 0) {
	  ncids = (int *) emalloc(numgrps * sizeof(int));
	  /* Get the list of group ids. */
	  NC_CHECK( nc_inq_grps(ncid, NULL, ncids) );
	  /* Call this function for each group. */
	  for (g = 0; g < numgrps; g++) {
	      init_types(ncids[g]);
	  }
	  free(ncids);
      }
   }
#endif /* USE_NETCDF4 */
}


/*
 * return 1 if varid identifies a coordinate variable
 * else return 0
 */
int
iscoordvar(int ncid, int varid)
{
    int ndims, ndims1;
    int dimid;
    int* dimids = 0;
    ncdim_t *dims = 0;
#ifdef USE_NETCDF4
    int include_parents = 1;
#endif
    int is_coord = 0;		/* true if variable is a coordinate variable */
    char varname[NC_MAX_NAME];
    int varndims;

    do {	  /* be safe in case someone is currently adding
		   * dimensions */
#ifdef USE_NETCDF4
	NC_CHECK( nc_inq_dimids(ncid, &ndims, NULL, include_parents ) );
#else
	NC_CHECK( nc_inq_ndims(ncid, &ndims) );
#endif
	if (dims)
	    free(dims);
	dims = (ncdim_t *) emalloc((ndims + 1) * sizeof(ncdim_t));
	if (dimids)
	    free(dimids);
	dimids = (int *) emalloc((ndims + 1) * sizeof(int));
#ifdef USE_NETCDF4
	NC_CHECK( nc_inq_dimids(ncid, &ndims1, dimids, include_parents ) );
#else
	{
	    int i;
	    for(i = 0; i < ndims; i++) {
		dimids[i] = i;	/* for netCDF-3, dimids are 0, 1, ..., ndims-1 */
	    }
	    NC_CHECK( nc_inq_ndims(ncid, &ndims1) );
	}
#endif	/* USE_NETCDF4 */
    } while (ndims != ndims1);

    for (dimid = 0; dimid < ndims; dimid++) {
	NC_CHECK( nc_inq_dimname(ncid, dimids[dimid], dims[dimid].name) );
    }
    NC_CHECK( nc_inq_varname(ncid, varid, varname) );
    NC_CHECK( nc_inq_varndims(ncid, varid, &varndims) );

    for (dimid = 0; dimid < ndims; dimid++) {
	if (strcmp(dims[dimid].name, varname) == 0 && varndims == 1) {
	    is_coord = 1;
	    break;
	}
    }
    if(dims)
	free(dims);
    if(dimids)
	free(dimids);
    return is_coord;
}


/* Return true if user-defined type */
int
is_user_defined_type(nc_type type) {
    nctype_t *typeinfop = get_typeinfo(type);
    return (typeinfop->class > 0);
}


/*
 * Return name of type in user-allocated space, whether built-in
 * primitive type or user-defined type.  Note: name must have enough
 * space allocated to hold type name.
 */
void
get_type_name(int ncid, nc_type type, char *name)
{
#ifdef USE_NETCDF4
    if (is_user_defined_type(type)) {
	NC_CHECK(nc_inq_user_type(ncid, type, name, NULL, NULL, NULL, NULL));
    } else {
	strncpy(name, prim_type_name(type), NC_MAX_NAME + 1);
    }
#else
    strncpy(name, prim_type_name(type), NC_MAX_NAME + 1);
#endif /* USE_NETCDF4 */
}

/*
 * Print type name with CDL escapes for special characters.  locid is
 * the id of the group in which the type is referenced, which is
 * needed to determine whether an absolute type name must be printed.
 * If the type is defined in the referenced group or in some ancestor
 * group, only the simple type name is printed.  If the type is
 * defined in some other non-ancestor group, an absolute path for the
 * typename is printed instead.
 */
void
print_type_name(int locid, int typeid) {
    char *ename;
#ifdef USE_NETCDF4
    char name[NC_MAX_NAME+1];
    int type_inherited = 0;
    int curlocid;		/* group we are searching in */
    int parent_groupid = locid;
    int ntypes;
    int stat;
#endif

    assert(typeid > 0 && typeid <= max_type);
    ename = escaped_name(nctypes[typeid]->name);
#ifdef USE_NETCDF4
    if(is_user_defined_type(typeid)) {
	/* determine if type is inherited, that is if defined in this
	 * group or any ancestor group */
	name[NC_MAX_NAME] = '\0';
	strncpy(name,nctypes[typeid]->name,NC_MAX_NAME);
	do {
	    curlocid = parent_groupid;
	    NC_CHECK( nc_inq_typeids(curlocid, &ntypes, NULL) );
	    if(ntypes > 0) {
		int *typeids = (int *) emalloc((ntypes + 1) * sizeof(int));
		int i;
		NC_CHECK( nc_inq_typeids(curlocid, &ntypes, typeids) );
		for(i = 0; i < ntypes; i++) {
		    char curname[NC_MAX_NAME];
		    NC_CHECK( nc_inq_type(curlocid, typeids[i], curname, NULL) );
		    if(strncmp(name, curname, NC_MAX_NAME) == 0) {
			type_inherited = 1;
			break;
		    }
		}
		free(typeids);
		if(type_inherited)
		    break;
	    }
	    stat = nc_inq_grp_parent(curlocid, &parent_groupid);
	} while (stat != NC_ENOGRP && stat != NC_ENOTNC4);

	if (type_inherited == 0) {
	    char *gname = nctypes[typeid]->grps;
	    print_name(gname);
	    fputs("/", stdout);
	}
    }
#endif	/*  USE_NETCDF4 */
    fputs(ename, stdout);
    free(ename);
}

/* Allocate and initialize table of unlimited dimensions for ncid, for
 * use by is_unlim_dim() function.  If ncid is a subgroup of a netCDF
 * dataset, the table will still be initialized for the whole dataset
 * in which the subgroup resides. */
#ifdef USE_NETCDF4
static int
init_is_unlim(int ncid, int **is_unlim_p)
{
    int num_grps;	 /* total number of groups */
    int num_dims = 0;    /* total number of dimensions in all groups */
    int num_undims = 0;  /* total number of unlimited dimensions in all groups */
    int *grpids = NULL;	 /* temporary list of all grpids */
    int igrp;
    int grpid;

    /* if ncid is not root group, find its ancestor root group id */
    int status = nc_inq_grp_parent(ncid, &grpid);
    while(status == NC_NOERR && grpid != ncid) {
	ncid = grpid;
	status = nc_inq_grp_parent(ncid, &grpid);
    }
    if (status != NC_ENOGRP)
	return NC_EBADGRPID;
    /* Now ncid is root group.  Get total number of groups and their ids */
    NC_CHECK( nc_inq_grps_full(ncid, &num_grps, NULL) );
    grpids = emalloc((num_grps + 1) * sizeof(int));
    NC_CHECK( nc_inq_grps_full(ncid, &num_grps, grpids) );
#define DONT_INCLUDE_PARENTS 0
    /* Get all dimensions in groups and info about which ones are unlimited */
    for(igrp = 0; igrp < num_grps; igrp++) {
	int ndims;
	grpid = grpids[igrp];
	NC_CHECK( nc_inq_dimids(grpid, &ndims, NULL, DONT_INCLUDE_PARENTS) );
	num_dims += ndims;
    }
    *is_unlim_p = emalloc((num_dims + 1) * sizeof(int));
    for(igrp = 0; igrp < num_grps; igrp++) {
	int ndims, idim, *dimids, nundims;
	grpid = grpids[igrp];
	NC_CHECK( nc_inq_dimids(grpid, &ndims, NULL, DONT_INCLUDE_PARENTS) );
	dimids = emalloc((ndims + 1) * sizeof(int));
	NC_CHECK( nc_inq_dimids(grpid, &ndims, dimids, DONT_INCLUDE_PARENTS) );
	/* mark all dims in this group as fixed-size */
	for(idim = 0; idim < ndims; idim++) {
	    (*is_unlim_p)[dimids[idim]] = 0;
	}
	NC_CHECK( nc_inq_unlimdims(grpid, &nundims, dimids) );
	assert(nundims <= ndims);
	/* mark the subset of dims in this group that are unlimited */
	for(idim = 0; idim < nundims; idim++) {
	    (*is_unlim_p)[dimids[idim]] = 1;
	    num_undims++;
	}
	if(dimids)
	    free(dimids);
    }
    free(grpids);
    return NC_NOERR;
}
#endif	/*  USE_NETCDF4 */

/* TODO: make list of these arrays for multiple open datasets, such as
 * the idnode_t lists above.  For now, we just have one of these, for
 * the unique input dataset for this invocation of ncdump. */

#define UNLIM_NOT_INITIALIZED (-1)

/* Is dimid the dimension ID of an unlimited dimension?  */
bool_t
is_unlim_dim(int ncid, int dimid) {
    bool_t result;		/* 0 if fixed, 1 if unlimited size */
    static int for_ncid = UNLIM_NOT_INITIALIZED; /* ensure only ever called for one ncid */
#ifdef USE_NETCDF4
    static int *is_unlim = NULL; /* gets allocated by init_is_unlim() */
    if(for_ncid == UNLIM_NOT_INITIALIZED) {
      NC_CHECK(init_is_unlim(ncid, &is_unlim));
      for_ncid = ncid;
    }
    assert(is_unlim);
    result = is_unlim[dimid];	/* 0 if fixed, 1 if unlimited size */
#else
    static int unlimdimid;
    if(for_ncid == UNLIM_NOT_INITIALIZED) {
	NC_CHECK( nc_inq_unlimdim(ncid, &unlimdimid) );
	for_ncid = ncid;
    }
    result = (dimid == unlimdimid) ;
#endif 	 /* USE_NETCDF4 */
    return result;
}

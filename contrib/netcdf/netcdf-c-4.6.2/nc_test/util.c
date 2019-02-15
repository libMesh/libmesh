/*********************************************************************
 *   Copyright 1996, UCAR/Unidata
 *   See netcdf/COPYRIGHT file for copying and redistribution conditions.
 *   $Id: util.c 2792 2014-10-27 06:02:59Z wkliao $
 *********************************************************************/

#include <math.h> /* floor() */
#include "tests.h"

void
print_nok(int nok)
{
    if (verbose || nfails > 0)
        print("\n");
    print(" %d good comparisons. ", nok);
}


/* Is value within external type range? */
int
inRange(const double value, const nc_type xtype)
{
  switch (xtype) {
  case NC_CHAR:   return value >= X_CHAR_MIN   && value <= X_CHAR_MAX;
  case NC_BYTE:   return value >= X_BYTE_MIN   && value <= X_BYTE_MAX;
  case NC_SHORT:  return value >= X_SHORT_MIN  && value <= X_SHORT_MAX;
  case NC_INT:    return value >= X_INT_MIN    && value <= X_INT_MAX;
  case NC_FLOAT:  return value >= X_FLOAT_MIN  && value <= X_FLOAT_MAX;
  case NC_DOUBLE: return value >= X_DOUBLE_MIN && value <= X_DOUBLE_MAX;
  case NC_UBYTE:  return value >= 0            && value <= X_UCHAR_MAX;
  case NC_USHORT: return value >= 0            && value <= X_USHORT_MAX;
  case NC_UINT:   return value >= 0            && value <= X_UINT_MAX;
  case NC_INT64:  return value >= X_INT64_MIN  && value <= X_INT64_MAX;
  case NC_UINT64: return value >= 0            && value <= X_UINT64_MAX;
  default:  assert(0);
    return(0);
  }
}

static int
inRange_uchar(const int     cdf_format,
              const double  value,
              const nc_type xtype)
{
    /* check value of type xtype if within uchar range */

    if (cdf_format < NC_FORMAT_CDF5 && xtype == NC_BYTE) {
        /* netCDF specification make a special case for type conversion between
         * uchar and scahr: do not check for range error. See
         * http://www.unidata.ucar.edu/software/netcdf/docs/data_type.html#type_conversion
         */
        return(value >= 0 && value <= 255);
        /* this is to ensure value is within the range of uchar */
    }
    /* else */
    return inRange(value, xtype);
}

static int
inRange_float(const double value, const nc_type xtype)
{
    double min, max;

    switch (xtype) {
	case NC_CHAR:   min = X_CHAR_MIN;   max = X_CHAR_MAX;  break;
	case NC_BYTE:   min = X_BYTE_MIN;   max = X_BYTE_MAX;  break;
	case NC_SHORT:  min = X_SHORT_MIN;  max = X_SHORT_MAX; break;
	case NC_INT:    min = X_INT_MIN;    max = X_INT_MAX;   break;
	case NC_FLOAT:
		if(FLT_MAX < X_FLOAT_MAX) {
			min = (-FLT_MAX);
			max = FLT_MAX;
		} else {
			min = X_FLOAT_MIN;
			max = X_FLOAT_MAX;
		}
		break;
	case NC_DOUBLE:
		if(FLT_MAX < X_DOUBLE_MAX) {
			min = (-FLT_MAX);
			max = FLT_MAX;
		} else {
			min = X_DOUBLE_MIN;
			max = X_DOUBLE_MAX;
		}
		break;
        case NC_UBYTE:  min = 0;            max = X_UCHAR_MAX;  break;
        case NC_USHORT: min = 0;            max = X_USHORT_MAX; break;
        case NC_UINT:   min = 0;            max = X_UINT_MAX;   break;
        case NC_INT64:  min = X_INT64_MIN;  max = X_INT64_MAX;  break;
        case NC_UINT64: min = 0;            max = X_UINT64_MAX; break;
	default:  assert(0);
    }
    if(!( value >= min && value <= max)) {
#if 0	/* DEBUG */
	if(xtype == NC_FLOAT) {
	fprintf(stderr, "\n");
	fprintf(stderr, "min   % .17e\n", min);
	fprintf(stderr, "value % .17e\n", value);
	fprintf(stderr, "max   % .17e\n", max);
	}
#endif
	return 0;
    }
#if FLT_MANT_DIG != DBL_MANT_DIG
    /* else */
    {
	const float fvalue = (float)value;
	return fvalue >= min && fvalue <= max;
    }
#else
    return 1;
#endif
}

/* wrapper for inRange to handle special NC_BYTE/uchar adjustment */
/* this function checks whether "value" to be casted to type "itype" is
 * within the range of external "xtype".
 */
int
inRange3(const int    cdf_format,
         const double value,
         const nc_type xtype,
         const nct_itype itype)
{
    /* netCDF specification make a special case for type conversion between
     * uchar and NC_BYTE: do not check for range error. See
     * http://www.unidata.ucar.edu/software/netcdf/docs/data_type.html#type_conversion
     * The _uchar and _schar functions were introduced in netCDF-3 to eliminate
     * an ambiguity, and support both signed and unsigned byte data. In
     * netCDF-2, whether the external NC_BYTE type represented signed or
     * unsigned values was left up to the user. In netcdf-3, we treat NC_BYTE
     * as signed for the purposes of conversion to short, int, long, float, or
     * double. (Of course, no conversion takes place when the internal type is
     * signed char.) In the _uchar functions, we treat NC_BYTE as if it were
     * unsigned. Thus, no NC_ERANGE error can occur converting between NC_BYTE
     * and unsigned char.
     */
    switch (itype) {
        case NCT_UCHAR:
	    return inRange_uchar(cdf_format, value, xtype);
        case NCT_FLOAT:
	    return inRange_float(value, xtype);
        default:
	    break;
    }
    return inRange(value, xtype);
}


/*
 *  Does x == y, where one is internal and other external (netCDF)?
 *  Use tolerant comparison based on IEEE FLT_EPSILON or DBL_EPSILON.
 */
int
equal(const double x,
      const double y,
      nc_type xtype, 	/* external data type */
      nct_itype itype)
{
    const double flt_epsilon = 1.19209290E-07;
    const double dbl_epsilon = 2.2204460492503131E-16;
    double epsilon;

    epsilon = xtype == NC_FLOAT ||
              itype == NCT_FLOAT ? flt_epsilon : dbl_epsilon;

    if (xtype == NC_CHAR && itype == NCT_TEXT) {
        /* because in-memory data type char can be signed or unsigned,
         * type cast the value from external NC_CHAR before the comparison
         */
        char x2 = (char) x;
        char y2 = (char) y;
        return ABS(x2-y2) <= epsilon * MAX( ABS(x2), ABS(y2));
    }

    return ABS(x-y) <= epsilon * MAX( ABS(x), ABS(y));
}

/* this function is for the APIs without itype, i.e. xtype == itype */
int
equal2(const double x,
       const double y,
       nc_type xtype)    /* external data type */
{
    const double flt_epsilon = 1.19209290E-07;
    const double dbl_epsilon = 2.2204460492503131E-16;
    double epsilon;

    epsilon = xtype == NC_FLOAT ? flt_epsilon : dbl_epsilon;

    if (xtype == NC_CHAR) {
        /* because in-memory data type char can be signed or unsigned,
         * type cast the value from external NC_CHAR before the comparison
         */
        char x2 = (char) x;
        char y2 = (char) y;
        return ABS(x2-y2) <= epsilon * MAX( ABS(x2), ABS(y2));
    }

    return ABS(x-y) <= epsilon * MAX( ABS(x), ABS(y));
}

/* Test whether two int vectors are equal. If so return 1, else 0  */
int
int_vec_eq(const int *v1, const int *v2, const int n)
{
    int i;
    for (i= 0; i < n && v1[i] == v2[i]; i++)
	;
    return i == n;
}


/*
 *  Generate random integer from 0 to n-1
 *  Like throwing an n-sided dice marked 0, 1, 2, ..., n-1
 */
size_t roll( size_t n )
{
    size_t r;

    do
	/*
	 * Compute a pseudo-random value between 0.0 and 1.0, multiply
	 * it by n-1, and then find the nearest integer.
	 *
	 * We don't use RAND_MAX here because not all compilation
	 * environments define it (e.g. gcc(1) under SunOS 4.1.4).
	 */
	r = (size_t)(((rand() % 32768) / 32767.0) * (n - 1) + 0.5);
    while (r >= n);

    return r;
}


/*
 *      Convert number to mixed base
 *
 *      E.g. to convert 41 inches to yards, feet and inches:
 *      size_t base[] = {1, 3, 12};
 *      size_t result[3];
 *      status = toMixedBase(41, 3, base, result);
 *
 *      Author: Harvey Davies, Unidata/UCAR, Boulder, Colorado
 */
int
toMixedBase(
    size_t number,        /* number to be converted to mixed base */
    int length,
    const size_t base[],        /* dimensioned [length], base[0] ignored */
    size_t result[])      /* dimensioned [length] */
{
    int i;

    if (length > 0) {
	for (i = length - 1; i > 0; i--) {
	    if (base[i] == 0) return 1;
	    result[i] = number % base[i];
	    number = number / base[i];
	}
        result[0] = number;
    }
    return 0;
}

/*
 *      Convert number from mixed base
 *
 *      E.g. to convert 1 yard, 0 feet, 5 inches to inches:
 *      size_t number[] = {1, 0, 5};
 *      size_t base[] = {1, 3, 12};
 *      inches = fromMixedBase(3, number, base);
 *
 *      Author: Harvey Davies, Unidata/UCAR, Boulder, Colorado
 */
size_t
fromMixedBase(int    length,
              size_t number[],      /* dimensioned [length] */
              size_t base[])        /* dimensioned [length], base[0] ignored */
{
    size_t i;
    size_t result = 0;

    for (i = 1; i < length; i++) {
        result += number[i-1];
        result *= base[i];
    }
    if (length > 0)
        result += number[i-1];
    return result;
}


/* Convert any nc_type to double */
int nc2dbl ( const nc_type xtype, const void *p, double *result)
{
    if ( ! p ) return 2;
    if ( ! result ) return 3;
    switch (xtype) {
        case NC_CHAR:   *result = *((char *)          p); break;
        case NC_BYTE:   *result = *((signed char *)    p); break;
        case NC_UBYTE:  *result = *((unsigned char *)  p); break;
        case NC_SHORT:  *result = *((short *)          p); break;
        case NC_USHORT: *result = *((unsigned short *) p); break;
        case NC_INT:
#if INT_MAX >= X_INT_MAX
		*result = *((int *) p); break;
#else
		*result = *((long *) p); break;
#endif
        case NC_UINT:
#if UINT_MAX >= X_UINT_MAX
            *result = *((unsigned int *) p); break;
#else
            *result = *((unsigned long *) p); break;
#endif
        case NC_FLOAT:  *result = *((float *)              p); break;
        case NC_DOUBLE: *result = *((double *)             p); break;
        case NC_INT64:  *result = *((long long *)          p); break;
        case NC_UINT64: *result = *((unsigned long long *) p); break;
        default: return 1;
    }
    return 0;
}


/* Convert double to any nc_type */
int dbl2nc ( const double d, const nc_type xtype, void *p)
{
    double r;   /* rounded value */

    if (p == NULL) return 1;
    switch (xtype) {
        case NC_CHAR:
            r = floor(0.5+d);
            /* d is obtained from hash() which may be set to X_CHAR_MIN (0)
             * or X_CHAR_MAX (255). When in-memory data type char is signed
             * (i.e. ranged from -128 to 127), we should still allow a type
             * cast a unsigned value > 127 to a signed char without
             * reporting it as a range error.
             */
            if ( r < X_CHAR_MIN || r > X_CHAR_MAX ) return 2;
            *((signed char*) p) = (signed char)r;
            break;
        case NC_BYTE:
            r = floor(0.5+d);
            if ( r < schar_min  ||  r > schar_max )  return 2;
            *((signed char *) p) = (signed char)r;
            break;
        case NC_UBYTE:
            r = floor(0.5+d);
            if ( r < 0.0  ||  r > uchar_max )  return 2;
            *((unsigned char *) p) = (unsigned char)r;
            break;
        case NC_SHORT:
            r = floor(0.5+d);
            if ( r < short_min  ||  r > short_max )  return 2;
            *((short  *) p) = (short)r;
            break;
        case NC_USHORT:
            r = floor(0.5+d);
            if ( r < 0.0  ||  r > ushort_max )  return 2;
            *((unsigned short *) p) = (unsigned short)r;
            break;
        case NC_INT:
            r = floor(0.5+d);
            if ( r < long_min  ||  r > long_max )  return 2;
#if INT_MAX >= X_INT_MAX
            *((int   *) p) = (int)r;
#else
            *((long   *) p) = (long)r;
#endif
            break;
        case NC_UINT:
            r = floor(0.5+d);
            if ( r < 0.0  ||  r > uint_max )  return 2;
#if UINT_MAX >= X_UINT_MAX
            *((unsigned int  *) p) = (unsigned int)r;
#else
            *((unsigned long *) p) = (unsigned long)r;
#endif
            break;
        case NC_FLOAT:
            if ( fabs(d) > float_max )  return 2;
            *((float  *) p) = (float)d;
            break;
        case NC_DOUBLE:
            *((double *) p) = (double)d;
            break;
        case NC_INT64:
            r = floor(0.5+d);
            if ( r < int64_min  ||  r > int64_max )  return 2;
            *((long long *) p) = (long long)r;
            break;
        case NC_UINT64:
            r = floor(0.5+d);
            if ( r < 0.0  ||  r > uint64_max )  return 2;
            *((unsigned long long *) p) = (unsigned long long)r;
            break;
        default:
            return 1;
    }
    return 0;
}

#define FUZZ (1.19209290E-07)

#ifdef USE_EXTREME_NUMBERS
/* Generate data values as function of type, rank (-1 for attribute), index */
double
hash( const nc_type xtype, const int rank, const size_t *index )
{
    double base;
    double result;
    int  d;       /* index of dimension */

	/* If vector then elements 0 & 1 are min & max. Elements 2 & 3 are */
	/* just < min & > max (except for NC_CHAR & NC_DOUBLE) */
    if (abs(rank) == 1 && index[0] <= 3) {
	switch (index[0]) {
	    case 0:
		switch (xtype) {
		    case NC_CHAR:   return X_CHAR_MIN;
		    case NC_BYTE:   return X_BYTE_MIN;
		    case NC_SHORT:  return X_SHORT_MIN;
		    case NC_INT:    return X_INT_MIN;
		    case NC_FLOAT:  return X_FLOAT_MIN;
		    case NC_DOUBLE: return X_DOUBLE_MIN;
                    case NC_UBYTE:  return 0;
                    case NC_USHORT: return 0;
                    case NC_UINT:   return 0;
                    case NC_INT64:  return X_INT_MIN - 128.0; /* slight smaller
                                                                 than INT_MIN */
                    case NC_UINT64: return 0;
		    default:  assert(0);
		}
	    case 1:
		switch (xtype) {
		    case NC_CHAR:   return X_CHAR_MAX;
		    case NC_BYTE:   return X_BYTE_MAX;
		    case NC_SHORT:  return X_SHORT_MAX;
		    case NC_INT:    return X_INT_MAX;
		    case NC_FLOAT:  return X_FLOAT_MAX;
		    case NC_DOUBLE: return X_DOUBLE_MAX;
                    case NC_UBYTE:  return X_UCHAR_MAX;
                    case NC_USHORT: return X_USHORT_MAX;
                    case NC_UINT:   return X_UINT_MAX;
                    case NC_INT64:  return X_INT_MAX + 128.0;
                                    /* slightly bigger than INT_MAX */
                    case NC_UINT64: return X_UINT_MAX + 128.0;
                                    /* slightly bigger than UINT_MAX */
		    default:  assert(0);
		}
	    case 2:
		switch (xtype) {
		    case NC_CHAR:   return 'A';
		    case NC_BYTE:   return X_BYTE_MIN-1.0;
		    case NC_SHORT:  return X_SHORT_MIN-1.0;
		    case NC_INT:    return X_INT_MIN-1.0;
		    case NC_FLOAT:  return X_FLOAT_MIN * (1.0 + FUZZ);
		    case NC_DOUBLE: return -1.0;
                    case NC_UBYTE:  return -1.0;
                    case NC_USHORT: return -1.0;
                    case NC_UINT:   return -1.0;
                    case NC_INT64:  return -1.0;  /* skip test */
                    case NC_UINT64: return -1.0;
		    default:  assert(0);
		}
	    case 3:
		switch (xtype) {
		    case NC_CHAR:   return 'Z';
		    case NC_BYTE:   return X_BYTE_MAX+1.0;
		    case NC_SHORT:  return X_SHORT_MAX+1.0;
		    case NC_INT:    return X_INT_MAX+1.0;
		    case NC_FLOAT:  return X_FLOAT_MAX * (1.0 + FUZZ);
		    case NC_DOUBLE: return 1.0;
                    case NC_UBYTE:  return X_UCHAR_MAX +1.0;
                    case NC_USHORT: return X_USHORT_MAX+1.0;
                    case NC_UINT:   return X_UINT_MAX  +1.0;
                    case NC_INT64:  return 1.0;    /* skip test */
                    case NC_UINT64: return 1.0;    /* skip test */
		    default:  assert(0);
		}
	}
    } else {
	switch (xtype) {
	    case NC_CHAR:   base =   2; break;
	    case NC_BYTE:   base =  -2; break;
	    case NC_SHORT:  base =  -5; break;
	    case NC_INT:    base = -20; break;
	    case NC_FLOAT:  base =  -9; break;
	    case NC_DOUBLE: base = -10; break;

            /* not sure what right values are */
            case NC_UBYTE:   base =   2;  break;
            case NC_USHORT:  base =   5;  break;
            case NC_UINT:    base =  20;  break;
            case NC_INT64:   base = -20;  break;
            case NC_UINT64:  base =  20;  break;
	    default:  assert(0);
	}
	result = rank < 0 ? base * 7 : base * (rank + 1);
	for (d = 0; d < abs(rank); d++)
	    result = base * (result + index[d]);
    }
    return result;
}
#else /* USE_EXTREME_NUMBERS */
#define SANE_SHORT 3333
#define SANE_INT 2222
#define SANE_FLOAT 300.0
#define SANE_DOUBLE 1000.0

/* Generate data values as function of type, rank (-1 for attribute), index */
double
hash( const nc_type xtype, const int rank, const size_t *index )
{
    double base;
    double result;
    int  d;       /* index of dimension */

	/* If vector then elements 0 & 1 are min & max. Elements 2 & 3 are */
	/* just < min & > max (except for NC_CHAR & NC_DOUBLE) */
    if (abs(rank) == 1 && index[0] <= 3) {
	switch (index[0]) {
	    case 0:
		switch (xtype) {
		    case NC_CHAR:   return X_CHAR_MIN;
		    case NC_BYTE:   return X_BYTE_MIN;
		    case NC_SHORT:  return SANE_SHORT;
		    case NC_INT:   return SANE_INT;
		    case NC_FLOAT:  return SANE_FLOAT;
		    case NC_DOUBLE: return SANE_DOUBLE;
                    case NC_UBYTE:  return 0;
                    case NC_USHORT: return 0;
                    case NC_UINT:   return 0;
                    case NC_INT64:  return X_INT_MIN - 128.0; /* slight smaller
                                                                 than INT_MIN */
                    case NC_UINT64: return 0;
		    default:  assert(0);
		}
	    case 1:
		switch (xtype) {
		    case NC_CHAR:   return X_CHAR_MAX;
		    case NC_BYTE:   return X_BYTE_MAX;
		    case NC_SHORT:  return SANE_SHORT;
		    case NC_INT:   return SANE_INT;
		    case NC_FLOAT:  return SANE_FLOAT;
		    case NC_DOUBLE: return SANE_DOUBLE;
                    case NC_UBYTE:  return X_UCHAR_MAX;
                    case NC_USHORT: return X_USHORT_MAX;
                    case NC_UINT:   return X_UINT_MAX;
                    case NC_INT64:  return X_INT_MAX + 128.0;
                                    /* slightly bigger than INT_MAX */
                    case NC_UINT64: return X_UINT_MAX + 128.0;
                                    /* slightly bigger than UINT_MAX */
		    default:  assert(0);
		}
	    case 2:
		switch (xtype) {
		    case NC_CHAR:   return 'A';
		    case NC_BYTE:   return X_BYTE_MIN-1.0;
		    case NC_SHORT:  return SANE_SHORT-1.0;
		    case NC_INT:   return SANE_INT-1.0;
		    case NC_FLOAT:  return SANE_FLOAT * (1.0 + FUZZ);
		    case NC_DOUBLE: return -1.0;
                    case NC_UBYTE:  return -1.0;
                    case NC_USHORT: return -1.0;
                    case NC_UINT:   return -1.0;
                    case NC_INT64:  return -1.0;  /* skip test */
                    case NC_UINT64: return -1.0;
		    default:  assert(0);
		}
	    case 3:
		switch (xtype) {
		    case NC_CHAR:   return 'Z';
		    case NC_BYTE:   return X_BYTE_MAX+1.0;
		    case NC_SHORT:  return SANE_SHORT+1.0;
		    case NC_INT:   return SANE_INT+1.0;
		    case NC_FLOAT:  return SANE_FLOAT * (1.0 + FUZZ);
		    case NC_DOUBLE: return 1.0;
                    case NC_UBYTE:  return X_UCHAR_MAX +1.0;
                    case NC_USHORT: return X_USHORT_MAX+1.0;
                    case NC_UINT:   return X_UINT_MAX  +1.0;
                    case NC_INT64:  return 1.0;    /* skip test */
                    case NC_UINT64: return 1.0;    /* skip test */
		    default:  assert(0);
		}
	}
    } else {
	switch (xtype) {
	    case NC_CHAR: base = 2; break;
	    case NC_BYTE: base = -2; break;
	    case NC_SHORT: base = -5; break;
	    case NC_INT: base = -20; break;
	    case NC_FLOAT: base = -9; break;
	    case NC_DOUBLE: base = -10; break;

            /* not sure what right values are */
            case NC_UBYTE:   base =   2;  break;
            case NC_USHORT:  base =   5;  break;
            case NC_UINT:    base =  20;  break;
            case NC_INT64:   base = -20;  break;
            case NC_UINT64:  base =  20;  break;
	    default:  assert(0);
	}
	result = rank < 0 ? base * 7 : base * (rank + 1);
	for (d = 0; d < abs(rank); d++)
	    result = base * (result + index[d]);
    }
    return result;
}
#endif
/* wrapper for hash to handle special NC_BYTE/uchar adjustment */
double
hash4(const int        cdf_format,
      const nc_type    xtype,
      const int        rank,
      const size_t    *index,
      const nct_itype  itype)
{
    double result;

    result = hash( xtype, rank, index );

    /* netCDF specification make a special case for type conversion between
     * uchar and NC_BYTE: do not check for range error. See
     * http://www.unidata.ucar.edu/software/netcdf/docs/data_type.html#type_conversion
     * The _uchar and _schar functions were introduced in netCDF-3 to eliminate
     * an ambiguity, and support both signed and unsigned byte data. In
     * netCDF-2, whether the external NC_BYTE type represented signed or
     * unsigned values was left up to the user. In netcdf-3, we treat NC_BYTE
     * as signed for the purposes of conversion to short, int, long, float, or
     * double. (Of course, no conversion takes place when the internal type is
     * signed char.) In the _uchar functions, we treat NC_BYTE as if it were
     * unsigned. Thus, no NC_ERANGE error can occur converting between NC_BYTE
     * and unsigned char.
     */
    if (cdf_format < NC_FORMAT_CDF5 &&
        itype == NCT_UCHAR && xtype == NC_BYTE &&
        result >= -128 && result < 0)
	result += 256;

    return result;
}

static nc_type
char2type(char letter) {
    switch (letter) {
        case 'c': return NC_CHAR;
        case 'b': return NC_BYTE;
        case 's': return NC_SHORT;
        case 'i': return NC_INT;
        case 'f': return NC_FLOAT;
        case 'd': return NC_DOUBLE;
        case 'y': return NC_UBYTE;
        case 't': return NC_USHORT;
        case 'u': return NC_UINT;
        case 'x': return NC_INT64;
        case 'z': return NC_UINT64;
        default:  assert(0);
    }
    return NC_CHAR;  /* Just to keep compiler happy */
}


static void
init_dims(const char *digit)
{
	int dimid;			/* index of dimension */
	for (dimid = 0; dimid < NDIMS; dimid++)
	{
		dim_len[dimid] = dimid == 0 ? NRECS : dimid;
		dim_name[dimid][0] = 'D';
		dim_name[dimid][1] = digit[dimid];
		dim_name[dimid][2] = '\0';
	}
}

static void
init_gatts(const char *type_letter)
{
	int attid;
	for (attid = 0; attid < numGatts; attid++)
	{
		gatt_name[attid][0] = 'G';
		gatt_name[attid][1] = type_letter[attid];
		gatt_name[attid][2] = '\0';
		gatt_len[attid] = 1 + attid;
		gatt_type[attid] = char2type (type_letter[attid]);
	}
}

static size_t
product(int nn, const size_t *sp)
{
	size_t result = 1;
	while(nn-- > 0)
		result *= *sp++;
	return result;
}

/*
   define global variables:
   dim_name, dim_len,
   var_name, var_type, var_rank, var_shape, var_natts, var_dimid, var_nels
   att_name, gatt_name, att_type, gatt_type, att_len, gatt_len
 */
void
init_gvars (void)
{
	const size_t max_dim_len[MAX_RANK] = {
		MAX_DIM_LEN +1,
		MAX_DIM_LEN,
		MAX_DIM_LEN
	};
    const char type_letter[] = "cbsifdytuxz";
    /* c:char, b:byte, s:short, i:int, f:float, d:double, y:ubyte, t:ushort,
     * u:uint, x:int64, z:uint64
     */
	const char digit[] = "r123456789";

	int rank;
	int vn;			/* var number */
	int xtype;		/* index of type */
	int an;			/* attribute number */

	assert(sizeof(max_dim_len)/sizeof(max_dim_len[0]) >= MAX_RANK);

	init_dims(digit);

	for (vn=0; vn<numVars; vn++)
	    memset(var_name[vn], 0, 2+MAX_RANK);

	for (rank = 0, vn = 0, xtype = 0, an = 0;  rank <= MAX_RANK; rank++)
	{
			/* number variables of a type and rank */
		const size_t nvars = product(rank, max_dim_len);
		size_t jj;

		for (jj = 0; jj < nvars; jj++)
		{
				/* number types of this shape */
			const int ntypes = rank < 2 ? numTypes : 1;

			int tc;
			for (tc = 0; tc < ntypes;
			     tc++, vn++, xtype = (xtype + 1) % numTypes)
			{
				size_t tmp[MAX_RANK];

				var_name[vn][0] = type_letter[xtype];
				var_type[vn] = char2type (type_letter[xtype]);
				var_rank[vn] = rank;
				var_natts[vn] = rank == 0 ? vn % (MAX_NATTS + 1) : 0;
				{
					int ac;
					for (ac = 0; ac < var_natts[vn]; ac++, an++)
					{
						att_name[vn][ac][0] = type_letter[an % numTypes];
						att_name[vn][ac][1] = '\0';
						att_len[vn][ac] = an;
						att_type[vn][ac] = char2type (type_letter[an % numTypes]);
					}
				} /* ac block */
#ifndef NDEBUG
				assert(toMixedBase (jj, rank, max_dim_len, tmp) == 0);
#else
				(void) toMixedBase (jj, rank, max_dim_len, tmp);
#endif
				{
					int dn; /* dimension number */
					for (dn = 0; dn < rank; dn++)
						var_dimid[vn][dn] = (int)tmp[dn];
					for (dn = 0, var_nels[vn] = 1; dn < rank; dn++)
					{
						var_dimid[vn][dn] += dn > 0;
						assert (var_dimid[vn][dn] <= 9);
						var_name[vn][dn + 1] = digit[var_dimid[vn][dn]];
						var_shape[vn][dn] = var_dimid[vn][dn] ?
							var_dimid[vn][dn] : NRECS;
						var_nels[vn] *= var_shape[vn][dn];
					}
				} /* dn block */
			}
		}
	}

	init_gatts(type_letter);
}


/* define dims defined by global variables */
void
def_dims(int ncid)
{
    int  err;             /* status */
    int  i;
    int  dimid;		/* dimension id */

    for (i = 0; i < NDIMS; i++) {
	err = nc_def_dim(ncid, dim_name[i], i==0 ? NC_UNLIMITED : dim_len[i],
	    &dimid);
	IF (err) error("nc_def_dim: %s", nc_strerror(err));
    }
}


/* define vars defined by global variables */
void
def_vars(int ncid)
{
    int  err;             /* status */
    int  i;
    int var_id;

    for (i = 0; i < numVars; i++) {
	err = nc_def_var(ncid, var_name[i], var_type[i], var_rank[i],
	    var_dimid[i], &var_id);
	IF (err) error("nc_def_var: %s", nc_strerror(err));
    }
}


/* put attributes defined by global variables */
void
put_atts(int ncid)
{
    int  err;             /* status */
    int  i;
    size_t  k;
    int  j;		/* index of attribute */
    int  allInRange;
    double att[MAX_NELS];
    char catt[MAX_NELS];

    for (i = -1; i < numVars; i++) {
	for (j = 0; j < NATTS(i); j++) {
	    if (ATT_TYPE(i,j) == NC_CHAR) {
		for (k = 0; k < ATT_LEN(i,j); k++) {
                    catt[k] = (char) hash(ATT_TYPE(i,j), -1, &k);
		}
		err = nc_put_att_text(ncid, i, ATT_NAME(i,j),
		    ATT_LEN(i,j), catt);
		IF (err)
		    error("nc_put_att_text: %s", nc_strerror(err));
	    } else {
		for (allInRange = 1, k = 0; k < ATT_LEN(i,j); k++) {
		    att[k] = hash(ATT_TYPE(i,j), -1, &k);
		    allInRange = allInRange && inRange(att[k], ATT_TYPE(i,j));
		}
		err = nc_put_att_double(ncid, i, ATT_NAME(i,j),
		    ATT_TYPE(i,j), ATT_LEN(i,j), att);
                if (allInRange) {
                    IF (err)
                        error("nc_put_att_double: %s", nc_strerror(err));
                } else {
                    IF (err != NC_ERANGE)
			error("type-conversion range error: status = %d", err);
                }
	    }
        }
    }
}

/* put variables defined by global variables */
void
put_vars(int ncid)
{
    size_t start[MAX_RANK];
    size_t index[MAX_RANK];
    int  err;             /* status */
    int  i;
    size_t  j;
    double value[MAX_NELS];
    char text[MAX_NELS];
    int  allInRange;

    for (j = 0; j < MAX_RANK; j++)
	start[j] = 0;
    for (i = 0; i < numVars; i++) {
	for (allInRange = 1, j = 0; j < var_nels[i]; j++) {
	    err = toMixedBase(j, var_rank[i], var_shape[i], index);
	    IF (err) error("toMixedBase");
	    if (var_name[i][0] == 'c') {
		text[j] = (char) hash(var_type[i], var_rank[i], index);
	    } else {
		value[j]  = hash(var_type[i], var_rank[i], index);
		allInRange = allInRange && inRange(value[j], var_type[i]);
	    }
	}
	if (var_name[i][0] == 'c') {
	    err = nc_put_vara_text(ncid, i, start, var_shape[i], text);
	    IF (err)
		error("nc_put_vara_text: %s", nc_strerror(err));
	} else {
	    err = nc_put_vara_double(ncid, i, start, var_shape[i], value);
	    if (allInRange) {
		IF (err)
		    error("nc_put_vara_double: %s", nc_strerror(err));
	    } else {
		IF (err != NC_ERANGE)
		    error("type-conversion range error: status = %d", err);
	    }
	}
    }
}


/* Create & write all of specified file using global variables */
void
write_file(char *filename)
{
    int  ncid; /* netCDF id */
    int  err;  /* status */
    err = file_create(filename, NC_CLOBBER, &ncid);
    IF (err)
	error("nc_create: %s", nc_strerror(err));

    def_dims(ncid);
    def_vars(ncid);
    put_atts(ncid);
    err = nc_enddef(ncid);
    IF (err)
	error("nc_enddef: %s", nc_strerror(err));

#ifdef USE_PNETCDF
    err = nc_var_par_access(ncid, NC_GLOBAL, NC_COLLECTIVE);
    IF (err) error("nc_var_par_access: %s", nc_strerror(err));
#endif

    put_vars(ncid);

    err = nc_close (ncid);
    IF (err)
	error("nc_close: %s", nc_strerror(err));
}


/*
 * check dimensions of specified file have expected name & length
 */
void
check_dims(int  ncid)
{
    char name[NC_MAX_NAME];
    size_t length;
    int  i;
    int  err;           /* status */

    for (i = 0; i < NDIMS; i++) {
	err = nc_inq_dim(ncid, i, name, &length);
	IF (err)
	    error("nc_inq_dim: %s", nc_strerror(err));
	IF (strcmp(name, dim_name[i]) != 0)
	    error("Unexpected name of dimension %d: '%s', expected: '%s'", i, name, dim_name[i]);
	IF (length != dim_len[i])
	    error("Unexpected length %d of dimension %d, expected %zu", length, i, dim_len[i]);
    }
}


/*
 * check variables of specified file have expected name, type, shape & values
 */
void
check_vars(int  ncid)
{
    size_t index[MAX_RANK];
    char  text, name[NC_MAX_NAME];
    int  i, err;		/* status */
    size_t  j;
    int nok = 0;      /* count of valid comparisons */
    int isChar, ndims, dimids[MAX_RANK];
    double value, expect;
    nc_type xtype;
    size_t length;

    for (i = 0; i < numVars; i++) {
        isChar = var_type[i] == NC_CHAR;
	err = nc_inq_var(ncid, i, name, &xtype, &ndims, dimids, NULL);
	IF (err)
	    error("nc_inq_var: %s", nc_strerror(err));
	IF (strcmp(name, var_name[i]) != 0)
	    error("Unexpected var_name");
	IF (xtype != var_type[i])
	    error("Unexpected type");
	IF (ndims != var_rank[i])
	    error("Unexpected rank");
	for (j = 0; j < ndims; j++) {
	    err = nc_inq_dim(ncid, dimids[j], 0, &length);
	    IF (err)
		error("nc_inq_dim: %s", nc_strerror(err));
	    IF (length != var_shape[i][j])
		error("Unexpected shape");
	}
	for (j = 0; j < var_nels[i]; j++) {
	    err = toMixedBase(j, var_rank[i], var_shape[i], index);
	    IF (err)
		error("error in toMixedBase 2");
	    expect = hash( var_type[i], var_rank[i], index );
	    if (isChar) {
          	err = nc_get_var1_text(ncid, i, index, &text);
            IF (err)
		    error("nc_get_var1_text: %s", nc_strerror(err));
            IF (text != (char)expect) {
              error("Var %s [%lu] value read %hhd not that expected %g ",
                  var_name[i], j, text, expect);
		    print_n_size_t(var_rank[i], index);
		} else {
		    nok++;
		}
	    } else {
		err = nc_get_var1_double(ncid, i, index, &value);
		if (inRange(expect,var_type[i])) {
		    IF (err) {
			error("nc_get_var1_double: %s", nc_strerror(err));
		    } else {
			IF (!equal(value,expect,var_type[i], NCT_DOUBLE)) {
			    error("Var %s [%lu] value read %g  not that expected %g ",
				  var_name[i], j, value, expect);
			    print_n_size_t(var_rank[i], index);
			} else {
			    nok++;
			}
		    }
		}
	    }
	}
    }
    print_nok(nok);
}


/*
 * check attributes of specified file have expected name, type, length & values
 */
void
check_atts(int  ncid)
{
    int  err;		/* status */
    int  i;
    int  j;
    size_t  k;
    nc_type xtype;
    char name[NC_MAX_NAME];
    size_t length;
    char text[MAX_NELS];
    double value[MAX_NELS];
    double expect;
    int nok = 0;      /* count of valid comparisons */

    for (i = -1; i < numVars; i++) {
	for (j = 0; j < NATTS(i); j++) {
            err = nc_inq_attname(ncid, i, j, name);
            IF (err)
                error("nc_inq_attname: %s", nc_strerror(err));
            IF (strcmp(name, ATT_NAME(i,j)) != 0)
                error("nc_inq_attname: unexpected name");
	    err = nc_inq_att(ncid, i, name, &xtype, &length);
	    IF (err)
		error("nc_inq_att: %s", nc_strerror(err));
	    IF (xtype != ATT_TYPE(i,j))
		error("nc_inq_att: unexpected type");
	    IF (length != ATT_LEN(i,j))
		error("nc_inq_att: unexpected length");
	    if (xtype == NC_CHAR) {
		err = nc_get_att_text(ncid, i, name, text);
		IF (err)
		    error("nc_get_att_text: %s", nc_strerror(err));
		for (k = 0; k < ATT_LEN(i,j); k++) {
		    expect = hash(xtype, -1, &k);
		    IF (text[k] != (char)expect) {
			error("nc_get_att_text: unexpected value");
            	    } else {
              		nok++;
            	    }
		}
	    } else {
		err = nc_get_att_double(ncid, i, name, value);
		for (k = 0; k < ATT_LEN(i,j); k++) {
		    expect = hash(xtype, -1, &k);
		    if (inRange(expect,ATT_TYPE(i,j))) {
			IF (err)
			    error("nc_get_att_double: %s", nc_strerror(err));
			IF (!equal(value[k], expect, ATT_TYPE(i,j), NCT_DOUBLE)) {
			    error("Att value read not that expected");
			} else {
			    nok++;
			}
		    }
		}
	    }
	}
    }
    print_nok(nok);
}


/* Check file (dims, vars, atts) corresponds to global variables */
void
check_file(char *filename)
{
    int  ncid;		/* netCDF id */
    int  err;		/* status */

    err = file_open(filename, NC_NOWRITE, &ncid);
    IF (err) {
        error("nc_open: %s", nc_strerror(err));
    } else {
	check_dims(ncid);
	check_vars(ncid);
	check_atts(ncid);
	err = nc_close (ncid);
	IF (err)
	    error("nc_close: %s", nc_strerror(err));
    }
}

/* TODO: Maybe this function belongs in the netcdf library. */
const char *
s_nc_type(nc_type xtype)
{
	switch((int)xtype){
        case NC_CHAR:   return "NC_CHAR";
        case NC_BYTE:   return "NC_BYTE";
        case NC_UBYTE:  return "NC_UBYTE";
        case NC_SHORT:  return "NC_SHORT";
        case NC_USHORT: return "NC_USHORT";
        case NC_INT:    return "NC_INT";
        case NC_UINT:   return "NC_UINT";
        case NC_FLOAT:  return "NC_FLOAT";
        case NC_DOUBLE: return "NC_DOUBLE";
        case NC_INT64:  return "NC_INT64";
        case NC_UINT64: return "NC_UINT64";
	}
	return "";
}


int file_create(const char *filename, int cmode, int *ncid)
{
    int err;

#ifdef USE_PNETCDF
    /* get the default file format */
    int default_format;
    nc_set_default_format(NC_FORMAT_CLASSIC, &default_format);
    /* set it back to the default */
    nc_set_default_format(default_format, NULL);

    if (default_format == NC_FORMAT_CLASSIC ||
        default_format == NC_FORMAT_64BIT_OFFSET ||
        default_format == NC_FORMAT_64BIT_DATA)
        err = nc_create_par(filename, cmode, MPI_COMM_WORLD, MPI_INFO_NULL, ncid);
    else
#endif
        err = nc_create(filename, cmode, ncid);

    return err;
}

int file__create(const char *filename,
                 int         cmode,
                 size_t      initialsz,
                 size_t     *bufrsizehintp,
                 int        *ncid)
{
    int err;

#ifdef USE_PNETCDF
    /* get the default file format */
    int default_format;
    err = nc_set_default_format(NC_FORMAT_CLASSIC, &default_format);
    /* set it back to the default */
    err = nc_set_default_format(default_format, NULL);

    if (default_format == NC_FORMAT_CLASSIC ||
        default_format == NC_FORMAT_64BIT_OFFSET ||
        default_format == NC_FORMAT_64BIT_DATA)
        err = nc_create_par(filename, cmode, MPI_COMM_WORLD, MPI_INFO_NULL, ncid);
    else
#endif
        err = nc__create(filename, cmode, initialsz, bufrsizehintp, ncid);

    return err;
}

int file_open(const char *filename, int omode, int *ncid)
{
    int err;

#ifdef USE_PNETCDF
    /* get the default file format */
    int default_format;
    err = nc_set_default_format(NC_FORMAT_CLASSIC, &default_format);
    /* set it back to the default */
    err = nc_set_default_format(default_format, NULL);

    if (default_format == NC_FORMAT_CLASSIC ||
        default_format == NC_FORMAT_64BIT_OFFSET ||
        default_format == NC_FORMAT_64BIT_DATA)
        err = nc_open_par(filename, omode, MPI_COMM_WORLD, MPI_INFO_NULL, ncid);
    else
#endif
        err = nc_open(filename, omode, ncid);

    return err;
}


#ifdef USE_PNETCDF
#include <pnetcdf.h>  /* to include PnetCDF error codes */
#endif

char* nc_err_code_name(int err)
{
    static char unknown_str[32];

    if (err > 0) { /* system error */
        const char *cp = (const char *) strerror(err);
        if (cp == NULL)
            sprintf(unknown_str,"Unknown error code %d",err);
        else
            sprintf(unknown_str,"Error code %d (%s)",err,cp);
        return unknown_str;
    }

    switch (err) {
        case (NC_NOERR):			return "NC_NOERR";
        case (NC_EBADID):			return "NC_EBADID";
        case (NC_ENFILE):			return "NC_ENFILE";
        case (NC_EEXIST):			return "NC_EEXIST";
        case (NC_EINVAL):			return "NC_EINVAL";
        case (NC_EPERM):			return "NC_EPERM";
        case (NC_ENOTINDEFINE):			return "NC_ENOTINDEFINE";
        case (NC_EINDEFINE):			return "NC_EINDEFINE";
        case (NC_EINVALCOORDS):			return "NC_EINVALCOORDS";
        case (NC_EMAXDIMS):			return "NC_EMAXDIMS"; /* not enforced after 4.5.0 */
        case (NC_ENAMEINUSE):			return "NC_ENAMEINUSE";
        case (NC_ENOTATT):			return "NC_ENOTATT";
        case (NC_EMAXATTS):			return "NC_EMAXATTS"; /* not enforced after 4.5.0 */
        case (NC_EBADTYPE):			return "NC_EBADTYPE";
        case (NC_EBADDIM):			return "NC_EBADDIM";
        case (NC_EUNLIMPOS):			return "NC_EUNLIMPOS";
        case (NC_EMAXVARS):			return "NC_EMAXVARS"; /* not enforced after 4.5.0 */
        case (NC_ENOTVAR):			return "NC_ENOTVAR";
        case (NC_EGLOBAL):			return "NC_EGLOBAL";
        case (NC_ENOTNC):			return "NC_ENOTNC";
        case (NC_ESTS):				return "NC_ESTS";
        case (NC_EMAXNAME):			return "NC_EMAXNAME";
        case (NC_EUNLIMIT):			return "NC_EUNLIMIT";
        case (NC_ENORECVARS):			return "NC_ENORECVARS";
        case (NC_ECHAR):			return "NC_ECHAR";
        case (NC_EEDGE):			return "NC_EEDGE";
        case (NC_ESTRIDE):			return "NC_ESTRIDE";
        case (NC_EBADNAME):			return "NC_EBADNAME";
        case (NC_ERANGE):			return "NC_ERANGE";
        case (NC_ENOMEM):			return "NC_ENOMEM";
        case (NC_EVARSIZE):			return "NC_EVARSIZE";
        case (NC_EDIMSIZE):			return "NC_EDIMSIZE";
        case (NC_ETRUNC):			return "NC_ETRUNC";
        case (NC_EAXISTYPE):			return "NC_EAXISTYPE";
        case (NC_EDAP):				return "NC_EDAP";
        case (NC_ECURL):			return "NC_ECURL";
        case (NC_EIO):				return "NC_EIO";
        case (NC_ENODATA):			return "NC_ENODATA";
        case (NC_EDAPSVC):			return "NC_EDAPSVC";
        case (NC_EDAS):				return "NC_EDAS";
        case (NC_EDDS):				return "NC_EDDS";
        case (NC_EDATADDS):			return "NC_EDATADDS";
        case (NC_EDAPURL):			return "NC_EDAPURL";
        case (NC_EDAPCONSTRAINT):		return "NC_EDAPCONSTRAINT";
        case (NC_ETRANSLATION):			return "NC_ETRANSLATION";
        case (NC_EACCESS):			return "NC_EACCESS";
        case (NC_EAUTH):			return "NC_EAUTH";
        case (NC_ENOTFOUND):			return "NC_ENOTFOUND";
        case (NC_ECANTREMOVE):			return "NC_ECANTREMOVE";
        case (NC_EINTERNAL):			return "NC_EINTERNAL";
        case (NC_EPNETCDF):			return "NC_EPNETCDF";
        case (NC_EHDFERR):			return "NC_EHDFERR";
        case (NC_ECANTREAD):			return "NC_ECANTREAD";
        case (NC_ECANTWRITE):			return "NC_ECANTWRITE";
        case (NC_ECANTCREATE):			return "NC_ECANTCREATE";
        case (NC_EFILEMETA):			return "NC_EFILEMETA";
        case (NC_EDIMMETA):			return "NC_EDIMMETA";
        case (NC_EATTMETA):			return "NC_EATTMETA";
        case (NC_EVARMETA):			return "NC_EVARMETA";
        case (NC_ENOCOMPOUND):			return "NC_ENOCOMPOUND";
        case (NC_EATTEXISTS):			return "NC_EATTEXISTS";
        case (NC_ENOTNC4):			return "NC_ENOTNC4";
        case (NC_ESTRICTNC3):			return "NC_ESTRICTNC3";
        case (NC_ENOTNC3):			return "NC_ENOTNC3";
        case (NC_ENOPAR):			return "NC_ENOPAR";
        case (NC_EPARINIT):			return "NC_EPARINIT";
        case (NC_EBADGRPID):			return "NC_EBADGRPID";
        case (NC_EBADTYPID):			return "NC_EBADTYPID";
        case (NC_ETYPDEFINED):			return "NC_ETYPDEFINED";
        case (NC_EBADFIELD):			return "NC_EBADFIELD";
        case (NC_EBADCLASS):			return "NC_EBADCLASS";
        case (NC_EMAPTYPE):			return "NC_EMAPTYPE";
        case (NC_ELATEFILL):			return "NC_ELATEFILL";
        case (NC_ELATEDEF):			return "NC_ELATEDEF";
        case (NC_EDIMSCALE):			return "NC_EDIMSCALE";
        case (NC_ENOGRP):			return "NC_ENOGRP";
        case (NC_ESTORAGE):			return "NC_ESTORAGE";
        case (NC_EBADCHUNK):			return "NC_EBADCHUNK";
        case (NC_ENOTBUILT):			return "NC_ENOTBUILT";
        case (NC_EDISKLESS):			return "NC_EDISKLESS";
        case (NC_ECANTEXTEND):			return "NC_ECANTEXTEND";
        case (NC_EMPI):				return "NC_EMPI";
        case (NC_ENULLPAD):			return "NC_NULLPAD";
        case (NC_EINMEMORY):			return "NC_EINMEMORY";
        // case (NC_EURL):				return "NC_EURL";
        // case (NC_ECONSTRAINT):			return "NC_ECONSTRAINT";
#ifdef USE_PNETCDF
        case (NC_ESMALL):			return "NC_ESMALL";
        case (NC_ENOTINDEP):			return "NC_ENOTINDEP";
        case (NC_EINDEP):			return "NC_EINDEP";
        case (NC_EFILE):			return "NC_EFILE";
        case (NC_EREAD):			return "NC_EREAD";
        case (NC_EWRITE):			return "NC_EWRITE";
        case (NC_EOFILE):			return "NC_EOFILE";
        case (NC_EMULTITYPES):			return "NC_EMULTITYPES";
        case (NC_EIOMISMATCH):			return "NC_EIOMISMATCH";
        case (NC_ENEGATIVECNT):			return "NC_ENEGATIVECNT";
        case (NC_EUNSPTETYPE):			return "NC_EUNSPTETYPE";
        case (NC_EINVAL_REQUEST):		return "NC_EINVAL_REQUEST";
        case (NC_EAINT_TOO_SMALL):		return "NC_EAINT_TOO_SMALL";
        case (NC_ENOENT):			return "NC_ENOENT";
#ifdef NC_EMULTIDEFINE
        case (NC_EMULTIDEFINE):			return "NC_EMULTIDEFINE";
#endif
#if PNETCDF_VERSION_MAJOR>=1 && PNETCDF_VERSION_MINOR>=3
        case (NC_ENOTSUPPORT):			return "NC_ENOTSUPPORT";
        case (NC_ENULLBUF):			return "NC_ENULLBUF";
        case (NC_EPREVATTACHBUF):		return "NC_EPREVATTACHBUF";
        case (NC_ENULLABUF):			return "NC_ENULLABUF";
        case (NC_EPENDINGBPUT):			return "NC_EPENDINGBPUT";
        case (NC_EINSUFFBUF):			return "NC_EINSUFFBUF";
#endif
#if PNETCDF_VERSION_MAJOR>=1 && PNETCDF_VERSION_MINOR>=4
        case (NC_EINTOVERFLOW):			return "NC_EINTOVERFLOW";
        case (NC_EMULTIDEFINE_OMODE):		return "NC_EMULTIDEFINE_OMODE";
        case (NC_EMULTIDEFINE_DIM_NUM):		return "NC_EMULTIDEFINE_DIM_NUM";
        case (NC_EMULTIDEFINE_DIM_SIZE):	return "NC_EMULTIDEFINE_DIM_SIZE";
        case (NC_EMULTIDEFINE_DIM_NAME):	return "NC_EMULTIDEFINE_DIM_NAME";
        case (NC_EMULTIDEFINE_VAR_NUM):		return "NC_EMULTIDEFINE_VAR_NUM";
        case (NC_EMULTIDEFINE_VAR_NAME):	return "NC_EMULTIDEFINE_VAR_NAME";
        case (NC_EMULTIDEFINE_VAR_NDIMS):	return "NC_EMULTIDEFINE_VAR_NDIMS";
        case (NC_EMULTIDEFINE_VAR_DIMIDS):	return "NC_EMULTIDEFINE_VAR_DIMIDS";
        case (NC_EMULTIDEFINE_VAR_TYPE):	return "NC_EMULTIDEFINE_VAR_TYPE";
        case (NC_EMULTIDEFINE_VAR_LEN):		return "NC_EMULTIDEFINE_VAR_LEN";
        case (NC_EMULTIDEFINE_NUMRECS):		return "NC_EMULTIDEFINE_NUMRECS";
        case (NC_EMULTIDEFINE_VAR_BEGIN):	return "NC_EMULTIDEFINE_VAR_BEGIN";
        case (NC_EMULTIDEFINE_ATTR_NUM):	return "NC_EMULTIDEFINE_ATTR_NUM";
        case (NC_EMULTIDEFINE_ATTR_SIZE):	return "NC_EMULTIDEFINE_ATTR_SIZE";
        case (NC_EMULTIDEFINE_ATTR_NAME):	return "NC_EMULTIDEFINE_ATTR_NAME";
        case (NC_EMULTIDEFINE_ATTR_TYPE):	return "NC_EMULTIDEFINE_ATTR_TYPE";
        case (NC_EMULTIDEFINE_ATTR_LEN):	return "NC_EMULTIDEFINE_ATTR_LEN";
        case (NC_EMULTIDEFINE_ATTR_VAL):	return "NC_EMULTIDEFINE_ATTR_VAL";
#endif
#if PNETCDF_VERSION_MAJOR>=1 && PNETCDF_VERSION_MINOR>=5
        case (NC_ENOTENABLED):			return "NC_ENOTENABLED";
        case (NC_EBAD_FILE):			return "NC_EBAD_FILE";
        case (NC_ENO_SPACE):			return "NC_ENO_SPACE";
        case (NC_EQUOTA):			return "NC_EQUOTA";
        case (NC_EMULTIDEFINE_FNC_ARGS):	return "NC_EMULTIDEFINE_FNC_ARGS";
#endif
#if PNETCDF_VERSION_MAJOR>=1 && PNETCDF_VERSION_MINOR>=6
        case (NC_EINVAL_CMODE):			return "NC_EINVAL_CMODE";
        case (NC_ENULLSTART):			return "NC_ENULLSTART";
        case (NC_ENULLCOUNT):			return "NC_ENULLCOUNT";
        case (NC_ETYPESIZE_MISMATCH):		return "NC_ETYPESIZE_MISMATCH";
        case (NC_ETYPESIZE):			return "NC_ETYPESIZE";
        case (NC_ETYPE_MISMATCH):		return "NC_ETYPE_MISMATCH";
        case (NC_ESTRICTCDF2):			return "NC_ESTRICTCDF2";
#endif
#if PNETCDF_VERSION_MAJOR>=1 && PNETCDF_VERSION_MINOR>=7
        case (NC_ENOTRECVAR):			return "NC_ENOTRECVAR";
        case (NC_ENOTFILL):			return "NC_ENOTFILL";
        case (NC_EMULTIDEFINE_FILL_MODE):	return "NC_EMULTIDEFINE_FILL_MODE";
        case (NC_EMULTIDEFINE_VAR_FILL_MODE):	return "NC_EMULTIDEFINE_VAR_FILL_MODE";
        case (NC_EMULTIDEFINE_VAR_FILL_VALUE):	return "NC_EMULTIDEFINE_VAR_FILL_VALUE";
#endif
#if PNETCDF_VERSION_MAJOR>=1 && PNETCDF_VERSION_MINOR>=8
        case (NC_EPENDING):			return "NC_EPENDING";
        case (NC_EINVAL_OMODE):			return "NC_EINVAL_OMODE";
        case (NC_EMULTIDEFINE_CMODE):		return "NC_EMULTIDEFINE_CMODE";
#endif
#endif
        default:
              sprintf(unknown_str,"Unknown code %d",err);
    }
    return unknown_str;
}


int
test_nc_against_pnetcdf(void)
{
    int format;

    nc_set_default_format(NC_FORMAT_CLASSIC, &format);
    nc_set_default_format(format, NULL); /* restore default */
    if (format == NC_FORMAT_NETCDF4 || format == NC_FORMAT_NETCDF4_CLASSIC)
        return 1; /* skip test for netcdf4 formats */

#ifdef USE_PNETCDF
    int  ncid; /* netCDF id */
    int  err;  /* status */

    /* Using netCDF library to create file */
    err = nc_create(scratch, NC_CLOBBER, &ncid);
    IF (err != NC_NOERR) error("nc_create: %s", nc_strerror(err));
    def_dims(ncid);
    def_vars(ncid);
    put_atts(ncid);
    err = nc_enddef(ncid);
    IF (err != NC_NOERR) error("nc_enddef: %s", nc_strerror(err));
    put_vars(ncid);
    err = nc_close (ncid);
    IF (err != NC_NOERR) error("nc_close: %s", nc_strerror(err));

    /* Using PnetCDF library to check file */
    err = nc_open_par(scratch, NC_NOWRITE, MPI_COMM_WORLD, MPI_INFO_NULL, &ncid);
    IF (err != NC_NOERR) error("nc_open_par: %s", nc_strerror(err));
    check_dims(ncid);
    check_vars(ncid);
    check_atts(ncid);
    err = nc_close (ncid);
    IF (err != NC_NOERR) error("nc_close: %s", nc_strerror(err));

    /* Using PnetCDF library to create file */
    err = nc_create_par(scratch, 0, MPI_COMM_WORLD, MPI_INFO_NULL, &ncid);
    IF (err != NC_NOERR) error("nc_create_par: %s", nc_strerror(err));
    def_dims(ncid);
    def_vars(ncid);
    put_atts(ncid);
    err = nc_enddef(ncid);
    IF (err != NC_NOERR) error("nc_enddef: %s", nc_strerror(err));
    put_vars(ncid);
    err = nc_close (ncid);
    IF (err != NC_NOERR) error("nc_close: %s", nc_strerror(err));

    /* Using NetCDF library to check file */
    err = nc_open(scratch, NC_NOWRITE, &ncid);
    IF (err != NC_NOERR) error("nc_open: %s", nc_strerror(err));
    check_dims(ncid);
    check_vars(ncid);
    check_atts(ncid);
    err = nc_close (ncid);
    IF (err != NC_NOERR) error("nc_close: %s", nc_strerror(err));
    err = nc_delete(scratch);
    IF (err != NC_NOERR) error("remove of %s failed", scratch);
#endif
    return 1;
}

// The libMesh Finite Element Library.
// Copyright (C) 2002-2022 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

#ifndef LIBMESH_FLOAT128_SHIMS_H
#define LIBMESH_FLOAT128_SHIMS_H

// The library configuration options
#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_DEFAULT_QUADRUPLE_PRECISION

# include <boost/multiprecision/float128.hpp>

// Boost doesn't add float128 overloads to namespace std.  But, we
// expect to be able to explicitly specify std::foo(Real) in trailing
// return types, where a using-declaration + ADL isn't an option.  So
// we add overloads ourselves.
namespace std
{
#define LIBMESH_FLOAT128_UNARY(funcname) \
inline boost::multiprecision::float128 funcname \
  (const boost::multiprecision::float128 in) \
{ \
  return boost::multiprecision::funcname(in); \
}

#define LIBMESH_FLOAT128_MATH_BOOL(funcname) \
inline bool funcname \
  (const boost::multiprecision::float128 in) \
{ \
  return boost::math::funcname(in); \
}

#define LIBMESH_FLOAT128_BINARY(funcname) \
inline boost::multiprecision::float128 funcname \
  (const boost::multiprecision::float128 in1, \
   const boost::multiprecision::float128 in2) \
{ \
  return boost::multiprecision::funcname(in1, in2); \
}

LIBMESH_FLOAT128_UNARY(sqrt)
LIBMESH_FLOAT128_UNARY(exp)
LIBMESH_FLOAT128_UNARY(log)
LIBMESH_FLOAT128_UNARY(log10)
LIBMESH_FLOAT128_UNARY(sin)
LIBMESH_FLOAT128_UNARY(cos)
LIBMESH_FLOAT128_UNARY(tan)
LIBMESH_FLOAT128_UNARY(asin)
LIBMESH_FLOAT128_UNARY(acos)
LIBMESH_FLOAT128_UNARY(atan)
LIBMESH_FLOAT128_UNARY(sinh)
LIBMESH_FLOAT128_UNARY(cosh)
LIBMESH_FLOAT128_UNARY(tanh)
LIBMESH_FLOAT128_UNARY(abs)
LIBMESH_FLOAT128_UNARY(fabs)
LIBMESH_FLOAT128_UNARY(ceil)
LIBMESH_FLOAT128_UNARY(floor)
// LIBMESH_FLOAT128_UNARY(norm)

inline boost::multiprecision::float128 norm
  (const boost::multiprecision::float128 in)
{
  return in * in;
}

inline boost::multiprecision::float128 real
  (const boost::multiprecision::float128 in)
{
  return in;
}

inline boost::multiprecision::float128 imag
  (const boost::multiprecision::float128 /*in*/)
{
  return 0;
}

template <>
struct plus<boost::multiprecision::float128>
{
  boost::multiprecision::float128 operator ()
    (const boost::multiprecision::float128 a,
     const boost::multiprecision::float128 b)
  {
    return a + b;
  }
};

template <>
struct multiplies<boost::multiprecision::float128>
{
  boost::multiprecision::float128 operator ()
    (const boost::multiprecision::float128 a,
     const boost::multiprecision::float128 b)
  {
    return a * b;
  }
};


LIBMESH_FLOAT128_MATH_BOOL(isinf)
LIBMESH_FLOAT128_MATH_BOOL(isnan)

// Shimming modf by hand since Boost didn't add a shim until 1.62 and
// I'd like to still support systems with older Boost in /usr/include/
inline boost::multiprecision::float128 modf
  (const boost::multiprecision::float128 in,
   boost::multiprecision::float128 * intpart)
{
#ifdef BOOST_MP_USE_QUAD
  return __modfq(in.backend().value(), &intpart->backend().value());
#elif defined(BOOST_MP_USE_FLOAT128)
  return modfq(in.backend().value(), &intpart->backend().value());
#endif
}

LIBMESH_FLOAT128_BINARY(pow)
LIBMESH_FLOAT128_BINARY(atan2)
// LIBMESH_FLOAT128_BINARY(max)
// LIBMESH_FLOAT128_BINARY(min)
LIBMESH_FLOAT128_BINARY(fmod)

inline boost::multiprecision::float128 pow
  (const boost::multiprecision::float128 in1,
   const int in2)
{
  return boost::multiprecision::pow(in1, in2);
}

// Boost float128 leaves a lot of C++11 math undefined?  So we'll just
// add shims as we need them, for maximum compatibility with older
// Boost versions.

// Stuff that was defined as far back as I've tested:
LIBMESH_FLOAT128_UNARY(trunc)
LIBMESH_FLOAT128_UNARY(round)

// log1p was added in Boost 1.63
#if BOOST_VERSION > 106300
LIBMESH_FLOAT128_UNARY(log1p)
#endif

// This doesn't take Real->Real:
inline long long llround
  (const boost::multiprecision::float128 in)
{
  return boost::multiprecision::llround(in);
}

// Stuff that wasn't, that we don't need yet:
// LIBMESH_FLOAT128_UNARY(exp2)
// LIBMESH_FLOAT128_UNARY(expm1)
// LIBMESH_FLOAT128_UNARY(log2)
// LIBMESH_FLOAT128_UNARY(cbrt)
// LIBMESH_FLOAT128_UNARY(asinh)
// LIBMESH_FLOAT128_UNARY(acosh)
// LIBMESH_FLOAT128_UNARY(atanh)
// LIBMESH_FLOAT128_UNARY(erf)
// LIBMESH_FLOAT128_UNARY(erfc)
// LIBMESH_FLOAT128_UNARY(nearbyint)
// LIBMESH_FLOAT128_UNARY(rint)

// LIBMESH_FLOAT128_BINARY(remainder)
// LIBMESH_FLOAT128_BINARY(fmax)
// LIBMESH_FLOAT128_BINARY(fmin)
// LIBMESH_FLOAT128_BINARY(fdim)
// LIBMESH_FLOAT128_BINARY(hypot)

}

#endif // LIBMESH_DEFAULT_QUADRUPLE_PRECISION

#endif // LIBMESH_FLOAT128_SHIMS_H

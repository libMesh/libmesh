// The libMesh Finite Element Library.
// Copyright (C) 2002-2019 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

LIBMESH_FLOAT128_MATH_BOOL(isinf)
LIBMESH_FLOAT128_MATH_BOOL(isnan)

inline boost::multiprecision::float128 modf
  (const boost::multiprecision::float128 in,
   boost::multiprecision::float128 * intpart)
{
  return boost::math::modf(in, intpart);
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

// Boost leaves a lot of C++11 math undefined??

// LIBMESH_FLOAT128_UNARY(exp2)
// LIBMESH_FLOAT128_UNARY(expm1)
// LIBMESH_FLOAT128_UNARY(log2)
// LIBMESH_FLOAT128_UNARY(log1p)
// LIBMESH_FLOAT128_UNARY(cbrt)
// LIBMESH_FLOAT128_UNARY(asinh)
// LIBMESH_FLOAT128_UNARY(acosh)
// LIBMESH_FLOAT128_UNARY(atanh)
// LIBMESH_FLOAT128_UNARY(erf)
// LIBMESH_FLOAT128_UNARY(erfc)
LIBMESH_FLOAT128_UNARY(trunc)
LIBMESH_FLOAT128_UNARY(round)
// LIBMESH_FLOAT128_UNARY(nearbyint)
// LIBMESH_FLOAT128_UNARY(rint)

// LIBMESH_FLOAT128_BINARY(remainder)
// LIBMESH_FLOAT128_BINARY(fmax)
// LIBMESH_FLOAT128_BINARY(fmin)
// LIBMESH_FLOAT128_BINARY(fdim)
// LIBMESH_FLOAT128_BINARY(hypot)

}

#endif // LIBMESH_FLOAT128_SHIMS_H

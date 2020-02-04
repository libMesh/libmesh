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


#ifndef LIBMESH_COMPARE_TYPES_H
#define LIBMESH_COMPARE_TYPES_H

// System includes
#include <complex>

#include "libmesh_config.h"

#ifdef LIBMESH_DEFAULT_QUADRUPLE_PRECISION
#include "libmesh_common.h" // for Real
#endif

namespace libMesh
{

// Copy of boost's enable_if_c

namespace boostcopy {
template <bool B, class T = void>
struct enable_if_c {
  typedef T type;
};

template <class T>
struct enable_if_c<false, T> {};
}


// TypesEqual takes two types as parameters.
// If they are the exact same type, then TypesEqual::value is the boolean true,
// otherwise TypesEqual::value is the boolean false.
template <typename T1, typename T2>
struct TypesEqual {
  static const bool value = false;
};

template <typename T>
struct TypesEqual<T,T> {
  static const bool value = true;
};


// Complete list of scalar classes, needed for disambiguation
template <typename T>
struct ScalarTraits {
  static const bool value = false;
};

#define ScalarTraits_true(type)                                 \
  template<>                                                    \
  struct ScalarTraits<type> { static const bool value = true; }

ScalarTraits_true(char);
ScalarTraits_true(short);
ScalarTraits_true(int);
ScalarTraits_true(long);
ScalarTraits_true(unsigned char);
ScalarTraits_true(unsigned short);
ScalarTraits_true(unsigned int);
ScalarTraits_true(unsigned long);
ScalarTraits_true(float);
ScalarTraits_true(double);
ScalarTraits_true(long double);
#ifdef LIBMESH_DEFAULT_QUADRUPLE_PRECISION
ScalarTraits_true(Real);
#endif

template<typename T>
struct ScalarTraits<std::complex<T>> { static const bool value = ScalarTraits<T>::value; };



// Operators using different but compatible types need a return value
// based on whichever type the other can be upconverted into.  For
// instance, the proper return type for
// TypeVector<float>::operator*(double) is TypeVector<double>.  In
// general, an operation using types S and T should return a value
// based on CompareTypes<S,T>::supertype

template<typename S, typename T>
struct CompareTypes {
  typedef void supertype;
};

template<typename T>
struct CompareTypes<T, T> {
  typedef T supertype;
};

template<typename T>
struct CompareTypes<T, std::complex<T>> {
  typedef std::complex<T> supertype;
};

template<typename T>
struct CompareTypes<std::complex<T>, T> {
  typedef std::complex<T> supertype;
};

// There's got to be some magic template way to do these better - but the best
// thing on the net requires a bunch of Alexandrescu's code and doesn't work
// with older compilers

#define CompareTypes_super(a,b,super)           \
  template<>                                    \
  struct CompareTypes<a, b> {                   \
    typedef super supertype;                    \
  }

#define SUPERTYPE(mysub,mysuper)                                        \
  CompareTypes_super(mysub, mysuper, mysuper);                          \
  CompareTypes_super(mysuper, mysub, mysuper);                          \
  CompareTypes_super(std::complex<mysub>, mysuper, std::complex<mysuper>); \
  CompareTypes_super(mysuper, std::complex<mysub>, std::complex<mysuper>); \
  CompareTypes_super(mysub, std::complex<mysuper>, std::complex<mysuper>); \
  CompareTypes_super(std::complex<mysuper>, mysub, std::complex<mysuper>); \
  CompareTypes_super(std::complex<mysub>, std::complex<mysuper>, std::complex<mysuper>); \
  CompareTypes_super(std::complex<mysuper>, std::complex<mysub>, std::complex<mysuper>)

SUPERTYPE(unsigned char, short);
SUPERTYPE(unsigned char, int);
SUPERTYPE(unsigned char, float);
SUPERTYPE(unsigned char, double);
SUPERTYPE(unsigned char, long double);
SUPERTYPE(unsigned short, int);
SUPERTYPE(unsigned short, float);
SUPERTYPE(unsigned short, double);
SUPERTYPE(unsigned short, long double);
SUPERTYPE(unsigned int, float);
SUPERTYPE(unsigned int, double);
SUPERTYPE(unsigned int, long double);
SUPERTYPE(char, short);
SUPERTYPE(char, int);
SUPERTYPE(char, float);
SUPERTYPE(char, double);
SUPERTYPE(char, long double);
SUPERTYPE(short, int);
SUPERTYPE(short, float);
SUPERTYPE(short, double);
SUPERTYPE(short, long double);
SUPERTYPE(int, float);
SUPERTYPE(int, double);
SUPERTYPE(int, long double);
SUPERTYPE(float, double);
SUPERTYPE(float, long double);
SUPERTYPE(double, long double);

#ifdef LIBMESH_DEFAULT_QUADRUPLE_PRECISION
SUPERTYPE(unsigned char, Real);
SUPERTYPE(unsigned short, Real);
SUPERTYPE(unsigned int, Real);
SUPERTYPE(char, Real);
SUPERTYPE(short, Real);
SUPERTYPE(int, Real);
SUPERTYPE(float, Real);
SUPERTYPE(double, Real);
SUPERTYPE(long double, Real);
#endif

#undef CompareTypes_super
#undef SUPERTYPE

template <typename T>
struct GetScalarType
{
  typedef typename std::enable_if<ScalarTraits<T>::value, T>::type type;
};

template <template <typename> class Wrapper, typename T>
struct GetScalarType<Wrapper<T>>
{
  typedef typename std::enable_if<ScalarTraits<T>::value, T>::type type;
};

// gcc can't tell which of the following is the most specialized?  Weak.
/*
  template<typename S, typename T>
  struct CompareTypes<std::complex<S>, std::complex<T>> {
  typedef std::complex<typename CompareTypes<S, T>::supertype> supertype;
  };

  template<typename S, typename T>
  struct CompareTypes<std::complex<S>, T> {
  typedef std::complex<typename CompareTypes<S, T>::supertype> supertype;
  };

  template<typename S, typename T>
  struct CompareTypes<S, std::complex<T>> {
  typedef std::complex<typename CompareTypes<S, T>::supertype> supertype;
  };
*/

} // namespace libMesh

#ifdef LIBMESH_HAVE_METAPHYSICL
#define LIBMESH_DUAL_NUMBER_COMPARE_TYPES

namespace MetaPhysicL
{
template <typename, typename>
class DualNumber;
} // namespace MetaPhysicL

namespace libMesh
{
template <typename T, typename T2, typename D>
struct CompareTypes<T, MetaPhysicL::DualNumber<T2, D>>
{
  typedef MetaPhysicL::DualNumber<typename CompareTypes<T, T2>::supertype,
                     typename D::template rebind<typename CompareTypes<T, T2>::supertype>::other>
      supertype;
};
template <typename T, typename D, typename T2>
struct CompareTypes<MetaPhysicL::DualNumber<T, D>, T2>
{
  typedef MetaPhysicL::DualNumber<typename CompareTypes<T, T2>::supertype,
                     typename D::template rebind<typename CompareTypes<T, T2>::supertype>::other>
      supertype;
};
template <typename T, typename D, typename T2, typename D2>
struct CompareTypes<MetaPhysicL::DualNumber<T, D>, MetaPhysicL::DualNumber<T2, D2>>
{
  typedef MetaPhysicL::DualNumber<typename CompareTypes<T, T2>::supertype,
                     typename D::template rebind<typename CompareTypes<T, T2>::supertype>::other>
      supertype;
};
template <typename T, typename D>
struct CompareTypes<MetaPhysicL::DualNumber<T, D>, MetaPhysicL::DualNumber<T, D>>
{
  typedef MetaPhysicL::DualNumber<T, D> supertype;
};
template <typename T, typename T2>
struct CompareTypes<T, MetaPhysicL::DualNumber<T2, T2>>
{
  typedef MetaPhysicL::DualNumber<typename CompareTypes<T, T2>::supertype,
                     typename CompareTypes<T, T2>::supertype>
      supertype;
};
template <typename T, typename T2>
struct CompareTypes<MetaPhysicL::DualNumber<T, T>, T2>
{
  typedef MetaPhysicL::DualNumber<typename CompareTypes<T, T2>::supertype,
                     typename CompareTypes<T, T2>::supertype>
      supertype;
};
template <typename T, typename T2>
struct CompareTypes<MetaPhysicL::DualNumber<T, T>, MetaPhysicL::DualNumber<T2, T2>>
{
  typedef MetaPhysicL::DualNumber<typename CompareTypes<T, T2>::supertype,
                     typename CompareTypes<T, T2>::supertype>
      supertype;
};
template <typename T, typename D>
struct ScalarTraits<MetaPhysicL::DualNumber<T, D>>
{
  static const bool value = ScalarTraits<T>::value;
};
} // namespace libMesh

#endif // LIBMESH_HAVE_METAPHYSICL

#endif // LIBMESH_COMPARE_TYPES_H

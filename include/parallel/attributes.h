// The libMesh Finite Element Library.
// Copyright (C) 2002-2018 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


#ifndef LIBMESH_ATTRIBUTES_H
#define LIBMESH_ATTRIBUTES_H

// C++ includes
#include <limits>
#include <set>
#include <vector>


namespace libMesh
{

namespace Parallel
{
//-------------------------------------------------------------------

/*
 * The unspecialized class gives default, lowest-common-denominator
 * attributes, for values which can't be used with Parallel min/max.
 * Specialized classes can set this to true, and should define
 * the lowest and highest values possible for the type.
 */
template<typename T>
struct Attributes
{
  static const bool has_min_max = false;
  static void set_lowest(T &) {}
  static void set_highest(T &) {}
};

// ------------------------------------------------------------
// Declare Attributes specializations for C++ built-in types

#define LIBMESH_INT_TYPE(cxxtype)                                       \
  template<>                                                            \
  struct Attributes<cxxtype>                                            \
  {                                                                     \
    static const bool has_min_max = true;                               \
    static void set_lowest(cxxtype & x) { x = std::numeric_limits<cxxtype>::min(); } \
    static void set_highest(cxxtype & x) { x = std::numeric_limits<cxxtype>::max(); } \
  }

#define LIBMESH_FLOAT_TYPE(cxxtype)                                     \
  template<>                                                            \
  struct Attributes<cxxtype>                                            \
  {                                                                     \
    static const bool has_min_max = true;                               \
    static void set_lowest(cxxtype & x) { x = -std::numeric_limits<cxxtype>::infinity(); } \
    static void set_highest(cxxtype & x) { x = std::numeric_limits<cxxtype>::infinity(); } \
  }

#define LIBMESH_CONTAINER_TYPE(cxxtype)                                 \
  struct Attributes<cxxtype>                                         \
  {                                                                     \
    static const bool has_min_max = Attributes<T>::has_min_max;         \
    static void set_lowest(cxxtype & x) {                            \
      for (auto & val : x)                                              \
        Attributes<T>::set_lowest(val); }                               \
    static void set_highest(cxxtype & x) {                           \
      for (auto & val : x)                                              \
        Attributes<T>::set_highest(val); }                              \
  }


LIBMESH_INT_TYPE(char);
LIBMESH_INT_TYPE(signed char);
LIBMESH_INT_TYPE(unsigned char);
LIBMESH_INT_TYPE(short int);
LIBMESH_INT_TYPE(unsigned short int);
LIBMESH_INT_TYPE(int);
LIBMESH_INT_TYPE(unsigned int);
LIBMESH_INT_TYPE(long);
LIBMESH_INT_TYPE(long long);
LIBMESH_INT_TYPE(unsigned long);
LIBMESH_INT_TYPE(unsigned long long);

LIBMESH_FLOAT_TYPE(float);
LIBMESH_FLOAT_TYPE(double);
LIBMESH_FLOAT_TYPE(long double);

#define LIBMESH_ATTRIBUTES_COMMA ,

template <typename T, typename C, typename A>
LIBMESH_CONTAINER_TYPE(std::set<T LIBMESH_ATTRIBUTES_COMMA C LIBMESH_ATTRIBUTES_COMMA A>);

template <typename T, typename A>
LIBMESH_CONTAINER_TYPE(std::vector<T LIBMESH_ATTRIBUTES_COMMA A>);

} // namespace Parallel

} // namespace libMesh

#endif // LIBMESH_ATTRIBUTES_H

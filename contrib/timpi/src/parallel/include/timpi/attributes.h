// The TIMPI Message-Passing Parallelism Library.
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


#ifndef TIMPI_ATTRIBUTES_H
#define TIMPI_ATTRIBUTES_H

// TIMPI includes
#include "timpi/timpi_config.h"

// C++ includes
#include <limits>
#include <set>
#include <vector>


namespace TIMPI
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

#define TIMPI_INT_TYPE(cxxtype)                                                    \
  template<>                                                                          \
  struct Attributes<cxxtype>                                                          \
  {                                                                                   \
    static const bool has_min_max = true;                                             \
    static void set_lowest(cxxtype & x) { x = std::numeric_limits<cxxtype>::min(); }  \
    static void set_highest(cxxtype & x) { x = std::numeric_limits<cxxtype>::max(); } \
  }

#define TIMPI_FLOAT_TYPE(cxxtype)                                                       \
  template<>                                                                               \
  struct Attributes<cxxtype>                                                               \
  {                                                                                        \
    static const bool has_min_max = true;                                                  \
    static void set_lowest(cxxtype & x) { x = -std::numeric_limits<cxxtype>::infinity(); } \
    static void set_highest(cxxtype & x) { x = std::numeric_limits<cxxtype>::infinity(); } \
  }

#define TIMPI_CONTAINER_TYPE(cxxtype)                        \
  struct Attributes<cxxtype>                                    \
  {                                                             \
    static const bool has_min_max = Attributes<T>::has_min_max; \
    static void set_lowest(cxxtype & x) {                       \
      for (auto & val : x)                                      \
        Attributes<T>::set_lowest(val); }                       \
    static void set_highest(cxxtype & x) {                      \
      for (auto & val : x)                                      \
        Attributes<T>::set_highest(val); }                      \
  }


TIMPI_INT_TYPE(char);
TIMPI_INT_TYPE(signed char);
TIMPI_INT_TYPE(unsigned char);
TIMPI_INT_TYPE(short int);
TIMPI_INT_TYPE(unsigned short int);
TIMPI_INT_TYPE(int);
TIMPI_INT_TYPE(unsigned int);
TIMPI_INT_TYPE(long);
TIMPI_INT_TYPE(long long);
TIMPI_INT_TYPE(unsigned long);
TIMPI_INT_TYPE(unsigned long long);

TIMPI_FLOAT_TYPE(float);
TIMPI_FLOAT_TYPE(double);
TIMPI_FLOAT_TYPE(long double);
#ifdef TIMPI_DEFAULT_QUADRUPLE_PRECISION
TIMPI_FLOAT_TYPE(Real);
#endif

#define TIMPI_ATTRIBUTES_COMMA ,

template <typename T, typename C, typename A>
TIMPI_CONTAINER_TYPE(std::set<T TIMPI_ATTRIBUTES_COMMA C TIMPI_ATTRIBUTES_COMMA A>);

template <typename T, typename A>
TIMPI_CONTAINER_TYPE(std::vector<T TIMPI_ATTRIBUTES_COMMA A>);

} // namespace TIMPI

#endif // TIMPI_ATTRIBUTES_H

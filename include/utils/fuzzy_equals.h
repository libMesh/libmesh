// The libMesh Finite Element Library.
// Copyright (C) 2002-2025 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#ifndef LIBMESH_FUZZY_EQUAL_H
#define LIBMESH_FUZZY_EQUAL_H

#include "libmesh/libmesh_common.h"
#include "libmesh/tensor_tools.h"
#include "libmesh/compare_types.h"
#include "libmesh/int_range.h"

#ifdef LIBMESH_HAVE_METAPHYSICL
#include "metaphysicl/raw_type.h"
#endif

namespace libMesh
{
/**
 * Computes the L1 norm
 */
template <typename T, typename std::enable_if<ScalarTraits<T>::value, int>::type = 0>
auto
l1_norm(const T & var)
{
  return std::abs(var);
}

/**
 * Computes the L1 norm of the diff between \p var1 and \p var2
 */
template <typename T,
          typename T2,
          typename std::enable_if<ScalarTraits<T>::value && ScalarTraits<T2>::value, int>::type = 0>
auto
l1_norm_diff(const T & var1, const T2 & var2)
{
  return l1_norm(var1 - var2);
}

#ifdef LIBMESH_HAVE_METAPHYSICL
/**
 * Function to check whether two variables are equal within an absolute tolerance
 * @param var1 The first variable to be checked
 * @param var2 The second variable to be checked
 * @param tol The tolerance to be used
 * @return true if var1 and var2 are equal within tol
 */
template <typename T, typename T2>
bool
absolute_fuzzy_equals(const T & var1, const T2 & var2, const Real tol = TOLERANCE * TOLERANCE)
{
  return MetaPhysicL::raw_value(l1_norm_diff(var1, var2)) <= tol;
}

/**
 * Function to check whether two variables are equal within a relative tolerance
 * @param var1 The first variable to be checked
 * @param var2 The second variable to be checked
 * @param tol The relative tolerance to be used
 * @return true if var1 and var2 are equal within relative tol
 */
template <typename T, typename T2>
bool
relative_fuzzy_equals(const T & var1, const T2 & var2, const Real tol = TOLERANCE * TOLERANCE)
{
  return absolute_fuzzy_equals(
      var1,
      var2,
      tol * (MetaPhysicL::raw_value(l1_norm(var1)) + MetaPhysicL::raw_value(l1_norm(var2))));
}
#else
/**
 * Function to check whether two variables are equal within an absolute tolerance
 * @param var1 The first variable to be checked
 * @param var2 The second variable to be checked
 * @param tol The tolerance to be used
 * @return true if var1 and var2 are equal within tol
 */
template <typename T, typename T2>
bool
absolute_fuzzy_equals(const T & var1, const T2 & var2, const Real tol = TOLERANCE * TOLERANCE)
{
  return l1_norm_diff(var1, var2) <= tol;
}

/**
 * Function to check whether two variables are equal within a relative tolerance
 * @param var1 The first variable to be checked
 * @param var2 The second variable to be checked
 * @param tol The relative tolerance to be used
 * @return true if var1 and var2 are equal within relative tol
 */
template <typename T, typename T2>
bool
relative_fuzzy_equals(const T & var1, const T2 & var2, const Real tol = TOLERANCE * TOLERANCE)
{
  return absolute_fuzzy_equals(var1, var2, tol * (l1_norm(var1) + l1_norm(var2)));
}
#endif

} // namespace libMesh

#endif // LIBMESH_FUZZY_EQUAL_H

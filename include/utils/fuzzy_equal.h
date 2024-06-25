// The libMesh Finite Element Library.
// Copyright (C) 2002-2024 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

namespace libMesh
{
#ifdef LIBMESH_HAVE_METAPHYSICL

#include "metaphysicl/raw_type.h"

/**
 * Function to check whether two variables are equal within an absolute tolerance
 * @param var1 The first variable to be checked
 * @param var2 The second variable to be checked
 * @param tol The tolerance to be used
 * @return true if var1 and var2 are equal within tol
 */
template <typename T,
          typename T2,
          typename T3 = T,
          typename std::enable_if<ScalarTraits<T>::value && ScalarTraits<T2>::value &&
                                      ScalarTraits<T3>::value,
                                  int>::type = 0>
bool
absolute_fuzzy_equal(const T & var1, const T2 & var2, const T3 & tol = TOLERANCE * TOLERANCE)
{
  return (std::abs(MetaPhysicL::raw_value(var1) - MetaPhysicL::raw_value(var2)) <=
          MetaPhysicL::raw_value(tol));
}

/**
 * Function to check whether two variables are equal within a relative tolerance
 * @param var1 The first variable to be checked
 * @param var2 The second variable to be checked
 * @param tol The relative tolerance to be used
 * @return true if var1 and var2 are equal within relative tol
 */
template <typename T, typename T2, typename T3 = Real>
bool
relative_fuzzy_equal(const T & var1, const T2 & var2, const T3 & tol = TOLERANCE * TOLERANCE)
{
  if constexpr (ScalarTraits<T>::value || TensorTools::MathWrapperTraits<T>::value)
  {
    static_assert(TensorTools::TensorTraits<T>::rank == TensorTools::TensorTraits<T2>::rank,
                  "Mathematical types must be same for arguments to relativelyFuzzEqual");
    if constexpr (TensorTools::TensorTraits<T>::rank == 0)
      return absolute_fuzzy_equal(
          var1,
          var2,
          tol * (std::abs(MetaPhysicL::raw_value(var1)) + std::abs(MetaPhysicL::raw_value(var2))));
    else if constexpr (TensorTools::TensorTraits<T>::rank == 1)
    {
      for (const auto i : make_range(max_mesh_dim))
        if (!relative_fuzzy_equal(var1(i), var2(i), tol))
          return false;

      return true;
    }
    else if constexpr (TensorTools::TensorTraits<T>::rank == 2)
    {
      for (const auto i : make_range(max_mesh_dim))
        for (const auto j : make_range(max_mesh_dim))
          if (!relative_fuzzy_equal(var1(i, j), var2(i, j), tol))
            return false;

      return true;
    }
  }
  else
  {
    // We dare to dream
    mooseAssert(var1.size() == var2.size(), "These must be the same size");
    for (const auto i : index_range(var1))
      if (!relative_fuzzy_equal(var1(i), var2(i), tol))
        return false;

    return true;
  }
}

#else

/**
 * Function to check whether two variables are equal within an absolute tolerance
 * @param var1 The first variable to be checked
 * @param var2 The second variable to be checked
 * @param tol The tolerance to be used
 * @return true if var1 and var2 are equal within tol
 */
template <typename T,
          typename T2,
          typename T3 = T,
          typename std::enable_if<ScalarTraits<T>::value && ScalarTraits<T2>::value &&
                                      ScalarTraits<T3>::value,
                                  int>::type = 0>
bool
absolute_fuzzy_equal(const T & var1, const T2 & var2, const T3 & tol = TOLERANCE * TOLERANCE)
{
  return (std::abs(var1 - var2) <= tol);
}

/**
 * Function to check whether two variables are equal within a relative tolerance
 * @param var1 The first variable to be checked
 * @param var2 The second variable to be checked
 * @param tol The relative tolerance to be used
 * @return true if var1 and var2 are equal within relative tol
 */
template <typename T, typename T2, typename T3 = Real>
bool
relative_fuzzy_equal(const T & var1, const T2 & var2, const T3 & tol = TOLERANCE * TOLERANCE)
{
  if constexpr (ScalarTraits<T>::value || TensorTools::MathWrapperTraits<T>::value)
  {
    static_assert(TensorTools::TensorTraits<T>::rank == TensorTools::TensorTraits<T2>::rank,
                  "Mathematical types must be same for arguments to relativelyFuzzEqual");
    if constexpr (TensorTools::TensorTraits<T>::rank == 0)
      return absolute_fuzzy_equal(var1, var2, tol * (std::abs(var1) + std::abs(var2)));
    else if constexpr (TensorTools::TensorTraits<T>::rank == 1)
    {
      for (const auto i : make_range(max_mesh_dim))
        if (!relative_fuzzy_equal(var1(i), var2(i), tol))
          return false;

      return true;
    }
    else if constexpr (TensorTools::TensorTraits<T>::rank == 2)
    {
      for (const auto i : make_range(max_mesh_dim))
        for (const auto j : make_range(max_mesh_dim))
          if (!relative_fuzzy_equal(var1(i, j), var2(i, j), tol))
            return false;

      return true;
    }
  }
  else
  {
    // We dare to dream
    mooseAssert(var1.size() == var2.size(), "These must be the same size");
    for (const auto i : index_range(var1))
      if (!relative_fuzzy_equal(var1(i), var2(i), tol))
        return false;

    return true;
  }
}

#endif // !LIBMESH_HAVE_METAPHYSICL

} // namespace libMesh

#endif // LIBMESH_FUZZY_EQUAL_H

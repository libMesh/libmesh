// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_PARALLEL_CONVERSION_UTILS_H
#define LIBMESH_PARALLEL_CONVERSION_UTILS_H

// Local includes
#include "libmesh/libmesh_common.h"

#ifdef LIBMESH_HAVE_LIBHILBERT
#include "hilbert.h"
#endif

// C++ includes
#include <vector>

namespace libMesh
{



//--------------------------------------------------------------------------
namespace Parallel {
namespace Utils {

/**
 * Utility function that returns true if the vector v
 * is sorted, false otherwise.  O(N), the length of the
 * vector.  This is implemented solely because the std::is_sorted
 * appears to be an STL extension.
 */
template <typename KeyType>
inline
bool is_sorted (const std::vector<KeyType> & v)
{
  if (v.empty())
    return true;

  for (std::size_t i=1; i<v.size(); i++)
    if (v[i] < v[i-1])
      return false;

  return true;
}

/**
 * A utility function which converts whatever \p KeyType is to
 * a \p double for the histogram bounds
 */
template <typename KeyType>
inline
double to_double (const KeyType & k)
{
  return static_cast<double>(k);
}

/**
 * A utility to convert a \p double to some sort of \p KeyType, for
 * interpreting how histogram bounds relate to \p KeyType positions.
 *
 * This is a class to allow partial template specialization for the
 * std::pair case without adding a "dummy" variable.
 */
template <typename KeyType>
struct Convert {
  inline static
  KeyType to_key_type (const double f)
  {
    return static_cast<KeyType>(f);
  }
};

/**
 * A pseudoinverse for converting bounds back to pairs of key types.
 */
template <typename FirstKeyType, typename SecondKeyType>
struct Convert<std::pair<FirstKeyType, SecondKeyType> > {
  inline static
  std::pair<FirstKeyType,SecondKeyType> to_key_type (const double f)
  {
    return std::make_pair
      (Convert<FirstKeyType>::to_key_type(f),SecondKeyType());
  }
};



/**
 * A utility function for pairs of key types.  When finding bounds,
 * the second entry of the pair is effectively "rounded away".
 */
template <typename FirstKeyType, typename SecondKeyType>
inline
double to_double (const std::pair<FirstKeyType,SecondKeyType> &k)
{
  return to_double(k.first);
}


#ifdef LIBMESH_HAVE_LIBHILBERT

template <>
inline
double to_double (const Hilbert::HilbertIndices & bvt)
{
  return static_cast<double>(bvt.rack2);
}

template <>
struct Convert<Hilbert::HilbertIndices> {
  inline static
  Hilbert::HilbertIndices
  to_key_type (const double f)
  {
    Hilbert::HilbertIndices bvt;

    bvt.rack0 = 0;
    bvt.rack1 = 0;
    bvt.rack2 = static_cast<Hilbert::inttype>(f);

    return bvt;
  }
};
#endif // LIBMESH_HAVE_LIBHILBERT
}
}

} // namespace libMesh

#endif // LIBMESH_PARALLEL_CONVERSION_UTILS_H

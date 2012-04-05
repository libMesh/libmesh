// The libMesh Finite Element Library.
// Copyright (C) 2002-2012 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef __parallel_conversion_utils_h__
#define __parallel_conversion_utils_h__

#include <vector>

#include "libmesh_common.h"

#ifdef LIBMESH_HAVE_LIBHILBERT
#include "hilbert.h"
#endif

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
      bool is_sorted (const std::vector<KeyType>& v)
      {
	if (v.empty())
	  return true;

	for (unsigned int i=1; i<v.size(); i++)
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
      double to_double (const KeyType &k)
      {
	return static_cast<double>(k);
      }

    /**
     * A utility to convert a \p double to some
     * sort of \p KeyType, for interpreting how
     * histogram bounds relate to \p KeyType
     * positions.
     */
    template <typename KeyType>
      inline
      KeyType to_key_type (const double f)
      {
	return static_cast<KeyType>(f);
      }

#ifdef LIBMESH_HAVE_LIBHILBERT

    template <>
      inline
      double to_double (const Hilbert::HilbertIndices &bvt)
      {
	return static_cast<double>(bvt.rack2);
      }

    template <>
      inline
      Hilbert::HilbertIndices
      to_key_type (const double f)
      {
	Hilbert::HilbertIndices bvt;

	bvt.rack0 = 0;
	bvt.rack1 = 0;
	bvt.rack2 = static_cast<Hilbert::inttype>(f);

	return bvt;
      }
#endif // LIBMESH_HAVE_LIBHILBERT
  }
}

} // namespace libMesh

#endif // #define __parallel_conversion_utils_h__

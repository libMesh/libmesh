// The libMesh Finite Element Library.
// Copyright (C) 2002-2016 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


#ifndef LIBMESH_PARALLEL_BIN_SORTER_H
#define LIBMESH_PARALLEL_BIN_SORTER_H

// This class contains all the functionality for bin sorting
// Templated on the type of keys you will be sorting and the
// type of iterator you will be using.

// libMesh includes
#include "libmesh/libmesh_common.h" // cast_int
#include "libmesh/parallel_object.h"

// C++ includes
#include <vector>
#include <iterator>

namespace libMesh
{

namespace Parallel {

template <typename KeyType, typename IdxType=unsigned int>
/**
 * Perform a parallel sort using a bin-sort method.
 */
class BinSorter : public ParallelObject
{
  // the type of iterator we will be using is inferred from KeyType
  typedef typename std::vector<KeyType>::const_iterator IterType;

public:

  // Constructor
  explicit
  BinSorter (const Parallel::Communicator & comm,
             const std::vector<KeyType> & d);

  // The actual function which sorts the data into
  // nbins.  Currently based on the global min and
  // max which you must provide e.g. by using MPI.
  void binsort (const IdxType nbins,
                KeyType max,
                KeyType min);

  // Returns the size of bin b as an unsigned int.
  IdxType sizeof_bin (const IdxType bin) const;


private:

  const std::vector<KeyType> & data;
  std::vector<IterType>       bin_iters;   // Iterators to the bin boundaries
                                           //  in data
};



//--------------------------------------------------------------------------
template <typename KeyType, typename IdxType>
inline
IdxType BinSorter<KeyType,IdxType>::sizeof_bin (const IdxType bin) const
{
  libmesh_assert_less ((bin+1), bin_iters.size());

  // The size of the bin is defined by the distance between
  // its bounding iterators
  return cast_int<IdxType>
    (std::distance (bin_iters[bin], bin_iters[bin+1]));
}

}

} // namespace libMesh

#endif // LIBMESH_PARALLEL_BIN_SORTER_H

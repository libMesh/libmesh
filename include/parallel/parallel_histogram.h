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


#ifndef LIBMESH_PARALLEL_HISTOGRAM_H
#define LIBMESH_PARALLEL_HISTOGRAM_H

// This class contains all the functionality for bin sorting
// Templated on the type of keys you will be sorting and the
// type of iterator you will be using.


// C++ includes
#include "libmesh/libmesh_common.h" // for libmesh_assert()
#include "libmesh/parallel_object.h"

// Local includes
#include <vector>
#include <iterator>

namespace libMesh
{

namespace Parallel {

/**
 * Defines a histogram to be used in parallel in conjuction with
 * a \p BinSorter.
 */
template <typename KeyType, typename IdxType=unsigned int>
class Histogram : public ParallelObject
{
  // The type of iterator we will be using is inferred from KeyType
  typedef typename std::vector<KeyType>::const_iterator IterType;

public:

  // Constructor
  explicit
  Histogram (const Parallel::Communicator & comm,
             const std::vector<KeyType> & d);

  // The actual function which sorts the data into
  // nbins.  Currently based on the global min and
  // max which you must provide e.g. by using MPI.
  void make_histogram (const IdxType nbins,
                       KeyType max,
                       KeyType min);

  // Build the histogram across all processors and store the
  // result in the input vector \p hist
  void build_histogram ();

  // Return the raw histogram data to the user
  const std::vector<IdxType> & get_histogram() const;

  // The number of bins in the histogram
  IdxType n_bins () const;

  // Returns the size of local bin b
  IdxType local_bin_size (const IdxType bin) const;

  // Returns the size of global bin b
  // Requires that the user first call \p build_histogram()
  IdxType global_bin_size (const IdxType bin) const;

  // Returns the lower boundary of bin \p bin
  double lower_bound (const IdxType bin) const;

  // Returns the upper boundary of bin \p bin
  double upper_bound (const IdxType bin) const;


private:


  const std::vector<KeyType> & data;
  std::vector<IdxType> hist;        // The actual histogram
  std::vector<double>  bin_bounds;  // The boundary values of each bin
  std::vector<IterType> bin_iters;  // Iterators to the bin boundaries in data
};




//--------------------------------------------------------------------------
template <typename KeyType, typename IdxType>
inline
const std::vector<IdxType> & Histogram<KeyType,IdxType>::get_histogram () const
{
  return hist;
}



template <typename KeyType, typename IdxType>
inline
IdxType Histogram<KeyType,IdxType>::n_bins () const
{
  if (bin_iters.empty())
    return 0;

  return cast_int<IdxType>(bin_iters.size()-1);
}



template <typename KeyType, typename IdxType>
inline
IdxType Histogram<KeyType,IdxType>::local_bin_size (const IdxType bin) const
{
  libmesh_assert_less ((bin+1), bin_iters.size());

  // The number of entries in the bin (locally)
  return cast_int<IdxType>
    (std::distance (bin_iters[bin], bin_iters[bin+1]));
}



template <typename KeyType, typename IdxType>
inline
IdxType Histogram<KeyType,IdxType>::global_bin_size (const IdxType bin) const
{
  libmesh_assert_less (bin, hist.size());

  // The number of entries in the bin (globally)
  return hist[bin];
}



template <typename KeyType, typename IdxType>
inline
double Histogram<KeyType,IdxType>::lower_bound (const IdxType bin) const
{
  libmesh_assert_less ((bin+1), bin_bounds.size());

  return bin_bounds[bin];
}



template <typename KeyType, typename IdxType>
inline
double Histogram<KeyType,IdxType>::upper_bound (const IdxType bin) const
{
  libmesh_assert_less ((bin+1), bin_bounds.size());

  return bin_bounds[bin+1];
}

}

} // namespace libMesh

#endif // LIBMESH_PARALLEL_HISTOGRAM_H

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


// C++ includes
#include <algorithm>  // std::lower_bound

// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/parallel.h"
#include "libmesh/parallel_bin_sorter.h"
#include "libmesh/parallel_hilbert.h"
#include "libmesh/parallel_histogram.h"
#include "libmesh/parallel_conversion_utils.h"

namespace libMesh
{



namespace Parallel {

template <typename KeyType, typename IdxType>
BinSorter<KeyType,IdxType>::BinSorter (const Parallel::Communicator & comm_in,
                                       const std::vector<KeyType> & d) :
  ParallelObject(comm_in),
  data(d)
{
  // Assume (& libmesh_assert) we are working with a sorted range

  // Ah...  is_sorted is an STL extension!
  //libmesh_assert (std::is_sorted (data.begin(), data.end()));

  // Home-grown is_sorted
  libmesh_assert (Parallel::Utils::is_sorted (data));
}



template <typename KeyType, typename IdxType>
void BinSorter<KeyType,IdxType>::binsort (const IdxType nbins,
                                          KeyType max,
                                          KeyType min)
{
  libmesh_assert_less (min, max);

  // Build a histogram in parallel from our data.
  // Use this to create quasi-uniform bins.
  Parallel::Histogram<KeyType,IdxType> phist (this->comm(), data);
  phist.make_histogram (nbins*50, max, min);
  phist.build_histogram ();

  const std::vector<IdxType> & histogram =
    phist.get_histogram();


  // Now we will locate the bin boundaries so
  // that each bin is roughly equal size
  {
    // Find the total size of the data set
    IdxType local_data_size = cast_int<IdxType>(data.size());
    IdxType global_data_size = cast_int<IdxType>(local_data_size);

    this->comm().sum(global_data_size);

    std::vector<IdxType> target_bin_size (nbins, global_data_size / nbins);

    // Equally distribute the remainder
    for (IdxType i=0; i<(global_data_size % nbins); i++)
      ++target_bin_size[i];

    // Set the iterators corresponding to the bin boundaries
    {
      std::vector<double> bin_bounds (nbins+1);
      bin_iters.resize  (nbins+1, data.begin());

      // Set the minimum bin boundary iterator
      bin_iters[0]  = data.begin();
      bin_bounds[0] = Parallel::Utils::to_double(min);

      // The current location in the histogram
      IdxType current_histogram_bin = 0;

      // How much above (+) or below (-) we are from the
      // target size for the last bin.
      // Note that when delta is (-) we will
      // accept a slightly larger size for the next bin,
      // the goal being to keep the whole mess average
      int delta = 0;

      // Set the internal bin boundary iterators
      for (IdxType b=0; b<nbins; ++b)
        {
          // The size of bin b.  We want this to
          // be ~= target_bin_size[b]
          int current_bin_size = 0;

          // Step through the histogram until we have the
          // desired bin size
          while ((current_bin_size + histogram[current_histogram_bin] + delta) <= target_bin_size[b])
            {
              // Don't index out of the histogram!
              if ((current_histogram_bin+1) == phist.n_bins())
                break;

              current_bin_size += histogram[current_histogram_bin++];
            }

          delta += current_bin_size - target_bin_size[b];

          // Set the upper bound of the bin
          bin_bounds[b+1] = phist.upper_bound (current_histogram_bin);
          bin_iters[b+1] =
            std::lower_bound(bin_iters[b], data.end(),
                             Parallel::Utils::Convert<KeyType>::to_key_type(bin_bounds[b+1]));
        }

      // Just be sure the last boundary points to the right place
      bin_iters[nbins]  = data.end();
      // This gets destructed here anyway
      // bin_bounds[nbins] = Parallel::Utils::to_double(max);
    }
  }
}

}


// Explicitly instantiate for int, double
template class Parallel::BinSorter<int, unsigned int>;
template class Parallel::BinSorter<double, unsigned int>;
#ifdef LIBMESH_HAVE_LIBHILBERT
template class Parallel::BinSorter<Parallel::DofObjectKey, unsigned int>;
#endif

} // namespace libMesh

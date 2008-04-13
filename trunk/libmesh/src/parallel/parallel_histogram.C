// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2007  Benjamin S. Kirk, John W. Peterson
  
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
#include <iterator>   // std::distance(), std::advance()
#include <algorithm>  // std::swap
#include <iostream>   // std::cout

// Local includes
#include "parallel_histogram.h"
#ifdef HAVE_LIBHILBERT
#  include "hilbert.h"
#endif
#include "parallel.h"
#include "parallel_conversion_utils.h"



namespace Parallel {
template <typename KeyType>
Histogram<KeyType>::Histogram (const std::vector<KeyType>& d) :
  data(d)
{
  libmesh_assert (Parallel::Utils::is_sorted (data));
}



template <typename KeyType>
void Histogram<KeyType>::make_histogram (const unsigned int nbins,
					  KeyType max,
					  KeyType min)
{
  libmesh_assert (min < max);
  
  // The width of each bin.  Store this as a floating point value
  double bin_width = (Parallel::Utils::to_double(max)-
		      Parallel::Utils::to_double(min))/static_cast<double>(nbins);


  // The idea for 4 bins of size d is this:
  //
  //  0          1          2           3          4
  //  |----------|----------|-----------|----------|
  // min   0   min+d  1   min+2d  2  min+3d   3   max


  
  // Set the iterators corresponding to the boundaries
  // as defined above.  This takes nbins * O(log N) time. 
  bin_bounds.resize (nbins+1);
  bin_iters.resize  (nbins+1, data.begin());

  // Set the minimum bin boundary iterator
  bin_iters[0]  = data.begin();
  bin_bounds[0] = Parallel::Utils::to_double(min);
  
  // Set the internal bin boundary iterators
  for (unsigned int b=1; b<nbins; ++b)
    {
      bin_bounds[b] = Parallel::Utils::to_double(min) + bin_width * b;

      bin_iters[b]  = std::lower_bound (bin_iters[b-1], data.end(), 
					Parallel::Utils::to_key_type<KeyType>(bin_bounds[b]));
    }

  bin_iters[nbins]  = data.end();
  bin_bounds[nbins] = Parallel::Utils::to_double(max);
}  



template <typename KeyType>
void Histogram<KeyType>::build_histogram ()
{
  // Build a local histogram
  std::vector<unsigned int> local_hist (this->n_bins());

  for (unsigned int b=0; b<this->n_bins(); b++)
    local_hist[b] = this->local_bin_size(b);

  // Add all the local histograms to get the global histogram
  hist = local_hist;
  Parallel::sum(hist);
  
  // All done!  
}

}


// Explicitly instantiate for int, double
template class Parallel::Histogram<int>;
template class Parallel::Histogram<double>;
#ifdef HAVE_LIBHILBERT
template class Parallel::Histogram<Hilbert::HilbertIndices>;
#endif


// $Id: mesh_base.h 2378 2007-11-09 07:19:22Z roystgnr $

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
#include <assert.h>
#include <iterator>   // std::distance(), std::advance()
#include <algorithm>  // std::swap
#include <iostream>   // std::cout

// Local includes
#include "parallel_histogram.h"
#include "libmesh_common.h"
#ifdef HAVE_LIBHILBERT
#  include "hilbert.h"
#endif


//--------------------------------------------------------------------------
namespace {

  // Utility function that returns true if the vector v
  // is sorted, false otherwise.  O(N), the length of the
  // vector
  template <typename KeyType>
  bool is_sorted (const std::vector<KeyType>& v)
  {
    if (v.empty())
      return true;

    for (unsigned int i=1; i<v.size(); i++)
      if (v[i] < v[i-1])
	return false;

    return true;    
  }



  // A utility function which converts whatever KeyType is to 
  // a double for the histogram bounds
  template <typename KeyType>
  double to_double (const KeyType &k)
  {
    return static_cast<double>(k);
  }

  template <typename KeyType>
  KeyType to_key_type (const double f)
  {
    return static_cast<KeyType>(f);
  }

#ifdef HAVE_LIBHILBERT
  template <>
  double to_double (const Hilbert::BitVecType &bvt)
  {
    assert (bvt.rackCount() == 3);

    return static_cast<double>(bvt.racks()[2]);
  }

  template <>
  Hilbert::BitVecType 
  to_key_type (const double f)
  {
    Hilbert::BitVecType bvt;

    assert (bvt.rackCount() == 3);

    bvt.racks()[0] = 0;
    bvt.racks()[1] = 0;
    bvt.racks()[2] = f;

    return bvt;
  }
#endif // #ifdef HAVE_LIBHILBERT
}


namespace Parallel {
template <typename KeyType>
Histogram<KeyType>::Histogram (const std::vector<KeyType>& d) :
  data(d)
{
  assert (is_sorted (data));
}



template <typename KeyType>
void Histogram<KeyType>::make_histogram (const unsigned int nbins,
					  KeyType max,
					  KeyType min)
{
//   // If the user did not provide us with global maximum
//   // and minimum for the data we must compute them.
//   if (max == min)
//     {
//       KeyType
// 	local_min_max[2],
// 	global_min_max[2];

//       // negate the min so that we can use a single
//       // MPI call and take the max.
//       local_min_max[0] = -data.front();
//       local_min_max[1] =  data.back();

//       // Fix this!
//       assert (sizeof(KeyType) == sizeof(int));
      
//       MPI_Allreduce (&local_min_max[0],
// 		     &global_min_max[0],
// 		     2,
// 		     MPI_INT,
// 		     MPI_MAX,
// 		     libMesh::COMM_WORLD);

//       // re-negate the min
//       global_min_max[0] *= -1;

//       min = global_min_max[0];
//       max = global_min_max[1];

//       assert (min < max);
//       assert (min <= data.front());
//       assert (max >= data.back());
//     }
//   else
    {
      assert (min < max);
    }


    // The width of each bin.  Store this as a floating point value
    double bin_width = (to_double(max)-to_double(min))/static_cast<double>(nbins);


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
  bin_bounds[0] = to_double(min);
  
  // Set the internal bin boundary iterators
  for (unsigned int b=1; b<nbins; ++b)
    {
      bin_bounds[b] = to_double(min) + bin_width * b;

      bin_iters[b]  = std::lower_bound (bin_iters[b-1], data.end(), to_key_type<KeyType>(bin_bounds[b]));
    }

  bin_iters[nbins]  = data.end();
  bin_bounds[nbins] = to_double(max);
}  



template <typename KeyType>
void Histogram<KeyType>::build_histogram ()
{
  // Build a local histogram
  std::vector<unsigned int> local_hist (this->n_bins());

  for (unsigned int b=0; b<this->n_bins(); b++)
    local_hist[b] = this->local_bin_size(b);

  // Add all the local histograms to get the global histogram
  hist.resize (this->n_bins());
  
  MPI_Allreduce (&local_hist[0],
		 &hist[0],
		 this->n_bins(),
		 MPI_UNSIGNED,
		 MPI_SUM,
		 libMesh::COMM_WORLD);

  // All done!  
}

}


// Explicitly instantiate for int, double
template class Parallel::Histogram<int>;
template class Parallel::Histogram<double>;
#ifdef HAVE_LIBHILBERT
template class Parallel::Histogram<Hilbert::BitVecType>;
#endif


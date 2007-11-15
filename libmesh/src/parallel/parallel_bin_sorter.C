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
#include "libmesh_common.h"
#include "parallel_bin_sorter.h"
#include "parallel_histogram.h"
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
#endif // HAVE_LIBHILBERT
}

namespace Parallel {

template <typename KeyType>
BinSorter<KeyType>::BinSorter (const std::vector<KeyType>& d) :
  data(d)
{
  // Assume (& assert) we are working with a sorted range

  // Ah, suck...  is_sorted is an STL extension!
  //assert (std::is_sorted (data.begin(), data.end()));

  // Home-grown is_sorted
  assert (is_sorted (data));
}



template <typename KeyType>
void BinSorter<KeyType>::binsort (const unsigned int nbins,
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

  
  // Build a histogram in parallel from our data.
  // Use this to create quasi-uniform bins.
  Parallel::Histogram<KeyType> phist (data);
  phist.make_histogram (nbins*50, max, min);
  phist.build_histogram ();

  const std::vector<unsigned int>& histogram =
    phist.get_histogram();


  // Now we will locate the bin boundaries so
  // that each bin is roughly equal size
  {
    // Find the total size of the data set
    unsigned int local_data_size = data.size();
    unsigned int global_data_size;
    
    MPI_Allreduce (&local_data_size,
		   &global_data_size,
		   1,
		   MPI_UNSIGNED,
		   MPI_SUM,
		   libMesh::COMM_WORLD);
    
    // Set the target size of each bin
    std::vector<unsigned int> target_bin_size (nbins, global_data_size / nbins);
    
    // Equally distribute the remainder
    for (unsigned int i=0; i<(global_data_size % nbins); i++)
      ++target_bin_size[i];
    
    // Set the iterators corresponding to the bin boundaries
    {
      std::vector<double> bin_bounds (nbins+1);
      bin_iters.resize  (nbins+1, data.begin());
      
      // Set the minimum bin boundary iterator
      bin_iters[0]  = data.begin();
      bin_bounds[0] = to_double(min);
      
      // The current location in the histogram
      unsigned int current_histogram_bin = 0;

      // How much above (+) or below (-) we are from the
      // target size for the last bin.
      // Note that when delta is (-) we will
      // accept a slightly larger size for the next bin,
      // the goal being to keep the whole mess average
      int delta = 0;
      
      // Set the internal bin boundary iterators
      for (unsigned int b=0; b<nbins; ++b)
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
	  bin_iters[b+1]  = std::lower_bound(bin_iters[b], data.end(), to_key_type<KeyType>(bin_bounds[b+1]));
	}

      // Just be sure the last boundaries point to the right place
      bin_iters[nbins]  = data.end();
      bin_bounds[nbins] = to_double(max);
    }
  }
}
  
}


// Explicitly instantiate for int, double
template class Parallel::BinSorter<int>;
template class Parallel::BinSorter<double>;
#ifdef HAVE_LIBHILBERT
template class Parallel::BinSorter<Hilbert::BitVecType>;
#endif


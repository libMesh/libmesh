// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2008 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner
  
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


#ifndef __parallel_histogram_h__
#define __parallel_histogram_h__

// This class contains all the functionality for bin sorting
// Templated on the type of keys you will be sorting and the
// type of iterator you will be using.


#include <vector>
#include <iterator>

#include "libmesh_common.h" // for libmesh_assert()

namespace libMesh
{

namespace Parallel {

  /**
   * Defines a histogram to be used in parallel in conjuction with 
   * a \p BinSorter.
   */
template <typename KeyType>
class Histogram
{
  // The type of iterator we will be using is inferred from KeyType
  typedef typename std::vector<KeyType>::const_iterator IterType;
  
public:

  // Constructor
  Histogram (const std::vector<KeyType>& d);
  
  // The actual function which sorts the data into
  // nbins.  Currently based on the global min and
  // max which you must provide e.g. by using MPI.
  void make_histogram (const unsigned int nbins,
		       KeyType max,
		       KeyType min);

  // Build the histogram across all processors and store the
  // result in the input vector \p hist
  void build_histogram ();

  // Return the raw histogram data to the user
  const std::vector<unsigned int>& get_histogram() const;
  
  // The number of bins in the histogram
  unsigned int n_bins () const;
  
  // Returns the size of local bin b as an unsigned int.
  unsigned int local_bin_size (const unsigned int bin) const;  
  
  // Returns the size of global bin b as an unsigned int.
  // Requires that the user first call \p build_histogram()
  unsigned int global_bin_size (const unsigned int bin) const;  

  // Returns the lower boundary of bin \p bin
  double lower_bound (const unsigned int bin) const;

  // Returns the upper boundary of bin \p bin
  double upper_bound (const unsigned int bin) const;

  
private:

  
  const std::vector<KeyType>& data;
  std::vector<unsigned int>   hist;        // The actual histogram
  std::vector<double>         bin_bounds;  // The boundary values of each bin
  std::vector<IterType>       bin_iters;   // Iterators to the bin boundaries
                                           //  in data
};




//--------------------------------------------------------------------------
template <typename KeyType>
inline
const std::vector<unsigned int>& Histogram<KeyType>::get_histogram () const
{
  return hist;
}



template <typename KeyType>
inline
unsigned int Histogram<KeyType>::n_bins () const
{
  if (bin_iters.empty())
    return 0;

  return (bin_iters.size()-1);
}



template <typename KeyType>
inline
unsigned int Histogram<KeyType>::local_bin_size (const unsigned int bin) const
{
  libmesh_assert ((bin+1) < bin_iters.size());

  // The number of entries in the bin (locally)
  return std::distance (bin_iters[bin], bin_iters[bin+1]);
}



template <typename KeyType>
inline
unsigned int Histogram<KeyType>::global_bin_size (const unsigned int bin) const
{
  libmesh_assert (bin < hist.size());

  // The number of entries in the bin (globally)
  return hist[bin];
}



template <typename KeyType>
inline
double Histogram<KeyType>::lower_bound (const unsigned int bin) const
{
  libmesh_assert ((bin+1) < bin_bounds.size());

  return bin_bounds[bin];
}



template <typename KeyType>
inline
double Histogram<KeyType>::upper_bound (const unsigned int bin) const
{
  libmesh_assert ((bin+1) < bin_bounds.size());

  return bin_bounds[bin+1];
}

}

} // namespace libMesh

#endif


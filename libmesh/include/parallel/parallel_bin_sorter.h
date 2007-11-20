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


#ifndef __parallel_bin_sorter_h__
#define __parallel_bin_sorter_h__

// This class contains all the functionality for bin sorting
// Templated on the type of keys you will be sorting and the
// type of iterator you will be using.


#include <assert.h>
#include <vector>
#include <iterator>

namespace Parallel {

template <typename KeyType>
  /**
   * Perform a parallel sort using a bin-sort method.
   */
class BinSorter
{
  // The type of iterator we will be using is inferred from KeyType
  typedef typename std::vector<KeyType>::const_iterator IterType;
  
public:

  // Constructor
  BinSorter (const std::vector<KeyType>& d);
  
  // The actual function which sorts the data into
  // nbins.  Currently based on the global min and
  // max which you must provide e.g. by using MPI.
  void binsort (const unsigned int nbins,
		KeyType max,
		KeyType min);
  
  // Returns the size of bin b as an unsigned int.
  unsigned int sizeof_bin (const unsigned int bin) const;
  
    
private:
  
  const std::vector<KeyType>& data;
  std::vector<IterType>       bin_iters;   // Iterators to the bin boundaries
                                           //  in data
};



//--------------------------------------------------------------------------
template <typename KeyType>
inline
unsigned int BinSorter<KeyType>::sizeof_bin (const unsigned int bin) const
{
  assert ((bin+1) < bin_iters.size());

  // The size of the bin is defined by the distance between
  // its bounding iterators
  return std::distance (bin_iters[bin], bin_iters[bin+1]);
}

}
#endif


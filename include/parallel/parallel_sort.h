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


#ifndef LIBMESH_PARALLEL_SORT_H
#define LIBMESH_PARALLEL_SORT_H

// Local Includes
#include "libmesh/parallel.h"
#include "libmesh/libmesh_common.h"
#include "libmesh/parallel_object.h"

// C++ Includes
#include <vector>

namespace libMesh
{


namespace Parallel
{
/**
 * The parallel sorting method is templated on the
 * type of data which is to be sorted.  It may later
 * be templated on other things if we are ambitious.
 * This class knows about MPI, and knows how many
 * processors there are.  It is responsible for
 * transmitting data between the processors and
 * ensuring that the data is properly sorted between
 * all the processors.  We assume that a Sort
 * is instantiated on all processors.
 */
template <typename KeyType, typename IdxType=unsigned int>
class Sort : public ParallelObject
{
public:
  /**
   * Constructor takes the number of processors,
   * the processor id, and a reference to a vector of data
   * to be sorted.  This vector is sorted by
   * the constructor, therefore, construction of
   * a Sort object takes O(nlogn) time,
   * where n is the length of the vector.
   */
  Sort (const Parallel::Communicator & comm,
        std::vector<KeyType> & d);


  /**
   * This is the only method which needs to be
   * called by the user.  Its only responsibility
   * is to call three private methods in the correct
   * order.
   */
  void sort();

  /**
   * Return a constant reference to _my_bin.  This allows
   * us to do things like check if sorting was successful
   * by printing _my_bin.
   */
  const std::vector<KeyType> & bin();

private:

  /**
   * The number of processors to work with.
   */
  const processor_id_type _n_procs;

  /**
   * The identity of this processor.
   */
  const processor_id_type _proc_id;

  /**
   * Flag which lets you know if sorting is complete
   */
  bool _bin_is_sorted;

  /**
   * The raw, unsorted data which will need to
   * be sorted (in parallel) across all
   * processors.
   */
  std::vector<KeyType> & _data;

  /**
   * Vector which holds the size of each
   * bin on this processor.  It has
   * size equal to _n_procs.
   */
  std::vector<IdxType> _local_bin_sizes;

  /**
   * The bin which will eventually be held
   * by this processor.  It may be shorter or
   * longer than _data.  It will be dynamically
   * resized when it is needed.
   */
  std::vector<KeyType> _my_bin;

  /**
   * Sorts the local data into bins across all processors.
   * Right now it constructs a BenSorter<KeyType> object.
   * In the future this could be a template parameter.
   */
  void binsort ();

  /**
   * Communicates the bins from each processor to the
   * appropriate processor.  By the time this function
   * is finished, each processor will hold only its
   * own bin(s).
   */
  void communicate_bins();

  /**
   * After all the bins have been communicated, we can
   * sort our local bin.  This is nothing more than a
   * call to std::sort
   */
  void sort_local_bin();

};
}

} // namespace libMesh

#endif // LIBMESH_PARALLEL_SORT_H

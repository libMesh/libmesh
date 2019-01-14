// The libMesh Finite Element Library.
// Copyright (C) 2002-2019 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


// System Includes
#include <algorithm>
#include <iostream>

// Local Includes
#include "libmesh/parallel_sort.h"

#include "libmesh/libmesh_common.h"
#include "libmesh/parallel.h"
#include "libmesh/parallel_bin_sorter.h"
#include "libmesh/parallel_hilbert.h"
#include "libmesh/parallel_sync.h"

namespace libMesh
{


namespace Parallel {

// The Constructor sorts the local data using
// std::sort().  Therefore, the construction of
// a Parallel::Sort object takes O(n log n) time,
// where n is the length of _data.
template <typename KeyType, typename IdxType>
Sort<KeyType,IdxType>::Sort(const Parallel::Communicator & comm_in,
                            std::vector<KeyType> & d) :
  ParallelObject(comm_in),
  _n_procs(cast_int<processor_id_type>(comm_in.size())),
  _proc_id(cast_int<processor_id_type>(comm_in.rank())),
  _bin_is_sorted(false),
  _data(d)
{
  std::sort(_data.begin(), _data.end());

  // Allocate storage
  _local_bin_sizes.resize(_n_procs);
}



template <typename KeyType, typename IdxType>
void Sort<KeyType,IdxType>::sort()
{
  // Find the global data size.  The sorting
  // algorithms assume they have a range to
  // work with, so catch the degenerate cases here
  IdxType global_data_size = cast_int<IdxType>(_data.size());

  this->comm().sum (global_data_size);

  if (global_data_size < 2)
    {
      // the entire global range is either empty
      // or contains only one element
      _my_bin = _data;

      this->comm().allgather (static_cast<IdxType>(_my_bin.size()),
                              _local_bin_sizes);
    }
  else
    {
      if (this->n_processors() > 1)
        {
          this->binsort();
          this->communicate_bins();
        }
      else
        _my_bin = _data;

      this->sort_local_bin();
    }

  // Set sorted flag to true
  _bin_is_sorted = true;
}



template <typename KeyType, typename IdxType>
void Sort<KeyType,IdxType>::binsort()
{
  // Find the global min and max from all the
  // processors.
  std::vector<KeyType> global_min_max(2);

  // Insert the local min and max for this processor
  global_min_max[0] = -_data.front();
  global_min_max[1] =  _data.back();

  // Communicate to determine the global
  // min and max for all processors.
  this->comm().max(global_min_max);

  // Multiply the min by -1 to obtain the true min
  global_min_max[0] *= -1;

  // Bin-Sort based on the global min and max
  Parallel::BinSorter<KeyType> bs(this->comm(), _data);
  bs.binsort(_n_procs, global_min_max[1], global_min_max[0]);

  // Now save the local bin sizes in a vector so
  // we don't have to keep around the BinSorter.
  for (processor_id_type i=0; i<_n_procs; ++i)
    _local_bin_sizes[i] = bs.sizeof_bin(i);
}



#if defined(LIBMESH_HAVE_LIBHILBERT) && defined(LIBMESH_HAVE_MPI)
// Full specialization for HilbertIndices, there is a fair amount of
// code duplication here that could potentially be consolidated with the
// above method
template <>
void Sort<Parallel::DofObjectKey,unsigned int>::binsort()
{
  // Find the global min and max from all the
  // processors.  Do this using MPI_Allreduce.
  Parallel::DofObjectKey
    local_min,  local_max,
    global_min, global_max;

  if (_data.empty())
    {
#ifdef LIBMESH_ENABLE_UNIQUE_ID
      local_min.first.rack0 = local_min.first.rack1 = local_min.first.rack2 = static_cast<Hilbert::inttype>(-1);
      local_min.second = std::numeric_limits<unique_id_type>::max();
      local_max.first.rack0 = local_max.first.rack1 = local_max.first.rack2 = 0;
      local_max.second = 0;
#else
      local_min.rack0 = local_min.rack1 = local_min.rack2 = static_cast<Hilbert::inttype>(-1);
      local_max.rack0 = local_max.rack1 = local_max.rack2 = 0;
#endif
    }
  else
    {
      local_min = _data.front();
      local_max = _data.back();
    }

  MPI_Op hilbert_max, hilbert_min;

  MPI_Op_create       ((MPI_User_function*)dofobjectkey_max_op, true, &hilbert_max);
  MPI_Op_create       ((MPI_User_function*)dofobjectkey_min_op, true, &hilbert_min);

  // Communicate to determine the global
  // min and max for all processors.
  MPI_Allreduce(&local_min,
                &global_min,
                1,
                Parallel::StandardType<Parallel::DofObjectKey>(&local_min),
                hilbert_min,
                this->comm().get());

  MPI_Allreduce(&local_max,
                &global_max,
                1,
                Parallel::StandardType<Parallel::DofObjectKey>(&local_max),
                hilbert_max,
                this->comm().get());

  MPI_Op_free   (&hilbert_max);
  MPI_Op_free   (&hilbert_min);

  // Bin-Sort based on the global min and max
  Parallel::BinSorter<Parallel::DofObjectKey> bs(this->comm(),_data);
  bs.binsort(_n_procs, global_max, global_min);

  // Now save the local bin sizes in a vector so
  // we don't have to keep around the BinSorter.
  for (processor_id_type i=0; i<_n_procs; ++i)
    _local_bin_sizes[i] = bs.sizeof_bin(i);
}

#endif // #ifdef LIBMESH_HAVE_LIBHILBERT


template <typename KeyType, typename IdxType>
void Sort<KeyType,IdxType>::communicate_bins()
{
#ifdef LIBMESH_HAVE_MPI
  // Find each section of our data to send
  IdxType local_offset = 0;
  std::map<processor_id_type, std::vector<KeyType> > pushed_keys, received_keys;

  for (processor_id_type i=0; i != _n_procs; ++i)
    {
      IdxType next_offset = local_offset + _local_bin_sizes[i];
        if (_local_bin_sizes[i])
          {
            auto begin = _data.begin() + local_offset;
            auto end = _data.begin() + next_offset;
            pushed_keys[i].assign(begin, end);
          }

      local_offset = next_offset;
    }

  auto keys_action_functor =
    [& received_keys]
    (processor_id_type pid,
     const std::vector<KeyType> & keys)
    {
      received_keys[pid] = keys;
    };

  Parallel::push_parallel_vector_data
    (this->comm(), pushed_keys, keys_action_functor);

  std::size_t my_bin_size = 0;
  for (auto & p : received_keys)
    my_bin_size += p.second.size();

  _my_bin.clear();
  _my_bin.reserve(my_bin_size);

  for (auto & p : received_keys)
    _my_bin.insert(_my_bin.end(), p.second.begin(), p.second.end());

#ifdef DEBUG
  std::vector<IdxType> global_bin_sizes = _local_bin_sizes;

  this->comm().sum(global_bin_sizes);

  libmesh_assert_equal_to
    (global_bin_sizes[this->processor_id()], _my_bin.size());
#endif

#endif // LIBMESH_HAVE_MPI
}



template <typename KeyType, typename IdxType>
void Sort<KeyType,IdxType>::sort_local_bin()
{
  std::sort(_my_bin.begin(), _my_bin.end());
}



template <typename KeyType, typename IdxType>
const std::vector<KeyType> & Sort<KeyType,IdxType>::bin()
{
  if (!_bin_is_sorted)
    {
      libMesh::out << "Warning! Bin is not yet sorted!" << std::endl;
    }

  return _my_bin;
}

}



// Explicitly instantiate for int, double
template class Parallel::Sort<int, unsigned int>;
template class Parallel::Sort<double, unsigned int>;
#if defined(LIBMESH_HAVE_LIBHILBERT) && defined(LIBMESH_HAVE_MPI)
template class Parallel::Sort<Parallel::DofObjectKey, unsigned int>;
#endif

} // namespace libMesh

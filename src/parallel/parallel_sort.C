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


// System Includes
#include <algorithm>
#include <iostream>

// Local Includes
#include "libmesh/libmesh_common.h"
#include "libmesh/parallel.h"
#include "libmesh/parallel_hilbert.h"
#include "libmesh/parallel_sort.h"
#include "libmesh/parallel_bin_sorter.h"

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
  // Create storage for the global bin sizes.  This
  // is the number of keys which will be held in
  // each bin over all processors.
  std::vector<IdxType> global_bin_sizes = _local_bin_sizes;

  // Sum to find the total number of entries in each bin.
  this->comm().sum(global_bin_sizes);

  // Create a vector to temporarily hold the results of MPI_Gatherv
  // calls.  The vector dest  may be saved away to _my_bin depending on which
  // processor is being MPI_Gatherv'd.
  std::vector<KeyType> dest;

  IdxType local_offset = 0;

  for (processor_id_type i=0; i<_n_procs; ++i)
    {
      // Vector to receive the total bin size for each
      // processor.  Processor i's bin size will be
      // held in proc_bin_size[i]
      std::vector<int> proc_bin_size;

      // Find the number of contributions coming from each
      // processor for this bin.  Note: allgather combines
      // the MPI_Gather and MPI_Bcast operations into one.
      this->comm().allgather(static_cast<int>(_local_bin_sizes[i]),
                             proc_bin_size);

      // Compute the offsets into my_bin for each processor's
      // portion of the bin.  These are basically partial sums
      // of the proc_bin_size vector.
      std::vector<int> displacements(_n_procs);
      for (processor_id_type j=1; j<_n_procs; ++j)
        displacements[j] = proc_bin_size[j-1] + displacements[j-1];

      // Resize the destination buffer
      dest.resize (global_bin_sizes[i]);

      // Points to the beginning of the bin to be sent
      void * sendbuf = (_data.size() > local_offset) ? &_data[local_offset] : libmesh_nullptr;

      // Enough storage to hold all bin contributions
      void * recvbuf = (dest.empty()) ? libmesh_nullptr : &dest[0];

      // If the sendbuf is NULL, make sure we aren't claiming to send something.
      if (sendbuf == libmesh_nullptr && _local_bin_sizes[i] != 0)
        libmesh_error_msg("Error: invalid MPI_Gatherv call constructed!");

      KeyType example;

      MPI_Gatherv(sendbuf,
                  _local_bin_sizes[i],                       // How much data is in the bin being sent.
                  Parallel::StandardType<KeyType>(&example), // The data type we are sorting
                  recvbuf,
                  &proc_bin_size[0],          // How much is to be received from each processor
                  &displacements[0],          // Offsets into the receive buffer
                  Parallel::StandardType<KeyType>(&example), // The data type we are sorting
                  i,                                         // The root process (we do this once for each proc)
                  this->comm().get());

      // Copy the destination buffer if it
      // corresponds to the bin for this processor
      if (i == _proc_id)
        _my_bin = dest;

      // Increment the local offset counter
      local_offset += _local_bin_sizes[i];
    }
#endif // LIBMESH_HAVE_MPI
}



#if defined(LIBMESH_HAVE_LIBHILBERT) && defined(LIBMESH_HAVE_MPI)
// Full specialization for HilbertIndices, there is a fair amount of
// code duplication here that could potentially be consolidated with the
// above method
template <>
void Sort<Parallel::DofObjectKey,unsigned int>::communicate_bins()
{
  // Create storage for the global bin sizes.  This
  // is the number of keys which will be held in
  // each bin over all processors.
  std::vector<unsigned int> global_bin_sizes(_n_procs);

  libmesh_assert_equal_to (_local_bin_sizes.size(), global_bin_sizes.size());

  // Sum to find the total number of entries in each bin.
  // This is stored in global_bin_sizes.  Note, we
  // explicitly know that we are communicating MPI_UNSIGNED's here.
  MPI_Allreduce(&_local_bin_sizes[0],
                &global_bin_sizes[0],
                _n_procs,
                MPI_UNSIGNED,
                MPI_SUM,
                this->comm().get());

  // Create a vector to temporarily hold the results of MPI_Gatherv
  // calls.  The vector dest  may be saved away to _my_bin depending on which
  // processor is being MPI_Gatherv'd.
  std::vector<Parallel::DofObjectKey> dest;

  unsigned int local_offset = 0;

  for (unsigned int i=0; i<_n_procs; ++i)
    {
      // Vector to receive the total bin size for each
      // processor.  Processor i's bin size will be
      // held in proc_bin_size[i]
      std::vector<int> proc_bin_size(_n_procs);

      // Find the number of contributions coming from each
      // processor for this bin.  Note: Allgather combines
      // the MPI_Gather and MPI_Bcast operations into one.
      // Note: Here again we know that we are communicating
      // MPI_UNSIGNED's so there is no need to check the MPI_traits.
      MPI_Allgather(&_local_bin_sizes[i], // Source: # of entries on this proc in bin i
                    1,                    // Number of items to gather
                    MPI_UNSIGNED,
                    &proc_bin_size[0],    // Destination: Total # of entries in bin i
                    1,
                    MPI_INT,
                    this->comm().get());

      // Compute the offsets into my_bin for each processor's
      // portion of the bin.  These are basically partial sums
      // of the proc_bin_size vector.
      std::vector<int> displacements(_n_procs);
      for (unsigned int j=1; j<_n_procs; ++j)
        displacements[j] = proc_bin_size[j-1] + displacements[j-1];

      // Resize the destination buffer
      dest.resize (global_bin_sizes[i]);

      // Points to the beginning of the bin to be sent
      void * sendbuf = (_data.size() > local_offset) ? &_data[local_offset] : libmesh_nullptr;

      // Enough storage to hold all bin contributions
      void * recvbuf = (dest.empty()) ? libmesh_nullptr : &dest[0];

      // If the sendbuf is NULL, make sure we aren't claiming to send something.
      if (sendbuf == libmesh_nullptr && _local_bin_sizes[i] != 0)
        libmesh_error_msg("Error: invalid MPI_Gatherv call constructed!");

      Parallel::DofObjectKey example;

      MPI_Gatherv(sendbuf,
                  _local_bin_sizes[i],      // How much data is in the bin being sent.
                  Parallel::StandardType<Parallel::DofObjectKey>(&example), // The data type we are sorting
                  recvbuf,
                  &proc_bin_size[0], // How much is to be received from each processor
                  &displacements[0], // Offsets into the receive buffer
                  Parallel::StandardType<Parallel::DofObjectKey>(&example), // The data type we are sorting
                  i,                        // The root process (we do this once for each proc)
                  this->comm().get());

      // Copy the destination buffer if it
      // corresponds to the bin for this processor
      if (i == _proc_id)
        _my_bin = dest;

      // Increment the local offset counter
      local_offset += _local_bin_sizes[i];
    }
}

#endif // #if defined(LIBMESH_HAVE_LIBHILBERT) && defined(LIBMESH_HAVE_MPI)



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

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



#include "libmesh/libmesh_common.h"

// C++ includes
#include <cstdlib> // *must* precede <cmath> for proper std:abs() on PGI, Sun Studio CC
#include <cmath> // for std::abs
#include <limits> // std::numeric_limits<T>::min()

// Local Includes
#include "libmesh/distributed_vector.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dense_subvector.h"
#include "libmesh/parallel.h"
#include "libmesh/parallel_sync.h"
#include "libmesh/tensor_tools.h"
#include "libmesh/int_range.h"

namespace libMesh
{



//--------------------------------------------------------------------------
// DistributedVector methods
template <typename T>
T DistributedVector<T>::sum () const
{
  // This function must be run on all processors at once
  parallel_object_only();

  libmesh_assert (this->initialized());
  libmesh_assert_equal_to (_values.size(), _local_size);
  libmesh_assert_equal_to ((_last_local_index - _first_local_index), _local_size);

  T local_sum = 0.;

  for (auto & val : _values)
    local_sum += val;

  this->comm().sum(local_sum);

  return local_sum;
}



template <typename T>
Real DistributedVector<T>::l1_norm () const
{
  // This function must be run on all processors at once
  parallel_object_only();

  libmesh_assert (this->initialized());
  libmesh_assert_equal_to (_values.size(), _local_size);
  libmesh_assert_equal_to ((_last_local_index - _first_local_index), _local_size);

  Real local_l1 = 0.;

  for (auto & val : _values)
    local_l1 += std::abs(val);

  this->comm().sum(local_l1);

  return local_l1;
}



template <typename T>
Real DistributedVector<T>::l2_norm () const
{
  // This function must be run on all processors at once
  parallel_object_only();

  libmesh_assert (this->initialized());
  libmesh_assert_equal_to (_values.size(), _local_size);
  libmesh_assert_equal_to ((_last_local_index - _first_local_index), _local_size);

  Real local_l2 = 0.;

  for (auto & val : _values)
    local_l2 += TensorTools::norm_sq(val);

  this->comm().sum(local_l2);

  return std::sqrt(local_l2);
}



template <typename T>
Real DistributedVector<T>::linfty_norm () const
{
  // This function must be run on all processors at once
  parallel_object_only();

  libmesh_assert (this->initialized());
  libmesh_assert_equal_to (_values.size(), _local_size);
  libmesh_assert_equal_to ((_last_local_index - _first_local_index), _local_size);

  Real local_linfty = 0.;

  for (auto & val : _values)
    local_linfty  = std::max(local_linfty,
                             static_cast<Real>(std::abs(val))
                             ); // Note we static_cast so that both
                                // types are the same, as required
                                // by std::max

  this->comm().max(local_linfty);

  return local_linfty;
}



template <typename T>
NumericVector<T> & DistributedVector<T>::operator += (const NumericVector<T> & v)
{
  libmesh_assert (this->closed());
  libmesh_assert (this->initialized());
  libmesh_assert_equal_to (_values.size(), _local_size);
  libmesh_assert_equal_to ((_last_local_index - _first_local_index), _local_size);

  add(1., v);

  return *this;
}



template <typename T>
NumericVector<T> & DistributedVector<T>::operator -= (const NumericVector<T> & v)
{
  libmesh_assert (this->closed());
  libmesh_assert (this->initialized());
  libmesh_assert_equal_to (_values.size(), _local_size);
  libmesh_assert_equal_to ((_last_local_index - _first_local_index), _local_size);

  add(-1., v);

  return *this;
}



template <typename T>
NumericVector<T> & DistributedVector<T>::operator /= (const NumericVector<T> & v)
{
  libmesh_assert_equal_to(size(), v.size());

  const DistributedVector<T> & v_vec = cast_ref<const DistributedVector<T> &>(v);

  for (auto i : index_range(_values))
    _values[i] /= v_vec._values[i];

  return *this;
}




template <typename T>
void DistributedVector<T>::reciprocal()
{
  for (auto & val : _values)
    {
      // Don't divide by zero
      libmesh_assert_not_equal_to (val, T(0));

      val = 1. / val;
    }
}




template <typename T>
void DistributedVector<T>::conjugate()
{
  // Replace values by complex conjugate
  for (auto & val : _values)
    val = libmesh_conj(val);
}





template <typename T>
void DistributedVector<T>::add (const T v)
{
  libmesh_assert (this->initialized());
  libmesh_assert_equal_to (_values.size(), _local_size);
  libmesh_assert_equal_to ((_last_local_index - _first_local_index), _local_size);

  for (auto & val : _values)
    val += v;
}



template <typename T>
void DistributedVector<T>::add (const NumericVector<T> & v)
{
  libmesh_assert (this->initialized());
  libmesh_assert_equal_to (_values.size(), _local_size);
  libmesh_assert_equal_to ((_last_local_index - _first_local_index), _local_size);

  add (1., v);
}



template <typename T>
void DistributedVector<T>::add (const T a, const NumericVector<T> & v_in)
{
  libmesh_assert (this->initialized());
  libmesh_assert_equal_to (_values.size(), _local_size);
  libmesh_assert_equal_to ((_last_local_index - _first_local_index), _local_size);

  // Make sure the NumericVector passed in is really a DistributedVector
  const DistributedVector<T> * v = cast_ptr<const DistributedVector<T> *>(&v_in);
  if (!v)
    libmesh_error_msg("Cannot add different types of NumericVectors.");

  for (auto i : index_range(_values))
    _values[i] += a * v->_values[i];
}



template <typename T>
void DistributedVector<T>::scale (const T factor)
{
  libmesh_assert (this->initialized());
  libmesh_assert_equal_to (_values.size(), _local_size);
  libmesh_assert_equal_to ((_last_local_index - _first_local_index), _local_size);

  for (auto & val : _values)
    val *= factor;
}

template <typename T>
void DistributedVector<T>::abs()
{
  libmesh_assert (this->initialized());
  libmesh_assert_equal_to ((_last_local_index - _first_local_index), _local_size);

  for (auto & val : _values)
    val = std::abs(val);
}





template <typename T>
T DistributedVector<T>::dot (const NumericVector<T> & V) const
{
  // This function must be run on all processors at once
  parallel_object_only();

  // Make sure the NumericVector passed in is really a DistributedVector
  const DistributedVector<T> * v = cast_ptr<const DistributedVector<T> *>(&V);

  // Make sure that the two vectors are distributed in the same way.
  libmesh_assert_equal_to ( this->first_local_index(), v->first_local_index() );
  libmesh_assert_equal_to ( this->last_local_index(), v->last_local_index()  );

  // The result of dotting together the local parts of the vector.
  T local_dot = 0;

  for (auto i : index_range(_values))
    local_dot += this->_values[i] * v->_values[i];

  // The local dot products are now summed via MPI
  this->comm().sum(local_dot);

  return local_dot;
}



template <typename T>
NumericVector<T> &
DistributedVector<T>::operator = (const T s)
{
  libmesh_assert (this->initialized());
  libmesh_assert_equal_to (_values.size(), _local_size);
  libmesh_assert_equal_to ((_last_local_index - _first_local_index), _local_size);

  for (auto & val : _values)
    val = s;

  return *this;
}



template <typename T>
NumericVector<T> &
DistributedVector<T>::operator = (const NumericVector<T> & v_in)
{
  // Make sure the NumericVector passed in is really a DistributedVector
  const DistributedVector<T> * v = cast_ptr<const DistributedVector<T> *>(&v_in);

  *this = *v;

  return *this;
}



template <typename T>
DistributedVector<T> &
DistributedVector<T>::operator = (const DistributedVector<T> & v)
{
  this->_is_initialized    = v._is_initialized;
  this->_is_closed         = v._is_closed;

  _global_size       = v._global_size;
  _local_size        = v._local_size;
  _first_local_index = v._first_local_index;
  _last_local_index  = v._last_local_index;

  if (v.local_size() == this->local_size())
    _values = v._values;

  else
    libmesh_error_msg("v.local_size() = " << v.local_size() << " must be equal to this->local_size() = " << this->local_size());

  return *this;
}



template <typename T>
NumericVector<T> &
DistributedVector<T>::operator = (const std::vector<T> & v)
{
  libmesh_assert (this->initialized());
  libmesh_assert_equal_to (_values.size(), _local_size);
  libmesh_assert_equal_to ((_last_local_index - _first_local_index), _local_size);

  if (v.size() == local_size())
    _values = v;

  else if (v.size() == size())
    for (std::size_t i=first_local_index(); i<last_local_index(); i++)
      _values[i-first_local_index()] = v[i];

  else
    libmesh_error_msg("Incompatible sizes in DistributedVector::operator=");

  return *this;
}



template <typename T>
void DistributedVector<T>::localize (NumericVector<T> & v_local_in) const

{
  libmesh_assert (this->initialized());
  libmesh_assert_equal_to (_values.size(), _local_size);
  libmesh_assert_equal_to ((_last_local_index - _first_local_index), _local_size);

  DistributedVector<T> * v_local = cast_ptr<DistributedVector<T> *>(&v_local_in);

  v_local->_first_local_index = 0;

  v_local->_global_size =
    v_local->_local_size =
    v_local->_last_local_index = size();

  v_local->_is_initialized =
    v_local->_is_closed = true;

  // Call localize on the vector's values.  This will help
  // prevent code duplication
  localize (v_local->_values);

#ifndef LIBMESH_HAVE_MPI

  libmesh_assert_equal_to (local_size(), size());

#endif
}



template <typename T>
void DistributedVector<T>::localize (NumericVector<T> & v_local_in,
                                     const std::vector<numeric_index_type> &) const
{
  libmesh_assert (this->initialized());
  libmesh_assert_equal_to (_values.size(), _local_size);
  libmesh_assert_equal_to ((_last_local_index - _first_local_index), _local_size);

  // TODO: We don't yet support the send list; this is inefficient:
  localize (v_local_in);
}



template <typename T>
void DistributedVector<T>::localize (std::vector<T> & v_local,
                                     const std::vector<numeric_index_type> & indices) const
{
  // Resize v_local so there is enough room to hold all the local values.
  v_local.resize(indices.size());

  // We need to know who has the values we want, so get everyone's _local_size
  std::vector<numeric_index_type> local_sizes;
  this->comm().allgather (_local_size, local_sizes);

  // Make a vector of partial sums of local sizes
  std::vector<numeric_index_type> local_size_sums(this->n_processors());
  local_size_sums[0] = local_sizes[0];
  for (auto i : IntRange<numeric_index_type>(1, local_sizes.size()))
    local_size_sums[i] = local_size_sums[i-1] + local_sizes[i];

  // We now fill in 'requested_ids' based on the indices.  Also keep
  // track of the local index (in the indices vector) for each of
  // these, since we need that when unpacking.
  std::map<processor_id_type, std::vector<numeric_index_type>>
    requested_ids, local_requested_ids;

  // We'll use this typedef a couple of times below.
  typedef typename std::vector<numeric_index_type>::iterator iter_t;

  // For each index in indices, determine which processor it is on.
  // This is an O(N*log(p)) algorithm that uses std::upper_bound().
  // Note: upper_bound() returns an iterator to the first entry which is
  // greater than the given value.
  for (auto i : index_range(indices))
    {
      iter_t ub = std::upper_bound(local_size_sums.begin(),
                                   local_size_sums.end(),
                                   indices[i]);

      processor_id_type on_proc = cast_int<processor_id_type>
        (std::distance(local_size_sums.begin(), ub));

      requested_ids[on_proc].push_back(indices[i]);
      local_requested_ids[on_proc].push_back(i);
    }

  auto gather_functor =
    [this]
    (processor_id_type, const std::vector<dof_id_type> & ids,
     std::vector<T> & values)
    {
      // The first send/receive we did was for indices, the second one will be
      // for corresponding floating point values, so create storage for that now...
      const std::size_t ids_size = ids.size();
      values.resize(ids_size);

      for (std::size_t i=0; i != ids_size; i++)
        {
          // The index of the requested value
          const numeric_index_type requested_index = ids[i];

          // Transform into local numbering, and get requested value.
          values[i] = _values[requested_index - _first_local_index];
        }
    };

  auto action_functor =
    [& v_local, & local_requested_ids]
    (processor_id_type pid,
     const std::vector<dof_id_type> &,
     const std::vector<T> & values)
    {
      // Now write the received values to the appropriate place(s) in v_local
      for (auto i : index_range(values))
        {
          libmesh_assert(local_requested_ids.count(pid));
          libmesh_assert_less(i, local_requested_ids[pid].size());

          // Get the index in v_local where this value needs to be inserted.
          const numeric_index_type local_requested_index =
            local_requested_ids[pid][i];

          // Actually set the value in v_local
          v_local[local_requested_index] = values[i];
        }
    };

  const T * ex = nullptr;
  Parallel::pull_parallel_vector_data
    (this->comm(), requested_ids, gather_functor, action_functor, ex);
}



template <typename T>
void DistributedVector<T>::localize (const numeric_index_type first_local_idx,
                                     const numeric_index_type last_local_idx,
                                     const std::vector<numeric_index_type> & send_list)
{
  // Only good for serial vectors
  libmesh_assert_equal_to (this->size(), this->local_size());
  libmesh_assert_greater (last_local_idx, first_local_idx);
  libmesh_assert_less_equal (send_list.size(), this->size());
  libmesh_assert_less (last_local_idx, this->size());

  const numeric_index_type my_size       = this->size();
  const numeric_index_type my_local_size = (last_local_idx - first_local_idx + 1);

  // Don't bother for serial cases
  if ((first_local_idx == 0) &&
      (my_local_size == my_size))
    return;


  // Build a parallel vector, initialize it with the local
  // parts of (*this)
  DistributedVector<T> parallel_vec(this->comm());

  parallel_vec.init (my_size, my_local_size, true, PARALLEL);

  // Copy part of *this into the parallel_vec
  for (numeric_index_type i=first_local_idx; i<=last_local_idx; i++)
    parallel_vec._values[i-first_local_idx] = _values[i];

  // localize like normal
  parallel_vec.localize (*this, send_list);
}



template <typename T>
void DistributedVector<T>::localize (std::vector<T> & v_local) const
{
  // This function must be run on all processors at once
  parallel_object_only();

  libmesh_assert (this->initialized());
  libmesh_assert_equal_to (_values.size(), _local_size);
  libmesh_assert_equal_to ((_last_local_index - _first_local_index), _local_size);

  v_local = this->_values;

  this->comm().allgather (v_local);

#ifndef LIBMESH_HAVE_MPI
  libmesh_assert_equal_to (local_size(), size());
#endif
}



template <typename T>
void DistributedVector<T>::localize_to_one (std::vector<T> & v_local,
                                            const processor_id_type pid) const
{
  // This function must be run on all processors at once
  parallel_object_only();

  libmesh_assert (this->initialized());
  libmesh_assert_equal_to (_values.size(), _local_size);
  libmesh_assert_equal_to ((_last_local_index - _first_local_index), _local_size);

  v_local = this->_values;

  this->comm().gather (pid, v_local);

#ifndef LIBMESH_HAVE_MPI
  libmesh_assert_equal_to (local_size(), size());
#endif
}



template <typename T>
void DistributedVector<T>::pointwise_mult (const NumericVector<T> &,
                                           const NumericVector<T> &)
//void DistributedVector<T>::pointwise_mult (const NumericVector<T> & vec1,
//   const NumericVector<T> & vec2)
{
  libmesh_not_implemented();
}



//--------------------------------------------------------------
// Explicit instantiations
template class DistributedVector<Number>;

} // namespace libMesh

// The libMesh Finite Element Library.
// Copyright (C) 2002-2012 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/tensor_tools.h"

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

  for (numeric_index_type i=0; i<local_size(); i++)
    local_sum += _values[i];

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

  double local_l1 = 0.;

  for (numeric_index_type i=0; i<local_size(); i++)
    local_l1 += std::abs(_values[i]);

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

  double local_l2 = 0.;

  for (numeric_index_type i=0; i<local_size(); i++)
    local_l2 += TensorTools::norm_sq(_values[i]);

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

  for (numeric_index_type i=0; i<local_size(); i++)
    local_linfty  = std::max(local_linfty,
			     static_cast<Real>(std::abs(_values[i]))
			     ); // Note we static_cast so that both
                                // types are the same, as required
                                // by std::max

  this->comm().max(local_linfty);

  return local_linfty;
}



template <typename T>
NumericVector<T>& DistributedVector<T>::operator += (const NumericVector<T>& v)
{
  libmesh_assert (this->closed());
  libmesh_assert (this->initialized());
  libmesh_assert_equal_to (_values.size(), _local_size);
  libmesh_assert_equal_to ((_last_local_index - _first_local_index), _local_size);

  add(1., v);

  return *this;
}



template <typename T>
NumericVector<T>& DistributedVector<T>::operator -= (const NumericVector<T>& v)
{
  libmesh_assert (this->closed());
  libmesh_assert (this->initialized());
  libmesh_assert_equal_to (_values.size(), _local_size);
  libmesh_assert_equal_to ((_last_local_index - _first_local_index), _local_size);

  add(-1., v);

  return *this;
}



template <typename T>
void DistributedVector<T>::add_vector (const std::vector<T>& v,
				       const std::vector<numeric_index_type>& dof_indices)
{
  libmesh_assert (!v.empty());
  libmesh_assert_equal_to (v.size(), dof_indices.size());
  libmesh_assert (this->initialized());
  libmesh_assert_equal_to (_values.size(), _local_size);
  libmesh_assert_equal_to ((_last_local_index - _first_local_index), _local_size);

  for (std::size_t i=0; i<v.size(); i++)
    add (dof_indices[i], v[i]);
}



template <typename T>
void DistributedVector<T>::add_vector (const NumericVector<T>& V,
				       const std::vector<numeric_index_type>& dof_indices)
{
  libmesh_assert_equal_to (V.size(), dof_indices.size());
  libmesh_assert (this->initialized());
  libmesh_assert_equal_to (_values.size(), _local_size);
  libmesh_assert_equal_to ((_last_local_index - _first_local_index), _local_size);

  for (std::size_t i=0; i<V.size(); i++)
    add (dof_indices[i], V(i));
}



template <typename T>
void DistributedVector<T>::add_vector (const DenseVector<T>& V,
				       const std::vector<numeric_index_type>& dof_indices)
{
  libmesh_assert_equal_to (V.size(), dof_indices.size());
  libmesh_assert (this->initialized());
  libmesh_assert_equal_to (_values.size(), _local_size);
  libmesh_assert_equal_to ((_last_local_index - _first_local_index), _local_size);

  for (unsigned int i=0; i<V.size(); i++)
    add (dof_indices[i], V(i));
}




template <typename T>
void DistributedVector<T>::reciprocal()
{
  for (numeric_index_type i=0; i<local_size(); i++)
    {
      // Don't divide by zero
      libmesh_assert_not_equal_to (_values[i], T(0));

      _values[i] = 1. / _values[i];
    }
}




template <typename T>
void DistributedVector<T>::conjugate()
{
  for (numeric_index_type i=0; i<local_size(); i++)
    {
      // Replace values by complex conjugate
      _values[i] = libmesh_conj(_values[i]);
    }
}





template <typename T>
void DistributedVector<T>::add (const T v)
{
  libmesh_assert (this->initialized());
  libmesh_assert_equal_to (_values.size(), _local_size);
  libmesh_assert_equal_to ((_last_local_index - _first_local_index), _local_size);

  for (numeric_index_type i=0; i<local_size(); i++)
    _values[i] += v;
}



template <typename T>
void DistributedVector<T>::add (const NumericVector<T>& v)
{
  libmesh_assert (this->initialized());
  libmesh_assert_equal_to (_values.size(), _local_size);
  libmesh_assert_equal_to ((_last_local_index - _first_local_index), _local_size);

  add (1., v);
}



template <typename T>
void DistributedVector<T>::add (const T a, const NumericVector<T>& v)
{
  libmesh_assert (this->initialized());
  libmesh_assert_equal_to (_values.size(), _local_size);
  libmesh_assert_equal_to ((_last_local_index - _first_local_index), _local_size);

  add(a, v);
}



template <typename T>
void DistributedVector<T>::insert (const std::vector<T>& v,
				   const std::vector<numeric_index_type>& dof_indices)
{
  libmesh_assert (!v.empty());
  libmesh_assert_equal_to (v.size(), dof_indices.size());
  libmesh_assert (this->initialized());
  libmesh_assert_equal_to (_values.size(), _local_size);
  libmesh_assert_equal_to ((_last_local_index - _first_local_index), _local_size);

  for (std::size_t i=0; i<v.size(); i++)
    this->set (dof_indices[i], v[i]);
}



template <typename T>
void DistributedVector<T>::insert (const NumericVector<T>& V,
				   const std::vector<numeric_index_type>& dof_indices)
{
  libmesh_assert_equal_to (V.size(), dof_indices.size());
  libmesh_assert (this->initialized());
  libmesh_assert_equal_to (_values.size(), _local_size);
  libmesh_assert_equal_to ((_last_local_index - _first_local_index), _local_size);

  for (std::size_t i=0; i<V.size(); i++)
    this->set (dof_indices[i], V(i));
}



template <typename T>
void DistributedVector<T>::insert (const DenseVector<T>& V,
				   const std::vector<numeric_index_type>& dof_indices)
{
  libmesh_assert_equal_to (V.size(), dof_indices.size());
  libmesh_assert (this->initialized());
  libmesh_assert_equal_to (_values.size(), _local_size);
  libmesh_assert_equal_to ((_last_local_index - _first_local_index), _local_size);

  for (unsigned int i=0; i<V.size(); i++)
    this->set (dof_indices[i], V(i));
}



template <typename T>
void DistributedVector<T>::insert (const DenseSubVector<T>& V,
				   const std::vector<numeric_index_type>& dof_indices)
{
  libmesh_assert_equal_to (V.size(), dof_indices.size());
  libmesh_assert (this->initialized());
  libmesh_assert_equal_to (_values.size(), _local_size);
  libmesh_assert_equal_to ((_last_local_index - _first_local_index), _local_size);

  for (unsigned int i=0; i<V.size(); i++)
    this->set (dof_indices[i], V(i));
}



template <typename T>
void DistributedVector<T>::scale (const T factor)
{
  libmesh_assert (this->initialized());
  libmesh_assert_equal_to (_values.size(), _local_size);
  libmesh_assert_equal_to ((_last_local_index - _first_local_index), _local_size);

  for (std::size_t i=0; i<local_size(); i++)
    _values[i] *= factor;
}

template <typename T>
void DistributedVector<T>::abs()
{
  libmesh_assert (this->initialized());
  libmesh_assert_equal_to ((_last_local_index - _first_local_index), _local_size);

  for (std::size_t i=0; i<local_size(); i++)
    this->set(i,std::abs(_values[i]));
}





template <typename T>
T DistributedVector<T>::dot (const NumericVector<T>& V) const
{
  // This function must be run on all processors at once
  parallel_object_only();

  // Make sure the NumericVector passed in is really a DistributedVector
  const DistributedVector<T>* v = libmesh_cast_ptr<const DistributedVector<T>*>(&V);

  // Make sure that the two vectors are distributed in the same way.
  libmesh_assert_equal_to ( this->first_local_index(), v->first_local_index() );
  libmesh_assert_equal_to ( this->last_local_index(), v->last_local_index()  );

  // The result of dotting together the local parts of the vector.
  T local_dot = 0;

  for (std::size_t i=0; i<this->local_size(); i++)
    local_dot += this->_values[i] * v->_values[i];

  // The local dot products are now summed via MPI
  this->comm().sum(local_dot);

  return local_dot;
}



template <typename T>
NumericVector<T>&
DistributedVector<T>::operator = (const T s)
{
  libmesh_assert (this->initialized());
  libmesh_assert_equal_to (_values.size(), _local_size);
  libmesh_assert_equal_to ((_last_local_index - _first_local_index), _local_size);

  for (std::size_t i=0; i<local_size(); i++)
    _values[i] = s;

  return *this;
}



template <typename T>
NumericVector<T>&
DistributedVector<T>::operator = (const NumericVector<T>& v_in)
{
  // Make sure the NumericVector passed in is really a DistributedVector
  const DistributedVector<T>* v = libmesh_cast_ptr<const DistributedVector<T>*>(&v_in);

  *this = *v;

  return *this;
}



template <typename T>
DistributedVector<T>&
DistributedVector<T>::operator = (const DistributedVector<T>& v)
{
  this->_is_initialized    = v._is_initialized;
  this->_is_closed         = v._is_closed;

  _global_size       = v._global_size;
  _local_size        = v._local_size;
  _first_local_index = v._first_local_index;
  _last_local_index  = v._last_local_index;

  if (v.local_size() == this->local_size())
    {
      _values = v._values;
    }
  else
    {
      libmesh_error();
    }

  return *this;
}



template <typename T>
NumericVector<T>&
DistributedVector<T>::operator = (const std::vector<T>& v)
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
    {
      libmesh_error();
    }


  return *this;
}



template <typename T>
void DistributedVector<T>::localize (NumericVector<T>& v_local_in) const

{
  libmesh_assert (this->initialized());
  libmesh_assert_equal_to (_values.size(), _local_size);
  libmesh_assert_equal_to ((_last_local_index - _first_local_index), _local_size);

  DistributedVector<T>* v_local = libmesh_cast_ptr<DistributedVector<T>*>(&v_local_in);

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
void DistributedVector<T>::localize (NumericVector<T>& v_local_in,
				     const std::vector<numeric_index_type>&) const
{
  libmesh_assert (this->initialized());
  libmesh_assert_equal_to (_values.size(), _local_size);
  libmesh_assert_equal_to ((_last_local_index - _first_local_index), _local_size);

  // TODO: We don't yet support the send list; this is inefficient:
  localize (v_local_in);
}



template <typename T>
void DistributedVector<T>::localize (const numeric_index_type first_local_idx,
				     const numeric_index_type last_local_idx,
				     const std::vector<numeric_index_type>& send_list)
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
void DistributedVector<T>::localize (std::vector<T>& v_local) const
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
void DistributedVector<T>::localize_to_one (std::vector<T>& v_local,
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
void DistributedVector<T>::pointwise_mult (const NumericVector<T>&,
					   const NumericVector<T>&)
//void DistributedVector<T>::pointwise_mult (const NumericVector<T>& vec1,
//					   const NumericVector<T>& vec2)
{
  libmesh_not_implemented();
}



//--------------------------------------------------------------
// Explicit instantiations
template class DistributedVector<Number>;

} // namespace libMesh

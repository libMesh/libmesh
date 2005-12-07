// $Id: distributed_vector.C,v 1.27 2005-12-07 15:51:55 roystgnr Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2005  Benjamin S. Kirk, John W. Peterson
  
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



#include "libmesh_common.h"

// C++ includes
#include <cmath> // for std::abs

// Local Includes
#include "distributed_vector.h"
#include "dense_vector.h"




//--------------------------------------------------------------------------
// DistributedVector methods
template <typename T>
Real DistributedVector<T>::l1_norm () const
{
  assert (this->initialized());
  assert (_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);

  double local_l1 = 0.;

  for (unsigned int i=0; i<local_size(); i++)
    local_l1 += std::abs(_values[i]);
  
  double global_l1 = local_l1;


#ifdef HAVE_MPI

  MPI_Allreduce (&local_l1, &global_l1, 1,
		 MPI_DOUBLE, MPI_SUM, libMesh::COMM_WORLD);
  
#endif

  return global_l1;
}



template <typename T>
Real DistributedVector<T>::l2_norm () const
{
  assert (this->initialized());
  assert (_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);

  double local_l2 = 0.;
  
  for (unsigned int i=0; i<local_size(); i++)
    local_l2 += _values[i]*_values[i];
  
  double global_l2 = local_l2;


#ifdef HAVE_MPI

  MPI_Allreduce (&local_l2, &global_l2, 1,
		 MPI_DOUBLE, MPI_SUM, libMesh::COMM_WORLD);
  
#endif

  return global_l2;
}



template <typename T>
Real DistributedVector<T>::linfty_norm () const
{
  assert (this->initialized());
  assert (_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);

  double local_linfty = 0.;
  
  for (unsigned int i=0; i<local_size(); i++)
    local_linfty  = std::max(local_linfty,
			     static_cast<double>(std::abs(_values[i]))
			     ); // Note we static_cast so that both
                                // types are the same, as required
                                // by std::max
  
  double global_linfty = local_linfty;


#ifdef HAVE_MPI

  MPI_Allreduce (&local_linfty, &global_linfty, 1,
		 MPI_DOUBLE, MPI_MAX, libMesh::COMM_WORLD);
  
#endif

  return global_linfty;
}



template <typename T>
NumericVector<T>& DistributedVector<T>::operator += (const NumericVector<T>& v)
{
  assert (this->closed());
  assert (this->initialized());
  assert (_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);
  
  add(1., v);
  
  return *this;
}



template <typename T>
NumericVector<T>& DistributedVector<T>::operator -= (const NumericVector<T>& v)
{
  assert (this->closed());
  assert (this->initialized());
  assert (_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);
  
  add(-1., v);
  
  return *this;
}



template <typename T>
void DistributedVector<T>::add_vector (const std::vector<T>& v,
				       const std::vector<unsigned int>& dof_indices)
{
  assert (!v.empty());
  assert (v.size() == dof_indices.size());
  assert (this->initialized());
  assert (_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);
  
  for (unsigned int i=0; i<v.size(); i++)
    add (dof_indices[i], v[i]);
}



template <typename T>
void DistributedVector<T>::add_vector (const NumericVector<T>& V,
				       const std::vector<unsigned int>& dof_indices)
{
  assert (V.size() == dof_indices.size());
  assert (this->initialized());
  assert (_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);

  for (unsigned int i=0; i<V.size(); i++)
    add (dof_indices[i], V(i));
}



template <typename T>
void DistributedVector<T>::add_vector (const DenseVector<T>& V,
				       const std::vector<unsigned int>& dof_indices)
{
  assert (V.size() == dof_indices.size());
  assert (this->initialized());
  assert (_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);

  for (unsigned int i=0; i<V.size(); i++)
    add (dof_indices[i], V(i));
}



template <typename T>
void DistributedVector<T>::add (const T v)
{
  assert (this->initialized());
  assert (_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);

  for (unsigned int i=0; i<local_size(); i++)
    _values[i] += v;
}



template <typename T>
void DistributedVector<T>::add (const NumericVector<T>& v)
{
  assert (this->initialized());
  assert (_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);

  add (1., v);
}



template <typename T>
void DistributedVector<T>::add (const T a, const NumericVector<T>& v)
{
  assert (this->initialized());
  assert (_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);

  add(a, v);
}



template <typename T>
void DistributedVector<T>::insert (const std::vector<T>& v,
				   const std::vector<unsigned int>& dof_indices)
{
  assert (!v.empty());
  assert (v.size() == dof_indices.size());
  assert (this->initialized());
  assert (_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);
  
  for (unsigned int i=0; i<v.size(); i++)
    this->set (dof_indices[i], v[i]);
}



template <typename T>
void DistributedVector<T>::insert (const NumericVector<T>& V,
				   const std::vector<unsigned int>& dof_indices)
{
  assert (V.size() == dof_indices.size());
  assert (this->initialized());
  assert (_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);

  for (unsigned int i=0; i<V.size(); i++)
    this->set (dof_indices[i], V(i));
}



template <typename T>
void DistributedVector<T>::insert (const DenseVector<T>& V,
				   const std::vector<unsigned int>& dof_indices)
{
  assert (V.size() == dof_indices.size());
  assert (this->initialized());
  assert (_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);

  for (unsigned int i=0; i<V.size(); i++)
    this->set (dof_indices[i], V(i));
}



template <typename T>
void DistributedVector<T>::scale (const T factor)
{
  assert (this->initialized());
  assert (_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);

  for (unsigned int i=0; i<local_size(); i++)
    _values[i] *= factor;  
}



template <typename T>
Real DistributedVector<T>::dot (const NumericVector<T>& V) const
{
  // Make sure the NumericVector passed in is really a DistributedVector
  const DistributedVector<T>* v = dynamic_cast<const DistributedVector<T>*>(&V);
  assert (v != NULL);

  // Make sure that the two vectors are distributed in the same way.
  assert ( this->first_local_index() == v->first_local_index() );
  assert ( this->last_local_index()  == v->last_local_index()  );
  
  // The result of dotting together the local parts of the vector.
  Real local_dot = 0;

  for (unsigned int i=0; i<this->local_size(); i++)
    local_dot += this->_values[i] * v->_values[i];

  // The local dot products are now summed via MPI
  Real global_dot = local_dot;

#ifdef HAVE_MPI
  MPI_Allreduce(&local_dot,       // sendbuf 
		&global_dot,      // recvbuf
		1,                // count
		MPI_DOUBLE,       // MPI_Datatype
		MPI_SUM,          // MPI_Op
		libMesh::COMM_WORLD);
#endif

  return global_dot;
} 



template <typename T>
NumericVector<T>& 
DistributedVector<T>::operator = (const T s)
{
  assert (this->initialized());
  assert (_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);

  for (unsigned int i=0; i<local_size(); i++)
    _values[i] = s;
  
  return *this;
}



template <typename T>
NumericVector<T>&
DistributedVector<T>::operator = (const NumericVector<T>& v_in)
{
  const DistributedVector<T>* v = dynamic_cast<const DistributedVector<T>*>(&v_in);

  assert (v != NULL);
  
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
      error();
    }
  
  return *this;
}



template <typename T>
NumericVector<T>&
DistributedVector<T>::operator = (const std::vector<T>& v)
{
  assert (this->initialized());
  assert (_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);

  if (v.size() == local_size())
    _values = v;

  else if (v.size() == size())
    for (unsigned int i=first_local_index(); i<last_local_index(); i++)
      _values[i-first_local_index()] = v[i];

  else
    {
      error();
    }

  
  return *this;
}



template <typename T>
void DistributedVector<T>::localize (NumericVector<T>& v_local_in) const

{
  assert (this->initialized());
  assert (_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);

  DistributedVector<T>* v_local = dynamic_cast<DistributedVector<T>*>(&v_local_in);

  assert (v_local != NULL);

  v_local->_first_local_index = 0;
  
  v_local->_global_size =
    v_local->_local_size =
    v_local->_last_local_index = size();

  v_local->_is_initialized =
    v_local->_is_closed = true;
  
  // Call localize on the vector's values.  This will help
  // prevent code duplication
  localize (v_local->_values);    
  
#ifndef HAVE_MPI

  assert (local_size() == size());
  
#endif
}



template <typename T>
void DistributedVector<T>::localize (NumericVector<T>& v_local_in,
				     const std::vector<unsigned int>&) const
{
  assert (this->initialized());
  assert (_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);

  // We don't support the send list.  Call the less efficient localize(v_local_in)
  localize (v_local_in);
}



template <typename T>
void DistributedVector<T>::localize (const unsigned int first_local_idx,
				     const unsigned int last_local_idx,
				     const std::vector<unsigned int>& send_list)
{
  // Only good for serial vectors
  assert (this->size() == this->local_size());
  assert (last_local_idx > first_local_idx);
  assert (send_list.size() <= this->size());
  assert (last_local_idx < this->size());
  
  const unsigned int size       = this->size();
  const unsigned int local_size = (last_local_idx - first_local_idx + 1);

    // Don't bother for serial cases
  if ((first_local_idx == 0) &&
      (local_size == size))
    return;
  
	  
  // Build a parallel vector, initialize it with the local
  // parts of (*this)
  DistributedVector<T> parallel_vec;

  parallel_vec.init (size, local_size);

  // Copy part of *this into the parallel_vec
  for (unsigned int i=first_local_idx; i<=last_local_idx; i++)
    parallel_vec._values[i-first_local_idx] = _values[i];

  // localize like normal
  parallel_vec.localize (*this, send_list);    
}



template <>
void DistributedVector<float>::localize (std::vector<float>& v_local) const

{
  assert (this->initialized());
  assert (_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);

  v_local.resize(size());
  
  std::fill (v_local.begin(),
	     v_local.end(),
	     0.);

  for (unsigned int i=0; i<local_size(); i++)
    v_local[i+first_local_index()] = _values[i];

#ifdef HAVE_MPI

  MPI_Allreduce (&v_local[0], &v_local[0], v_local.size(),
		 MPI_FLOAT, MPI_SUM, libMesh::COMM_WORLD);
  
#else

  assert (local_size() == size());
  
#endif  
}



template <>
void DistributedVector<double>::localize (std::vector<double>& v_local) const

{
  assert (this->initialized());
  assert (_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);

  v_local.resize(size());
  
  std::fill (v_local.begin(),
	     v_local.end(),
	     0.);

  for (unsigned int i=0; i<local_size(); i++)
    v_local[i+first_local_index()] = _values[i];

#ifdef HAVE_MPI

  MPI_Allreduce (&v_local[0], &v_local[0], v_local.size(),
		 MPI_DOUBLE, MPI_SUM, libMesh::COMM_WORLD);
  
#else

  assert (local_size() == size());
  
#endif  
}



template <>
void DistributedVector<float>::localize_to_one (std::vector<float>& v_local,
						const unsigned int pid) const
{
  assert (this->initialized());
  assert (_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);

  v_local.resize(size());
  
  std::fill (v_local.begin(),
	     v_local.end(),
	     0.);

  for (unsigned int i=0; i<local_size(); i++)
    v_local[i+first_local_index()] = _values[i];

#ifdef HAVE_MPI

  MPI_Reduce (&v_local[0], &v_local[0], v_local.size(),
	      MPI_FLOAT, MPI_SUM, pid, libMesh::COMM_WORLD);
  
#else

  assert (this->local_size() == this->size());
  assert (pid == 0);
  
#endif  
}



template <>
void DistributedVector<double>::localize_to_one (std::vector<double>& v_local,
						 const unsigned int pid) const
{
  assert (this->initialized());
  assert (_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);

  v_local.resize(size());
  
  std::fill (v_local.begin(),
	     v_local.end(),
	     0.);

  for (unsigned int i=0; i<local_size(); i++)
    v_local[i+first_local_index()] = _values[i];

#ifdef HAVE_MPI

  MPI_Reduce (&v_local[0], &v_local[0], v_local.size(),
	      MPI_DOUBLE, MPI_SUM, pid, libMesh::COMM_WORLD);
  
#else

  assert (this->local_size() == this->size());
  assert (pid == 0);
  
#endif  
}


//--------------------------------------------------------------
// Explicit instantiations
template class DistributedVector<float>;
template class DistributedVector<double>;

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



// C++ includes
#include <algorithm> // for std::min
#include <limits>

// Local Includes
#include "libmesh/dense_subvector.h"
#include "libmesh/dense_vector.h"
#include "libmesh/eigen_sparse_vector.h"
#include "libmesh/eigen_sparse_matrix.h"


#ifdef LIBMESH_HAVE_EIGEN

namespace libMesh
{

template <typename T>
T EigenSparseVector<T>::sum () const
{
  libmesh_assert (this->closed());
  libmesh_assert (this->initialized());

  return _vec.sum();
}



template <typename T>
Real EigenSparseVector<T>::l1_norm () const
{
  libmesh_assert (this->closed());
  libmesh_assert (this->initialized());

  return _vec.lpNorm<1>();
}



template <typename T>
Real EigenSparseVector<T>::l2_norm () const
{
  libmesh_assert (this->closed());
  libmesh_assert (this->initialized());

  return _vec.lpNorm<2>();
}



template <typename T>
Real EigenSparseVector<T>::linfty_norm () const
{
  libmesh_assert (this->closed());
  libmesh_assert (this->initialized());

  return _vec.lpNorm<Eigen::Infinity>();
}



template <typename T>
NumericVector<T>& EigenSparseVector<T>::operator += (const NumericVector<T>& v_in)
{
  libmesh_assert (this->closed());

  const EigenSparseVector<T>& v = libmesh_cast_ref<const EigenSparseVector<T>&>(v_in);

  _vec += v._vec;

  return *this;
}




template <typename T>
NumericVector<T>& EigenSparseVector<T>::operator -= (const NumericVector<T>& v_in)
{
  libmesh_assert (this->closed());

  const EigenSparseVector<T>& v = libmesh_cast_ref<const EigenSparseVector<T>&>(v_in);

  _vec -= v._vec;

  return *this;
}



template <typename T>
void EigenSparseVector<T>::reciprocal()
{
#ifndef NDEBUG
  const numeric_index_type n = this->size();

  for (numeric_index_type i=0; i<n; i++)
    // Don't divide by zero!
    libmesh_assert_not_equal_to ((*this)(i), T(0));
#endif

  _vec = _vec.cwiseInverse();
}



template <typename T>
void EigenSparseVector<T>::conjugate()
{
  _vec = _vec.conjugate();
}



template <typename T>
void EigenSparseVector<T>::add (const T v)
{
  _vec += EigenSV::Constant(this->size(), v);

#ifndef NDEBUG
  this->_is_closed = false;
#endif
}




template <typename T>
void EigenSparseVector<T>::add (const NumericVector<T>& v_in)
{
  libmesh_assert (this->initialized());

  const EigenSparseVector<T>& v = libmesh_cast_ref<const EigenSparseVector<T>&>(v_in);

  _vec += v._vec;
}



template <typename T>
void EigenSparseVector<T>::add (const T a, const NumericVector<T>& v_in)
{
  libmesh_assert (this->initialized());

  const EigenSparseVector<T>& v = libmesh_cast_ref<const EigenSparseVector<T>&>(v_in);

  _vec += v._vec*a;
}



template <typename T>
void EigenSparseVector<T>::add_vector (const std::vector<T>& v,
				       const std::vector<numeric_index_type>& dof_indices)
{
  libmesh_assert (!v.empty());
  libmesh_assert_equal_to (v.size(), dof_indices.size());

  for (numeric_index_type i=0; i<v.size(); i++)
    this->add (dof_indices[i], v[i]);
}



template <typename T>
void EigenSparseVector<T>::add_vector (const NumericVector<T>& V,
				       const std::vector<numeric_index_type>& dof_indices)
{
  libmesh_assert_equal_to (V.size(), dof_indices.size());

  for (numeric_index_type i=0; i<V.size(); i++)
    this->add (dof_indices[i], V(i));
}



template <typename T>
void EigenSparseVector<T>::add_vector (const DenseVector<T>& V,
				       const std::vector<numeric_index_type>& dof_indices)
{
  libmesh_assert_equal_to (V.size(), dof_indices.size());

  for (unsigned int i=0; i<V.size(); i++)
    this->add (dof_indices[i], V(i));
}



template <typename T>
void EigenSparseVector<T>::insert (const std::vector<T>& v,
				   const std::vector<numeric_index_type>& dof_indices)
{
  libmesh_assert (!v.empty());
  libmesh_assert_equal_to (v.size(), dof_indices.size());

  for (numeric_index_type i=0; i<v.size(); i++)
    this->set (dof_indices[i], v[i]);
}



template <typename T>
void EigenSparseVector<T>::insert (const NumericVector<T>& V,
				   const std::vector<numeric_index_type>& dof_indices)
{
  libmesh_assert_equal_to (V.size(), dof_indices.size());

  for (numeric_index_type i=0; i<V.size(); i++)
   this->set (dof_indices[i], V(i));
}



template <typename T>
void EigenSparseVector<T>::insert (const DenseVector<T>& V,
				   const std::vector<numeric_index_type>& dof_indices)
{
  libmesh_assert_equal_to (V.size(), dof_indices.size());

  for (unsigned int i=0; i<V.size(); i++)
    this->set (dof_indices[i], V(i));
}



template <typename T>
void EigenSparseVector<T>::insert (const DenseSubVector<T>& V,
				   const std::vector<numeric_index_type>& dof_indices)
{
  libmesh_assert_equal_to (V.size(), dof_indices.size());

  for (unsigned int i=0; i<V.size(); i++)
    this->set (dof_indices[i], V(i));
}



template <typename T>
void EigenSparseVector<T>::add_vector (const NumericVector<T> &vec_in,
				       const SparseMatrix<T>  &mat_in)
{
  // Make sure the data passed in are really in Eigen types
  const EigenSparseVector<T>* vec = libmesh_cast_ptr<const EigenSparseVector<T>*>(&vec_in);
  const EigenSparseMatrix<T>* mat = libmesh_cast_ptr<const EigenSparseMatrix<T>*>(&mat_in);

  libmesh_assert(vec);
  libmesh_assert(mat);

  _vec += mat->_mat*vec->_vec;
}



template <typename T>
void EigenSparseVector<T>::add_vector_transpose (const NumericVector<T> &vec_in,
						 const SparseMatrix<T>  &mat_in)
{
  // Make sure the data passed in are really in Eigen types
  const EigenSparseVector<T>* vec = libmesh_cast_ptr<const EigenSparseVector<T>*>(&vec_in);
  const EigenSparseMatrix<T>* mat = libmesh_cast_ptr<const EigenSparseMatrix<T>*>(&mat_in);

  libmesh_assert(vec);
  libmesh_assert(mat);

  _vec += mat->_mat.transpose()*vec->_vec;
}



template <typename T>
void EigenSparseVector<T>::scale (const T factor)
{
  libmesh_assert (this->initialized());

  _vec *= factor;
}



template <typename T>
void EigenSparseVector<T>::abs()
{
  libmesh_assert (this->initialized());

  const numeric_index_type n = this->size();

  for (numeric_index_type i=0; i!=n; ++i)
    this->set(i,std::abs((*this)(i)));
}



template <typename T>
T EigenSparseVector<T>::dot (const NumericVector<T>& V) const
{
  libmesh_assert (this->initialized());

  // Make sure the NumericVector passed in is really a EigenSparseVector
  const EigenSparseVector<T>* v = libmesh_cast_ptr<const EigenSparseVector<T>*>(&V);
  libmesh_assert(v);

  return _vec.dot(v->_vec);
}



template <typename T>
NumericVector<T>&
EigenSparseVector<T>::operator = (const T s)
{
  libmesh_assert (this->initialized());
  libmesh_assert (this->closed());

  _vec.fill(s);

  return *this;
}



template <typename T>
NumericVector<T>&
EigenSparseVector<T>::operator = (const NumericVector<T>& v_in)
{
  // Make sure the NumericVector passed in is really a EigenSparseVector
  const EigenSparseVector<T>* v =
    libmesh_cast_ptr<const EigenSparseVector<T>*>(&v_in);

  libmesh_assert(v);

  *this = *v;

  return *this;
}



template <typename T>
EigenSparseVector<T>&
EigenSparseVector<T>::operator = (const EigenSparseVector<T>& v)
{
  libmesh_assert (this->initialized());
  libmesh_assert (v.closed());
  libmesh_assert_equal_to (this->size(), v.size());

  _vec = v._vec;

#ifndef NDEBUG
  this->_is_closed = true;
#endif

  return *this;
}



template <typename T>
NumericVector<T>&
EigenSparseVector<T>::operator = (const std::vector<T>& v)
{
  /**
   * Case 1:  The vector is the same size of
   * The global vector.  Only add the local components.
   */
  if (this->size() == v.size())
    for (numeric_index_type i=0; i<v.size(); i++)
      this->set (i, v[i]);

  else
    libmesh_error();

  return *this;
}


template <typename T>
void EigenSparseVector<T>::localize (NumericVector<T>& v_local_in) const
{
  // Make sure the NumericVector passed in is really a EigenSparseVector
  EigenSparseVector<T>* v_local =
    libmesh_cast_ptr<EigenSparseVector<T>*>(&v_local_in);

  libmesh_assert(v_local);

  *v_local = *this;
}



template <typename T>
void EigenSparseVector<T>::localize (NumericVector<T>& v_local_in,
				     const std::vector<numeric_index_type>& libmesh_dbg_var(send_list)) const
{
  // Make sure the NumericVector passed in is really a EigenSparseVector
  EigenSparseVector<T>* v_local =
    libmesh_cast_ptr<EigenSparseVector<T>*>(&v_local_in);

  libmesh_assert(v_local);
  libmesh_assert_less_equal (send_list.size(), v_local->size());

  *v_local = *this;
}



template <typename T>
void EigenSparseVector<T>::localize (const numeric_index_type libmesh_dbg_var(first_local_idx),
				     const numeric_index_type libmesh_dbg_var(last_local_idx),
				     const std::vector<numeric_index_type>& libmesh_dbg_var(send_list))
{
  libmesh_assert_equal_to (first_local_idx, 0);
  libmesh_assert_equal_to (last_local_idx+1, this->size());

  libmesh_assert_less_equal (send_list.size(), this->size());

#ifndef NDEBUG
  this->_is_closed = true;
#endif
}



template <typename T>
void EigenSparseVector<T>::localize (std::vector<T>& v_local) const

{
  v_local.resize(this->size());

  for (numeric_index_type i=0; i<v_local.size(); i++)
    v_local[i] = (*this)(i);
}



template <typename T>
void EigenSparseVector<T>::localize_to_one (std::vector<T>& v_local,
					    const processor_id_type libmesh_dbg_var(pid)) const
{
  libmesh_assert_equal_to (pid, 0);

  this->localize (v_local);
}



template <typename T>
void EigenSparseVector<T>::pointwise_mult (const NumericVector<T>& /*vec1*/,
					   const NumericVector<T>& /*vec2*/)
{
  libmesh_not_implemented();
}



template <typename T>
Real EigenSparseVector<T>::max() const
{
  libmesh_assert (this->initialized());
  if (!this->size())
    return -std::numeric_limits<Real>::max();

#ifdef LIBMESH_USE_COMPLEX_NUMBERS
  Real the_max = libmesh_real((*this)(0));

  const numeric_index_type n = this->size();

  for (numeric_index_type i=1; i<n; i++)
    the_max = std::max (the_max, libmesh_real((*this)(i)));

  return the_max;
#else
  return libmesh_real(_vec.maxCoeff());
#endif
}



template <typename T>
Real EigenSparseVector<T>::min () const
{
  libmesh_assert (this->initialized());
  if (!this->size())
    return std::numeric_limits<Real>::max();

#ifdef LIBMESH_USE_COMPLEX_NUMBERS
  Real the_min = libmesh_real((*this)(0));

  const numeric_index_type n = this->size();

  for (numeric_index_type i=1; i<n; i++)
    the_min = std::min (the_min, libmesh_real((*this)(i)));

  return the_min;
#else
  return libmesh_real(_vec.minCoeff());
#endif
}


//------------------------------------------------------------------
// Explicit instantiations
template class EigenSparseVector<Number>;

} // namespace libMesh


#endif // #ifdef LIBMESH_HAVE_EIGEN

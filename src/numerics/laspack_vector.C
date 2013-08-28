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
#include "libmesh/laspack_vector.h"
#include "libmesh/laspack_matrix.h"


#ifdef LIBMESH_HAVE_LASPACK

namespace libMesh
{

template <typename T>
T LaspackVector<T>::sum () const
{
  libmesh_assert (this->closed());

  T _sum = 0;

  const numeric_index_type n = this->size();

  for (numeric_index_type i=0; i!=n; ++i)
    _sum += (*this)(i);

  return _sum;
}



template <typename T>
Real LaspackVector<T>::l1_norm () const
{
  libmesh_assert (this->closed());

  return static_cast<Real>(l1Norm_V(const_cast<QVector*>(&_vec)));
}



template <typename T>
Real LaspackVector<T>::l2_norm () const
{
  libmesh_assert (this->closed());

  return static_cast<Real>(l2Norm_V(const_cast<QVector*>(&_vec)));
}



template <typename T>
Real LaspackVector<T>::linfty_norm () const
{
  libmesh_assert (this->closed());

  return static_cast<Real>(MaxNorm_V(const_cast<QVector*>(&_vec)));
}



template <typename T>
NumericVector<T>& LaspackVector<T>::operator += (const NumericVector<T>& v)
{
  libmesh_assert (this->closed());

  this->add(1., v);

  return *this;
}




template <typename T>
NumericVector<T>& LaspackVector<T>::operator -= (const NumericVector<T>& v)
{
  libmesh_assert (this->closed());

  this->add(-1., v);

  return *this;
}



template <typename T>
void LaspackVector<T>::reciprocal()
{
  const numeric_index_type n = this->size();

  for (numeric_index_type i=0; i<n; i++)
    {
      T v = (*this)(i);

      // Don't divide by zero!
      libmesh_assert_not_equal_to (v, T(0));

      this->set(i, 1. / v);
    }
}



template <typename T>
void LaspackVector<T>::conjugate()
{
  const numeric_index_type n = this->size();

  for (numeric_index_type i=0; i<n; i++)
    {
      T v = (*this)(i);

      this->set(i, libmesh_conj(v) );
    }
}


template <typename T>
void LaspackVector<T>::add (const T v)
{
  const numeric_index_type n = this->size();

  for (numeric_index_type i=0; i<n; i++)
    this->add (i, v);

#ifndef NDEBUG
  this->_is_closed = false;
#endif
}




template <typename T>
void LaspackVector<T>::add (const NumericVector<T>& v)
{
  this->add (1., v);
}



template <typename T>
void LaspackVector<T>::add (const T a, const NumericVector<T>& v_in)
{
  // Make sure the vector passed in is really a LaspackVector
  const LaspackVector* v = libmesh_cast_ptr<const LaspackVector*>(&v_in);

#ifndef NDEBUG
  const bool was_closed = this->_is_closed;
#endif

  libmesh_assert(v);
  libmesh_assert_equal_to (this->size(), v->size());

  for (numeric_index_type i=0; i<v->size(); i++)
    this->add (i, a*(*v)(i));

#ifndef NDEBUG
  this->_is_closed = was_closed;
#endif
}



template <typename T>
void LaspackVector<T>::add_vector (const std::vector<T>& v,
				   const std::vector<numeric_index_type>& dof_indices)
{
  libmesh_assert (!v.empty());
  libmesh_assert_equal_to (v.size(), dof_indices.size());

  for (numeric_index_type i=0; i<v.size(); i++)
    this->add (dof_indices[i], v[i]);
}



template <typename T>
void LaspackVector<T>::add_vector (const NumericVector<T>& V,
				   const std::vector<numeric_index_type>& dof_indices)
{
  libmesh_assert_equal_to (V.size(), dof_indices.size());

  for (numeric_index_type i=0; i<V.size(); i++)
    this->add (dof_indices[i], V(i));
}



template <typename T>
void LaspackVector<T>::add_vector (const DenseVector<T>& V,
				   const std::vector<numeric_index_type>& dof_indices)
{
  libmesh_assert_equal_to (V.size(), dof_indices.size());

  for (unsigned int i=0; i<V.size(); i++)
    this->add (dof_indices[i], V(i));
}



template <typename T>
void LaspackVector<T>::insert (const std::vector<T>& v,
			       const std::vector<numeric_index_type>& dof_indices)
{
  libmesh_assert (!v.empty());
  libmesh_assert_equal_to (v.size(), dof_indices.size());

  for (numeric_index_type i=0; i<v.size(); i++)
    this->set (dof_indices[i], v[i]);
}



template <typename T>
void LaspackVector<T>::insert (const NumericVector<T>& V,
			       const std::vector<numeric_index_type>& dof_indices)
{
  libmesh_assert_equal_to (V.size(), dof_indices.size());

  for (numeric_index_type i=0; i<V.size(); i++)
   this->set (dof_indices[i], V(i));
}



template <typename T>
void LaspackVector<T>::insert (const DenseVector<T>& V,
			       const std::vector<numeric_index_type>& dof_indices)
{
  libmesh_assert_equal_to (V.size(), dof_indices.size());

  for (unsigned int i=0; i<V.size(); i++)
    this->set (dof_indices[i], V(i));
}



template <typename T>
void LaspackVector<T>::insert (const DenseSubVector<T>& V,
			       const std::vector<numeric_index_type>& dof_indices)
{
  libmesh_assert_equal_to (V.size(), dof_indices.size());

  for (unsigned int i=0; i<V.size(); i++)
    this->set (dof_indices[i], V(i));
}



template <typename T>
void LaspackVector<T>::add_vector (const NumericVector<T> &vec_in,
				   const SparseMatrix<T> &mat_in)
{
  // Make sure the data passed in are really in Laspack types
  const LaspackVector<T>* vec = libmesh_cast_ptr<const LaspackVector<T>*>(&vec_in);
  const LaspackMatrix<T>* mat = libmesh_cast_ptr<const LaspackMatrix<T>*>(&mat_in);

  libmesh_assert(vec);
  libmesh_assert(mat);

  // += mat*vec
  AddAsgn_VV (&_vec, Mul_QV(const_cast<QMatrix*>(&mat->_QMat),
			    const_cast<QVector*>(&vec->_vec)));
}


template <typename T>
void LaspackVector<T>::add_vector_transpose (const NumericVector<T> &,
				             const SparseMatrix<T> &)
{
  libmesh_not_implemented();
}



template <typename T>
void LaspackVector<T>::scale (const T factor)
{
  libmesh_assert (this->initialized());

  Asgn_VV(&_vec, Mul_SV (factor, &_vec));
}

template <typename T>
void LaspackVector<T>::abs()
{
  libmesh_assert (this->initialized());

  const numeric_index_type n = this->size();

  for (numeric_index_type i=0; i!=n; ++i)
    this->set(i,std::abs((*this)(i)));
}

template <typename T>
T LaspackVector<T>::dot (const NumericVector<T>& V) const
{
  libmesh_assert (this->initialized());

  // Make sure the NumericVector passed in is really a LaspackVector
  const LaspackVector<T>* v = libmesh_cast_ptr<const LaspackVector<T>*>(&V);
  libmesh_assert(v);

  return Mul_VV (const_cast<QVector*>(&(this->_vec)),
		 const_cast<QVector*>(&(v->_vec)));
}



template <typename T>
NumericVector<T>&
LaspackVector<T>::operator = (const T s)
{
  libmesh_assert (this->initialized());
  libmesh_assert (this->closed());

  V_SetAllCmp (&_vec, s);

  return *this;
}



template <typename T>
NumericVector<T>&
LaspackVector<T>::operator = (const NumericVector<T>& v_in)
{
  // Make sure the NumericVector passed in is really a LaspackVector
  const LaspackVector<T>* v =
    libmesh_cast_ptr<const LaspackVector<T>*>(&v_in);

  libmesh_assert(v);

  *this = *v;

  return *this;
}



template <typename T>
LaspackVector<T>&
LaspackVector<T>::operator = (const LaspackVector<T>& v)
{
  libmesh_assert (this->initialized());
  libmesh_assert (v.closed());
  libmesh_assert_equal_to (this->size(), v.size());

  if (v.size() != 0)
    Asgn_VV (const_cast<QVector*>(&_vec),
	     const_cast<QVector*>(&v._vec)
	     );

#ifndef NDEBUG
  this->_is_closed = true;
#endif

  return *this;
}



template <typename T>
NumericVector<T>&
LaspackVector<T>::operator = (const std::vector<T>& v)
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
void LaspackVector<T>::localize (NumericVector<T>& v_local_in) const
{
  // Make sure the NumericVector passed in is really a LaspackVector
  LaspackVector<T>* v_local =
    libmesh_cast_ptr<LaspackVector<T>*>(&v_local_in);

  libmesh_assert(v_local);

  *v_local = *this;
}



template <typename T>
void LaspackVector<T>::localize (NumericVector<T>& v_local_in,
				 const std::vector<numeric_index_type>& libmesh_dbg_var(send_list)) const
{
  // Make sure the NumericVector passed in is really a LaspackVector
  LaspackVector<T>* v_local =
    libmesh_cast_ptr<LaspackVector<T>*>(&v_local_in);

  libmesh_assert(v_local);
  libmesh_assert_less_equal (send_list.size(), v_local->size());

  *v_local = *this;
}



template <typename T>
void LaspackVector<T>::localize (const numeric_index_type libmesh_dbg_var(first_local_idx),
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
void LaspackVector<T>::localize (std::vector<T>& v_local) const

{
  v_local.resize(this->size());

  for (numeric_index_type i=0; i<v_local.size(); i++)
    v_local[i] = (*this)(i);
}



template <typename T>
void LaspackVector<T>::localize_to_one (std::vector<T>& v_local,
					const processor_id_type libmesh_dbg_var(pid)) const
{
  libmesh_assert_equal_to (pid, 0);

  this->localize (v_local);
}



template <typename T>
void LaspackVector<T>::pointwise_mult (const NumericVector<T>& /*vec1*/,
				       const NumericVector<T>& /*vec2*/)
{
  libmesh_not_implemented();
}



template <typename T>
Real LaspackVector<T>::max() const
{
  libmesh_assert (this->initialized());
  if (!this->size())
    return -std::numeric_limits<Real>::max();

  Real the_max = libmesh_real((*this)(0));

  const numeric_index_type n = this->size();

  for (numeric_index_type i=1; i<n; i++)
    the_max = std::max (the_max, libmesh_real((*this)(i)));

  return the_max;
}



template <typename T>
Real LaspackVector<T>::min () const
{
  libmesh_assert (this->initialized());
  if (!this->size())
    return std::numeric_limits<Real>::max();

  Real the_min = libmesh_real((*this)(0));

  const numeric_index_type n = this->size();

  for (numeric_index_type i=1; i<n; i++)
    the_min = std::min (the_min, libmesh_real((*this)(i)));

  return the_min;
}


//------------------------------------------------------------------
// Explicit instantiations
template class LaspackVector<Number>;

} // namespace libMesh


#endif // #ifdef LIBMESH_HAVE_LASPACK

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



// C++ includes
#include <algorithm> // for std::min

// Local Includes
#include "laspack_vector.h"
#include "laspack_matrix.h"


#ifdef HAVE_LASPACK

template <typename T>
T LaspackVector<T>::sum () const
{
  assert (this->closed());

  T _sum = 0;

  const unsigned int n = this->size();
  
  for (unsigned int i=0; i!=n; ++i)
    _sum += (*this)(i);

  return _sum;
}



template <typename T>
Real LaspackVector<T>::l1_norm () const
{
  assert (this->closed());
  
  return static_cast<Real>(l1Norm_V(const_cast<QVector*>(&_vec)));
}



template <typename T>
Real LaspackVector<T>::l2_norm () const
{
  assert (this->closed());
  
  return static_cast<Real>(l2Norm_V(const_cast<QVector*>(&_vec)));
}



template <typename T>
Real LaspackVector<T>::linfty_norm () const
{
  assert (this->closed());
  
  return static_cast<Real>(MaxNorm_V(const_cast<QVector*>(&_vec)));
}



template <typename T> 
Real LaspackVector<T>::max() const
{
  assert (this->initialized());

  return static_cast<Real>(MaxNorm_V(const_cast<QVector*>(&_vec)));
}




template <typename T>
NumericVector<T>& LaspackVector<T>::operator += (const NumericVector<T>& v)
{
  assert (this->closed());
  
  this->add(1., v);
  
  return *this;
}




template <typename T>
NumericVector<T>& LaspackVector<T>::operator -= (const NumericVector<T>& v)
{
  assert (this->closed());
  
  this->add(-1., v);
  
  return *this;
}



template <typename T>
void LaspackVector<T>::add (const T v)
{
  const unsigned int n = this->size();
  
  for (unsigned int i=0; i<n; i++)
    this->add (i, v);
}




template <typename T>
void LaspackVector<T>::add (const NumericVector<T>& v)
{
  this->add (1., v);
}



template <typename T>
void LaspackVector<T>::add (const T a, const NumericVector<T>& v_in)
{
  const LaspackVector* v = dynamic_cast<const LaspackVector*>(&v_in);

  assert (v != NULL);
  assert (this->size() == v->size());

  for (unsigned int i=0; i<v->size(); i++)
    this->add (i, a*(*v)(i));
}



template <typename T> 
void LaspackVector<T>::add_vector (const std::vector<T>& v,
				   const std::vector<unsigned int>& dof_indices)
{
  assert (!v.empty());
  assert (v.size() == dof_indices.size());
  
  for (unsigned int i=0; i<v.size(); i++)
    this->add (dof_indices[i], v[i]);
}



template <typename T> 
void LaspackVector<T>::add_vector (const NumericVector<T>& V,
				   const std::vector<unsigned int>& dof_indices)
{
  assert (V.size() == dof_indices.size());

  for (unsigned int i=0; i<V.size(); i++)
    this->add (dof_indices[i], V(i));
}



template <typename T> 
void LaspackVector<T>::add_vector (const DenseVector<T>& V,
				   const std::vector<unsigned int>& dof_indices)
{
  assert (V.size() == dof_indices.size());

  for (unsigned int i=0; i<V.size(); i++)
    this->add (dof_indices[i], V(i));
}



template <typename T> 
void LaspackVector<T>::insert (const std::vector<T>& v,
			       const std::vector<unsigned int>& dof_indices)
{
  assert (!v.empty());
  assert (v.size() == dof_indices.size());
  
  for (unsigned int i=0; i<v.size(); i++)
    this->set (dof_indices[i], v[i]);
}



template <typename T> 
void LaspackVector<T>::insert (const NumericVector<T>& V,
			       const std::vector<unsigned int>& dof_indices)
{
  assert (V.size() == dof_indices.size());

  for (unsigned int i=0; i<V.size(); i++)
   this->set (dof_indices[i], V(i));
}



template <typename T> 
void LaspackVector<T>::insert (const DenseVector<T>& V,
			       const std::vector<unsigned int>& dof_indices)
{
  assert (V.size() == dof_indices.size());

  for (unsigned int i=0; i<V.size(); i++)
    this->set (dof_indices[i], V(i));
}



template <typename T> 
void LaspackVector<T>::add_vector (const NumericVector<T> &vec_in,
				   const SparseMatrix<T> &mat_in)
{
  // Convert from input types
  const LaspackVector<T>* vec = dynamic_cast<const LaspackVector<T>*>(&vec_in);
  const LaspackMatrix<T>* mat = dynamic_cast<const LaspackMatrix<T>*>(&mat_in);

  assert (vec != NULL);
  assert (mat != NULL);
  
  // += mat*vec
  AddAsgn_VV (&_vec, Mul_QV(const_cast<QMatrix*>(&mat->_QMat),
			    const_cast<QVector*>(&vec->_vec)));
}



template <typename T>
void LaspackVector<T>::scale (const T factor)
{
  assert (this->initialized());
  
  Mul_SV (factor, &_vec);
}


template <typename T>
T LaspackVector<T>::dot (const NumericVector<T>& V) const
{
  assert (this->initialized());

  // Make sure the NumericVector passed in is really a LasPackVector
  const LaspackVector<T>* v = dynamic_cast<const LaspackVector<T>*>(&V);
  assert (v != NULL);
  
  return Mul_VV (const_cast<QVector*>(&(this->_vec)),
		 const_cast<QVector*>(&(v->_vec)));
}



template <typename T>
NumericVector<T>& 
LaspackVector<T>::operator = (const T s)
{
  assert (this->initialized());

  V_SetAllCmp (&_vec, s);
  
  return *this;
}



template <typename T>
NumericVector<T>&
LaspackVector<T>::operator = (const NumericVector<T>& v_in)
{
  const LaspackVector<T>* v =
    dynamic_cast<const LaspackVector<T>*>(&v_in);

  assert (v != NULL);
  
  *this = *v;

  return *this;
}



template <typename T>
LaspackVector<T>&
LaspackVector<T>::operator = (const LaspackVector<T>& v)
{
  if (!this->initialized())
    this->init (v.size());

  if (this->size() != v.size())
    this->init (v.size());
  

  this->_is_closed = v._is_closed;

  if (v.size() != 0)    
    Asgn_VV (const_cast<QVector*>(&_vec),
	     const_cast<QVector*>(&v._vec)
	     );
  
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
    for (unsigned int i=0; i<v.size(); i++)
      this->set (i, v[i]);
  
  else
    error();

  return *this;
}


template <typename T>
void LaspackVector<T>::localize (NumericVector<T>& v_local_in) const
{
  LaspackVector<T>* v_local =
    dynamic_cast<LaspackVector<T>*>(&v_local_in);

  assert (v_local != NULL);

  *v_local = *this;
}



template <typename T>
void LaspackVector<T>::localize (NumericVector<T>& v_local_in,
				 const std::vector<unsigned int>& send_list) const
{
  LaspackVector<T>* v_local =
    dynamic_cast<LaspackVector<T>*>(&v_local_in);

  assert (v_local != NULL);
  assert (send_list.size() == v_local->size());

  *v_local = *this;
}



template <typename T>
void LaspackVector<T>::localize (const unsigned int first_local_idx,
				 const unsigned int last_local_idx,
				 const std::vector<unsigned int>& send_list)
{
  assert (first_local_idx  == 0);
  assert (last_local_idx+1 == this->size());
  
  assert (send_list.size() == this->size());
}



template <typename T>
void LaspackVector<T>::localize (std::vector<T>& v_local) const

{
  v_local.resize(this->size());

  for (unsigned int i=0; i<v_local.size(); i++)
    v_local[i] = (*this)(i);  
}



template <typename T>
void LaspackVector<T>::localize_to_one (std::vector<T>& v_local,
					const unsigned int pid) const
{
  assert (pid == 0);

  this->localize (v_local);
}


// Full specialization for Real datatypes
#ifdef USE_REAL_NUMBERS

template <> 
Real LaspackVector<Real>::min () const
{
  assert (this->initialized());

  Real min = (*this)(0);

  const unsigned int n = this->size();
  
  for (unsigned int i=1; i<n; i++)
    min = std::min (min, (*this)(i));

  return min;
}


#endif



// Full specialization for Complex datatypes
#ifdef USE_COMPLEX_NUMBERS

template <> 
Real LaspackVector<Complex>::min () const
{
  assert (this->initialized());

  Real min = (*this)(0).real();

  const unsigned int n = this->size();
  
  for (unsigned int i=1; i<n; i++)
    min = std::min (min, (*this)(i).real() );

  return min;
}

#endif



//------------------------------------------------------------------
// Explicit instantiations
template class LaspackVector<Number>;
 

#endif // #ifdef HAVE_LASPACK

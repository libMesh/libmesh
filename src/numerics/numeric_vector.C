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
#include <cstdlib> // *must* precede <cmath> for proper std:abs() on PGI, Sun Studio CC
#include <cmath> // for std::abs
#include <limits>

// Local Includes
#include "libmesh/numeric_vector.h"
#include "libmesh/distributed_vector.h"
#include "libmesh/laspack_vector.h"
#include "libmesh/eigen_sparse_vector.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/trilinos_epetra_vector.h"
#include "libmesh/shell_matrix.h"
#include "libmesh/tensor_tools.h"

namespace libMesh
{



//------------------------------------------------------------------
// NumericVector methods

// Full specialization for Real datatypes
template <typename T>
AutoPtr<NumericVector<T> >
NumericVector<T>::build(const SolverPackage solver_package, const Parallel::Communicator &comm)
{
  // Build the appropriate vector
  switch (solver_package)
    {


#ifdef LIBMESH_HAVE_LASPACK
    case LASPACK_SOLVERS:
      {
	AutoPtr<NumericVector<T> > ap(new LaspackVector<T>(AUTOMATIC, comm));
	return ap;
      }
#endif


#ifdef LIBMESH_HAVE_PETSC
    case PETSC_SOLVERS:
      {
	AutoPtr<NumericVector<T> > ap(new PetscVector<T>(AUTOMATIC, comm));
	return ap;
      }
#endif

/*
#ifdef LIBMESH_HAVE_TRILINOS
    case TRILINOS_SOLVERS:
      {
	AutoPtr<NumericVector<T> > ap(new EpetraVector<T>(comm));
	return ap;
      }
#endif


#ifdef LIBMESH_HAVE_EIGEN
    case EIGEN_SOLVERS:
      {
	AutoPtr<NumericVector<T> > ap(new EigenSparseVector<T>(comm));
	return ap;
      }
#endif


    default:
      AutoPtr<NumericVector<T> > ap(new DistributedVector<T>(comm));
      return ap;
*/
    }

  AutoPtr<NumericVector<T> > ap(NULL);
  return ap;
}


template <typename T>
int NumericVector<T>::compare (const NumericVector<T> &other_vector,
			       const Real threshold) const
{
  libmesh_assert (this->initialized());
  libmesh_assert (other_vector.initialized());
  libmesh_assert_equal_to (this->first_local_index(), other_vector.first_local_index());
  libmesh_assert_equal_to (this->last_local_index(), other_vector.last_local_index());

  int first_different_i = std::numeric_limits<int>::max();
  numeric_index_type i = first_local_index();

  do
    {
      if ( std::abs( (*this)(i) - other_vector(i) ) > threshold )
	first_different_i = i;
      else
	i++;
    }
  while (first_different_i==std::numeric_limits<int>::max()
         && i<last_local_index());

  // Find the correct first differing index in parallel
  this->communicator().min(first_different_i);

  if (first_different_i == std::numeric_limits<int>::max())
    return -1;

  return first_different_i;
}


template <typename T>
int NumericVector<T>::local_relative_compare (const NumericVector<T> &other_vector,
			                      const Real threshold) const
{
  libmesh_assert (this->initialized());
  libmesh_assert (other_vector.initialized());
  libmesh_assert_equal_to (this->first_local_index(), other_vector.first_local_index());
  libmesh_assert_equal_to (this->last_local_index(), other_vector.last_local_index());

  int first_different_i = std::numeric_limits<int>::max();
  numeric_index_type i = first_local_index();

  do
    {
      if ( std::abs( (*this)(i) - other_vector(i) ) > threshold *
           std::max(std::abs((*this)(i)), std::abs(other_vector(i))))
	first_different_i = i;
      else
	i++;
    }
  while (first_different_i==std::numeric_limits<int>::max()
         && i<last_local_index());

  // Find the correct first differing index in parallel
  this->communicator().min(first_different_i);

  if (first_different_i == std::numeric_limits<int>::max())
    return -1;

  return first_different_i;
}


template <typename T>
int NumericVector<T>::global_relative_compare (const NumericVector<T> &other_vector,
			                       const Real threshold) const
{
  libmesh_assert (this->initialized());
  libmesh_assert (other_vector.initialized());
  libmesh_assert_equal_to (this->first_local_index(), other_vector.first_local_index());
  libmesh_assert_equal_to (this->last_local_index(), other_vector.last_local_index());

  int first_different_i = std::numeric_limits<int>::max();
  numeric_index_type i = first_local_index();

  const Real my_norm = this->linfty_norm();
  const Real other_norm = other_vector.linfty_norm();
  const Real abs_threshold = std::max(my_norm, other_norm) * threshold;

  do
    {
      if ( std::abs( (*this)(i) - other_vector(i) ) > abs_threshold )
	first_different_i = i;
      else
	i++;
    }
  while (first_different_i==std::numeric_limits<int>::max()
         && i<last_local_index());

  // Find the correct first differing index in parallel
  this->communicator().min(first_different_i);

  if (first_different_i == std::numeric_limits<int>::max())
    return -1;

  return first_different_i;
}

/*
// Full specialization for float datatypes (DistributedVector wants this)

template <>
int NumericVector<float>::compare (const NumericVector<float> &other_vector,
				   const Real threshold) const
{
  libmesh_assert (this->initialized());
  libmesh_assert (other_vector.initialized());
  libmesh_assert_equal_to (this->first_local_index(), other_vector.first_local_index());
  libmesh_assert_equal_to (this->last_local_index(), other_vector.last_local_index());

  int rvalue     = -1;
  numeric_index_type i = first_local_index();

  do
    {
      if ( std::abs( (*this)(i) - other_vector(i) ) > threshold )
	rvalue = i;
      else
	i++;
    }
  while (rvalue==-1 && i<last_local_index());

  return rvalue;
}

// Full specialization for double datatypes
template <>
int NumericVector<double>::compare (const NumericVector<double> &other_vector,
				    const Real threshold) const
{
  libmesh_assert (this->initialized());
  libmesh_assert (other_vector.initialized());
  libmesh_assert_equal_to (this->first_local_index(), other_vector.first_local_index());
  libmesh_assert_equal_to (this->last_local_index(), other_vector.last_local_index());

  int rvalue     = -1;
  numeric_index_type i = first_local_index();

  do
    {
      if ( std::abs( (*this)(i) - other_vector(i) ) > threshold )
	rvalue = i;
      else
	i++;
    }
  while (rvalue==-1 && i<last_local_index());

  return rvalue;
}

#ifdef LIBMESH_DEFAULT_TRIPLE_PRECISION
// Full specialization for long double datatypes
template <>
int NumericVector<long double>::compare (const NumericVector<long double> &other_vector,
				         const Real threshold) const
{
  libmesh_assert (this->initialized());
  libmesh_assert (other_vector.initialized());
  libmesh_assert_equal_to (this->first_local_index(), other_vector.first_local_index());
  libmesh_assert_equal_to (this->last_local_index(), other_vector.last_local_index());

  int rvalue     = -1;
  numeric_index_type i = first_local_index();

  do
    {
      if ( std::abs( (*this)(i) - other_vector(i) ) > threshold )
	rvalue = i;
      else
	i++;
    }
  while (rvalue==-1 && i<last_local_index());

  return rvalue;
}
#endif


// Full specialization for Complex datatypes
template <>
int NumericVector<Complex>::compare (const NumericVector<Complex> &other_vector,
				     const Real threshold) const
{
  libmesh_assert (this->initialized());
  libmesh_assert (other_vector.initialized());
  libmesh_assert_equal_to (this->first_local_index(), other_vector.first_local_index());
  libmesh_assert_equal_to (this->last_local_index(), other_vector.last_local_index());

  int rvalue     = -1;
  numeric_index_type i = first_local_index();

  do
    {
      if (( std::abs( (*this)(i).real() - other_vector(i).real() ) > threshold ) ||
	  ( std::abs( (*this)(i).imag() - other_vector(i).imag() ) > threshold ))
	rvalue = i;
      else
	i++;
    }
  while (rvalue==-1 && i<this->last_local_index());

  return rvalue;
}
*/


template <class T>
Real NumericVector<T>::subset_l1_norm (const std::set<numeric_index_type> & indices) const
{
  const NumericVector<T> & v = *this;

  std::set<numeric_index_type>::const_iterator it = indices.begin();
  const std::set<numeric_index_type>::const_iterator it_end = indices.end();

  Real norm = 0;

  for(; it!=it_end; ++it)
    norm += std::abs(v(*it));

  this->communicator().sum(norm);

  return norm;
}

template <class T>
Real NumericVector<T>::subset_l2_norm (const std::set<numeric_index_type> & indices) const
{
  const NumericVector<T> & v = *this;

  std::set<numeric_index_type>::const_iterator it = indices.begin();
  const std::set<numeric_index_type>::const_iterator it_end = indices.end();

  Real norm = 0;

  for(; it!=it_end; ++it)
    norm += TensorTools::norm_sq(v(*it));

  this->communicator().sum(norm);

  return std::sqrt(norm);
}

template <class T>
Real NumericVector<T>::subset_linfty_norm (const std::set<numeric_index_type> & indices) const
{
  const NumericVector<T> & v = *this;

  std::set<numeric_index_type>::const_iterator it = indices.begin();
  const std::set<numeric_index_type>::const_iterator it_end = indices.end();

  Real norm = 0;

  for(; it!=it_end; ++it)
    {
      Real value = std::abs(v(*it));
      if(value > norm)
        norm = value;
    }

  this->communicator().max(norm);

  return norm;
}



template <typename T>
void NumericVector<T>::add_vector (const NumericVector<T>& v,
				   const ShellMatrix<T>& a)
{
  a.vector_mult_add(*this,v);
}



//------------------------------------------------------------------
// Explicit instantiations
template class NumericVector<Number>;

} // namespace libMesh

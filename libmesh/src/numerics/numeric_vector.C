// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2008 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner
  
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
#include <cmath> // for std::abs

// Local Includes
#include "numeric_vector.h"
#include "distributed_vector.h"
#include "laspack_vector.h"
#include "petsc_vector.h"
#include "trilinos_epetra_vector.h"
#include "shell_matrix.h"



//------------------------------------------------------------------
// NumericVector methods

// Full specialization for Real datatypes
template <typename T>
AutoPtr<NumericVector<T> >
NumericVector<T>::build(const SolverPackage solver_package)
{
  // Build the appropriate vector
  switch (solver_package)
    {


#ifdef LIBMESH_HAVE_LASPACK
    case LASPACK_SOLVERS:
      {
	AutoPtr<NumericVector<T> > ap(new LaspackVector<T>);
	return ap;
      }
#endif


#ifdef LIBMESH_HAVE_PETSC
    case PETSC_SOLVERS:
      {
	AutoPtr<NumericVector<T> > ap(new PetscVector<T>);
	return ap;
      }
#endif


#ifdef LIBMESH_HAVE_TRILINOS
    case TRILINOS_SOLVERS:
      {
	AutoPtr<NumericVector<T> > ap(new EpetraVector<T>);
	return ap;
      }
#endif


    default:
      AutoPtr<NumericVector<T> > ap(new DistributedVector<T>);
      return ap;
    }
    
  AutoPtr<NumericVector<T> > ap(NULL);
  return ap;    
}



// Full specialization for float datatypes (DistributedVector wants this)

template <>
int NumericVector<float>::compare (const NumericVector<float> &other_vector,
				   const Real threshold) const
{
  libmesh_assert (this->initialized());
  libmesh_assert (other_vector.initialized());
  libmesh_assert (this->first_local_index() == other_vector.first_local_index());
  libmesh_assert (this->last_local_index()  == other_vector.last_local_index());

  int rvalue     = -1;
  unsigned int i = first_local_index();

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
  libmesh_assert (this->first_local_index() == other_vector.first_local_index());
  libmesh_assert (this->last_local_index()  == other_vector.last_local_index());

  int rvalue     = -1;
  unsigned int i = first_local_index();

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

#ifdef TRIPLE_PRECISION
// Full specialization for long double datatypes
template <>
int NumericVector<long double>::compare (const NumericVector<long double> &other_vector,
				         const Real threshold) const
{
  libmesh_assert (this->initialized());
  libmesh_assert (other_vector.initialized());
  libmesh_assert (this->first_local_index() == other_vector.first_local_index());
  libmesh_assert (this->last_local_index()  == other_vector.last_local_index());

  int rvalue     = -1;
  unsigned int i = first_local_index();

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
  libmesh_assert (this->first_local_index() == other_vector.first_local_index());
  libmesh_assert (this->last_local_index()  == other_vector.last_local_index());

  int rvalue     = -1;
  unsigned int i = first_local_index();

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


template <class T>
Real NumericVector<T>::subset_l1_norm (const std::set<unsigned int> & indices)
{
  NumericVector<T> & v = *this;
  
  std::set<unsigned int>::iterator it = indices.begin();
  const std::set<unsigned int>::iterator it_end = indices.end();

  Real norm = 0;
  
  for(; it!=it_end; ++it)
    norm += std::abs(v(*it));

  Parallel::sum(norm);

  return norm;
}

template <class T>
Real NumericVector<T>::subset_l2_norm (const std::set<unsigned int> & indices)
{
  NumericVector<T> & v = *this;

  std::set<unsigned int>::iterator it = indices.begin();
  const std::set<unsigned int>::iterator it_end = indices.end();

  Real norm = 0;
  
  for(; it!=it_end; ++it)
    norm += v(*it) * v(*it);

  Parallel::sum(norm);

  return std::sqrt(norm);
}

template <class T>
Real NumericVector<T>::subset_linfty_norm (const std::set<unsigned int> & indices)
{
  NumericVector<T> & v = *this;

  std::set<unsigned int>::iterator it = indices.begin();
  const std::set<unsigned int>::iterator it_end = indices.end();

  Real norm = 0;
  
  for(; it!=it_end; ++it)
    {
      Real value = std::abs(v(*it));
      if(value > norm)
        norm = value;
    }

  Parallel::max(norm);

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

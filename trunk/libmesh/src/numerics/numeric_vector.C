// $Id: numeric_vector.C,v 1.7 2003-04-08 03:21:26 benkirk Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2002  Benjamin S. Kirk, John W. Peterson
  
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
#include <math.h>

// Local Includes
#include "numeric_vector.h"
#include "laspack_vector.h"
#include "petsc_vector.h"



//------------------------------------------------------------------
// NumericVector methods

// Full specialization for Real datatypes
template <typename T>
AutoPtr<NumericVector<T> >
NumericVector<T>::build(const SolverPackage solver_package_in)
{
  // Possibly overload the solver package based on
  // command-line arguments
  SolverPackage solver_package = solver_package_in;
  
#ifdef HAVE_PETSC
  if (libMesh::on_command_line ("--use-petsc"))
    solver_package = PETSC_SOLVERS;
#endif
  
#ifdef HAVE_LASPACK
  if (libMesh::on_command_line("--use-laspack"))
    solver_package = LASPACK_SOLVERS;
#endif




  // Build the appropriate vector
  switch (solver_package)
    {


#ifdef HAVE_LASPACK
    case LASPACK_SOLVERS:
      {
	AutoPtr<NumericVector<T> > ap(new LaspackVector<T>);
	return ap;
      }
#endif


#ifdef HAVE_PETSC
    case PETSC_SOLVERS:
      {
	AutoPtr<NumericVector<T> > ap(new PetscVector<T>);
	return ap;
      }
#endif

    default:
      std::cerr << "ERROR:  Unrecognized solver package: "
		<< solver_package
		<< std::endl;
      error();
    }
    
  AutoPtr<NumericVector<T> > ap(NULL);
  return ap;    
}



// Full specialization for float datatypes (DistributedVector wants this)
template <>
int NumericVector<float>::compare (const NumericVector<float> &other_vector,
				   const Real threshold) const
{
  assert (this->initialized());
  assert (other_vector.initialized());
  assert (this->first_local_index() == other_vector.first_local_index());
  assert (this->last_local_index()  == other_vector.last_local_index());

  int rvalue     = -1;
  unsigned int i = first_local_index();

  do
    {
      if ( fabs( (*this)(i) - other_vector(i) ) > threshold )
	rvalue = i;
      else
	i++;
    }
  while (rvalue==-1 && i<last_local_index());

  return rvalue;
}



// Full specialization for Real datatypes
template <>
int NumericVector<Real>::compare (const NumericVector<Real> &other_vector,
				  const Real threshold) const
{
  assert (this->initialized());
  assert (other_vector.initialized());
  assert (this->first_local_index() == other_vector.first_local_index());
  assert (this->last_local_index()  == other_vector.last_local_index());

  int rvalue     = -1;
  unsigned int i = first_local_index();

  do
    {
      if ( fabs( (*this)(i) - other_vector(i) ) > threshold )
	rvalue = i;
      else
	i++;
    }
  while (rvalue==-1 && i<last_local_index());

  return rvalue;
}



// Full specialization for Complex datatypes
template <>
int NumericVector<Complex>::compare (const NumericVector<Complex> &other_vector,
				     const Real threshold) const
{
  assert (this->initialized());
  assert (other_vector.initialized());
  assert (this->first_local_index() == other_vector.first_local_index());
  assert (this->last_local_index()  == other_vector.last_local_index());

  int rvalue     = -1;
  unsigned int i = first_local_index();

  do
    {
      if (( fabs( (*this)(i).real() - other_vector(i).real() ) > threshold ) ||
	  ( fabs( (*this)(i).imag() - other_vector(i).imag() ) > threshold ))
	rvalue = i;
      else
	i++;
    }
  while (rvalue==-1 && i<last_local_index());

  return rvalue;
}



//------------------------------------------------------------------
// Explicit instantiations
template class NumericVector<Number>;

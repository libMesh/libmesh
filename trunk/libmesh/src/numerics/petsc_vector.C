// $Id: petsc_vector.C,v 1.19 2003-05-28 22:03:15 benkirk Exp $

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

// Local Includes
#include "petsc_vector.h"
#include "petsc_matrix.h"

#ifdef HAVE_PETSC

#include "utility.h"
#include "dense_vector.h"




//-----------------------------------------------------------------------
// PetscVector members

// void PetscVector<T>::init (const NumericVector<T>& v, const bool fast)
// {
//   error();
  
//   init (v.local_size(), v.size(), fast);

//   vec = dynamic_cast<const PetscVector<T>&>(v).vec;
// }



template <typename T>
Real PetscVector<T>::l1_norm () const
{
  assert(this->closed());
  
  int ierr=0;
  double value=0.;
  
  ierr = VecNorm (vec, NORM_1, &value);
         CHKERRQ(ierr);
  
  return static_cast<Real>(value);
}



template <typename T>
Real PetscVector<T>::l2_norm () const
{
  assert(this->closed());
  
  int ierr=0;
  double value=0.;
  
  ierr = VecNorm (vec, NORM_2, &value);
         CHKERRQ(ierr);
  
  return static_cast<Real>(value);
}




template <typename T>
Real PetscVector<T>::linfty_norm () const
{
  assert(this->closed());
  
  int ierr=0;
  double value=0.;
  
  ierr = VecNorm (vec, NORM_INFINITY, &value);
         CHKERRQ(ierr);
  
  return static_cast<Real>(value);
}




template <typename T>
NumericVector<T>&
PetscVector<T>::operator += (const NumericVector<T>& v)
{
  assert(this->closed());
  
  this->add(1., v);
  
  return *this;
}



template <typename T>
NumericVector<T>&
PetscVector<T>::operator -= (const NumericVector<T>& v)
{
  assert(this->closed());
  
  this->add(-1., v);
  
  return *this;
}



template <typename T>
void PetscVector<T>::set (const unsigned int i, const T value)
{
  assert(i<size());
  
  int ierr=0;
  int i_val = static_cast<int>(i);
  PetscScalar petsc_value = static_cast<PetscScalar>(value);

  ierr = VecSetValues (vec, 1, &i_val, &petsc_value, INSERT_VALUES);
         CHKERRQ(ierr);

  return;
}



template <typename T>
void PetscVector<T>::add (const unsigned int i, const T value)
{
  assert(i<size());
  
  int ierr=0;
  int i_val = static_cast<int>(i);
  PetscScalar petsc_value = static_cast<PetscScalar>(value);

  ierr = VecSetValues (vec, 1, &i_val, &petsc_value, ADD_VALUES);
         CHKERRQ(ierr);

  return;
}



template <typename T>
void PetscVector<T>::add_vector (const std::vector<T>& v,
				 const std::vector<unsigned int>& dof_indices)
{
  assert (!v.empty());
  assert (v.size() == dof_indices.size());
  
  for (unsigned int i=0; i<v.size(); i++)
    this->add (dof_indices[i], v[i]);
}



template <typename T>
void PetscVector<T>::add_vector (const NumericVector<T>& V,
				 const std::vector<unsigned int>& dof_indices)
{
  assert (V.size() == dof_indices.size());

  for (unsigned int i=0; i<V.size(); i++)
    this->add (dof_indices[i], V(i));
}



template <typename T>
inline
void PetscVector<T>::add_vector (const NumericVector<T>& V_in,
				 SparseMatrix<T>& A_in)
{
  const PetscVector<T>& V = dynamic_cast<const PetscVector<T>&>(V_in);
  PetscMatrix<T>&       A = dynamic_cast<PetscMatrix<T>&>(A_in);
  int ierr=0;

  A.close();

  ierr = MatMultAdd(A.mat, V.vec, vec, vec);
         CHKERRQ(ierr); 
  
  return;
}



template <typename T>
void PetscVector<T>::add_vector (const DenseVector<T>& V,
				 const std::vector<unsigned int>& dof_indices)
{
  assert (V.size() == dof_indices.size());

  for (unsigned int i=0; i<V.size(); i++)
    add (dof_indices[i], V(i));
}



template <typename T>
void PetscVector<T>::add (const T v_in)
{
  int ierr=0;
  PetscScalar* values;
  const PetscScalar v = static_cast<PetscScalar>(v_in);  
  const int n   = static_cast<int>(this->local_size());
  const int fli = static_cast<int>(this->first_local_index());
  
  for (int i=0; i<n; i++)
    {
      ierr = VecGetArray (vec, &values);
  	     CHKERRQ(ierr);
      
      int ig = fli + i;      
      
      PetscScalar value = (values[ig] + v);
      
      ierr = VecRestoreArray (vec, &values);
  	     CHKERRQ(ierr);
      
      ierr = VecSetValues (vec, 1, &ig, &value, INSERT_VALUES);
 	     CHKERRQ(ierr); 
    }
}



template <typename T>
void PetscVector<T>::add (const NumericVector<T>& v)
{
  this->add (1., v);
}



template <typename T>
void PetscVector<T>::add (const T a_in, const NumericVector<T>& v_in)
{
  int ierr = 0;
  PetscScalar a = static_cast<PetscScalar>(a_in);

  const PetscVector<T>& v = dynamic_cast<const PetscVector<T>&>(v_in);
  
  assert(this->size() == v.size());
  
  ierr = VecAXPY(&a, v.vec, vec);
         CHKERRQ(ierr);
}



template <typename T>
void PetscVector<T>::scale (const T factor_in)
{
  int ierr = 0;
  PetscScalar factor = static_cast<PetscScalar>(factor_in);
  
  ierr = VecScale(&factor, vec);
         CHKERRQ(ierr);
}



template <typename T>
NumericVector<T>& 
PetscVector<T>::operator = (const T s_in)
{
  int ierr = 0;
  PetscScalar s = static_cast<PetscScalar>(s_in);

  if (this->size() != 0)
    {
      ierr = VecSet(&s, vec);
             CHKERRQ(ierr);
    }
  
  return *this;
}



template <typename T>
NumericVector<T>&
PetscVector<T>::operator = (const NumericVector<T>& v_in)
{
  const PetscVector<T>& v = dynamic_cast<const PetscVector<T>&>(v_in);

  *this = v;
  
  return *this;
}



template <typename T>
PetscVector<T>&
PetscVector<T>::operator = (const PetscVector<T>& v)
{
  if (v.initialized())
    {
      this->init (v.size(), v.local_size());
      _is_closed      = v._is_closed;
      _is_initialized = v._is_initialized;
  
      if (v.size() != 0)
	{
	  int ierr = 0;
	  
	  ierr = VecCopy (v.vec, vec);
       	         CHKERRQ(ierr);
	}
    }
  
  return *this;
}



template <typename T>
NumericVector<T>&
PetscVector<T>::operator = (const std::vector<T>& v)
{
  const unsigned int nl   = this->local_size();
  const unsigned int ioff = this->first_local_index();
  int ierr=0;
  PetscScalar* values;
      
  /**
   * Case 1:  The vector is the same size of
   * The global vector.  Only add the local components.
   */
  if (this->size() == v.size())
    {
      ierr = VecGetArray (vec, &values);
 	     CHKERRQ(ierr);

      for (unsigned int i=0; i<nl; i++)
	values[i] =  static_cast<PetscScalar>(v[i+ioff]);
      
      ierr = VecRestoreArray (vec, &values);
	     CHKERRQ(ierr);
    }

  /**
   * Case 2: The vector is the same size as our local
   * piece.  Insert directly to the local piece.
   */
  else
    {
      assert (this->local_size() == v.size());

      ierr = VecGetArray (vec, &values);
	     CHKERRQ(ierr);

      for (unsigned int i=0; i<nl; i++)
	values[i] = static_cast<PetscScalar>(v[i]);
      
      ierr = VecRestoreArray (vec, &values);
	     CHKERRQ(ierr);
    }

  return *this;
}



template <typename T>
void PetscVector<T>::localize (NumericVector<T>& v_local_in) const
{
  PetscVector<T>& v_local = dynamic_cast<PetscVector<T>&>(v_local_in);

  assert (v_local.local_size() == this->size());

  int ierr = 0;
  const int n = this->size();

  IS is;
  VecScatter scatter;

  // Create idx, idx[i] = i;
  std::vector<int> idx(n); Utility::iota (idx.begin(), idx.end(), 0);

  // Create the index set & scatter object
  ierr = ISCreateGeneral(PETSC_COMM_WORLD, n, &idx[0], &is);
         CHKERRQ(ierr);

  ierr = VecScatterCreate(vec,         is,
			  v_local.vec, is,
			  &scatter);
         CHKERRQ(ierr);

  // Perform the scatter
  ierr = VecScatterBegin(vec, v_local.vec, INSERT_VALUES,
			 SCATTER_FORWARD, scatter);
         CHKERRQ(ierr);
  
  ierr = VecScatterEnd  (vec, v_local.vec, INSERT_VALUES,
			 SCATTER_FORWARD, scatter);
         CHKERRQ(ierr);

  // Clean up
  ierr = ISDestroy (is);
         CHKERRQ(ierr);
  
  ierr = VecScatterDestroy(scatter);
         CHKERRQ(ierr);
}



template <typename T>
void PetscVector<T>::localize (NumericVector<T>& v_local_in,
			       const std::vector<unsigned int>& send_list) const
{
  PetscVector<T>& v_local = dynamic_cast<PetscVector<T>&>(v_local_in);

  assert (v_local.local_size() == this->size());
  assert (send_list.size()     <= v_local.size());
  
  int ierr=0;
  const int n_sl = send_list.size();

  IS is;
  VecScatter scatter;

  std::vector<int> idx(n_sl);
  
  for (int i=0; i<n_sl; i++)
    idx[i] = static_cast<int>(send_list[i]);
  
  // Create the index set & scatter object
  ierr = ISCreateGeneral(PETSC_COMM_WORLD, n_sl, &idx[0], &is);
         CHKERRQ(ierr);

  ierr = VecScatterCreate(vec,         is,
			  v_local.vec, is,
			  &scatter);
         CHKERRQ(ierr);

  
  // Perform the scatter
  ierr = VecScatterBegin(vec, v_local.vec, INSERT_VALUES,
			 SCATTER_FORWARD, scatter);
         CHKERRQ(ierr);
  
  ierr = VecScatterEnd  (vec, v_local.vec, INSERT_VALUES,
			 SCATTER_FORWARD, scatter);
         CHKERRQ(ierr);

  // Clean up
  ierr = ISDestroy (is);
         CHKERRQ(ierr);
  
  ierr = VecScatterDestroy(scatter);
         CHKERRQ(ierr);
}


template <typename T>
void PetscVector<T>::localize (const unsigned int first_local_idx,
			       const unsigned int last_local_idx,
			       const std::vector<unsigned int>& send_list)
{
  // Only good for serial vectors.
  assert (this->size() == this->local_size());
  assert (last_local_idx > first_local_idx);
  assert (send_list.size() <= this->size());
  assert (last_local_idx < this->size());
  
  const unsigned int size       = this->size();
  const unsigned int local_size = (last_local_idx - first_local_idx + 1);
  int ierr=0;
  PetscScalar *my_values, *their_values;

  
  // Don't bother for serial cases
  if ((first_local_idx == 0) &&
      (local_size == size))
    return;
  
	  
  // Build a parallel vector, initialize it with the local
  // parts of (*this)
  PetscVector<T> parallel_vec;

  parallel_vec.init (size, local_size);


  // Copy part of *this into the parallel_vec
  {
    ierr = VecGetArray (vec, &my_values);
           CHKERRQ(ierr);
  
    ierr = VecGetArray (parallel_vec.vec, &their_values);
           CHKERRQ(ierr);
  
    for (unsigned int i=first_local_idx; i<=last_local_idx; i++)
      their_values[i-first_local_idx] = my_values[i];

    ierr = VecRestoreArray (vec, &my_values);
           CHKERRQ(ierr);

    ierr = VecRestoreArray (parallel_vec.vec, &their_values);
          CHKERRQ(ierr);
  }

  // localize like normal
  parallel_vec.localize (*this, send_list);  
}



// Full specialization for Real datatypes
#ifdef USE_REAL_NUMBERS

template <>
void PetscVector<Real>::localize (std::vector<Real>& v_local) const
{
  int ierr=0;
  const int n  = this->size();
  const int nl = this->local_size();
  PetscScalar *values;

  
  v_local.resize(n);

  
  for (int i=0; i<n; i++)
    v_local[i] = 0.;
  
  // only one processor
  if (n == nl)
    {      
      ierr = VecGetArray (vec, &values);
	     CHKERRQ(ierr);

      for (int i=0; i<n; i++)
	v_local[i] = static_cast<Real>(values[i]);

      ierr = VecRestoreArray (vec, &values);
	     CHKERRQ(ierr);
    }

  // otherwise multiple processors
  else
    {
      unsigned int ioff = first_local_index();
      std::vector<Real> local_values(n, 0.);

      {
	ierr = VecGetArray (vec, &values);
	       CHKERRQ(ierr);
	
	for (int i=0; i<nl; i++)
	  local_values[i+ioff] = static_cast<Real>(values[i]);
	
	ierr = VecRestoreArray (vec, &values);
	       CHKERRQ(ierr);
      }

      if (sizeof(Real) == sizeof(double))
	MPI_Allreduce (&local_values[0], &v_local[0], n, MPI_DOUBLE, MPI_SUM,
		       PETSC_COMM_WORLD);
      
      else if (sizeof(Real) == sizeof(float))
	MPI_Allreduce (&local_values[0], &v_local[0], n, MPI_FLOAT, MPI_SUM,
		       PETSC_COMM_WORLD);

      else
	error();
    }  
}

#endif



// Full specialization for Complex datatypes
#ifdef USE_COMPLEX_NUMBERS

template <>
void PetscVector<Complex>::localize (std::vector<Complex>& v_local) const
{
  int ierr=0;
  const int n  = size();
  const int nl = local_size();
  PetscScalar *values;

  v_local.resize(n);

  
  for (int i=0; i<n; i++)
    v_local[i] = 0.;
  
  // only one processor
  if (n == nl)
    {      
      ierr = VecGetArray (vec, &values);
	     CHKERRQ(ierr);

      for (int i=0; i<n; i++)
	v_local[i] = static_cast<Complex>(values[i]);

      ierr = VecRestoreArray (vec, &values);
	     CHKERRQ(ierr);
    }

  // otherwise multiple processors
  else
    {
      unsigned int ioff = first_local_index();

      /* in here the local values are stored, acting as send buffer for MPI
       * initialize to zero, since we collect using MPI_SUM
       */
      std::vector<Real> real_local_values(n, 0.);
      std::vector<Real> imag_local_values(n, 0.);

      {
	ierr = VecGetArray (vec, &values);
	       CHKERRQ(ierr);
	
	// provide my local share to the real and imag buffers
	for (int i=0; i<nl; i++)
	  {
	    real_local_values[i+ioff] = static_cast<Complex>(values[i]).real();
	    imag_local_values[i+ioff] = static_cast<Complex>(values[i]).imag();
	  }

	ierr = VecRestoreArray (vec, &values);
	       CHKERRQ(ierr);
      }
   
      /* have buffers of the real and imaginary part of v_local.
       * Once MPI_Reduce() collected all the real and imaginary
       * parts in these std::vector<double>, the values can be 
       * copied to v_local
       */
      std::vector<Real> real_v_local(n);
      std::vector<Real> imag_v_local(n);

      // collect entries from other proc's in real_v_local, imag_v_local
      MPI_Allreduce (&real_local_values[0], &real_v_local[0], n, 
		     MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);	

      MPI_Allreduce (&imag_local_values[0], &imag_v_local[0], n, 
		     MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);	

      // copy real_v_local and imag_v_local to v_local
      for (int i=0; i<n; i++)
	v_local[i] = Complex(real_v_local[i], imag_v_local[i]);

    }
}

#endif



// Full specialization for Real datatypes
#ifdef USE_REAL_NUMBERS

template <>
void PetscVector<Real>::localize_to_one (std::vector<Real>& v_local,
					 const unsigned int pid) const
{
  int ierr=0;
  const int n  = size();
  const int nl = local_size();
  PetscScalar *values;

  
  v_local.resize(n);

  
  // only one processor
  if (n == nl)
    {      
      ierr = VecGetArray (vec, &values);
	     CHKERRQ(ierr);

      for (int i=0; i<n; i++)
	v_local[i] = static_cast<Real>(values[i]);

      ierr = VecRestoreArray (vec, &values);
	     CHKERRQ(ierr);
    }

  // otherwise multiple processors
  else
    {
      unsigned int ioff = this->first_local_index();
      std::vector<Real> local_values (n, 0.);
      
      {
	ierr = VecGetArray (vec, &values);
	       CHKERRQ(ierr);
	
	for (int i=0; i<nl; i++)
	  local_values[i+ioff] = static_cast<Real>(values[i]);
	
	ierr = VecRestoreArray (vec, &values);
	       CHKERRQ(ierr);
      }
      

      if (sizeof(Real) == sizeof(double))
	MPI_Reduce (&local_values[0], &v_local[0], n, MPI_DOUBLE, MPI_SUM,
		    pid, PETSC_COMM_WORLD);
      
      else if (sizeof(Real) == sizeof(float))
	MPI_Reduce (&local_values[0], &v_local[0], n, MPI_FLOAT, MPI_SUM,
		    pid, PETSC_COMM_WORLD);

      else
	error();
    }
}

#endif


// Full specialization for Complex datatypes
#ifdef USE_COMPLEX_NUMBERS

template <>
void PetscVector<Complex>::localize_to_one (std::vector<Complex>& v_local,
					    const unsigned int pid) const
{
  int ierr=0;
  const int n  = size();
  const int nl = local_size();
  PetscScalar *values;

  
  v_local.resize(n);

  
  for (int i=0; i<n; i++)
    v_local[i] = 0.;
  
  // only one processor
  if (n == nl)
    {      
      ierr = VecGetArray (vec, &values);
	     CHKERRQ(ierr);

      for (int i=0; i<n; i++)
	v_local[i] = static_cast<Complex>(values[i]);

      ierr = VecRestoreArray (vec, &values);
	     CHKERRQ(ierr);
    }

  // otherwise multiple processors
  else
    {
      unsigned int ioff = this->first_local_index();

      /* in here the local values are stored, acting as send buffer for MPI
       * initialize to zero, since we collect using MPI_SUM
       */
      std::vector<Real> real_local_values(n, 0.);
      std::vector<Real> imag_local_values(n, 0.);

      {
	ierr = VecGetArray (vec, &values);
	       CHKERRQ(ierr);
	
	// provide my local share to the real and imag buffers
	for (int i=0; i<nl; i++)
	  {
	    real_local_values[i+ioff] = static_cast<Complex>(values[i]).real();
	    imag_local_values[i+ioff] = static_cast<Complex>(values[i]).imag();
	  }

	ierr = VecRestoreArray (vec, &values);
	       CHKERRQ(ierr);
      }
   
      /* have buffers of the real and imaginary part of v_local.
       * Once MPI_Reduce() collected all the real and imaginary
       * parts in these std::vector<double>, the values can be 
       * copied to v_local
       */
      std::vector<Real> real_v_local(n);
      std::vector<Real> imag_v_local(n);

      // collect entries from other proc's in real_v_local, imag_v_local
      MPI_Reduce (&real_local_values[0], &real_v_local[0], n, 
		  MPI_DOUBLE, MPI_SUM,
		  pid, PETSC_COMM_WORLD);	

      MPI_Reduce (&imag_local_values[0], &imag_v_local[0], n, 
		  MPI_DOUBLE, MPI_SUM,
		  pid, PETSC_COMM_WORLD);	

      // copy real_v_local and imag_v_local to v_local
      for (int i=0; i<n; i++)
	v_local[i] = Complex(real_v_local[i], imag_v_local[i]);
    }  
}

#endif




//------------------------------------------------------------------
// Explicit instantiations
template class PetscVector<Number>;



#endif // #ifdef HAVE_PETSC

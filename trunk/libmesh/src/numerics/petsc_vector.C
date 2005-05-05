// $Id: petsc_vector.C,v 1.36 2005-05-05 20:20:49 benkirk Exp $

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
  
  ierr = VecNorm (_vec, NORM_1, &value);
         CHKERRABORT(PETSC_COMM_WORLD,ierr);
  
  return static_cast<Real>(value);
}



template <typename T>
Real PetscVector<T>::l2_norm () const
{
  assert(this->closed());
  
  int ierr=0;
  double value=0.;
  
  ierr = VecNorm (_vec, NORM_2, &value);
         CHKERRABORT(PETSC_COMM_WORLD,ierr);
  
  return static_cast<Real>(value);
}




template <typename T>
Real PetscVector<T>::linfty_norm () const
{
  assert(this->closed());
  
  int ierr=0;
  double value=0.;
  
  ierr = VecNorm (_vec, NORM_INFINITY, &value);
         CHKERRABORT(PETSC_COMM_WORLD,ierr);
  
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

  ierr = VecSetValues (_vec, 1, &i_val, &petsc_value, INSERT_VALUES);
         CHKERRABORT(PETSC_COMM_WORLD,ierr);
}



template <typename T>
void PetscVector<T>::add (const unsigned int i, const T value)
{
  assert(i<size());
  
  int ierr=0;
  int i_val = static_cast<int>(i);
  PetscScalar petsc_value = static_cast<PetscScalar>(value);

  ierr = VecSetValues (_vec, 1, &i_val, &petsc_value, ADD_VALUES);
         CHKERRABORT(PETSC_COMM_WORLD,ierr);
}



template <typename T>
void PetscVector<T>::add_vector (const std::vector<T>& v,
				 const std::vector<unsigned int>& dof_indices)
{
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
void PetscVector<T>::add_vector (const NumericVector<T>& V_in,
				 const SparseMatrix<T>& A_in)
{
  const PetscVector<T>* V = dynamic_cast<const PetscVector<T>*>(&V_in);
  const PetscMatrix<T>* A = dynamic_cast<const PetscMatrix<T>*>(&A_in);

  assert (V != NULL);
  assert (A != NULL);
  
  int ierr=0;

  A->close();

  // The const_cast<> is not elegant, but it is required since PETSc
  // is not const-correct.  
  ierr = MatMultAdd(const_cast<PetscMatrix<T>*>(A)->mat(), V->_vec, _vec, _vec);
         CHKERRABORT(PETSC_COMM_WORLD,ierr); 
}



template <typename T>
void PetscVector<T>::add_vector (const DenseVector<T>& V,
				 const std::vector<unsigned int>& dof_indices)
{
  assert (V.size() == dof_indices.size());

  for (unsigned int i=0; i<V.size(); i++)
    this->add (dof_indices[i], V(i));
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
      ierr = VecGetArray (_vec, &values);
  	     CHKERRABORT(PETSC_COMM_WORLD,ierr);
      
      int ig = fli + i;      
      
      PetscScalar value = (values[ig] + v);
      
      ierr = VecRestoreArray (_vec, &values);
  	     CHKERRABORT(PETSC_COMM_WORLD,ierr);
      
      ierr = VecSetValues (_vec, 1, &ig, &value, INSERT_VALUES);
 	     CHKERRABORT(PETSC_COMM_WORLD,ierr); 
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

  const PetscVector<T>* v = dynamic_cast<const PetscVector<T>*>(&v_in);

  assert (v != NULL);
  assert(this->size() == v->size());
  
// 2.2.x & earlier style
#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 2)
  
  ierr = VecAXPY(&a, v->_vec, _vec);
         CHKERRABORT(PETSC_COMM_WORLD,ierr);
	 
// 2.3.x & later style
#else
  
  ierr = VecAXPY(_vec, a, v->_vec);
         CHKERRABORT(PETSC_COMM_WORLD,ierr);
	 
#endif
}



template <typename T>
void PetscVector<T>::insert (const std::vector<T>& v,
			     const std::vector<unsigned int>& dof_indices)
{
  assert (v.size() == dof_indices.size());

  for (unsigned int i=0; i<v.size(); i++)
    this->set (dof_indices[i], v[i]);
}



template <typename T>
void PetscVector<T>::insert (const NumericVector<T>& V,
			     const std::vector<unsigned int>& dof_indices)
{
  assert (V.size() == dof_indices.size());

  for (unsigned int i=0; i<V.size(); i++)
    this->set (dof_indices[i], V(i));
}



template <typename T>
void PetscVector<T>::insert (const DenseVector<T>& V,
			     const std::vector<unsigned int>& dof_indices)
{
  assert (V.size() == dof_indices.size());

  for (unsigned int i=0; i<V.size(); i++)
    this->set (dof_indices[i], V(i));
}



template <typename T>
void PetscVector<T>::scale (const T factor_in)
{
  int ierr = 0;
  PetscScalar factor = static_cast<PetscScalar>(factor_in);
  
// 2.2.x & earlier style
#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 2)
  
  ierr = VecScale(&factor, _vec);
         CHKERRABORT(PETSC_COMM_WORLD,ierr);

// 2.3.x & later style	 
#else
  
  ierr = VecScale(_vec, factor);
         CHKERRABORT(PETSC_COMM_WORLD,ierr);

#endif
}



template <typename T>
NumericVector<T>& 
PetscVector<T>::operator = (const T s_in)
{
  int ierr = 0;
  PetscScalar s = static_cast<PetscScalar>(s_in);

  if (this->size() != 0)
    {
// 2.2.x & earlier style
#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 2)
  
      ierr = VecSet(&s, _vec);
             CHKERRABORT(PETSC_COMM_WORLD,ierr);
	     
// 2.3.x & later style	 
#else

      ierr = VecSet(_vec, s);
             CHKERRABORT(PETSC_COMM_WORLD,ierr);
	     
#endif
    }
  
  return *this;
}



template <typename T>
NumericVector<T>&
PetscVector<T>::operator = (const NumericVector<T>& v_in)
{
  const PetscVector<T>* v = dynamic_cast<const PetscVector<T>*>(&v_in);

  assert (v != NULL);
  
  *this = *v;
  
  return *this;
}



template <typename T>
PetscVector<T>&
PetscVector<T>::operator = (const PetscVector<T>& v)
{
  if (v.initialized())
    {
      this->init (v.size(), v.local_size());
      this->_is_closed      = v._is_closed;
      this->_is_initialized = v._is_initialized;
  
      if (v.size() != 0)
	{
	  int ierr = 0;

	  ierr = VecCopy (v._vec, this->_vec);
       	         CHKERRABORT(PETSC_COMM_WORLD,ierr);
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
      ierr = VecGetArray (_vec, &values);
 	     CHKERRABORT(PETSC_COMM_WORLD,ierr);

      for (unsigned int i=0; i<nl; i++)
	values[i] =  static_cast<PetscScalar>(v[i+ioff]);
      
      ierr = VecRestoreArray (_vec, &values);
	     CHKERRABORT(PETSC_COMM_WORLD,ierr);
    }

  /**
   * Case 2: The vector is the same size as our local
   * piece.  Insert directly to the local piece.
   */
  else
    {
      assert (this->local_size() == v.size());

      ierr = VecGetArray (_vec, &values);
	     CHKERRABORT(PETSC_COMM_WORLD,ierr);

      for (unsigned int i=0; i<nl; i++)
	values[i] = static_cast<PetscScalar>(v[i]);
      
      ierr = VecRestoreArray (_vec, &values);
	     CHKERRABORT(PETSC_COMM_WORLD,ierr);
    }

  return *this;
}



template <typename T>
void PetscVector<T>::localize (NumericVector<T>& v_local_in) const
{
  PetscVector<T>* v_local = dynamic_cast<PetscVector<T>*>(&v_local_in);

  assert (v_local != NULL);
  assert (v_local->local_size() == this->size());

  int ierr = 0;
  const int n = this->size();

  IS is;
  VecScatter scatter;

  // Create idx, idx[i] = i;
  std::vector<int> idx(n); Utility::iota (idx.begin(), idx.end(), 0);

  // Create the index set & scatter object
  ierr = ISCreateGeneral(PETSC_COMM_WORLD, n, &idx[0], &is);
         CHKERRABORT(PETSC_COMM_WORLD,ierr);

  ierr = VecScatterCreate(_vec,          is,
			  v_local->_vec, is,
			  &scatter);
         CHKERRABORT(PETSC_COMM_WORLD,ierr);

  // Perform the scatter
  ierr = VecScatterBegin(_vec, v_local->_vec, INSERT_VALUES,
			 SCATTER_FORWARD, scatter);
         CHKERRABORT(PETSC_COMM_WORLD,ierr);
  
  ierr = VecScatterEnd  (_vec, v_local->_vec, INSERT_VALUES,
			 SCATTER_FORWARD, scatter);
         CHKERRABORT(PETSC_COMM_WORLD,ierr);

  // Clean up
  ierr = ISDestroy (is);
         CHKERRABORT(PETSC_COMM_WORLD,ierr);
  
  ierr = VecScatterDestroy(scatter);
         CHKERRABORT(PETSC_COMM_WORLD,ierr);
}



template <typename T>
void PetscVector<T>::localize (NumericVector<T>& v_local_in,
			       const std::vector<unsigned int>& send_list) const
{
  PetscVector<T>* v_local = dynamic_cast<PetscVector<T>*>(&v_local_in);

  assert (v_local != NULL);
  assert (v_local->local_size() == this->size());
  assert (send_list.size()     <= v_local->size());
  
  int ierr=0;
  const int n_sl = send_list.size();

  IS is;
  VecScatter scatter;

  std::vector<int> idx(n_sl);
  
  for (int i=0; i<n_sl; i++)
    idx[i] = static_cast<int>(send_list[i]);
  
  // Create the index set & scatter object
  ierr = ISCreateGeneral(PETSC_COMM_WORLD, n_sl, &idx[0], &is);
         CHKERRABORT(PETSC_COMM_WORLD,ierr);

  ierr = VecScatterCreate(_vec,          is,
			  v_local->_vec, is,
			  &scatter);
         CHKERRABORT(PETSC_COMM_WORLD,ierr);

  
  // Perform the scatter
  ierr = VecScatterBegin(_vec, v_local->_vec, INSERT_VALUES,
			 SCATTER_FORWARD, scatter);
         CHKERRABORT(PETSC_COMM_WORLD,ierr);
  
  ierr = VecScatterEnd  (_vec, v_local->_vec, INSERT_VALUES,
			 SCATTER_FORWARD, scatter);
         CHKERRABORT(PETSC_COMM_WORLD,ierr);

  // Clean up
  ierr = ISDestroy (is);
         CHKERRABORT(PETSC_COMM_WORLD,ierr);
  
  ierr = VecScatterDestroy(scatter);
         CHKERRABORT(PETSC_COMM_WORLD,ierr);
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
    IS is;
    VecScatter scatter;

    // Create idx, idx[i] = i+first_local_idx;
    std::vector<int> idx(local_size);
    Utility::iota (idx.begin(), idx.end(), first_local_idx);

    // Create the index set & scatter object
    ierr = ISCreateGeneral(PETSC_COMM_WORLD, local_size, &idx[0], &is); 
           CHKERRABORT(PETSC_COMM_WORLD,ierr);

    ierr = VecScatterCreate(_vec,              is,
			    parallel_vec._vec, is,
			    &scatter);
           CHKERRABORT(PETSC_COMM_WORLD,ierr);

    // Perform the scatter
    ierr = VecScatterBegin(_vec, parallel_vec._vec, INSERT_VALUES,
			   SCATTER_FORWARD, scatter);
           CHKERRABORT(PETSC_COMM_WORLD,ierr);
  
    ierr = VecScatterEnd  (_vec, parallel_vec._vec, INSERT_VALUES,
			   SCATTER_FORWARD, scatter);
           CHKERRABORT(PETSC_COMM_WORLD,ierr);

    // Clean up
    ierr = ISDestroy (is);
           CHKERRABORT(PETSC_COMM_WORLD,ierr);
  
    ierr = VecScatterDestroy(scatter);
           CHKERRABORT(PETSC_COMM_WORLD,ierr);
  }

  // localize like normal
  parallel_vec.close();
  parallel_vec.localize (*this, send_list);
  this->close();
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
      ierr = VecGetArray (_vec, &values);
	     CHKERRABORT(PETSC_COMM_WORLD,ierr);

      for (int i=0; i<n; i++)
	v_local[i] = static_cast<Real>(values[i]);

      ierr = VecRestoreArray (_vec, &values);
	     CHKERRABORT(PETSC_COMM_WORLD,ierr);
    }

  // otherwise multiple processors
  else
    {
      unsigned int ioff = first_local_index();
      std::vector<Real> local_values(n, 0.);

      {
	ierr = VecGetArray (_vec, &values);
	       CHKERRABORT(PETSC_COMM_WORLD,ierr);
	
	for (int i=0; i<nl; i++)
	  local_values[i+ioff] = static_cast<Real>(values[i]);
	
	ierr = VecRestoreArray (_vec, &values);
	       CHKERRABORT(PETSC_COMM_WORLD,ierr);
      }

      MPI_Allreduce (&local_values[0], &v_local[0], n, MPI_REAL, MPI_SUM,
		     PETSC_COMM_WORLD);
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
      ierr = VecGetArray (_vec, &values);
	     CHKERRABORT(PETSC_COMM_WORLD,ierr);

      for (int i=0; i<n; i++)
	v_local[i] = static_cast<Complex>(values[i]);

      ierr = VecRestoreArray (_vec, &values);
	     CHKERRABORT(PETSC_COMM_WORLD,ierr);
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
	ierr = VecGetArray (_vec, &values);
	       CHKERRABORT(PETSC_COMM_WORLD,ierr);
	
	// provide my local share to the real and imag buffers
	for (int i=0; i<nl; i++)
	  {
	    real_local_values[i+ioff] = static_cast<Complex>(values[i]).real();
	    imag_local_values[i+ioff] = static_cast<Complex>(values[i]).imag();
	  }

	ierr = VecRestoreArray (_vec, &values);
	       CHKERRABORT(PETSC_COMM_WORLD,ierr);
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
      ierr = VecGetArray (_vec, &values);
	     CHKERRABORT(PETSC_COMM_WORLD,ierr);

      for (int i=0; i<n; i++)
	v_local[i] = static_cast<Real>(values[i]);

      ierr = VecRestoreArray (_vec, &values);
	     CHKERRABORT(PETSC_COMM_WORLD,ierr);
    }

  // otherwise multiple processors
  else
    {
      unsigned int ioff = this->first_local_index();
      std::vector<Real> local_values (n, 0.);
      
      {
	ierr = VecGetArray (_vec, &values);
	       CHKERRABORT(PETSC_COMM_WORLD,ierr);
	
	for (int i=0; i<nl; i++)
	  local_values[i+ioff] = static_cast<Real>(values[i]);
	
	ierr = VecRestoreArray (_vec, &values);
	       CHKERRABORT(PETSC_COMM_WORLD,ierr);
      }
      

      MPI_Reduce (&local_values[0], &v_local[0], n, MPI_REAL, MPI_SUM,
		  pid, PETSC_COMM_WORLD);
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
      ierr = VecGetArray (_vec, &values);
	     CHKERRABORT(PETSC_COMM_WORLD,ierr);

      for (int i=0; i<n; i++)
	v_local[i] = static_cast<Complex>(values[i]);

      ierr = VecRestoreArray (_vec, &values);
	     CHKERRABORT(PETSC_COMM_WORLD,ierr);
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
	ierr = VecGetArray (_vec, &values);
	       CHKERRABORT(PETSC_COMM_WORLD,ierr);
	
	// provide my local share to the real and imag buffers
	for (int i=0; i<nl; i++)
	  {
	    real_local_values[i+ioff] = static_cast<Complex>(values[i]).real();
	    imag_local_values[i+ioff] = static_cast<Complex>(values[i]).imag();
	  }

	ierr = VecRestoreArray (_vec, &values);
	       CHKERRABORT(PETSC_COMM_WORLD,ierr);
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



template <typename T>
void PetscVector<T>::print_matlab (const std::string name) const
{
  assert (this->initialized());
  assert (this->closed());
  
  int ierr=0; 
  PetscViewer petsc_viewer;


  ierr = PetscViewerCreate (PETSC_COMM_WORLD,
			    &petsc_viewer);
         CHKERRABORT(PETSC_COMM_WORLD,ierr);

  /**
   * Create an ASCII file containing the matrix
   * if a filename was provided.  
   */
  if (name != "NULL")
    {
      ierr = PetscViewerASCIIOpen( PETSC_COMM_WORLD,
				   name.c_str(),
				   &petsc_viewer);
             CHKERRABORT(PETSC_COMM_WORLD,ierr);
      
      ierr = PetscViewerSetFormat (petsc_viewer,
				   PETSC_VIEWER_ASCII_MATLAB);
             CHKERRABORT(PETSC_COMM_WORLD,ierr);
  
      ierr = VecView (_vec, petsc_viewer);
             CHKERRABORT(PETSC_COMM_WORLD,ierr);
    }

  /**
   * Otherwise the matrix will be dumped to the screen.
   */
  else
    {
      ierr = PetscViewerSetFormat (PETSC_VIEWER_STDOUT_WORLD,
				   PETSC_VIEWER_ASCII_MATLAB);
             CHKERRABORT(PETSC_COMM_WORLD,ierr);
  
      ierr = VecView (_vec, PETSC_VIEWER_STDOUT_WORLD);
             CHKERRABORT(PETSC_COMM_WORLD,ierr);
    }


  /**
   * Destroy the viewer.
   */
  ierr = PetscViewerDestroy (petsc_viewer);
         CHKERRABORT(PETSC_COMM_WORLD,ierr);
}





template <typename T>
void PetscVector<T>::create_subvector(NumericVector<T>& subvector,
				      const std::vector<unsigned int>& rows) const
{
  // PETSc data structures
  IS parent_is, subvector_is;
  VecScatter scatter;
  int ierr = 0;
  
  // Make sure the passed int subvector is really a PetscVector
  PetscVector<T>* petsc_subvector = dynamic_cast<PetscVector<T>*>(&subvector);
  assert(petsc_subvector != NULL);
  
  // If the petsc_subvector is already initialized, we assume that the
  // user has already allocated the *correct* amount of space for it.
  // If not, we use the appropriate PETSc routines to initialize it.
  if (!petsc_subvector->initialized())
    {
      // Initialize the petsc_subvector to have enough space to hold
      // the entries which will be scattered into it.  Note: such an
      // init() function (where we let PETSc decide the number of local
      // entries) is not currently offered by the PetscVector
      // class.  Should we differentiate here between sequential and
      // parallel vector creation based on libMesh::n_processors() ?
      ierr = VecCreateMPI(PETSC_COMM_WORLD,
			  PETSC_DECIDE,          // n_local
			  rows.size(),           // n_global
			  &(petsc_subvector->_vec)); CHKERRABORT(PETSC_COMM_WORLD,ierr);

      ierr = VecSetFromOptions (petsc_subvector->_vec); CHKERRABORT(PETSC_COMM_WORLD,ierr);

      // Mark the subvector as initialized
      petsc_subvector->_is_initialized = true;
    }
  
  // Use iota to fill an array with entries [0,1,2,3,4,...rows.size()]
  std::vector<int> idx(rows.size());
  Utility::iota (idx.begin(), idx.end(), 0);

  // Construct index sets
  ierr = ISCreateGeneral(PETSC_COMM_WORLD,
			 rows.size(),
			 (int*) &rows[0],
			 &parent_is); CHKERRABORT(PETSC_COMM_WORLD,ierr);

  ierr = ISCreateGeneral(PETSC_COMM_WORLD,
			 rows.size(),
			 (int*) &idx[0],
			 &subvector_is); CHKERRABORT(PETSC_COMM_WORLD,ierr);

  // Construct the scatter object
  ierr = VecScatterCreate(this->_vec,
			  parent_is,
			  petsc_subvector->_vec,
			  subvector_is,
			  &scatter); CHKERRABORT(PETSC_COMM_WORLD,ierr);

  // Actually perform the scatter
  ierr = VecScatterBegin(this->_vec,
			 petsc_subvector->_vec,
			 INSERT_VALUES,
			 SCATTER_FORWARD,
			 scatter); CHKERRABORT(PETSC_COMM_WORLD,ierr);

  ierr = VecScatterEnd(this->_vec,
		       petsc_subvector->_vec,
		       INSERT_VALUES,
		       SCATTER_FORWARD,
		       scatter); CHKERRABORT(PETSC_COMM_WORLD,ierr);
  
  // Clean up 
  ierr = ISDestroy(parent_is);       CHKERRABORT(PETSC_COMM_WORLD,ierr);
  ierr = ISDestroy(subvector_is);    CHKERRABORT(PETSC_COMM_WORLD,ierr);
  ierr = VecScatterDestroy(scatter); CHKERRABORT(PETSC_COMM_WORLD,ierr); 

}




//------------------------------------------------------------------
// Explicit instantiations
template class PetscVector<Number>;



#endif // #ifdef HAVE_PETSC

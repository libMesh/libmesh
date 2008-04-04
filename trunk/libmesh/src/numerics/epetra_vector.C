// $Id: epetra_vector.C 2606 2008-01-23 20:21:47Z roystgnr $

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

// Local Includes
#include "epetra_vector.h"

#ifdef HAVE_TRILINOS

#include "parallel.h"
#include "utility.h"
#include "dense_vector.h"
#include "epetra_macro.h"

//-----------------------------------------------------------------------
// EpetraVector members

// void EpetraVector<T>::init (const NumericVector<T>& v, const bool fast)
// {
//   error();
  
//   init (v.local_size(), v.size(), fast);

//   vec = dynamic_cast<const EpetraVector<T>&>(v).vec;
// }

template <typename T>
T EpetraVector<T>::sum () const
{
  assert(this->closed());
  
  int ierr=0;
  EpetraScalar value=0.;
  
  ierr = VecSum (_vec, &value);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);
  
  return static_cast<T>(value);
}


template <typename T>
Real EpetraVector<T>::l1_norm () const
{
  assert(this->closed());
  
  int ierr=0;
  EpetraReal value=0.;
  
  ierr = VecNorm (_vec, NORM_1, &value);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);
  
  return static_cast<Real>(value);
}



template <typename T>
Real EpetraVector<T>::l2_norm () const
{
  assert(this->closed());
  
  int ierr=0;
  EpetraReal value=0.;
  
  ierr = VecNorm (_vec, NORM_2, &value);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);
  
  return static_cast<Real>(value);
}




template <typename T>
Real EpetraVector<T>::linfty_norm () const
{
  assert(this->closed());
  
  int ierr=0;
  EpetraReal value=0.;
  
  ierr = VecNorm (_vec, NORM_INFINITY, &value);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);
  
  return static_cast<Real>(value);
}




template <typename T>
NumericVector<T>&
EpetraVector<T>::operator += (const NumericVector<T>& v)
{
  assert(this->closed());
  
  this->add(1., v);
  
  return *this;
}



template <typename T>
NumericVector<T>&
EpetraVector<T>::operator -= (const NumericVector<T>& v)
{
  assert(this->closed());
  
  this->add(-1., v);
  
  return *this;
}



template <typename T>
void EpetraVector<T>::set (const unsigned int i, const T value)
{
  assert(i<size());
  
  int ierr=0;
  int i_val = static_cast<int>(i);
  EpetraScalar epetra_value = static_cast<EpetraScalar>(value);

  ierr = VecSetValues (_vec, 1, &i_val, &epetra_value, INSERT_VALUES);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);

  this->_is_closed = false;
}



template <typename T>
void EpetraVector<T>::add (const unsigned int i, const T value)
{
  assert(i<size());
  
  int ierr=0;
  int i_val = static_cast<int>(i);
  EpetraScalar epetra_value = static_cast<EpetraScalar>(value);

  ierr = VecSetValues (_vec, 1, &i_val, &epetra_value, ADD_VALUES);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);

  this->_is_closed = false;
}



template <typename T>
void EpetraVector<T>::add_vector (const std::vector<T>& v,
				 const std::vector<unsigned int>& dof_indices)
{
  assert (v.size() == dof_indices.size());

  for (unsigned int i=0; i<v.size(); i++)
    this->add (dof_indices[i], v[i]);
}



template <typename T>
void EpetraVector<T>::add_vector (const NumericVector<T>& V,
				 const std::vector<unsigned int>& dof_indices)
{
  assert (V.size() == dof_indices.size());

  for (unsigned int i=0; i<V.size(); i++)
    this->add (dof_indices[i], V(i));
}



template <typename T>
void EpetraVector<T>::add_vector (const NumericVector<T>& V_in,
				 const SparseMatrix<T>& A_in)
{
  const EpetraVector<T>* V = dynamic_cast<const EpetraVector<T>*>(&V_in);
  const EpetraMatrix<T>* A = dynamic_cast<const EpetraMatrix<T>*>(&A_in);

  assert (V != NULL);
  assert (A != NULL);
  
  int ierr=0;

  A->close();

  // The const_cast<> is not elegant, but it is required since Epetra
  // is not const-correct.  
  ierr = MatMultAdd(const_cast<EpetraMatrix<T>*>(A)->mat(), V->_vec, _vec, _vec);
         CHKERRABORT(libMesh::COMM_WORLD,ierr); 
}



template <typename T>
void EpetraVector<T>::add_vector (const DenseVector<T>& V,
				 const std::vector<unsigned int>& dof_indices)
{
  assert (V.size() == dof_indices.size());

  for (unsigned int i=0; i<V.size(); i++)
    this->add (dof_indices[i], V(i));
}



template <typename T>
void EpetraVector<T>::add (const T v_in)
{
  int ierr=0;
  EpetraScalar* values;
  const EpetraScalar v = static_cast<EpetraScalar>(v_in);  
  const int n   = static_cast<int>(this->local_size());
  const int fli = static_cast<int>(this->first_local_index());
  
  for (int i=0; i<n; i++)
    {
      ierr = VecGetArray (_vec, &values);
  	     CHKERRABORT(libMesh::COMM_WORLD,ierr);
      
      int ig = fli + i;      
      
      EpetraScalar value = (values[ig] + v);
      
      ierr = VecRestoreArray (_vec, &values);
  	     CHKERRABORT(libMesh::COMM_WORLD,ierr);
      
      ierr = VecSetValues (_vec, 1, &ig, &value, INSERT_VALUES);
 	     CHKERRABORT(libMesh::COMM_WORLD,ierr); 
    }

  this->_is_closed = false;
}



template <typename T>
void EpetraVector<T>::add (const NumericVector<T>& v)
{
  this->add (1., v);
}



template <typename T>
void EpetraVector<T>::add (const T a_in, const NumericVector<T>& v_in)
{
  int ierr = 0;
  EpetraScalar a = static_cast<EpetraScalar>(a_in);

  const EpetraVector<T>* v = dynamic_cast<const EpetraVector<T>*>(&v_in);

  assert (v != NULL);
  assert(this->size() == v->size());
  

#if EPETRA_VERSION_LESS_THAN(2,3,0)
	 
  // 2.2.x & earlier style
  ierr = VecAXPY(&a, v->_vec, _vec);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);

#else
	 
  // 2.3.x & later style
  ierr = VecAXPY(_vec, a, v->_vec);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);
	 
#endif
}



template <typename T>
void EpetraVector<T>::insert (const std::vector<T>& v,
			     const std::vector<unsigned int>& dof_indices)
{
  assert (v.size() == dof_indices.size());

  for (unsigned int i=0; i<v.size(); i++)
    this->set (dof_indices[i], v[i]);
}



template <typename T>
void EpetraVector<T>::insert (const NumericVector<T>& V,
			     const std::vector<unsigned int>& dof_indices)
{
  assert (V.size() == dof_indices.size());

  for (unsigned int i=0; i<V.size(); i++)
    this->set (dof_indices[i], V(i));
}



template <typename T>
void EpetraVector<T>::insert (const DenseVector<T>& V,
			     const std::vector<unsigned int>& dof_indices)
{
  assert (V.size() == dof_indices.size());

  for (unsigned int i=0; i<V.size(); i++)
    this->set (dof_indices[i], V(i));
}



template <typename T>
void EpetraVector<T>::scale (const T factor_in)
{
  int ierr = 0;
  EpetraScalar factor = static_cast<EpetraScalar>(factor_in);
  
#if EPETRA_VERSION_LESS_THAN(2,3,0)

  // 2.2.x & earlier style
  ierr = VecScale(&factor, _vec);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);

#else
  
  // 2.3.x & later style	 
  ierr = VecScale(_vec, factor);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);

#endif
}




template <typename T>
T EpetraVector<T>::dot (const NumericVector<T>& V) const
{
  // Error flag
  int ierr = 0;
  
  // Return value
  EpetraScalar value=0.;
  
  // Make sure the NumericVector passed in is really a EpetraVector
  const EpetraVector<T>* v = dynamic_cast<const EpetraVector<T>*>(&V);
  assert (v != NULL);

  // 2.3.x (at least) style.  Untested for previous versions.
  ierr = VecDot(this->_vec, v->_vec, &value);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);

  return static_cast<T>(value);
}




template <typename T>
NumericVector<T>& 
EpetraVector<T>::operator = (const T s_in)
{
  assert(this->closed());

  int ierr = 0;
  EpetraScalar s = static_cast<EpetraScalar>(s_in);

  if (this->size() != 0)
    {
#if EPETRA_VERSION_LESS_THAN(2,3,0)
      
  // 2.2.x & earlier style
  ierr = VecSet(&s, _vec);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);
	     
#else

  // 2.3.x & later style	 
  ierr = VecSet(_vec, s);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);
	     
#endif
    }
  
  return *this;
}



template <typename T>
NumericVector<T>&
EpetraVector<T>::operator = (const NumericVector<T>& v_in)
{
  const EpetraVector<T>* v = dynamic_cast<const EpetraVector<T>*>(&v_in);

  assert (v != NULL);
  
  *this = *v;
  
  return *this;
}



template <typename T>
EpetraVector<T>&
EpetraVector<T>::operator = (const EpetraVector<T>& v)
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
       	         CHKERRABORT(libMesh::COMM_WORLD,ierr);
	}
    }
  
  return *this;
}



template <typename T>
NumericVector<T>&
EpetraVector<T>::operator = (const std::vector<T>& v)
{
  const unsigned int nl   = this->local_size();
  const unsigned int ioff = this->first_local_index();
  int ierr=0;
  EpetraScalar* values;
      
  /**
   * Case 1:  The vector is the same size of
   * The global vector.  Only add the local components.
   */
  if (this->size() == v.size())
    {
      ierr = VecGetArray (_vec, &values);
 	     CHKERRABORT(libMesh::COMM_WORLD,ierr);

      for (unsigned int i=0; i<nl; i++)
	values[i] =  static_cast<EpetraScalar>(v[i+ioff]);
      
      ierr = VecRestoreArray (_vec, &values);
	     CHKERRABORT(libMesh::COMM_WORLD,ierr);
    }

  /**
   * Case 2: The vector is the same size as our local
   * piece.  Insert directly to the local piece.
   */
  else
    {
      assert (this->local_size() == v.size());

      ierr = VecGetArray (_vec, &values);
	     CHKERRABORT(libMesh::COMM_WORLD,ierr);

      for (unsigned int i=0; i<nl; i++)
	values[i] = static_cast<EpetraScalar>(v[i]);
      
      ierr = VecRestoreArray (_vec, &values);
	     CHKERRABORT(libMesh::COMM_WORLD,ierr);
    }

  return *this;
}



template <typename T>
void EpetraVector<T>::localize (NumericVector<T>& v_local_in) const
{
  EpetraVector<T>* v_local = dynamic_cast<EpetraVector<T>*>(&v_local_in);

  assert (v_local != NULL);
  assert (v_local->local_size() == this->size());

  int ierr = 0;
  const int n = this->size();

  IS is;
  VecScatter scatter;

  // Create idx, idx[i] = i;
  std::vector<int> idx(n); Utility::iota (idx.begin(), idx.end(), 0);

  // Create the index set & scatter object
  ierr = ISCreateGeneral(libMesh::COMM_WORLD, n, &idx[0], &is);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);

  ierr = VecScatterCreate(_vec,          is,
			  v_local->_vec, is,
			  &scatter);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);

  // Perform the scatter
#if EPETRA_VERSION_LESS_THAN(2,3,3)
	 
  ierr = VecScatterBegin(_vec, v_local->_vec, INSERT_VALUES,
			 SCATTER_FORWARD, scatter);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);
  
  ierr = VecScatterEnd  (_vec, v_local->_vec, INSERT_VALUES,
			 SCATTER_FORWARD, scatter);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);
#else
  // API argument order change in Epetra 2.3.3
  ierr = VecScatterBegin(scatter, _vec, v_local->_vec,
			 INSERT_VALUES, SCATTER_FORWARD);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);
  
  ierr = VecScatterEnd  (scatter, _vec, v_local->_vec,
                         INSERT_VALUES, SCATTER_FORWARD);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);
#endif

  // Clean up
  ierr = ISDestroy (is);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);
  
  ierr = VecScatterDestroy(scatter);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);
}



template <typename T>
void EpetraVector<T>::localize (NumericVector<T>& v_local_in,
			       const std::vector<unsigned int>& send_list) const
{
  EpetraVector<T>* v_local = dynamic_cast<EpetraVector<T>*>(&v_local_in);

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
  if (idx.empty())
    ierr = ISCreateGeneral(libMesh::COMM_WORLD, n_sl, EPETRA_NULL, &is);
  else
    ierr = ISCreateGeneral(libMesh::COMM_WORLD, n_sl, &idx[0], &is);
           CHKERRABORT(libMesh::COMM_WORLD,ierr);

  ierr = VecScatterCreate(_vec,          is,
			  v_local->_vec, is,
			  &scatter);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);

  
  // Perform the scatter
#if EPETRA_VERSION_LESS_THAN(2,3,3)
	 
  ierr = VecScatterBegin(_vec, v_local->_vec, INSERT_VALUES,
			 SCATTER_FORWARD, scatter);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);
  
  ierr = VecScatterEnd  (_vec, v_local->_vec, INSERT_VALUES,
			 SCATTER_FORWARD, scatter);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);

#else
	 
  // API argument order change in Epetra 2.3.3
  ierr = VecScatterBegin(scatter, _vec, v_local->_vec,
                         INSERT_VALUES, SCATTER_FORWARD);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);
  
  ierr = VecScatterEnd  (scatter, _vec, v_local->_vec,
                         INSERT_VALUES, SCATTER_FORWARD);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);

#endif
	 

  // Clean up
  ierr = ISDestroy (is);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);
  
  ierr = VecScatterDestroy(scatter);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);
}


template <typename T>
void EpetraVector<T>::localize (const unsigned int first_local_idx,
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
  EpetraVector<T> parallel_vec;

  parallel_vec.init (size, local_size);


  // Copy part of *this into the parallel_vec
  {
    IS is;
    VecScatter scatter;

    // Create idx, idx[i] = i+first_local_idx;
    std::vector<int> idx(local_size);
    Utility::iota (idx.begin(), idx.end(), first_local_idx);

    // Create the index set & scatter object
    ierr = ISCreateGeneral(libMesh::COMM_WORLD, local_size, &idx[0], &is); 
           CHKERRABORT(libMesh::COMM_WORLD,ierr);

    ierr = VecScatterCreate(_vec,              is,
			    parallel_vec._vec, is,
			    &scatter);
           CHKERRABORT(libMesh::COMM_WORLD,ierr);

    // Perform the scatter
#if EPETRA_VERSION_LESS_THAN(2,3,3)

    ierr = VecScatterBegin(_vec, parallel_vec._vec, INSERT_VALUES,
			   SCATTER_FORWARD, scatter);
           CHKERRABORT(libMesh::COMM_WORLD,ierr);
  
    ierr = VecScatterEnd  (_vec, parallel_vec._vec, INSERT_VALUES,
			   SCATTER_FORWARD, scatter);
           CHKERRABORT(libMesh::COMM_WORLD,ierr);

#else
	   
      // API argument order change in Epetra 2.3.3
    ierr = VecScatterBegin(scatter, _vec, parallel_vec._vec,
			   INSERT_VALUES, SCATTER_FORWARD);
           CHKERRABORT(libMesh::COMM_WORLD,ierr);
  
    ierr = VecScatterEnd  (scatter, _vec, parallel_vec._vec,
			   INSERT_VALUES, SCATTER_FORWARD);
           CHKERRABORT(libMesh::COMM_WORLD,ierr);
	   
#endif

    // Clean up
    ierr = ISDestroy (is);
           CHKERRABORT(libMesh::COMM_WORLD,ierr);
  
    ierr = VecScatterDestroy(scatter);
           CHKERRABORT(libMesh::COMM_WORLD,ierr);
  }

  // localize like normal
  parallel_vec.close();
  parallel_vec.localize (*this, send_list);
  this->close();
}



template <typename T>
void EpetraVector<T>::localize (std::vector<T>& v_local) const
{
  // This function must be run on all processors at once
  parallel_only();

  int ierr=0;
  const int n = this->size();
  const int nl = this->local_size();
  EpetraScalar *values;

  v_local.clear();
  v_local.resize(n, 0.);

  ierr = VecGetArray (_vec, &values);
	 CHKERRABORT(libMesh::COMM_WORLD,ierr);

  unsigned int ioff = first_local_index();

  for (int i=0; i<nl; i++)
    v_local[i+ioff] = static_cast<T>(values[i]);

  ierr = VecRestoreArray (_vec, &values);
	 CHKERRABORT(libMesh::COMM_WORLD,ierr);

  Parallel::sum(v_local);
}



// Full specialization for Real datatypes
#ifdef USE_REAL_NUMBERS

template <>
void EpetraVector<Real>::localize_to_one (std::vector<Real>& v_local,
					 const unsigned int pid) const
{
  int ierr=0;
  const int n  = size();
  const int nl = local_size();
  EpetraScalar *values;

  
  v_local.resize(n);

  
  // only one processor
  if (n == nl)
    {      
      ierr = VecGetArray (_vec, &values);
	     CHKERRABORT(libMesh::COMM_WORLD,ierr);

      for (int i=0; i<n; i++)
	v_local[i] = static_cast<Real>(values[i]);

      ierr = VecRestoreArray (_vec, &values);
	     CHKERRABORT(libMesh::COMM_WORLD,ierr);
    }

  // otherwise multiple processors
  else
    {
      unsigned int ioff = this->first_local_index();
      std::vector<Real> local_values (n, 0.);
      
      {
	ierr = VecGetArray (_vec, &values);
	       CHKERRABORT(libMesh::COMM_WORLD,ierr);
	
	for (int i=0; i<nl; i++)
	  local_values[i+ioff] = static_cast<Real>(values[i]);
	
	ierr = VecRestoreArray (_vec, &values);
	       CHKERRABORT(libMesh::COMM_WORLD,ierr);
      }
      

      MPI_Reduce (&local_values[0], &v_local[0], n, MPI_REAL, MPI_SUM,
		  pid, libMesh::COMM_WORLD);
    }
}

#endif


// Full specialization for Complex datatypes
#ifdef USE_COMPLEX_NUMBERS

template <>
void EpetraVector<Complex>::localize_to_one (std::vector<Complex>& v_local,
					    const unsigned int pid) const
{
  int ierr=0;
  const int n  = size();
  const int nl = local_size();
  EpetraScalar *values;

  
  v_local.resize(n);

  
  for (int i=0; i<n; i++)
    v_local[i] = 0.;
  
  // only one processor
  if (n == nl)
    {      
      ierr = VecGetArray (_vec, &values);
	     CHKERRABORT(libMesh::COMM_WORLD,ierr);

      for (int i=0; i<n; i++)
	v_local[i] = static_cast<Complex>(values[i]);

      ierr = VecRestoreArray (_vec, &values);
	     CHKERRABORT(libMesh::COMM_WORLD,ierr);
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
	       CHKERRABORT(libMesh::COMM_WORLD,ierr);
	
	// provide my local share to the real and imag buffers
	for (int i=0; i<nl; i++)
	  {
	    real_local_values[i+ioff] = static_cast<Complex>(values[i]).real();
	    imag_local_values[i+ioff] = static_cast<Complex>(values[i]).imag();
	  }

	ierr = VecRestoreArray (_vec, &values);
	       CHKERRABORT(libMesh::COMM_WORLD,ierr);
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
		  pid, libMesh::COMM_WORLD);	

      MPI_Reduce (&imag_local_values[0], &imag_v_local[0], n, 
		  MPI_DOUBLE, MPI_SUM,
		  pid, libMesh::COMM_WORLD);	

      // copy real_v_local and imag_v_local to v_local
      for (int i=0; i<n; i++)
	v_local[i] = Complex(real_v_local[i], imag_v_local[i]);
    }  
}

#endif



template <typename T>
void EpetraVector<T>::print_matlab (const std::string name) const
{
  assert (this->initialized());
  assert (this->closed());
  
  int ierr=0; 
  EpetraViewer epetra_viewer;


  ierr = EpetraViewerCreate (libMesh::COMM_WORLD,
			    &epetra_viewer);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);

  /**
   * Create an ASCII file containing the matrix
   * if a filename was provided.  
   */
  if (name != "NULL")
    {
      ierr = EpetraViewerASCIIOpen( libMesh::COMM_WORLD,
				   name.c_str(),
				   &epetra_viewer);
             CHKERRABORT(libMesh::COMM_WORLD,ierr);
      
      ierr = EpetraViewerSetFormat (epetra_viewer,
				   EPETRA_VIEWER_ASCII_MATLAB);
             CHKERRABORT(libMesh::COMM_WORLD,ierr);
  
      ierr = VecView (_vec, epetra_viewer);
             CHKERRABORT(libMesh::COMM_WORLD,ierr);
    }

  /**
   * Otherwise the matrix will be dumped to the screen.
   */
  else
    {
      ierr = EpetraViewerSetFormat (EPETRA_VIEWER_STDOUT_WORLD,
				   EPETRA_VIEWER_ASCII_MATLAB);
             CHKERRABORT(libMesh::COMM_WORLD,ierr);
  
      ierr = VecView (_vec, EPETRA_VIEWER_STDOUT_WORLD);
             CHKERRABORT(libMesh::COMM_WORLD,ierr);
    }


  /**
   * Destroy the viewer.
   */
  ierr = EpetraViewerDestroy (epetra_viewer);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);
}





template <typename T>
void EpetraVector<T>::create_subvector(NumericVector<T>& subvector,
				      const std::vector<unsigned int>& rows) const
{
  // Epetra data structures
  IS parent_is, subvector_is;
  VecScatter scatter;
  int ierr = 0;
  
  // Make sure the passed int subvector is really a EpetraVector
  EpetraVector<T>* epetra_subvector = dynamic_cast<EpetraVector<T>*>(&subvector);
  assert(epetra_subvector != NULL);
  
  // If the epetra_subvector is already initialized, we assume that the
  // user has already allocated the *correct* amount of space for it.
  // If not, we use the appropriate Epetra routines to initialize it.
  if (!epetra_subvector->initialized())
    {
      // Initialize the epetra_subvector to have enough space to hold
      // the entries which will be scattered into it.  Note: such an
      // init() function (where we let Epetra decide the number of local
      // entries) is not currently offered by the EpetraVector
      // class.  Should we differentiate here between sequential and
      // parallel vector creation based on libMesh::n_processors() ?
      ierr = VecCreateMPI(libMesh::COMM_WORLD,
			  EPETRA_DECIDE,          // n_local
			  rows.size(),           // n_global
			  &(epetra_subvector->_vec)); CHKERRABORT(libMesh::COMM_WORLD,ierr);

      ierr = VecSetFromOptions (epetra_subvector->_vec); CHKERRABORT(libMesh::COMM_WORLD,ierr);

      // Mark the subvector as initialized
      epetra_subvector->_is_initialized = true;
    }
  
  // Use iota to fill an array with entries [0,1,2,3,4,...rows.size()]
  std::vector<int> idx(rows.size());
  Utility::iota (idx.begin(), idx.end(), 0);

  // Construct index sets
  ierr = ISCreateGeneral(libMesh::COMM_WORLD,
			 rows.size(),
			 (int*) &rows[0],
			 &parent_is); CHKERRABORT(libMesh::COMM_WORLD,ierr);

  ierr = ISCreateGeneral(libMesh::COMM_WORLD,
			 rows.size(),
			 (int*) &idx[0],
			 &subvector_is); CHKERRABORT(libMesh::COMM_WORLD,ierr);

  // Construct the scatter object
  ierr = VecScatterCreate(this->_vec,
			  parent_is,
			  epetra_subvector->_vec,
			  subvector_is,
			  &scatter); CHKERRABORT(libMesh::COMM_WORLD,ierr);

  // Actually perform the scatter
#if EPETRA_VERSION_LESS_THAN(2,3,3)
  ierr = VecScatterBegin(this->_vec,
			 epetra_subvector->_vec,
			 INSERT_VALUES,
			 SCATTER_FORWARD,
			 scatter); CHKERRABORT(libMesh::COMM_WORLD,ierr);

  ierr = VecScatterEnd(this->_vec,
		       epetra_subvector->_vec,
		       INSERT_VALUES,
		       SCATTER_FORWARD,
		       scatter); CHKERRABORT(libMesh::COMM_WORLD,ierr);
#else
  // API argument order change in Epetra 2.3.3
  ierr = VecScatterBegin(scatter,
			 this->_vec,
			 epetra_subvector->_vec,
			 INSERT_VALUES,
			 SCATTER_FORWARD); CHKERRABORT(libMesh::COMM_WORLD,ierr);

  ierr = VecScatterEnd(scatter,
		       this->_vec,
		       epetra_subvector->_vec,
		       INSERT_VALUES,
		       SCATTER_FORWARD); CHKERRABORT(libMesh::COMM_WORLD,ierr);
#endif
  
  // Clean up 
  ierr = ISDestroy(parent_is);       CHKERRABORT(libMesh::COMM_WORLD,ierr);
  ierr = ISDestroy(subvector_is);    CHKERRABORT(libMesh::COMM_WORLD,ierr);
  ierr = VecScatterDestroy(scatter); CHKERRABORT(libMesh::COMM_WORLD,ierr); 

}




//------------------------------------------------------------------
// Explicit instantiations
template class EpetraVector<Number>;



#endif // #ifdef HAVE_EPETRA

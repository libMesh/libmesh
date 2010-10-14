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

// Local Includes
#include "petsc_vector.h"
#include "petsc_matrix.h"

#ifdef LIBMESH_HAVE_PETSC

#include "dense_subvector.h"
#include "dense_vector.h"
#include "parallel.h"
#include "petsc_macro.h"
#include "utility.h"

namespace libMesh
{



//-----------------------------------------------------------------------
// PetscVector members

// void PetscVector<T>::init (const NumericVector<T>& v, const bool fast)
// {
//   libmesh_error();
  
//   init (v.local_size(), v.size(), fast);

//   vec = libmesh_cast_ref<const PetscVector<T>&>(v).vec;
// }

template <typename T>
T PetscVector<T>::sum () const
{
  this->_restore_array();
  libmesh_assert(this->closed());
  
  int ierr=0;
  PetscScalar value=0.;
  
  ierr = VecSum (_vec, &value);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);
  
  return static_cast<T>(value);
}


template <typename T>
Real PetscVector<T>::l1_norm () const
{
  this->_restore_array();
  libmesh_assert(this->closed());
  
  int ierr=0;
  PetscReal value=0.;
  
  ierr = VecNorm (_vec, NORM_1, &value);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);
  
  return static_cast<Real>(value);
}



template <typename T>
Real PetscVector<T>::l2_norm () const
{
  this->_restore_array();
  libmesh_assert(this->closed());
  
  int ierr=0;
  PetscReal value=0.;
  
  ierr = VecNorm (_vec, NORM_2, &value);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);
  
  return static_cast<Real>(value);
}




template <typename T>
Real PetscVector<T>::linfty_norm () const
{
  this->_restore_array();
  libmesh_assert(this->closed());
  
  int ierr=0;
  PetscReal value=0.;
  
  ierr = VecNorm (_vec, NORM_INFINITY, &value);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);
  
  return static_cast<Real>(value);
}




template <typename T>
NumericVector<T>&
PetscVector<T>::operator += (const NumericVector<T>& v)
{
  this->_restore_array();
  libmesh_assert(this->closed());
  
  this->add(1., v);
  
  return *this;
}



template <typename T>
NumericVector<T>&
PetscVector<T>::operator -= (const NumericVector<T>& v)
{
  this->_restore_array();
  libmesh_assert(this->closed());
  
  this->add(-1., v);
  
  return *this;
}



template <typename T>
void PetscVector<T>::set (const unsigned int i, const T value)
{
  this->_restore_array();
  libmesh_assert(i<size());
  
  int ierr=0;
  int i_val = static_cast<int>(i);
  PetscScalar petsc_value = static_cast<PetscScalar>(value);

  ierr = VecSetValues (_vec, 1, &i_val, &petsc_value, INSERT_VALUES);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);

  this->_is_closed = false;
}



template <typename T>
void PetscVector<T>::add (const unsigned int i, const T value)
{
  this->_restore_array();
  libmesh_assert(i<size());
  
  int ierr=0;
  int i_val = static_cast<int>(i);
  PetscScalar petsc_value = static_cast<PetscScalar>(value);

  ierr = VecSetValues (_vec, 1, &i_val, &petsc_value, ADD_VALUES);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);

  this->_is_closed = false;
}



template <typename T>
void PetscVector<T>::add_vector (const std::vector<T>& v,
				 const std::vector<unsigned int>& dof_indices)
{
  this->_restore_array();
  libmesh_assert (v.size() == dof_indices.size());

  for (unsigned int i=0; i<v.size(); i++)
    this->add (dof_indices[i], v[i]);
}



template <typename T>
void PetscVector<T>::add_vector (const NumericVector<T>& V,
				 const std::vector<unsigned int>& dof_indices)
{
  libmesh_assert (V.size() == dof_indices.size());

  for (unsigned int i=0; i<V.size(); i++)
    this->add (dof_indices[i], V(i));
}



template <typename T>
void PetscVector<T>::add_vector (const NumericVector<T>& V_in,
				 const SparseMatrix<T>& A_in)
{
  this->_restore_array();
  // Make sure the data passed in are really of Petsc types
  const PetscVector<T>* V = libmesh_cast_ptr<const PetscVector<T>*>(&V_in);
  const PetscMatrix<T>* A = libmesh_cast_ptr<const PetscMatrix<T>*>(&A_in);

  int ierr=0;

  A->close();

  // The const_cast<> is not elegant, but it is required since PETSc
  // is not const-correct.  
  ierr = MatMultAdd(const_cast<PetscMatrix<T>*>(A)->mat(), V->_vec, _vec, _vec);
         CHKERRABORT(libMesh::COMM_WORLD,ierr); 
}



template <typename T>
void PetscVector<T>::add_vector (const DenseVector<T>& V,
				 const std::vector<unsigned int>& dof_indices)
{
  libmesh_assert (V.size() == dof_indices.size());

  for (unsigned int i=0; i<V.size(); i++)
    this->add (dof_indices[i], V(i));
}



template <typename T>
void PetscVector<T>::add (const T v_in)
{
  this->_restore_array();

  int ierr=0;
  PetscScalar* values;
  const PetscScalar v = static_cast<PetscScalar>(v_in);  

  if(this->type() != GHOSTED)
    {
      const int n   = static_cast<int>(this->local_size());
      const int fli = static_cast<int>(this->first_local_index());
      
      for (int i=0; i<n; i++)
	{
	  ierr = VecGetArray (_vec, &values);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);
	  
	  int ig = fli + i;      
	  
	  PetscScalar value = (values[i] + v);
	  
	  ierr = VecRestoreArray (_vec, &values);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);
	  
	  ierr = VecSetValues (_vec, 1, &ig, &value, INSERT_VALUES);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr); 
	}
    }
  else
    {
      /* Vectors that include ghost values require a special
	 handling.  */
      Vec loc_vec;
      ierr = VecGhostGetLocalForm (_vec,&loc_vec);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);

      int n=0;
      ierr = VecGetSize(loc_vec, &n);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
      
      for (int i=0; i<n; i++)
	{
	  ierr = VecGetArray (loc_vec, &values);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);
	  
	  PetscScalar value = (values[i] + v);
	  
	  ierr = VecRestoreArray (loc_vec, &values);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);
	  
	  ierr = VecSetValues (loc_vec, 1, &i, &value, INSERT_VALUES);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr); 
	}

      ierr = VecGhostRestoreLocalForm (_vec,&loc_vec);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
    }

  this->_is_closed = false;
}



template <typename T>
void PetscVector<T>::add (const NumericVector<T>& v)
{
  this->add (1., v);
}



template <typename T>
void PetscVector<T>::add (const T a_in, const NumericVector<T>& v_in)
{
  this->_restore_array();

  int ierr = 0;
  PetscScalar a = static_cast<PetscScalar>(a_in);

  // Make sure the NumericVector passed in is really a PetscVector
  const PetscVector<T>* v = libmesh_cast_ptr<const PetscVector<T>*>(&v_in);
  v->_restore_array();

  libmesh_assert(this->size() == v->size());
  
  if(this->type() != GHOSTED)
    {
#if PETSC_VERSION_LESS_THAN(2,3,0)
      // 2.2.x & earlier style
      ierr = VecAXPY(&a, v->_vec, _vec);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
#else
      // 2.3.x & later style
      ierr = VecAXPY(_vec, a, v->_vec);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
#endif
    }
  else
    {
      Vec loc_vec;
      Vec v_loc_vec;
      ierr = VecGhostGetLocalForm (_vec,&loc_vec);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
      ierr = VecGhostGetLocalForm (v->_vec,&v_loc_vec);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);

#if PETSC_VERSION_LESS_THAN(2,3,0)
      // 2.2.x & earlier style
      ierr = VecAXPY(&a, v_loc_vec, loc_vec);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
#else
      // 2.3.x & later style
      ierr = VecAXPY(loc_vec, a, v_loc_vec);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
#endif

      ierr = VecGhostRestoreLocalForm (v->_vec,&v_loc_vec);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
      ierr = VecGhostRestoreLocalForm (_vec,&loc_vec);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
    }
}



template <typename T>
void PetscVector<T>::insert (const std::vector<T>& v,
			     const std::vector<unsigned int>& dof_indices)
{
  libmesh_assert (v.size() == dof_indices.size());

  for (unsigned int i=0; i<v.size(); i++)
    this->set (dof_indices[i], v[i]);
}



template <typename T>
void PetscVector<T>::insert (const NumericVector<T>& V,
			     const std::vector<unsigned int>& dof_indices)
{
  libmesh_assert (V.size() == dof_indices.size());

  for (unsigned int i=0; i<V.size(); i++)
    this->set (dof_indices[i], V(i));
}



template <typename T>
void PetscVector<T>::insert (const DenseVector<T>& V,
			     const std::vector<unsigned int>& dof_indices)
{
  libmesh_assert (V.size() == dof_indices.size());

  for (unsigned int i=0; i<V.size(); i++)
    this->set (dof_indices[i], V(i));
}



template <typename T>
void PetscVector<T>::insert (const DenseSubVector<T>& V,
			     const std::vector<unsigned int>& dof_indices)
{
  libmesh_assert (V.size() == dof_indices.size());

  for (unsigned int i=0; i<V.size(); i++)
    this->set (dof_indices[i], V(i));
}



template <typename T>
void PetscVector<T>::scale (const T factor_in)
{
  this->_restore_array();

  int ierr = 0;
  PetscScalar factor = static_cast<PetscScalar>(factor_in);
  
  if(this->type() != GHOSTED)
    {
#if PETSC_VERSION_LESS_THAN(2,3,0)
      // 2.2.x & earlier style
      ierr = VecScale(&factor, _vec);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
#else
      // 2.3.x & later style	 
      ierr = VecScale(_vec, factor);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
#endif
    }
  else
    {
      Vec loc_vec;
      ierr = VecGhostGetLocalForm (_vec,&loc_vec);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);

#if PETSC_VERSION_LESS_THAN(2,3,0)
      // 2.2.x & earlier style
      ierr = VecScale(&factor, loc_vec);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
#else
      // 2.3.x & later style	 
      ierr = VecScale(loc_vec, factor);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
#endif

      ierr = VecGhostRestoreLocalForm (_vec,&loc_vec);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
    }
}

template <typename T>
void PetscVector<T>::abs()
{
  this->_restore_array();

  int ierr = 0;

  if(this->type() != GHOSTED)
    {
      ierr = VecAbs(_vec);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);  
    }
  else
    {
      Vec loc_vec;
      ierr = VecGhostGetLocalForm (_vec,&loc_vec);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);

      ierr = VecAbs(loc_vec);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);  

      ierr = VecGhostRestoreLocalForm (_vec,&loc_vec);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
    }
}

template <typename T>
T PetscVector<T>::dot (const NumericVector<T>& V) const
{
  this->_restore_array();

  // Error flag
  int ierr = 0;
  
  // Return value
  PetscScalar value=0.;
  
  // Make sure the NumericVector passed in is really a PetscVector
  const PetscVector<T>* v = libmesh_cast_ptr<const PetscVector<T>*>(&V);

  // 2.3.x (at least) style.  Untested for previous versions.
  ierr = VecDot(this->_vec, v->_vec, &value);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);

  return static_cast<T>(value);
}




template <typename T>
NumericVector<T>& 
PetscVector<T>::operator = (const T s_in)
{
  this->_restore_array();
  libmesh_assert(this->closed());

  int ierr = 0;
  PetscScalar s = static_cast<PetscScalar>(s_in);

  if (this->size() != 0)
    {
      if(this->type() != GHOSTED)
	{
#if PETSC_VERSION_LESS_THAN(2,3,0)
	  // 2.2.x & earlier style
	  ierr = VecSet(&s, _vec);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);
#else
	  // 2.3.x & later style	 
	  ierr = VecSet(_vec, s);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);
#endif
	}
      else
	{
	  Vec loc_vec;
	  ierr = VecGhostGetLocalForm (_vec,&loc_vec);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);

#if PETSC_VERSION_LESS_THAN(2,3,0)
	  // 2.2.x & earlier style
	  ierr = VecSet(&s, loc_vec);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);
#else
	  // 2.3.x & later style	 
	  ierr = VecSet(loc_vec, s);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);
#endif

	  ierr = VecGhostRestoreLocalForm (_vec,&loc_vec);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);
	}
    }
  
  return *this;
}



template <typename T>
NumericVector<T>&
PetscVector<T>::operator = (const NumericVector<T>& v_in)
{
  // Make sure the NumericVector passed in is really a PetscVector
  const PetscVector<T>* v = libmesh_cast_ptr<const PetscVector<T>*>(&v_in);

  *this = *v;
  
  return *this;
}



template <typename T>
PetscVector<T>&
PetscVector<T>::operator = (const PetscVector<T>& v)
{
  this->_restore_array();
  v._restore_array();

  libmesh_assert (this->size() == v.size());
  libmesh_assert (this->local_size() == v.local_size());
  libmesh_assert (v.closed());
  this->_is_closed = true;

  int ierr = 0;

  if((this->type()==PARALLEL) && (v.type()==GHOSTED))
    {
      /* Allow assignment of a ghosted to a parallel vector since this
	 causes no difficulty.  See discussion in libmesh-devel of
	 June 24, 2010.  */
      ierr = VecCopy (v._vec, this->_vec);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
    }
  else
    {
      /* In all other cases, we assert that both vectors are of equal
	 type.  */
      libmesh_assert (this->_type == v._type);
      libmesh_assert (this->_global_to_local_map ==
		      v._global_to_local_map);
      
      if (v.size() != 0)
	{
	  if(this->type() != GHOSTED)
	    {
	      ierr = VecCopy (v._vec, this->_vec);
	      CHKERRABORT(libMesh::COMM_WORLD,ierr);
	    }
	  else
	    {
	      Vec loc_vec;
	      Vec v_loc_vec;
	      ierr = VecGhostGetLocalForm (_vec,&loc_vec);
	      CHKERRABORT(libMesh::COMM_WORLD,ierr);
	      ierr = VecGhostGetLocalForm (v._vec,&v_loc_vec);
	      CHKERRABORT(libMesh::COMM_WORLD,ierr);
	      
	      ierr = VecCopy (v_loc_vec, loc_vec);
	      CHKERRABORT(libMesh::COMM_WORLD,ierr);
	      
	      ierr = VecGhostRestoreLocalForm (v._vec,&v_loc_vec);
	      CHKERRABORT(libMesh::COMM_WORLD,ierr);
	      ierr = VecGhostRestoreLocalForm (_vec,&loc_vec);
	      CHKERRABORT(libMesh::COMM_WORLD,ierr);
	    }
	}
    }
  
  return *this;
}



template <typename T>
NumericVector<T>&
PetscVector<T>::operator = (const std::vector<T>& v)
{
  this->_restore_array();

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
 	     CHKERRABORT(libMesh::COMM_WORLD,ierr);

      for (unsigned int i=0; i<nl; i++)
	values[i] =  static_cast<PetscScalar>(v[i+ioff]);
      
      ierr = VecRestoreArray (_vec, &values);
	     CHKERRABORT(libMesh::COMM_WORLD,ierr);
    }

  /**
   * Case 2: The vector is the same size as our local
   * piece.  Insert directly to the local piece.
   */
  else
    {
      libmesh_assert (this->local_size() == v.size());

      ierr = VecGetArray (_vec, &values);
	     CHKERRABORT(libMesh::COMM_WORLD,ierr);

      for (unsigned int i=0; i<nl; i++)
	values[i] = static_cast<PetscScalar>(v[i]);
      
      ierr = VecRestoreArray (_vec, &values);
	     CHKERRABORT(libMesh::COMM_WORLD,ierr);
    }

  // Make sure ghost dofs are up to date
  if (this->type() == GHOSTED)
    this->close();

  return *this;
}



template <typename T>
void PetscVector<T>::localize (NumericVector<T>& v_local_in) const
{
  this->_restore_array();

  // Make sure the NumericVector passed in is really a PetscVector
  PetscVector<T>* v_local = libmesh_cast_ptr<PetscVector<T>*>(&v_local_in);

  libmesh_assert (v_local != NULL);
  libmesh_assert (v_local->size() == this->size());

  int ierr = 0;
  const int n = this->size();

  IS is;
  VecScatter scatter;

  // Create idx, idx[i] = i;
  std::vector<int> idx(n); Utility::iota (idx.begin(), idx.end(), 0);

  // Create the index set & scatter object
  ierr = ISCreateGeneral(libMesh::COMM_WORLD, n, &idx[0], PETSC_USE_POINTER, &is);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);

  ierr = VecScatterCreate(_vec,          is,
			  v_local->_vec, is,
			  &scatter);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);

  // Perform the scatter
#if PETSC_VERSION_LESS_THAN(2,3,3)
	 
  ierr = VecScatterBegin(_vec, v_local->_vec, INSERT_VALUES,
			 SCATTER_FORWARD, scatter);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);
  
  ierr = VecScatterEnd  (_vec, v_local->_vec, INSERT_VALUES,
			 SCATTER_FORWARD, scatter);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);
#else
  // API argument order change in PETSc 2.3.3
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

  // Make sure ghost dofs are up to date
  if (v_local->type() == GHOSTED)
    v_local->close();
}



template <typename T>
void PetscVector<T>::localize (NumericVector<T>& v_local_in,
			       const std::vector<unsigned int>& send_list) const
{
  this->_restore_array();

  // Make sure the NumericVector passed in is really a PetscVector
  PetscVector<T>* v_local = libmesh_cast_ptr<PetscVector<T>*>(&v_local_in);

  libmesh_assert (v_local != NULL);
  libmesh_assert (v_local->size() == this->size());
  libmesh_assert (send_list.size()     <= v_local->size());
  
  int ierr=0;
  const unsigned int n_sl = send_list.size();

  IS is;
  VecScatter scatter;

  std::vector<int> idx(n_sl + this->local_size());
  
  for (unsigned int i=0; i<n_sl; i++)
    idx[i] = static_cast<int>(send_list[i]);
  for (unsigned int i = 0; i != this->local_size(); ++i)
    idx[n_sl+i] = i + this->first_local_index();
  
  // Create the index set & scatter object
  if (idx.empty())
    ierr = ISCreateGeneral(libMesh::COMM_WORLD,
                           n_sl+this->local_size(), PETSC_NULL, PETSC_USE_POINTER, &is);
  else
    ierr = ISCreateGeneral(libMesh::COMM_WORLD,
			   n_sl+this->local_size(), &idx[0], PETSC_USE_POINTER, &is);
           CHKERRABORT(libMesh::COMM_WORLD,ierr);

  ierr = VecScatterCreate(_vec,          is,
			  v_local->_vec, is,
			  &scatter);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);

  
  // Perform the scatter
#if PETSC_VERSION_LESS_THAN(2,3,3)
	 
  ierr = VecScatterBegin(_vec, v_local->_vec, INSERT_VALUES,
			 SCATTER_FORWARD, scatter);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);
  
  ierr = VecScatterEnd  (_vec, v_local->_vec, INSERT_VALUES,
			 SCATTER_FORWARD, scatter);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);

#else
	 
  // API argument order change in PETSc 2.3.3
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

  // Make sure ghost dofs are up to date
  if (v_local->type() == GHOSTED)
    v_local->close();
}


template <typename T>
void PetscVector<T>::localize (const unsigned int first_local_idx,
			       const unsigned int last_local_idx,
			       const std::vector<unsigned int>& send_list)
{
  this->_restore_array();

  // Only good for serial vectors.
  // libmesh_assert (this->size() == this->local_size());
  libmesh_assert (last_local_idx >= first_local_idx);
  libmesh_assert (send_list.size() <= this->size());
  libmesh_assert (last_local_idx < this->size());
  
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

  parallel_vec.init (size, local_size, true, PARALLEL);


  // Copy part of *this into the parallel_vec
  {
    IS is;
    VecScatter scatter;

    // Create idx, idx[i] = i+first_local_idx;
    std::vector<int> idx(local_size);
    Utility::iota (idx.begin(), idx.end(), first_local_idx);

    // Create the index set & scatter object
    ierr = ISCreateGeneral(libMesh::COMM_WORLD, local_size, &idx[0], PETSC_USE_POINTER, &is);
           CHKERRABORT(libMesh::COMM_WORLD,ierr);

    ierr = VecScatterCreate(_vec,              is,
			    parallel_vec._vec, is,
			    &scatter);
           CHKERRABORT(libMesh::COMM_WORLD,ierr);

    // Perform the scatter
#if PETSC_VERSION_LESS_THAN(2,3,3)

    ierr = VecScatterBegin(_vec, parallel_vec._vec, INSERT_VALUES,
			   SCATTER_FORWARD, scatter);
           CHKERRABORT(libMesh::COMM_WORLD,ierr);
  
    ierr = VecScatterEnd  (_vec, parallel_vec._vec, INSERT_VALUES,
			   SCATTER_FORWARD, scatter);
           CHKERRABORT(libMesh::COMM_WORLD,ierr);

#else
	   
      // API argument order change in PETSc 2.3.3
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
void PetscVector<T>::localize (std::vector<T>& v_local) const
{
  this->_restore_array();

  // This function must be run on all processors at once
  parallel_only();

  int ierr=0;
  const int n = this->size();
  const int nl = this->local_size();
  PetscScalar *values;

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
#ifdef LIBMESH_USE_REAL_NUMBERS

template <>
void PetscVector<Real>::localize_to_one (std::vector<Real>& v_local,
					 const unsigned int pid) const
{
  this->_restore_array();

  int ierr=0;
  const int n  = size();
  const int nl = local_size();
  PetscScalar *values;

  
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
#ifdef LIBMESH_USE_COMPLEX_NUMBERS

template <>
void PetscVector<Complex>::localize_to_one (std::vector<Complex>& v_local,
					    const unsigned int pid) const
{
  this->_restore_array();

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
void PetscVector<T>::pointwise_mult (const NumericVector<T>& vec1,
				     const NumericVector<T>& vec2)
{
  this->_restore_array();

  int ierr = 0;

  // Convert arguments to PetscVector*.
  const PetscVector<T>* vec1_petsc = libmesh_cast_ptr<const PetscVector<T>*>(&vec1);
  const PetscVector<T>* vec2_petsc = libmesh_cast_ptr<const PetscVector<T>*>(&vec2);

  // Call PETSc function.

#if PETSC_VERSION_LESS_THAN(2,3,1)

  libMesh::out << "This method has been developed with PETSc 2.3.1.  "
	        << "No one has made it backwards compatible with older "
	        << "versions of PETSc so far; however, it might work "
	        << "without any change with some older version." << std::endl;
  libmesh_error();

#else
  
  if(this->type() != GHOSTED)
    {
      ierr = VecPointwiseMult(this->vec(),
			      const_cast<PetscVector<T>*>(vec1_petsc)->vec(),
			      const_cast<PetscVector<T>*>(vec2_petsc)->vec());
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
    }
  else
    {
	  Vec loc_vec;
	  Vec v1_loc_vec;
	  Vec v2_loc_vec;
	  ierr = VecGhostGetLocalForm (_vec,&loc_vec);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);
	  ierr = VecGhostGetLocalForm (const_cast<PetscVector<T>*>(vec1_petsc)->vec(),&v1_loc_vec);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);
	  ierr = VecGhostGetLocalForm (const_cast<PetscVector<T>*>(vec2_petsc)->vec(),&v2_loc_vec);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);

	  ierr = VecPointwiseMult(loc_vec,v1_loc_vec,v2_loc_vec);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);

	  ierr = VecGhostRestoreLocalForm (const_cast<PetscVector<T>*>(vec1_petsc)->vec(),&v1_loc_vec);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);
	  ierr = VecGhostRestoreLocalForm (const_cast<PetscVector<T>*>(vec2_petsc)->vec(),&v2_loc_vec);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);
	  ierr = VecGhostRestoreLocalForm (_vec,&loc_vec);
	  CHKERRABORT(libMesh::COMM_WORLD,ierr);
    }
      
#endif

}



template <typename T>
void PetscVector<T>::print_matlab (const std::string name) const
{
  this->_restore_array();
  libmesh_assert (this->closed());
  
  int ierr=0; 
  PetscViewer petsc_viewer;


  ierr = PetscViewerCreate (libMesh::COMM_WORLD,
			    &petsc_viewer);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);

  /**
   * Create an ASCII file containing the matrix
   * if a filename was provided.  
   */
  if (name != "NULL")
    {
      ierr = PetscViewerASCIIOpen( libMesh::COMM_WORLD,
				   name.c_str(),
				   &petsc_viewer);
             CHKERRABORT(libMesh::COMM_WORLD,ierr);
      
      ierr = PetscViewerSetFormat (petsc_viewer,
				   PETSC_VIEWER_ASCII_MATLAB);
             CHKERRABORT(libMesh::COMM_WORLD,ierr);
  
      ierr = VecView (_vec, petsc_viewer);
             CHKERRABORT(libMesh::COMM_WORLD,ierr);
    }

  /**
   * Otherwise the matrix will be dumped to the screen.
   */
  else
    {
      ierr = PetscViewerSetFormat (PETSC_VIEWER_STDOUT_WORLD,
				   PETSC_VIEWER_ASCII_MATLAB);
             CHKERRABORT(libMesh::COMM_WORLD,ierr);
  
      ierr = VecView (_vec, PETSC_VIEWER_STDOUT_WORLD);
             CHKERRABORT(libMesh::COMM_WORLD,ierr);
    }


  /**
   * Destroy the viewer.
   */
  ierr = PetscViewerDestroy (petsc_viewer);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);
}





template <typename T>
void PetscVector<T>::create_subvector(NumericVector<T>& subvector,
				      const std::vector<unsigned int>& rows) const
{
  this->_restore_array();

  // PETSc data structures
  IS parent_is, subvector_is;
  VecScatter scatter;
  int ierr = 0;
  
  // Make sure the passed in subvector is really a PetscVector
  PetscVector<T>* petsc_subvector = libmesh_cast_ptr<PetscVector<T>*>(&subvector);
  
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
      ierr = VecCreateMPI(libMesh::COMM_WORLD,
			  PETSC_DECIDE,          // n_local
			  rows.size(),           // n_global
			  &(petsc_subvector->_vec)); CHKERRABORT(libMesh::COMM_WORLD,ierr);

      ierr = VecSetFromOptions (petsc_subvector->_vec); CHKERRABORT(libMesh::COMM_WORLD,ierr);

      // Mark the subvector as initialized
      petsc_subvector->_is_initialized = true;
    }
  else
    {
      petsc_subvector->_restore_array();
    }
  
  // Use iota to fill an array with entries [0,1,2,3,4,...rows.size()]
  std::vector<int> idx(rows.size());
  Utility::iota (idx.begin(), idx.end(), 0);

  // Construct index sets
  ierr = ISCreateGeneral(libMesh::COMM_WORLD,
			 rows.size(),
			 (int*) &rows[0],
			 PETSC_USE_POINTER,
			 &parent_is); CHKERRABORT(libMesh::COMM_WORLD,ierr);

  ierr = ISCreateGeneral(libMesh::COMM_WORLD,
			 rows.size(),
			 (int*) &idx[0],
			 PETSC_USE_POINTER,
			 &subvector_is); CHKERRABORT(libMesh::COMM_WORLD,ierr);

  // Construct the scatter object
  ierr = VecScatterCreate(this->_vec,
			  parent_is,
			  petsc_subvector->_vec,
			  subvector_is,
			  &scatter); CHKERRABORT(libMesh::COMM_WORLD,ierr);

  // Actually perform the scatter
#if PETSC_VERSION_LESS_THAN(2,3,3)
  ierr = VecScatterBegin(this->_vec,
			 petsc_subvector->_vec,
			 INSERT_VALUES,
			 SCATTER_FORWARD,
			 scatter); CHKERRABORT(libMesh::COMM_WORLD,ierr);

  ierr = VecScatterEnd(this->_vec,
		       petsc_subvector->_vec,
		       INSERT_VALUES,
		       SCATTER_FORWARD,
		       scatter); CHKERRABORT(libMesh::COMM_WORLD,ierr);
#else
  // API argument order change in PETSc 2.3.3
  ierr = VecScatterBegin(scatter,
			 this->_vec,
			 petsc_subvector->_vec,
			 INSERT_VALUES,
			 SCATTER_FORWARD); CHKERRABORT(libMesh::COMM_WORLD,ierr);

  ierr = VecScatterEnd(scatter,
		       this->_vec,
		       petsc_subvector->_vec,
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
template class PetscVector<Number>;

} // namespace libMesh



#endif // #ifdef LIBMESH_HAVE_PETSC

// $Id: petsc_vector.C,v 1.10 2003-02-20 04:59:58 benkirk Exp $

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



#ifdef HAVE_PETSC



// void PetscVector<Tp>::init (const NumericVector<Tp>& v, const bool fast)
// {
//   error();
  
//   init (v.local_size(), v.size(), fast);

//   vec = reinterpret_cast<const PetscVector<Tp>&>(v).vec;
// }



template <typename Tp>
Real PetscVector<Tp>::l1_norm () const
{
  assert(closed());
  
  int ierr=0;
  double value=0.;
  
  ierr = VecNorm (vec, NORM_1, &value);
  
  return static_cast<Real>(value);
}



template <typename Tp>
Real PetscVector<Tp>::l2_norm () const
{
  assert(closed());
  
  int ierr=0;
  double value=0.;
  
  ierr = VecNorm (vec, NORM_2, &value);
  
  return static_cast<Real>(value);
}




template <typename Tp>
Real PetscVector<Tp>::linfty_norm () const
{
  assert(closed());
  
  int ierr=0;
  double value=0.;
  
  ierr = VecNorm (vec, NORM_INFINITY, &value);
  
  return static_cast<Real>(value);
}




template <typename Tp>
NumericVector<Tp>&
PetscVector<Tp>::operator += (const NumericVector<Tp>& v)
{
  assert(closed());
  
  add(1., v);
  
  return *this;
}



template <typename Tp>
NumericVector<Tp>&
PetscVector<Tp>::operator -= (const NumericVector<Tp>& v)
{
  assert(closed());
  
  add(-1., v);
  
  return *this;
}



template <typename Tp>
void PetscVector<Tp>::set (const unsigned int i, const Tp value)
{
  assert(i<size());
  
  int ierr=0;
  int i_val = static_cast<int>(i);
  PetscScalar petsc_value = static_cast<PetscScalar>(value);

  ierr = VecSetValues (vec, 1, &i_val, &petsc_value, INSERT_VALUES); CHKERRQ(ierr);

  return;
}



template <typename Tp>
void PetscVector<Tp>::add (const unsigned int i, const Tp value)
{
  assert(i<size());
  
  int ierr=0;
  int i_val = static_cast<int>(i);
  PetscScalar petsc_value = static_cast<PetscScalar>(value);

  ierr = VecSetValues (vec, 1, &i_val, &petsc_value, ADD_VALUES); CHKERRQ(ierr);

  return;
}



template <typename Tp>
void PetscVector<Tp>::add_vector (const std::vector<Tp>& v,
				  const std::vector<unsigned int>& dof_indices)
{
  assert (!v.empty());
  assert (v.size() == dof_indices.size());
  
  for (unsigned int i=0; i<v.size(); i++)
    add (dof_indices[i], v[i]);
}



template <typename Tp>
void PetscVector<Tp>::add_vector (const NumericVector<Tp>& V,
				  const std::vector<unsigned int>& dof_indices)
{
  assert (V.size() == dof_indices.size());

  for (unsigned int i=0; i<V.size(); i++)
    add (dof_indices[i], V(i));
}



template <typename Tp>
void PetscVector<Tp>::add (const Tp v)
{
  int ierr=0;
  PetscScalar* values;

  for (int i=0; i<static_cast<int>(local_size()); i++)
    {
      int ig = first_local_index()+i;
      ierr = VecGetArray (vec, &values); CHKERRQ(ierr);
      PetscScalar value = (values[ig] + static_cast<PetscScalar>(v));
      ierr = VecRestoreArray (vec, &values); CHKERRQ(ierr);
      ierr = VecSetValues (vec, 1, &ig, &value, INSERT_VALUES); CHKERRQ(ierr); 
    }
}



template <typename Tp>
void PetscVector<Tp>::add (const NumericVector<Tp>& v)
{
  add (1., v);
}



template <typename Tp>
void PetscVector<Tp>::add (const Tp a, const NumericVector<Tp>& v_in)
{
  int ierr=0;
  PetscScalar petsc_a=static_cast<PetscScalar>(a);

  const PetscVector<Tp>& v = reinterpret_cast<const PetscVector<Tp>&>(v_in);
  
  assert(size() == v.size());
  
  ierr = VecAXPY(&petsc_a, v.vec, vec); CHKERRQ(ierr);
}



template <typename Tp>
void PetscVector<Tp>::scale (const Tp factor)
{
  int ierr=0;
  PetscScalar petsc_factor = static_cast<PetscScalar>(factor);
  
  ierr = VecScale(&petsc_factor, vec); CHKERRQ(ierr);
}



template <typename Tp>
NumericVector<Tp>& 
PetscVector<Tp>::operator = (const Tp s)
{
  int ierr=0;
  PetscScalar petsc_s=static_cast<PetscScalar>(s);

  if (size() != 0)
    {
      ierr = VecSet(&petsc_s, vec); CHKERRQ(ierr);
    }
  
  return *this;
}



template <typename Tp>
NumericVector<Tp>&
PetscVector<Tp>::operator = (const NumericVector<Tp>& v_in)
{
  int ierr=0;
  
  const PetscVector<Tp>& v = reinterpret_cast<const PetscVector<Tp>&>(v_in);

  assert (size() == v.size());

  if (size() != 0)
    {
      ierr = VecCopy (v.vec, vec); CHKERRQ(ierr);
    }
  
  return *this;
}



template <typename Tp>
PetscVector<Tp>&
PetscVector<Tp>::operator = (const PetscVector<Tp>& v)
{
  int ierr=0;

  assert (size() == v.size());

  if (size() != 0)
    {
      ierr = VecCopy (v.vec, vec); CHKERRQ(ierr);
    }
  
  return *this;
}



template <typename Tp>
NumericVector<Tp>&
PetscVector<Tp>::operator = (const std::vector<Tp>& v)
{
  const unsigned int nl   = local_size();
  const unsigned int ioff = first_local_index();
  int ierr=0;
  PetscScalar* values;
      
  /**
   * Case 1:  The vector is the same size of
   * The global vector.  Only add the local components.
   */
  if (size() == v.size())
    {
      ierr = VecGetArray (vec, &values); CHKERRQ(ierr);

      for (unsigned int i=0; i<nl; i++)
	values[i] =  static_cast<PetscScalar>(v[i+ioff]);
      
      ierr = VecRestoreArray (vec, &values); CHKERRQ(ierr);
    }

  /**
   * Case 2: The vector is the same size as our local
   * piece.  Insert directly to the local piece.
   */
  else
    {
      assert (local_size() == v.size());

      ierr = VecGetArray (vec, &values); CHKERRQ(ierr);

      for (unsigned int i=0; i<nl; i++)
	values[i] = static_cast<PetscScalar>(v[i]);
      
      ierr = VecRestoreArray (vec, &values); CHKERRQ(ierr);
    }

  return *this;
}



template <typename Tp>
void PetscVector<Tp>::localize (NumericVector<Tp>& v_local_in) const
{
  const PetscVector<Tp>& v_local = reinterpret_cast<const PetscVector<Tp>&>(v_local_in);

  assert (v_local.local_size() == size());

  int ierr=0;
  const int n  = size();

  IS is;
  VecScatter scatter;

  std::vector<int> idx(n);

  for (int i=0; i<n; i++)
    idx[i] = i;

  ierr = ISCreateGeneral(PETSC_COMM_WORLD, n, &idx[0], &is);  CHKERRQ(ierr);

  ierr = VecScatterCreate(vec,         is,
			  v_local.vec, is,
			  &scatter);                      CHKERRQ(ierr);

  ierr = VecScatterBegin(vec, v_local.vec, INSERT_VALUES,
			 SCATTER_FORWARD, scatter);       CHKERRQ(ierr);
  
  ierr = VecScatterEnd  (vec, v_local.vec, INSERT_VALUES,
			 SCATTER_FORWARD, scatter);       CHKERRQ(ierr);

  ierr = ISDestroy (is);              CHKERRQ(ierr);
  ierr = VecScatterDestroy(scatter);  CHKERRQ(ierr);
}



template <typename Tp>
void PetscVector<Tp>::localize (NumericVector<Tp>& v_local_in,
				const std::vector<unsigned int>& send_list) const
{
  const PetscVector<Tp>& v_local = reinterpret_cast<const PetscVector<Tp>&>(v_local_in);

  assert (v_local.local_size() == size());
  assert (send_list.size() <= v_local.size());
  
  int ierr=0;
  const int n_sl = send_list.size();
  //  const int nl = local_size();
  //  PetscScalar *values;

  IS is;
  VecScatter scatter;

  std::vector<int> idx(n_sl);

  for (int i=0; i<n_sl; i++)
    idx[i] = static_cast<int>(send_list[i]);

  ierr = ISCreateGeneral(PETSC_COMM_WORLD, n_sl, &idx[0], &is); CHKERRQ(ierr);

  ierr = VecScatterCreate(vec,         is,
			  v_local.vec, is,
			  &scatter);                        CHKERRQ(ierr);

  ierr = VecScatterBegin(vec, v_local.vec, INSERT_VALUES,
			 SCATTER_FORWARD, scatter);         CHKERRQ(ierr);
  
  ierr = VecScatterEnd  (vec, v_local.vec, INSERT_VALUES,
			 SCATTER_FORWARD, scatter);         CHKERRQ(ierr);

  ierr = ISDestroy (is);              CHKERRQ(ierr);
  ierr = VecScatterDestroy(scatter);  CHKERRQ(ierr);
}



// Full specialization for Real datatypes
#ifdef USE_REAL_NUMBERS

template <>
void PetscVector<Real>::localize (std::vector<Real>& v_local) const
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
      ierr = VecGetArray (vec, &values); CHKERRQ(ierr);

      for (int i=0; i<n; i++)
	v_local[i] = static_cast<Real>(values[i]);

      ierr = VecRestoreArray (vec, &values); CHKERRQ(ierr);
    }

  // otherwise multiple processors
  else
    {
      unsigned int ioff = first_local_index();
      std::vector<Real> local_values(n, 0.);

      {
	ierr = VecGetArray (vec, &values); CHKERRQ(ierr);
	
	for (int i=0; i<nl; i++)
	  local_values[i+ioff] = static_cast<Real>(values[i]);
	
	ierr = VecRestoreArray (vec, &values); CHKERRQ(ierr);
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
  //TODO:[DD/BSK] Will this work in parallel?
  
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
      ierr = VecGetArray (vec, &values); CHKERRQ(ierr);

      for (int i=0; i<n; i++)
	v_local[i] = static_cast<Complex>(values[i]);

      ierr = VecRestoreArray (vec, &values); CHKERRQ(ierr);
    }

  // otherwise multiple processors
  else
    {
      unsigned int ioff = first_local_index();
      std::vector<Complex> local_values(n, 0.);

      {
	ierr = VecGetArray (vec, &values); CHKERRQ(ierr);
	
	for (int i=0; i<nl; i++)
	  local_values[i+ioff] = static_cast<Complex>(values[i]);
	
	ierr = VecRestoreArray (vec, &values); CHKERRQ(ierr);
      }
      
      MPI_Allreduce (&local_values[0], &v_local[0], n, MPIU_SCALAR, MPI_SUM,
		     PETSC_COMM_WORLD);	
    }  
}

#endif



/*
  ifdef USE_COMPLEX_NUMBERS
  {
  //TODO:[DD] localize may be done in a better way...
  // There is no MPI_COMPLEX in C. For now, use PETSc 
	    
  // have an appropriately sized PetscVector<Tp> handy
  PetscVector<Tp> pv(n);
	  
  // localize ourselves to this vector
  localize(&pv);

  // copy data to v_local
  for (int i=0; i<n; i++)
  v_local[i] =  static_cast<Tp>(pv(i));

  // note that the destructor of pv calls clear()
  }
  else
  endif
*/





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

  
  for (int i=0; i<n; i++)
    v_local[i] = 0.;
  
  // only one processor
  if (n == nl)
    {      
      ierr = VecGetArray (vec, &values); CHKERRQ(ierr);

      for (int i=0; i<n; i++)
	v_local[i] = static_cast<Real>(values[i]);

      ierr = VecRestoreArray (vec, &values); CHKERRQ(ierr);
    }

  // otherwise multiple processors
  else
    {
      unsigned int ioff = first_local_index();
      std::vector<Real> local_values(n, 0.);

      {
	ierr = VecGetArray (vec, &values); CHKERRQ(ierr);
	
	for (int i=0; i<nl; i++)
	  local_values[i+ioff] = static_cast<Real>(values[i]);
	
	ierr = VecRestoreArray (vec, &values); CHKERRQ(ierr);
      }
      
      if (sizeof(Real) == sizeof(double))     
	MPI_Reduce (&local_values[0], &v_local[0], n, MPI_DOUBLE, pid, MPI_SUM,
		    PETSC_COMM_WORLD);

      else if (sizeof(Real) == sizeof(float))     
	MPI_Reduce (&local_values[0], &v_local[0], n, MPI_FLOAT, pid, MPI_SUM,
		    PETSC_COMM_WORLD);

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
  //TODO:[DD/BSK] Will this work in parallel?
  
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
      ierr = VecGetArray (vec, &values); CHKERRQ(ierr);

      for (int i=0; i<n; i++)
	v_local[i] = static_cast<Complex>(values[i]);

      ierr = VecRestoreArray (vec, &values); CHKERRQ(ierr);
    }

  // otherwise multiple processors
  else
    {
      unsigned int ioff = first_local_index();
      std::vector<Complex> local_values(n, 0.);

      {
	ierr = VecGetArray (vec, &values); CHKERRQ(ierr);
	
	for (int i=0; i<nl; i++)
	  local_values[i+ioff] = static_cast<Complex>(values[i]);
	
	ierr = VecRestoreArray (vec, &values); CHKERRQ(ierr);
      }
      
      MPI_Reduce (&local_values[0], &v_local[0], n, MPIU_SCALAR, pid, MPI_SUM,
		  PETSC_COMM_WORLD);	
    }  
}

#endif


/*
  ifdef USE_COMPLEX_NUMBERS
  {
  //TODO:[DD] localize_to_one may be done in a better way...

  int my_proc_id = 0;
  MPI_Comm_rank (PETSC_COMM_WORLD, &my_proc_id);

  if (my_proc_id == proc_id)
  {
  // have an appropriately sized PetscVector<Tp> handy
  PetscVector<Tp> pv(n);
	  
  // localize ourselves to this vector
  localize(&pv);

  // copy data to v_local
  for (int i=0; i<n; i++)
  v_local[i] =  static_cast<Tp>(pv(i));

  }
  }
  else
  error();
  endif
*/



//------------------------------------------------------------------
// Explicit instantiations
template class PetscVector<Number>;



#endif // #ifdef HAVE_PETSC

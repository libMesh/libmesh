//    $Id: petsc_vector.C,v 1.6 2003-02-03 03:51:50 ddreyer Exp $

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



#include "mesh_common.h"

#ifdef HAVE_PETSC


// C++ includes



// Local Includes
#include "petsc_vector.h"




/*
PetscVector::PetscVector (const PetscVector& v) :
  is_closed(false),
  initialized(false)
{
  error();
  int ierr=0;
  
  ierr = VecDuplicate(v.vec, &vec); CHKERRQ(ierr);

  ierr = VecCopy(v.vec, vec); CHKERRQ(ierr);
};
*/



void PetscVector::init (const PetscVector& v, const bool fast)
{
  init (v.local_size(), v.size(), fast);

  vec = v.vec;
};



Real PetscVector::l1_norm () const
{
  assert(is_closed);
  
  int ierr=0;
  double value=0.;
  
  ierr = VecNorm (vec, NORM_1, &value);
  
  return static_cast<Real>(value);
};



Real PetscVector::l2_norm () const
{
  assert(is_closed);
  
  int ierr=0;
  double value=0.;
  
  ierr = VecNorm (vec, NORM_2, &value);
  
  return static_cast<Real>(value);
};




Real PetscVector::linfty_norm () const
{
  assert(is_closed);
  
  int ierr=0;
  double value=0.;
  
  ierr = VecNorm (vec, NORM_INFINITY, &value);
  
  return static_cast<Real>(value);
};




PetscVector& PetscVector::operator += (const PetscVector& v)
{
  assert(is_closed);
  
  add(1., v);
  
  return *this;
};




PetscVector& PetscVector::operator -= (const PetscVector& v)
{
  assert(is_closed);
  
  add(-1., v);
  
  return *this;
};




void PetscVector::set (const unsigned int i, const Complex value)
{
  assert(i<size());
  
  int ierr=0;
  int i_val = static_cast<int>(i);
  PetscScalar petsc_value = static_cast<PetscScalar>(value);

  ierr = VecSetValues (vec, 1, &i_val, &petsc_value, INSERT_VALUES); CHKERRQ(ierr);

  return;
};




void PetscVector::add (const unsigned int i, const Complex value)
{
  assert(i<size());
  
  int ierr=0;
  int i_val = static_cast<int>(i);
  PetscScalar petsc_value = static_cast<PetscScalar>(value);

  ierr = VecSetValues (vec, 1, &i_val, &petsc_value, ADD_VALUES); CHKERRQ(ierr);

  return;
};




void PetscVector::add_vector (const std::vector<Complex>& v,
			      const std::vector<unsigned int>& dof_indices)
{
  assert (!v.empty());
  assert (v.size() == dof_indices.size());
  
  for (unsigned int i=0; i<v.size(); i++)
    add (dof_indices[i], v[i]);
};



void PetscVector::add_petsc_vector (const PetscVector& V,
				    const std::vector<unsigned int>& dof_indices)
{
  assert (V.size() == dof_indices.size());

  for (unsigned int i=0; i<V.size(); i++)
    add (dof_indices[i], V(i));
};


void PetscVector::add (const Complex v)
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
    };
};




void PetscVector::add (const PetscVector& v)
{
  add (1., v);
};




void PetscVector::add (const Complex a, const PetscVector& v)
{
  int ierr=0;
  PetscScalar petsc_a=static_cast<PetscScalar>(a);
  
  assert(size() == v.size());
  
  ierr = VecAXPY(&petsc_a, v.vec, vec); CHKERRQ(ierr);
};



void PetscVector::scale (const Complex factor)
{
  int ierr=0;
  PetscScalar petsc_factor = static_cast<PetscScalar>(factor);
  
  ierr = VecScale(&petsc_factor, vec); CHKERRQ(ierr);
};



PetscVector& 
PetscVector::operator = (const Complex s)
{
  int ierr=0;
  PetscScalar petsc_s=static_cast<PetscScalar>(s);

  if (size() != 0)
    {
      ierr = VecSet(&petsc_s, vec); CHKERRQ(ierr);
    }
  
  return *this;
};



PetscVector&
PetscVector::operator = (const PetscVector& v)
{
  int ierr=0;

  assert (size() == v.size());

  if (size() != 0)
    {
      ierr = VecCopy (v.vec, vec); CHKERRQ(ierr);
    }
  
  return *this;
};



PetscVector&
PetscVector::operator = (const std::vector<Complex>& v)
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
};



void PetscVector::localize (PetscVector& v_local) const

{
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
};



void PetscVector::localize (PetscVector& v_local,
			    const std::vector<unsigned int>& send_list) const
{
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
};



void PetscVector::localize (std::vector<Complex>& v_local) const

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
      
      if (sizeof(Real) == sizeof(double))     
      {
	if (sizeof(Real) == sizeof(Complex))
	  MPI_Allreduce (&local_values[0], &v_local[0], n, MPI_DOUBLE, MPI_SUM,
		         PETSC_COMM_WORLD);
	else
	  MPI_Allreduce (&local_values[0], &v_local[0], n, MPIU_SCALAR, MPI_SUM,
		         PETSC_COMM_WORLD);	
      }

      else if (sizeof(Real) == sizeof(float))     

      {
	if (sizeof(Real) == sizeof(Complex))
  	  MPI_Allreduce (&local_values[0], &v_local[0], n, MPI_FLOAT, MPI_SUM,
		         PETSC_COMM_WORLD);
	else
	  MPI_Allreduce (&local_values[0], &v_local[0], n, MPIU_SCALAR, MPI_SUM,
		         PETSC_COMM_WORLD);
      }

      else
	error();
    }  
};



/*
ifdef USE_COMPLEX_NUMBERS
        {
          //TODO:[DD] localize may be done in a better way...
          // There is no MPI_COMPLEX in C. For now, use PETSc 
	    
	  // have an appropriately sized PetscVector handy
	  PetscVector pv(n);
	  
	  // localize ourselves to this vector
	  localize(&pv);

	  // copy data to v_local
	  for (int i=0; i<n; i++)
	      v_local[i] =  static_cast<Complex>(pv(i));

	  // note that the destructor of pv calls clear()
	};
else
endif
*/


void PetscVector::localize_to_one (std::vector<Complex>& v_local,
				   const unsigned int pid) const

{
  int ierr=0, proc_id=static_cast<int>(pid);
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


      if (sizeof(Real) == sizeof(double))     
      {
	if (sizeof(Real) == sizeof(Complex))
	  MPI_Reduce (&local_values[0], &v_local[0], n, MPI_DOUBLE, MPI_SUM,
		      proc_id, PETSC_COMM_WORLD);
	else
	  MPI_Reduce (&local_values[0], &v_local[0], n, MPIU_SCALAR, MPI_SUM,
		      proc_id, PETSC_COMM_WORLD);	
      }

      else if (sizeof(Real) == sizeof(float))     

      {
	if (sizeof(Real) == sizeof(Complex))
  	  MPI_Reduce (&local_values[0], &v_local[0], n, MPI_FLOAT, MPI_SUM,
		      proc_id, PETSC_COMM_WORLD);
	else
	{
	  std::cout << "Error: complex<float> not supported by PETSc." << std::endl;
	  error();
	};
      }

      else
	error();

    }

};


/*
ifdef USE_COMPLEX_NUMBERS
        {
          //TODO:[DD] localize_to_one may be done in a better way...

	  int my_proc_id = 0;
	  MPI_Comm_rank (PETSC_COMM_WORLD, &my_proc_id);

	  if (my_proc_id == proc_id)
	    {
	      // have an appropriately sized PetscVector handy
	      PetscVector pv(n);
	  
	      // localize ourselves to this vector
	      localize(&pv);

	      // copy data to v_local
	      for (int i=0; i<n; i++)
		  v_local[i] =  static_cast<Complex>(pv(i));

	    };
	};
else
        error();
endif
*/


#endif

// $Id: system_base.C,v 1.11 2003-03-21 15:29:28 ddreyer Exp $

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
#include <algorithm>

// Local includes
#include "mesh.h"
#include "libmesh.h"
#include "system_base.h"

// typedef
typedef std::map<std::string, SparseMatrix<Number>* >::iterator        other_matrices_iterator;
typedef std::map<std::string, SparseMatrix<Number>* >::const_iterator  other_matrices_const_iterator;
typedef std::map<std::string, NumericVector<Number>* >::iterator       other_vectors_iterator;
typedef std::map<std::string, NumericVector<Number>* >::const_iterator other_vectors_const_iterator;


// ------------------------------------------------------------
// SystemBase implementation

SystemBase::SystemBase (Mesh& mesh,
			const std::string&  name,
			const unsigned int  number,
			const SolverPackage solver_package) :
  
  solution                (NumericVector<Number>::build(solver_package)),
  rhs                     (NumericVector<Number>::build(solver_package)),
  matrix                  (SparseMatrix<Number>::build (solver_package)),
  linear_solver_interface (LinearSolverInterface<Number>::build(solver_package)),
  _sys_name               (name),
  _sys_number             (number),
  _can_add_matrices       (true),
  _can_add_vectors        (true),
  _solver_package         (solver_package),
  _dof_map                (number),
  _mesh                   (mesh)
{
}



SystemBase::~SystemBase ()
{
  SystemBase::clear ();

  assert (!libMesh::closed());
}


void SystemBase::clear ()
{
  _var_names.clear ();

  _var_type.clear ();
  
  _var_num.clear ();
  
  _dof_map.clear ();
  
  solution->clear ();

  rhs->clear ();

  matrix->clear ();

  // clear other matrices, if existent
  {
    other_matrices_iterator pos = _other_matrices.begin();
    while (pos != _other_matrices.end())
      {
        // clear matrix, delete pointer, erase entry in map & get next entry
        pos->second->clear ();
	delete pos->second;
	_other_matrices.erase(pos);
	pos = _other_matrices.begin();
      }

    _can_add_matrices=true;
  }

  // clear other vectors, if existent
  {
    other_vectors_iterator pos = _other_vectors.begin();
    while (pos != _other_vectors.end())
      {
        pos->second->clear();
	delete pos->second;
	_other_vectors.erase(pos);
	pos = _other_vectors.begin();
      }

    _can_add_vectors=true;
  }
}




void SystemBase::init ()
{
  assert (_mesh.is_prepared());
  
  _dof_map.distribute_dofs (_mesh);

#ifdef ENABLE_AMR

  _dof_map.create_dof_constraints(_mesh);

#endif

  solution->init (n_dofs(), n_local_dofs());

  rhs->init      (n_dofs(), n_local_dofs());

  // from now on, no chance to add additional vectors
  _can_add_vectors=false;

  // initialize & zero other vectors, if necessary
  for(other_vectors_iterator pos = _other_vectors.begin();
      pos != _other_vectors.end(); ++pos)
    {
      pos->second->init (n_dofs(), n_local_dofs());
      pos->second->zero ();
    }
}




void SystemBase::assemble ()
{
  // no chance to add other matrices
  _can_add_matrices=false;
  
  // initialize the matrix if not 
  // already initialized
  if (!matrix->initialized())
    {
      // Tell the matrix and the dof map
      // about each other.  This is necessary
      // because the _dof_map needs to share
      // information with the matrix during the
      // sparsity computation phase
      _dof_map.attach_matrix (*matrix.get());
//      matrix->attach_dof_map (_dof_map);

      // Also tell the additional matrices about
      // the dof map, and vice versa
      for(other_matrices_iterator pos = _other_matrices.begin(); 
	  pos != _other_matrices.end(); ++pos)
	if (!pos->second->initialized())
	  {
	    _dof_map.attach_other_matrix (*pos);
//	    pos->second->attach_dof_map (_dof_map);
	  }
	else
	  {
	    std::cerr << "ERROR: Something wrong: Mayor matrix is not initialized, "
		      << std::endl  << " but additional matrix "
		      << pos->first << " is!" << std::endl;
	    error();
	  }  



      // Compute the sparsity pattern for the current
      // mesh and DOF distribution.  This also updates
      // additional matrices, \p DofMap knows them
      _dof_map.compute_sparsity (_mesh);


      // Initialize the matrix conformal to this
      // sparsity pattern, though normally, _dof_map.compute_sparsity()
      // already did this
      matrix->init ();


      // Initialize additional matrices conformal
      // to this sparsity pattern, though normally, _dof_map.compute_sparsity()
      // already did this
      for(other_matrices_iterator pos = _other_matrices.begin(); 
	  pos != _other_matrices.end(); ++pos)
	    pos->second->init ();
	    
    }
  else
    {
      // Better check whether the other matrices are also initialized
      for(other_matrices_const_iterator pos = _other_matrices.begin(); 
	  pos != _other_matrices.end(); ++pos)
	if (!pos->second->initialized())
	  {
	    // theoretically, with the proper access tests in get_matrix(),
	    // this cannot happen.  But, who knows?
	    std::cerr << "ERROR: Mayor matrix is initialized, but additional matrix "
		      << pos->first
		      << " is not!"
		      << std::endl;
	    error();
	  }
    }


  // Clear the matrix and right-hand side.
  matrix->zero ();
  rhs->zero    ();

  // Clear the additional matrices
  for(other_matrices_iterator pos = _other_matrices.begin(); 
      pos != _other_matrices.end(); ++pos)
        pos->second->zero ();

  // Now everything is set up and ready for matrix assembly
}



bool SystemBase::compare (const SystemBase& other_system, 
			  const Real threshold,
			  const bool verbose) const
{
  // we do not care for matrices, but for vectors
  assert (!_can_add_vectors);
  assert (!other_system._can_add_vectors);

  if (verbose)
    {
      std::cout << "  Systems \"" << _sys_name << "\"" << std::endl;
      std::cout << "   comparing matrices not supported." << std::endl;
      std::cout << "   comparing names...";
    }

  // compare the name: 0 means identical
  const int name_result = _sys_name.compare(other_system.name());
  if (verbose)
    {
      if (name_result == 0)
	std::cout << " identical." << std::endl;
      else
	std::cout << "  names not identical." << std::endl;
      std::cout << "   comparing solution vector...";
    }


  // compare the solution: -1 means identical
  const int solu_result = solution->compare (*other_system.solution.get(),
					     threshold);

  if (verbose)
    {
      if (solu_result == -1)
	std::cout << " identical up to threshold." << std::endl;
      else
	std::cout << "  first difference occured at index = " 
		  << solu_result << "." << std::endl;
    }


  // safety check, whether we handle at least the same number
  // of vectors
  std::vector<int> ov_result;

  if (this->n_additional_vectors() != other_system.n_additional_vectors())
    {
      if (verbose)
        {
	  std::cout << "   Fatal difference. This system handles " 
		    << this->n_additional_vectors() << " add'l vectors," << std::endl
		    << "   while the other system handles "
		    << other_system.n_additional_vectors() 
		    << " add'l vectors." << std::endl
		    << "   Aborting comparison." << std::endl;
	}
      return false;
    }
  else if (this->n_additional_vectors() == 0)
    {
      // there are no additional vectors...
      ov_result.clear ();
    }
  else
    {
      // compare other vectors
      for(other_vectors_const_iterator pos = _other_vectors.begin();
	  pos != _other_vectors.end(); ++pos)
        {
	  if (verbose)
	      std::cout << "   comparing vector \""
			<< pos->first << "\" ...";

	  // assume they have the same name
	  const NumericVector<Number>& other_system_vector = 
	      other_system.get_vector(pos->first);

	  ov_result.push_back(pos->second->compare (other_system_vector,
						    threshold));

	  if (verbose)
	    {
	      if (solu_result == -1)
		std::cout << " identical up to threshold." << std::endl;
	      else
		std::cout << " first difference occured at" << std::endl
			  << "  index = " << solu_result << "." << std::endl;
	    }

	}

    } // finished comparing additional vectors


  bool overall_result;
       
  // sum up the results
  if ((name_result==0) && (solu_result==-1))
    {
      if (ov_result.size()==0)
	overall_result = true;
      else
        {
	  bool ov_identical;
	  unsigned int n    = 0;
	  do
	    {
	      ov_identical = (ov_result[n]==-1);
	      n++;
	    }
	  while (ov_identical && n<ov_result.size());
	  overall_result = ov_identical;
	}
    }
  else
    overall_result = false;

  if (verbose)
    {
      std::cout << "   finished comparisons, ";
      if (overall_result)
	std::cout << "found no differences." << std::endl << std::endl;
      else 
	std::cout << "found differences." << std::endl << std::endl;
    }
	  
  return overall_result;
}



void SystemBase::update_global_solution (std::vector<Number>& global_soln) const
{
  global_soln.resize (solution->size());

  solution->localize  (global_soln);
}



void SystemBase::update_global_solution (std::vector<Number>& global_soln,
					 const unsigned int dest_proc) const
{
  global_soln.resize       (solution->size());

  solution->localize_to_one (global_soln, dest_proc);
}




void SystemBase::add_matrix (const std::string& mat_name)
{
  // only add matrices _before_ assembly...
  if (!_can_add_matrices)
    {
      std::cerr << "ERROR: Too late.  Cannot add matrices to the system"
		<< std::endl
		<< " any more.  You should have done this earlier."
		<< std::endl;
      error();
    }

  // Make sure the matrix isn't there already
  if (_other_matrices.count(mat_name))
    {
      std::cerr << "ERROR: matrix "
		<< mat_name
		<< " has already been added to this system!"
		<< std::endl;
      error();
    }
  
  // build the matrix, add it to the map
  SparseMatrix<Number>* buf(SparseMatrix<Number>::build(_solver_package).release());
  _other_matrices[mat_name] = buf;

  // assemble() adds the matrix to the _dof_map
}




const SparseMatrix<Number> &  SystemBase::get_matrix(const std::string& mat_name) const
{
  // only enable access _after_ the matrices are properly initialized
  if (_can_add_matrices)
    {
      std::cerr << "ERROR: Too early.  Access to additional matrices granted only when "
		<< std::endl
		<< "these are already properly initialized."
		<< std::endl;
      error();
    }

   // Make sure the matrix exists
  other_matrices_const_iterator pos = _other_matrices.find(mat_name);  
  if (pos == _other_matrices.end())
    {
      std::cerr << "ERROR: matrix "
		<< mat_name
		<< " does not exist in this system!"
		<< std::endl;      
      error();
    }
  
  return *(pos->second);
}




SparseMatrix<Number> &  SystemBase::get_matrix(const std::string& mat_name)
{
  // only enable access _after_ the matrices are properly initialized
  if (_can_add_matrices)
    {
      std::cerr << "ERROR: Too early.  Access to additional matrices granted only when "
		<< std::endl
		<< "these are already properly initialized."
		<< std::endl;
      error();
    }

   // Make sure the matrix exists
  other_matrices_const_iterator pos = _other_matrices.find(mat_name);  
  if (pos == _other_matrices.end())
    {
      std::cerr << "ERROR: matrix "
		<< mat_name
		<< " does not exist in this system!"
		<< std::endl;      
      error();
    }
  
  return *(pos->second);
}




void SystemBase::add_vector (const std::string& vec_name)
{
  // only add vectors before initializing...
  if (!_can_add_vectors)
    {
      std::cerr << "ERROR: Too late.  Cannot add vectors to the system"
		<< std::endl
		<< " any more.  You should have done this earlier."
		<< std::endl;
      error();
    }

  // Make sure the vector isn't there already
  if (_other_vectors.count(vec_name))
    {
      std::cerr << "ERROR: vector "
		<< vec_name
		<< " has already been added for this system!"
		<< std::endl;
      error();
    }
  
  // Add the vector to the map
  NumericVector<Number>* buf(NumericVector<Number>::build(_solver_package).release());
  _other_vectors[vec_name] = buf;
}




const NumericVector<Number> &  SystemBase::get_vector(const std::string& vec_name) const
{
  // only enable access after the vectors are properly initialized
  if (_can_add_vectors)
    {
      std::cerr << "ERROR: Too early.  Access to vectors granted only when "
		<< std::endl
		<< "these are already properly initialized."
		<< std::endl;
      error();
    }

  // Make sure the vector exists
  other_vectors_const_iterator pos = _other_vectors.find(vec_name);
  
  if (pos == _other_vectors.end())
    {
      std::cerr << "ERROR: vector "
		<< vec_name
		<< " does not exist in this system!"
		<< std::endl;      
      error();
    }
  
  return *(pos->second);
}



NumericVector<Number> &  SystemBase::get_vector(const std::string& vec_name)
{
  // only enable access after the vectors are properly initialized
  if (_can_add_vectors)
    {
      std::cerr << "ERROR: Too early.  Access to vectors granted only when "
		<< std::endl
		<< "these are already properly initialized."
		<< std::endl;
      error();
    }

  // Make sure the vector exists
  other_vectors_const_iterator pos = _other_vectors.find(vec_name);
  
  if (pos == _other_vectors.end())
    {
      std::cerr << "ERROR: vector "
		<< vec_name
		<< " does not exist in this system!"
		<< std::endl;      
      error();
    }
  
  return *(pos->second);
}



// void SystemBase::zero_vectors()
// {
//   // do this only when we cannot add vectors anymore, since
//   // then the vectors are already initialized
//   if (_can_add_vectors)
//     {
//       std::cerr << "ERROR: Too early.  Can only zero the vectors when"
// 		<< std::endl
// 		<< "these are already properly initialized."
// 		<< std::endl;
//       error();
//     }

//   other_vectors_iterator pos = _other_vectors.begin();
//   if (pos == _other_vectors.end())
//     {
//       std::cerr << "ERROR: No additional vectors to zero!"
// 		<< std::endl;      
//       error(); 
//     }
//   else
//     // zero the other vectors
//     for(; pos != _other_vectors.end(); ++pos)
// 	pos->second->zero ();
// }



void SystemBase::add_variable (const std::string& var,
			       const FEType& type)
{  
  // Make sure the variable isn't there already
  if (_var_num.count(var))
    {
      std::cerr << "ERROR: variable "
		<< var
		<< " has already been added for this system!"
		<< std::endl;
      error();
    }
  
  // Add the variable to the list
  _var_names.push_back (var);
  _var_type[var]  = type;
  _var_num[var]   = (n_vars()-1);

  // Add the variable to the _dof_map
  _dof_map.add_variable (type);
}



void SystemBase::add_variable (const std::string& var,
			       const Order order)
{
  FEType fe_type(order, LAGRANGE);
  
  add_variable(var, fe_type);
}



unsigned short int SystemBase::variable_number (const std::string& var) const
{
  // Make sure the variable exists
  std::map<std::string, unsigned short int>::const_iterator
    pos = _var_num.find(var);
  
  if (pos == _var_num.end())
    {
      std::cerr << "ERROR: variable "
		<< var
		<< " does not exist in this system!"
		<< std::endl;      
      error();
    }
  
  return pos->second;
}


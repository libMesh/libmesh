// $Id: system_base.C,v 1.1 2003-02-12 02:03:49 ddreyer Exp $

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
#include "mesh_config.h"  /*  doxygen needs this for #ifdef ENABLE_AMR ??? */
#include "mesh.h"
#include "system_base.h"


// ------------------------------------------------------------
// SystemBase implementation

SystemBase::SystemBase (const Mesh& mesh,
			const std::string& name,
			const SolverPackage solver_package) :
  solution                (NumericVector::build(solver_package)),
  rhs                     (NumericVector::build(solver_package)),
  matrix                  (SparseMatrix::build (solver_package)),
  linear_solver_interface (LinearSolverInterface::build(solver_package)),
  _sys_name               (name),
  _dof_map                (mesh),
  _mesh                   (mesh)
{
};



SystemBase::~SystemBase ()
{
  SystemBase::clear ();
};


void SystemBase::clear ()
{
  _var_names.clear ();

  _var_type.clear ();
  
  _var_num.clear ();
  
  _dof_map.clear ();
  
  solution->clear ();

  rhs->clear ();

  matrix->clear ();
};




void SystemBase::init ()
{
  _dof_map.distribute_dofs (_mesh.processor_id());

#ifdef ENABLE_AMR

  _dof_map.create_dof_constraints();

#endif

  solution->init (n_dofs(), n_local_dofs());

  rhs->init      (n_dofs(), n_local_dofs());
};




void SystemBase::assemble ()
{
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
      matrix->attach_dof_map (_dof_map);

      // Compute the sparisty pattern for the current
      // mesh and DOF distribution
      _dof_map.compute_sparsity (_mesh.processor_id());

      // Initialize the matrix conformal to this
      // sparsity pattern
      matrix->init ();
    };

  // Clear the matrix and right-hand side.
  matrix->zero ();
  rhs->zero    ();

  // Now everything is set up and ready for matrix assembly
};





void SystemBase::update_global_solution (std::vector<Complex>& global_soln) const
{
  global_soln.resize (solution->size());

  solution->localize  (global_soln);
};



void SystemBase::update_global_solution (std::vector<Complex>& global_soln,
					 const unsigned int dest_proc) const
{
  global_soln.resize       (solution->size());

  solution->localize_to_one (global_soln, dest_proc);
};




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
  _dof_map.add_component (type);
};



void SystemBase::add_variable (const std::string& var,
			       const Order order)
{
  FEType fe_type(order, LAGRANGE);
  
  add_variable(var, fe_type);
};



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
};

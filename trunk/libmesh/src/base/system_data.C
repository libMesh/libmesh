// $Id: system_data.C,v 1.9 2003-02-10 03:55:51 benkirk Exp $

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
#include "mesh_config.h"
#include "petsc_interface.h"
#include "system_data.h"
#include "equation_systems.h"


// ------------------------------------------------------------
// SystemBase implementation

SystemBase::SystemBase (EquationSystems& es,
			const std::string& name) :
  dof_map(es.get_mesh()),
  init_system_fptr(NULL),
  assemble_fptr(NULL),
  sys_name(name),
  equation_systems(es),
  mesh(es.get_mesh())
{
};



SystemBase::~SystemBase ()
{
};



void GeneralSystem::clear ()
{
  var_names.clear ();

  var_type.clear ();
  
  var_num.clear ();
  
  dof_map.clear ();

  //init_system_fptr = assemble_fptr = NULL;
  
#ifdef HAVE_PETSC

  solution.clear ();

  rhs.clear ();

  matrix.clear ();
  
  current_local_solution.clear ();

  old_local_solution.clear ();
  
  older_local_solution.clear ();
  
#endif
};



void GeneralSystem::update ()
{
#ifndef HAVE_PETSC

  std::cerr << "ERROR: this feature requires Petsc support!"
	    << std::endl;
  error();

#else
  
  assert (current_local_solution.local_size() == solution.size());

  const std::vector<unsigned int>& send_list = dof_map.get_send_list ();

  assert (!send_list.empty());

  assert (send_list.size() <= solution.size());
  
  solution.localize (current_local_solution, send_list);
  
#endif
};



void SystemBase::update_global_solution (std::vector<Complex>& global_soln) const
{
#ifndef HAVE_PETSC

  std::cerr << "ERROR: this feature requires Petsc support!"
	    << std::endl;
  error();

#else
  
  global_soln.resize (solution.size());

  solution.localize  (global_soln);

#endif
};



void SystemBase::update_global_solution (std::vector<Complex>& global_soln,
					 const unsigned int dest_proc) const
{
#ifndef HAVE_PETSC

  std::cerr << "ERROR: this feature requires Petsc support!"
	    << std::endl;
  error();

#else
  
  global_soln.resize       (solution.size());

  solution.localize_to_one (global_soln, dest_proc);

#endif
};



void GeneralSystem::init ()
{
  dof_map.distribute_dofs (mesh.processor_id());

#ifdef ENABLE_AMR

  dof_map.create_dof_constraints();

#endif

#ifdef HAVE_PETSC
  
  solution.init (n_dofs(), n_local_dofs());

  rhs.init      (n_dofs(), n_local_dofs());

  current_local_solution.init (n_dofs ());

  old_local_solution.init     (n_dofs ());
  
  older_local_solution.init   (n_dofs ());
  
#endif

  // Possibly call a user-supplied initialization
  // method.
  if (init_system_fptr != NULL)
    {
      init_system_fptr (equation_systems, name());

      update();

#ifdef HAVE_PETSC
      old_local_solution   = current_local_solution;
      older_local_solution = current_local_solution; 
#endif
    };
};



void SystemBase::add_variable (const std::string& var,
			       const FEType& type)
{  
  // Make sure the variable isn't there already
  if (var_num.count(var))
    {
      std::cerr << "ERROR: variable "
		<< var
		<< " has already been added for this system!"
		<< std::endl;
      error();
    }
  
  // Add the variable to the list
  var_names.push_back (var);
  var_type[var]  = type;
  var_num[var]   = (n_vars()-1);

  // Add the variable to the dof_map
  dof_map.add_component (type);
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
    pos = var_num.find(var);
  
  if (pos == var_num.end())
    {
      std::cerr << "ERROR: variable "
		<< var
		<< " does not exist in this system!"
		<< std::endl;      
      error();
    }
  
  return pos->second;
};



#ifdef HAVE_PETSC

void GeneralSystem::assemble ()
{
  assert (assemble_fptr != NULL);

  // initialize the matrix if not 
  // already initialized
  if (!matrix.initialized())
    {
      // Tell the matrix and the dof map
      // about each other.  This is necessary
      // because the dof_map needs to share
      // information with the matrix during the
      // sparsity computation phase
      dof_map.attach_matrix (matrix);
      matrix.attach_dof_map (dof_map);

      // Compute the sparisty pattern for the current
      // mesh and DOF distribution
      dof_map.compute_sparsity (mesh.processor_id());

      // Initialize the matrix conformal to this
      // sparsity pattern
      matrix.init ();
    };

  // Clear the matrix and right-hand side.
  matrix.zero ();
  rhs.zero    ();

  // Call the user-specified matrix assembly function
  assemble_fptr (equation_systems, name());

  //matrix.print ();
  //rhs.print    ();
};



std::pair<unsigned int, Real>
GeneralSystem::solve ()
{
  // Assemble the linear system
  assemble (); 

  // Get the user-specifiied linear solver tolerance
  const Real tol            =
    equation_systems.parameter("linear solver tolerance");

  // Get the user-specified maximum # of linear solver iterations
  const unsigned int maxits =
    static_cast<unsigned int>(equation_systems.parameter("linear solver maximum iterations"));

  // Solve the linear system
  const std::pair<unsigned int, Real> rval = 
    petsc_interface.solve (matrix, solution, rhs, tol, maxits);

  // Update the local solution to reflect the new values
  update ();
  
  return rval; 
};

#endif



void SystemBase::attach_init_function(void fptr(EquationSystems& es,
						const std::string& name))
{
  assert (fptr != NULL);
  
  init_system_fptr = fptr;
};



void SystemBase::attach_assemble_function(void fptr(EquationSystems& es,
						    const std::string& name))
{
  assert (fptr != NULL);
  
  assemble_fptr = fptr;  
};

// $Id: system_data.C,v 1.2 2003-01-20 16:31:32 jwpeterson Exp $

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
#include <algorithm>

// Local includes
#include "petsc_interface.h"
#include "system_data.h"
#include "equation_systems.h"


// ------------------------------------------------------------
// SystemData implementation

SystemData::SystemData (EquationSystems& es,
			const std::string& name) :
  dof_map(es.get_mesh()),
  init_system(NULL),
  assemble(NULL),
  sys_name(name),
  equation_systems(es),
  mesh(es.get_mesh())
{
};



SystemData::~SystemData ()
{
};



void SystemData::clear ()
{
  var_names.clear ();

  var_type.clear ();
  
  var_num.clear ();
  
  dof_map.clear ();

  //init_system = assemble = NULL;
  
#ifdef HAVE_PETSC

  solution.clear ();

  rhs.clear ();

  matrix.clear ();
  
  current_local_solution.clear ();

  old_local_solution.clear ();
  
  older_local_solution.clear ();
  
#endif
};



void SystemData::update ()
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



void SystemData::update_global_solution (std::vector<number>& global_soln) const
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



void SystemData::update_global_solution (std::vector<number>& global_soln,
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



void SystemData::init ()
{
  dof_map.distribute_dofs (mesh.processor_id());

#ifdef ENABLE_AMR

  dof_map.create_dof_constraints();

#endif

#ifdef HAVE_PETSC
  
  solution.init (n_dofs(), n_local_dofs());

  rhs.init (n_dofs(), n_local_dofs());

  current_local_solution.init (n_dofs ());

  old_local_solution.init (n_dofs ());
  
  older_local_solution.init (n_dofs ());
  
#endif

  if (init_system != NULL)
    {
      init_system (equation_systems, name());

      update();

#ifdef HAVE_PETSC
      old_local_solution   = current_local_solution;
      older_local_solution = current_local_solution; 
#endif
    };
};



void SystemData::add_variable (const std::string& var,
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



void SystemData::add_variable (const std::string& var,
			       const Order order)
{
  FEType fe_type(order, LAGRANGE);
  
  add_variable(var, fe_type);
};



unsigned short int SystemData::variable_number (const std::string& var) const
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

std::pair<unsigned int, real>
SystemData::solve ()
{
  assert (assemble != NULL);

  // initialize the matrix if not 
  // already initialized
  if (!matrix.initialized())
    {
      dof_map.compute_sparsity (mesh.processor_id());
      
      matrix.init (dof_map);
    };

  matrix.zero ();
  rhs.zero    ();
  
  assemble (equation_systems, name()); 
  
  const real tol            =
    equation_systems.parameter("linear solver tolerance");
  
  const unsigned int maxits =
    static_cast<unsigned int>(equation_systems.parameter("linear solver maximum iterations"));

  const std::pair<unsigned int, real> rval = 
    petsc_interface.solve (matrix, solution, rhs, tol, maxits);

  // Update the local solution to reflect the new values
  update ();
  
  return rval; 
};

#endif



void SystemData::attach_init_function(void fptr(EquationSystems& es,
						const std::string& name))
{
  assert (fptr != NULL);
  
  init_system = fptr;
};



void SystemData::attach_assemble_function(void fptr(EquationSystems& es,
						    const std::string& name))
{
  assert (fptr != NULL);
  
  assemble = fptr;  
};

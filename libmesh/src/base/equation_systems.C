// $Id: equation_systems.C,v 1.37 2003-05-22 18:31:19 ddreyer Exp $

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


// System includes
#include <sstream>

// Local Includes
#include "fe_interface.h"
#include "libmesh.h"
#include "mesh.h"
#include "system_base.h"
#include "frequency_system.h"
#include "newmark_system.h"
#include "steady_system.h"
#include "transient_system.h"
// Include the systems before this one to avoid
// overlapping forward declarations.
#include "equation_systems.h"


// Forward Declarations




// ------------------------------------------------------------
// EquationSystems class implementation
EquationSystems::EquationSystems (const Mesh& m) :
  _mesh(m)
{
  // Set default parameters
  this->set_parameter("linear solver tolerance")          = 1.e-12;
  this->set_parameter("linear solver maximum iterations") = 5000;
}



EquationSystems::~EquationSystems ()
{
  this->clear ();
}



void EquationSystems::clear ()
{
  // Clear any user-supplied additional data
  data_map.clear ();

  // Clear any flags
  _flags.clear ();

  // Clear any parameters
  _parameters.clear ();

  // clear the systems.  We must delete them
  // since we newed them!
  {
    std::map<std::string, SystemBase*>::iterator
      pos = _systems.begin();
    
    for (; pos != _systems.end(); ++pos)
      {
	delete pos->second; pos->second = NULL;
      }
    
    _systems.clear ();
  }
}



void EquationSystems::init ()
{
   const unsigned int n_sys = this->n_systems();

 assert (n_sys != 0);

  /**
   * Tell all the \p DofObject entities how many systems
   * there are.
   */
  {
    // All the nodes
    const_node_iterator       node_it  (_mesh.nodes_begin());
    const const_node_iterator node_end (_mesh.nodes_end());
    
    for ( ; node_it != node_end; ++node_it)
      (*node_it)->set_n_systems(n_sys);
    
    // All the elements
    const_elem_iterator       elem_it (_mesh.elements_begin());
    const const_elem_iterator elem_end(_mesh.elements_end());
    
    for ( ; elem_it != elem_end; ++elem_it)
      (*elem_it)->set_n_systems(n_sys);
  }

  std::map<std::string, SystemBase*>::iterator
    pos = _systems.begin();
  
  for (; pos != _systems.end();  ++pos)
    {
      /**
       * Initialize the system.
       */
      pos->second->init();
    }
}



void EquationSystems::reinit ()
{
 const unsigned int n_sys = this->n_systems();

 assert (n_sys != 0);

  /**
   * Tell all the \p DofObject entities how many systems
   * there are.
   */
  {
    // All the nodes
    const_node_iterator       node_it  (_mesh.nodes_begin());
    const const_node_iterator node_end (_mesh.nodes_end());
    
    for ( ; node_it != node_end; ++node_it)
      (*node_it)->set_n_systems(n_sys);
    
    // All the elements
    const_elem_iterator       elem_it (_mesh.elements_begin());
    const const_elem_iterator elem_end(_mesh.elements_end());
    
    for ( ; elem_it != elem_end; ++elem_it)
      (*elem_it)->set_n_systems(n_sys);
  }

  std::map<std::string, SystemBase*>::iterator
    pos = _systems.begin();
  
  for (; pos != _systems.end();  ++pos)
    {
      /**
       * Initialize the system.
       */
      pos->second->reinit();
    }
}



void EquationSystems::add_system (const std::string& sys_type,
				  const std::string& name)
{
  // Build a Newmark system
  if      (sys_type == "Newmark")
    this->add_system<NewmarkSystem> (name);

  // Build a steady system
  else if (sys_type == "Steady")
    this->add_system<SteadySystem> (name);

  // build a transient system
  else if (sys_type == "Transient")
    this->add_system<TransientSystem> (name);

#if defined(USE_COMPLEX_NUMBERS)
  // build a frequency system
  else if (sys_type == "Frequency")
    this->add_system<FrequencySystem> (name);
#endif

  else
    {
      std::cerr << "ERROR: Unknown system type: "
		<< sys_type
		<< std::endl;
    }
}



template <typename T_sys>
void EquationSystems::add_system (const std::string& name)
{
  if (!_systems.count(name))
    {
      const unsigned int num = this->n_systems();

      _systems.insert (std::pair<std::string, SystemBase*>(name,
							   new T_sys(*this,
								     name,
								     num)
							   )
		       );
      
    }
  else
    {
      std::cerr << "ERROR: There was already a system"
		<< " named " << name
		<< std::endl;

      error();
    }

  
  /**
   * Tell all the \p DofObject entities to add a system.
   */
  {
    // All the nodes
    const_node_iterator       node_it  (_mesh.nodes_begin());
    const const_node_iterator node_end (_mesh.nodes_end());
    
    for ( ; node_it != node_end; ++node_it)
      (*node_it)->add_system();
    
    // All the elements
    const_elem_iterator       elem_it (_mesh.elements_begin());
    const const_elem_iterator elem_end(_mesh.elements_end());
    
    for ( ; elem_it != elem_end; ++elem_it)
      (*elem_it)->add_system();
  }
}



void EquationSystems::delete_system (const std::string& name)
{
  if (!_systems.count(name))
    {
      std::cerr << "ERROR: no system named "
                << name  << std::endl;

      error();
    }
  
  delete _systems[name];
  
  _systems.erase (name);
}



const SystemBase& EquationSystems::get_system (const std::string& name) const
{
  return this->get_system<SystemBase> (name);
}



SystemBase& EquationSystems::get_system (const std::string& name)
{
  return this->get_system<SystemBase> (name);
}




SystemBase & EquationSystems::operator () (const std::string& name)
{
  std::map<std::string, SystemBase*>::iterator
    pos = _systems.find(name);
   
  if (pos == _systems.end())
    {
      std::cerr << "ERROR: system "
                << name
                << " not found!"
                << std::endl;
 
      error();
    }
 
  return *(pos->second);
}
  
 
 
const SystemBase & EquationSystems::operator () (const std::string& name) const
{
  std::map<std::string, SystemBase*>::const_iterator
    pos = _systems.find(name);
   
  if (pos == _systems.end())
    {
      std::cerr << "ERROR: system "
                << name
                << " not found!"
                << std::endl;
 
      error();
    }
 
  return *(pos->second);
}
  
 
 
SystemBase & EquationSystems::operator () (const unsigned int num)
{
  assert (num < this->n_systems());
 
  std::map<std::string, SystemBase*>::iterator
    pos = _systems.begin();
   
  // New code
#if (__GNUC__ == 2)
  std::advance (pos, static_cast<int>(num));
#else
  std::advance (pos, num);
#endif
  
  return *(pos->second);
}



const SystemBase & EquationSystems::operator ()  (const unsigned int num) const
{
  assert (num < this->n_systems());

  std::map<std::string, SystemBase*>::const_iterator
    pos = _systems.begin();

  // New code
#if (__GNUC__ == 2)
  std::advance (pos, static_cast<int>(num));
#else
  std::advance (pos, num);
#endif

  return *(pos->second);
}

 
 
void EquationSystems::build_variable_names (std::vector<std::string>& var_names) const
{
  assert (n_systems());
  
  var_names.resize (n_vars());

  unsigned int var_num=0;

  std::map<std::string, SystemBase*>::const_iterator
    pos = _systems.begin();

  for (; pos != _systems.end(); ++pos)
    for (unsigned int vn=0; vn<pos->second->n_vars(); vn++)
      var_names[var_num++] = pos->second->variable_name(vn);       
}



void EquationSystems::build_solution_vector (std::vector<Number>& soln,
						 std::string& system_name,
						 std::string& variable_name) const
{
  error();

  // Get a reference to the named system
  const SystemBase& system = this->get_system(system_name);

  // Get the number associated with the variable_name we are passed
  const unsigned short int variable_num = system.variable_number(variable_name);

  // Get the dimension of the current mesh
  const unsigned int dim = _mesh.mesh_dimension();

  // If we're on processor 0, allocate enough memory to hold the solution.
  // Since we're only looking at one variable, there will be one solution value
  // for each node in the mesh.
  if (_mesh.processor_id() == 0)
    soln.resize(_mesh.n_nodes());

  // Vector to hold the global solution from all processors
  std::vector<Number> sys_soln;
  
  // Update the global solution from all processors
  system.update_global_solution (sys_soln, 0);
  
  // Temporary vector to store the solution on an individual element.
  std::vector<Number>       elem_soln;   

  // The FE solution interpolated to the nodes
  std::vector<Number>       nodal_soln;  

  // The DOF indices for the element
  std::vector<unsigned int> dof_indices; 

  // Determine the finite/infinite element type used in this system 
  const FEType& fe_type    = system.variable_type(variable_num);

  // Define iterators to iterate over all the elements of the mesh
  const_active_elem_iterator       it (_mesh.elements_begin());
  const const_active_elem_iterator end(_mesh.elements_end());

  // Loop over elements
  for ( ; it != end; ++it)
    {
      // Convenient shortcut to the element pointer
      const Elem* elem = *it;

      // Fill the dof_indices vector for this variable
      system.get_dof_map().dof_indices(elem,
				       dof_indices,
				       variable_num);

      // Resize the element solution vector to fit the
      // dof_indices for this element.
      elem_soln.resize(dof_indices.size());

      // Transfer the system solution to the element
      // solution by mapping it through the dof_indices vector.
      for (unsigned int i=0; i<dof_indices.size(); i++)
	elem_soln[i] = sys_soln[dof_indices[i]];

      // Using the FE interface, compute the nodal_soln
      // for the current elemnt type given the elem_soln
      FEInterface::nodal_soln (dim,
			       fe_type,
			       elem,
			       elem_soln,
			       nodal_soln);

      // Sanity check -- make sure that there are the same number
      // of entries in the nodal_soln as there are nodes in the
      // element!
      assert (nodal_soln.size() == elem->n_nodes());

      // Copy the nodal solution over into the correct place in
      // the global soln vector which will be returned to the user.
      for (unsigned int n=0; n<elem->n_nodes(); n++)
	soln[elem->node(n)] = nodal_soln[n];
    }
}




void EquationSystems::build_solution_vector (std::vector<Number>& soln) const
{
  assert (this->n_systems());

  const unsigned int dim = _mesh.mesh_dimension();
  const unsigned int nn  = _mesh.n_nodes();
  const unsigned int nv  = this->n_vars();

  // Only if we are on processor zero, allocate the storage
  // to hold (number_of_nodes)*(number_of_variables) entries.
  if (_mesh.processor_id() == 0)
    soln.resize(nn*nv);

  std::vector<Number> sys_soln; 
  
  unsigned int var_num=0;

  // For each system in this EquationSystems object,
  // update the global solution and if we are on processor 0,
  // loop over the elements and build the nodal solution
  // from the element solution.  Then insert this nodal solution
  // into the vector passed to build_solution_vector.
  std::map<std::string, SystemBase*>::const_iterator
    pos = _systems.begin();

  for (; pos != _systems.end(); ++pos)
    {  
      const SystemBase& system  = *(pos->second);
      const unsigned int nv_sys = system.n_vars();
      
      system.update_global_solution (sys_soln, 0);
      
      if (_mesh.processor_id() == 0)
	{
	  std::vector<Number>       elem_soln;   // The finite element solution
	  std::vector<Number>       nodal_soln;  // The FE solution interpolated to the nodes
	  std::vector<unsigned int> dof_indices; // The DOF indices for the finite element 
	      
	  for (unsigned int var=0; var<nv_sys; var++)
	    {
	      const FEType& fe_type    = system.variable_type(var);

	      const_active_elem_iterator       it (_mesh.elements_begin());
	      const const_active_elem_iterator end(_mesh.elements_end());

	      for ( ; it != end; ++it)
		{
		  const Elem* elem = *it;
		  system.get_dof_map().dof_indices (elem, dof_indices, var);
		  
		  elem_soln.resize(dof_indices.size());
		  
		  for (unsigned int i=0; i<dof_indices.size(); i++)
		    elem_soln[i] = sys_soln[dof_indices[i]];
		  
		  FEInterface::nodal_soln (dim,
					   fe_type,
					   elem,
					   elem_soln,
					   nodal_soln);

#ifdef ENABLE_INFINITE_ELEMENTS
		  // infinite elements should be skipped...
		  if (!elem->infinite())
#endif
		    { 
		      assert (nodal_soln.size() == elem->n_nodes());
		  
		      for (unsigned int n=0; n<elem->n_nodes(); n++)
			soln[nv*(elem->node(n)) + (var + var_num)] = nodal_soln[n];
		    }

		}
	    }	 
	}

      var_num += nv_sys;
    }
}



bool EquationSystems::compare (const EquationSystems& other_es, 
			       const Real threshold,
			       const bool verbose) const
{
  // safety check, whether we handle at least the same number
  // of systems
  std::vector<bool> os_result;

  if (this->n_systems() != other_es.n_systems())
    {
      if (verbose)
        {
	  std::cout << "  Fatal difference. This system handles " 
		    << this->n_systems() << " systems," << std::endl
		    << "  while the other system handles "
		    << other_es.n_systems() 
		    << " systems." << std::endl
		    << "  Aborting comparison." << std::endl;
	}
      return false;
    }
  else
    {
      // start comparing each system
      std::map<std::string, SystemBase*>::const_iterator pos=_systems.begin();
      
      for (; pos != _systems.end(); ++pos)
        {
	  const std::string& sys_name = pos->first;
	  const SystemBase&  system        = *(pos->second);
      
	  // get the other system
	  const SystemBase& other_system   = other_es.get_system (sys_name);

	  os_result.push_back (system.compare (other_system, threshold, verbose));

	}

    }
  

  // sum up the results
  if (os_result.size()==0)
    return true;
  else
    {
      bool os_identical;
      unsigned int n = 0;
	  do
	    {
	      os_identical = os_result[n];
	      n++;
	    }
	  while (os_identical && n<os_result.size());
	  return os_identical;
	}
}



std::string EquationSystems::get_info () const
{
  std::ostringstream out;
  
  out << " EquationSystems" << std::endl
      << "  n_systems()=" << this->n_systems() << std::endl;

  // Print the info for the individual systems
  std::map<std::string, SystemBase*>::const_iterator it=_systems.begin();

  for (; it != _systems.end(); ++it)
    out << it->second->get_info();
    
  
  // Possibly print the flags
  if (!_flags.empty())
    {  
      out << "  n_flags()=" << this->n_flags() << std::endl;
      out << "   Flags:" << std::endl;
      
      for (std::set<std::string>::const_iterator flag = _flags.begin();
	   flag != _flags.end(); ++flag)
	out << "    "
	    << "\""
	    << *flag
	    << "\""
	    << std::endl;
    }

  
  // Possibly print the parameters  
  if (!_parameters.empty())
    {  
      out << "  n_parameters()=" << this->n_parameters() << std::endl;
      out << "   Parameters:" << std::endl;
      
      for (std::map<std::string, Real>::const_iterator
	     param = _parameters.begin(); param != _parameters.end();
	   ++param)
	out << "    "
	    << "\""
	    << param->first
	    << "\""
	    << "="
	    << param->second
	    << std::endl;
    }
  
  return out.str();
}



void EquationSystems::print_info () const
{
  std::cout << this->get_info()
	    << std::endl;
}



unsigned int EquationSystems::n_vars () const
{
  if (_systems.empty())
    return 0;

  unsigned int tot=0;

  std::map<std::string, SystemBase*>::const_iterator pos = _systems.begin();
  
  for (; pos != _systems.end(); ++pos)
    tot += pos->second->n_vars();

  return tot;
}



unsigned int EquationSystems::n_dofs () const
{
  if (_systems.empty())
    return 0;

  unsigned int tot=0;

  std::map<std::string, SystemBase*>::const_iterator pos = _systems.begin();
  
  for (; pos != _systems.end(); ++pos)
    tot += pos->second->n_dofs();

  return tot;      
}





// void* & EquationSystems::additional_data (const std::string& name)
// {
//   // Check for the entry already.  If it is there return the pointer,
//   // but make sure it isn't NULL
//   if (_additional_data.count(name) != 0)
//     {
//       assert (_additional_data[name] != NULL);

//       return _additional_data[name];
//     }

//   return _additional_data[name];
// }



// void EquationSystems::unset_additional_data (const std::string& name)
// {
//   // Look for an entry matching name
//   std::map<std::string, void*>::iterator
//     pos = _additional_data.find(name);

//   // Remove it if an entry was found
//   if (pos != _additional_data.end())
//     _additional_data.erase(pos);
// }



bool EquationSystems::flag (const std::string& fl) const
{
  return (_flags.count(fl) != 0);
}



void EquationSystems::set_flag (const std::string& fl)
{
  _flags.insert (fl);
}



void EquationSystems::unset_flag (const std::string& fl)
{
  // Look for the flag in the database
  std::set<std::string>::iterator
    pos = _flags.find(fl);

  // If the flag was found remove it
  if (pos != _flags.end())
    _flags.erase(pos);
}



Real EquationSystems::parameter (const std::string& id) const
{
  return this->parameter<Real> (id);
}



Real & EquationSystems::set_parameter (const std::string& id)
{
#ifdef DEBUG
  /*
  // Make sure the parameter isn't already set
  if (parameters.count(id))
    std::cerr << "WARNING: parameter "
	      << "\""
	      << id
	      << "\""
	      << " is already set!"
	      << std::endl;
  */
#endif
  
  // Insert the parameter/value pair into the database
  return _parameters[id];
}



void EquationSystems::unset_parameter (const std::string& id)
{
  // Look for the id in the database
  std::map<std::string, Real>::iterator
    pos = _parameters.find(id);
  
  // If the parameter was found remove it
  if (pos != _parameters.end())
    _parameters.erase(pos);
}



bool EquationSystems::parameter_exists (const std::string& id) const
{
  // Look for the id in the database
  std::map<std::string, Real>::const_iterator
    pos = _parameters.find(id);

  return (pos != _parameters.end());
}


//--------------------------------------------------------------------------
// Instantiations of templated functions

// Some systems only work with complex numbers
#ifdef USE_COMPLEX_NUMBERS
template void EquationSystems::add_system<FrequencySystem> (const std::string&);
#endif

template void EquationSystems::add_system<NewmarkSystem>   (const std::string&);
template void EquationSystems::add_system<SteadySystem>    (const std::string&);
template void EquationSystems::add_system<TransientSystem> (const std::string&);

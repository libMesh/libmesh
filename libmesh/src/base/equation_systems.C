// $Id: equation_systems.C,v 1.31 2003-04-09 19:26:57 ddreyer Exp $

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
#include "libmesh.h"
#include "fe_interface.h"
#include "linear_solver_interface.h"
#include "equation_systems.h"
#include "general_system.h"
#include "frequency_system.h"
#include "thin_system.h"
#include "newmark_system.h"
#include "equation_systems_macro.h"

// Forward Declarations




// ------------------------------------------------------------
// EquationSystems<T_sys> class implementation
template <typename T_sys>
EquationSystems<T_sys>::EquationSystems (Mesh& m,
					 const SolverPackage sp) :
  EquationSystemsBase(m, sp)
{
  // Default parameters
  this->set_parameter("linear solver tolerance")          = 1.e-12;
  this->set_parameter("linear solver maximum iterations") = 5000;
}



template <typename T_sys>
EquationSystems<T_sys>::~EquationSystems ()
{
  this->clear();

  assert (!libMesh::closed());
}



template <typename T_sys>
void EquationSystems<T_sys>::clear ()
{
  typename std::map<std::string, T_sys*>::iterator pos = _systems.begin();
  
  for (; pos != _systems.end(); ++pos)
    delete pos->second;
  
 _systems.clear ();

  EquationSystemsBase::clear();
}



template <typename T_sys>
void EquationSystems<T_sys>::init ()
{
 const unsigned int n_sys = this->n_systems();

 assert (n_sys != 0);

 EquationSystemsBase::init();

  /**
   * Tell all the \p DofObject entities how many systems
   * there are.
   */
  {
    // All the nodes
    node_iterator       node_it  (_mesh.nodes_begin());
    const node_iterator node_end (_mesh.nodes_end());
    
    for ( ; node_it != node_end; ++node_it)
      (*node_it)->set_n_systems(n_sys);
    
    // All the elements
    elem_iterator       elem_it (_mesh.elements_begin());
    const elem_iterator elem_end(_mesh.elements_end());
    
    for ( ; elem_it != elem_end; ++elem_it)
      (*elem_it)->set_n_systems(n_sys);
  }

  typename std::map<std::string, T_sys*>::iterator sys = _systems.begin();
  
  for (; sys != _systems.end();  ++sys)
    {
      /**
       * Initialize the system.
       */
      sys->second->init();
    }
}



template <typename T_sys>
void EquationSystems<T_sys>::add_system (const std::string& name)
{
  if (!_systems.count(name))
    {
      const unsigned int num = this->n_systems();

      _systems.insert (std::pair<std::string, T_sys*>(name,
						      new T_sys(*this,
								name,
								num,
								_solver_package)
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
    node_iterator       node_it  (_mesh.nodes_begin());
    const node_iterator node_end (_mesh.nodes_end());
    
    for ( ; node_it != node_end; ++node_it)
      (*node_it)->add_system();
    
    // All the elements
    elem_iterator       elem_it (_mesh.elements_begin());
    const elem_iterator elem_end(_mesh.elements_end());
    
    for ( ; elem_it != elem_end; ++elem_it)
      (*elem_it)->add_system();
  }
}



template <typename T_sys>
void EquationSystems<T_sys>::delete_system (const std::string& name)
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



template <typename T_sys>
unsigned int EquationSystems<T_sys>::n_vars () const
{
  if (_systems.empty())
    return 0;

  unsigned int tot=0;

  typename std::map<std::string, T_sys*>::const_iterator pos = _systems.begin();
  
  for (; pos != _systems.end(); ++pos)
    tot += pos->second->n_vars();

  return tot;
}



template <typename T_sys>
unsigned int EquationSystems<T_sys>::n_dofs () const
{
  if (_systems.empty())
    return 0;

  unsigned int tot=0;

  typename std::map<std::string, T_sys*>::const_iterator pos = _systems.begin();
  
  for (; pos != _systems.end(); ++pos)
    tot += pos->second->n_dofs();

  return tot;      
}




template <typename T_sys>
T_sys & EquationSystems<T_sys>::operator () (const std::string& name)
{
  typename std::map<std::string, T_sys*>::iterator pos = _systems.find(name);
  
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




template <typename T_sys>
const T_sys & EquationSystems<T_sys>::operator () (const std::string& name) const
{
  typename std::map<std::string, T_sys*>::const_iterator pos = _systems.find(name);
  
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



template <typename T_sys>
const std::string & EquationSystems<T_sys>::name (const unsigned int num) const
{
  assert (num < this->n_systems());

  typename std::map<std::string, T_sys*>::const_iterator pos = _systems.begin();

  // New code
#if (__GNUC__ == 2)
  std::advance (pos, static_cast<int>(num));
#else
  std::advance (pos, num);
#endif

  return pos->first;
}



template <typename T_sys>
T_sys & EquationSystems<T_sys>::operator () (const unsigned int num)
{
  assert (num < this->n_systems());

  typename std::map<std::string, T_sys*>::iterator pos = _systems.begin();
  
  // New code
#if (__GNUC__ == 2)
  std::advance (pos, static_cast<int>(num));
#else
  std::advance (pos, num);
#endif

  return *(pos->second);
}



template <typename T_sys>
const T_sys & EquationSystems<T_sys>::operator ()  (const unsigned int num) const
{
  assert (num < this->n_systems());

  typename std::map<std::string, T_sys*>::const_iterator pos = _systems.begin();
  
    // New code
#if (__GNUC__ == 2)
  std::advance (pos, static_cast<int>(num));
#else
  std::advance (pos, num);
#endif
  
  return *(pos->second);
}




template <typename T_sys>
void EquationSystems<T_sys>::build_variable_names (std::vector<std::string>& var_names)
{
  assert (n_systems());
  
  var_names.resize (n_vars());

  unsigned int var_num=0;
  
  for (unsigned int sys=0; sys<this->n_systems(); sys++)
    for (unsigned int vn=0; vn < (*this)(sys).n_vars(); vn++)
      var_names[var_num++] = (*this)(sys).variable_name(vn);	   
}





template <typename T_sys>
void EquationSystems<T_sys>::build_solution_vector (std::vector<Number>& soln,
						    std::string& system_name,
						    std::string& variable_name)
{
  error();

  // Get a reference to the named system
  const T_sys& system = (*this)(system_name);

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
  active_elem_iterator       it (_mesh.elements_begin());
  const active_elem_iterator end(_mesh.elements_end());

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




template <typename T_sys>
void EquationSystems<T_sys>::build_solution_vector (std::vector<Number>& soln)
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
  for (unsigned int sys=0; sys<this->n_systems(); sys++)
    {
      const T_sys& system       = (*this)(sys);	      
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

	      active_elem_iterator       it (_mesh.elements_begin());
	      const active_elem_iterator end(_mesh.elements_end());

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

 		  if (nodal_soln.size() == elem->n_nodes())
		    for (unsigned int n=0; n<elem->n_nodes(); n++)
 		      soln[nv*(elem->node(n)) + (var + var_num)] =
 		        nodal_soln[n];
// OLD CODE		  
// 		  assert (nodal_soln.size() == elem->n_nodes());
		  assert (nodal_soln.size() == elem->n_nodes());
		  
		  for (unsigned int n=0; n<elem->n_nodes(); n++)
		    soln[nv*(elem->node(n)) + (var + var_num)] = nodal_soln[n];
		}
	    }	 
	}

      var_num += nv_sys;
    }
}



template <typename T_sys>
bool EquationSystems<T_sys>::compare (const EquationSystems<T_sys>& other_es, 
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
      typename std::map<std::string, T_sys*>::const_iterator pos=_systems.begin();
      
      for (; pos != _systems.end(); ++pos)
        {
	  const std::string& sys_name = pos->first;
	  const T_sys&  system        = *(pos->second);
      
	  // get the other system
	  const T_sys& other_system   = other_es (sys_name);

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





template <typename T_sys>
std::string EquationSystems<T_sys>::get_info () const
{
  std::ostringstream out;
  
  out << " EquationSystems<" << T_sys::system_type() << ">:" << std::endl
      << "  n_systems()=" << this->n_systems() << std::endl;

  typename std::map<std::string, T_sys*>::const_iterator it=_systems.begin();

  for (; it != _systems.end(); ++it)
    out << it->second->get_info();

  return out.str();
}






//--------------------------------------------------------------
// Explicit instantiations using the macro from equation_systems_macro.h

INSTANTIATE_EQUATION_SYSTEMS;

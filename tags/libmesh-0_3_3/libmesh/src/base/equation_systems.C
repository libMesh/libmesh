// $Id: equation_systems.C,v 1.26 2003-03-21 15:29:11 ddreyer Exp $

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

// Forward Declarations




// ------------------------------------------------------------
// EquationSystems<T_sys> class implementation
template <typename T_sys>
EquationSystems<T_sys>::EquationSystems (Mesh& m,
					 const SolverPackage sp) :
  EquationSystemsBase(m, sp)
{
  // remember the type of system as flag, useful when writing equation systems
  this->set_flag(T_sys::system_type());
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
void EquationSystems<T_sys>::build_solution_vector (std::vector<Number>& soln)
{
  assert (this->n_systems());

  const unsigned int dim = _mesh.mesh_dimension();
  const unsigned int nn  = _mesh.n_nodes();
  const unsigned int nv  = this->n_vars();

  if (_mesh.processor_id() == 0)
    soln.resize(nn*nv);

  std::vector<Number> sys_soln; 
  
  unsigned int var_num=0;

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
		  
		  FEInterface::nodal_soln (dim, fe_type, elem,
					   elem_soln, nodal_soln);
		  
		  assert (nodal_soln.size() == elem->n_nodes());
		  
		  for (unsigned int n=0; n<elem->n_nodes(); n++)
		    soln[nv*(elem->node(n)) + (var + var_num)] =
		      nodal_soln[n];
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
    {
      const std::string& sys_name = it->first;
      const T_sys&  system        = *(it->second);
      
      out << "   System \"" << sys_name << "\"" << std::endl
	  << "    Variables=";
      for (unsigned int vn=0; vn<system.n_vars(); vn++)
	out << "\"" << system.variable_name(vn) << "\" ";
     
      out << std::endl;

#ifndef ENABLE_INFINITE_ELEMENTS
      out << "    Finite Element Types=";
      for (unsigned int vn=0; vn<system.n_vars(); vn++)
      {
	out << "\"" << system.get_dof_map().variable_type(vn).family << "\" ";
      }
#else
      out << "    Finite Element Types=";
      for (unsigned int vn=0; vn<system.n_vars(); vn++)
      {
	out << "\"" << system.get_dof_map().variable_type(vn).family << "\", ";
	out << "\"" << system.get_dof_map().variable_type(vn).radial_family << "\" ";
      }

      out << std::endl << "    Infinite Element Mapping=";
      for (unsigned int vn=0; vn<system.n_vars(); vn++)
      {
	out << "\"" << system.get_dof_map().variable_type(vn).inf_map << "\" ";
      }
#endif      

      out << std::endl;
      
      out << "    Approximation Orders=";
      for (unsigned int vn=0; vn<system.n_vars(); vn++)
      {
#ifndef ENABLE_INFINITE_ELEMENTS
	out << "\"" << system.get_dof_map().variable_type(vn).order << "\" ";
#else
	out << "\"" << system.get_dof_map().variable_type(vn).order << "\", ";
	out << "\"" << system.get_dof_map().variable_type(vn).radial_order << "\" ";
#endif
      }

      out << std::endl;
      
      out << "    n_dofs()="             << system.n_dofs()             << std::endl;
      out << "    n_local_dofs()="       << system.n_local_dofs()       << std::endl;
#ifdef ENABLE_AMR
      out << "    n_constrained_dofs()=" << system.n_constrained_dofs() << std::endl;
#endif
    }
  
  return out.str();
}






//--------------------------------------------------------------
// Explicit instantiations
template class EquationSystems<GeneralSystem>;
template class EquationSystems<ThinSystem>;

#if defined(USE_COMPLEX_NUMBERS) 
template class EquationSystems<FrequencySystem>;
#endif

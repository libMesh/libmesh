// $Id: equation_systems.C,v 1.18 2003-02-17 05:33:09 jwpeterson Exp $

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
//#include "enum_solver_package.h"
#include "libmesh.h"
#include "fe_interface.h"
//#include "petsc_vector.h"
//#include "petsc_matrix.h"
//#include "petsc_interface.h"
#include "linear_solver_interface.h"
#include "equation_systems.h"
#include "general_system.h"

// Forward Declarations




// ------------------------------------------------------------
// EquationSystem class implementation
EquationSystems::EquationSystems (Mesh& m,
				  const SolverPackage sp) :
  _mesh(m),
  _solver_package(sp)
{
  // Default parameters
  set_parameter("linear solver tolerance")          = 1.e-12;
  set_parameter("linear solver maximum iterations") = 5000;
}



EquationSystems::~EquationSystems ()
{
  clear();

  assert (!libMesh::closed());
}



void EquationSystems::clear ()
{
  for (std::map<std::string, GeneralSystem*>::iterator
	 pos = _systems.begin(); pos != _systems.end();
       ++pos)
    delete pos->second;
  
  _systems.clear ();

  _flags.clear ();

  _parameters.clear ();
}



void EquationSystems::init ()
{
  const unsigned int n_sys = n_systems();
  
  assert (n_sys != 0);

//   /**
//    * Tell all the \p DofObject entities how many systems
//    * there are.
//    */
//   {
//     // All the nodes
//     node_iterator       node_it  (_mesh.nodes_begin());
//     const node_iterator node_end (_mesh.nodes_end());
    
//     for ( ; node_it != node_end; ++node_it)
//       (*node_it)->set_n_systems(n_sys);
    
//     // All the elements
//     elem_iterator       elem_it (_mesh.elements_begin());
//     const elem_iterator elem_end(_mesh.elements_end());
    
//     for ( ; elem_it != elem_end; ++elem_it)
//       (*elem_it)->set_n_systems(n_sys);
//   }


  for (std::map<std::string, GeneralSystem*>::iterator
	 sys = _systems.begin(); sys != _systems.end();
       ++sys)
    {
      /**
       * Initialize the system.
       */
      sys->second->init();
    }
}



void EquationSystems::add_system (const std::string& name)
{
  if (!_systems.count(name))
    {
      const unsigned int num = n_systems();
      
      _systems.insert (std::pair<std::string,
		                 GeneralSystem*>(name,
						 new GeneralSystem(*this,
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



unsigned int EquationSystems::n_vars () const
{
  if (_systems.empty())
    return 0;

  unsigned int tot=0;
  
  for (std::map<std::string, GeneralSystem*>::const_iterator
	 pos = _systems.begin(); pos != _systems.end(); ++pos)
    tot += pos->second->n_vars();

  return tot;      
}



unsigned int EquationSystems::n_dofs () const
{
  if (_systems.empty())
    return 0;

  unsigned int tot=0;
  
  for (std::map<std::string, GeneralSystem*>::const_iterator
	 pos = _systems.begin(); pos != _systems.end(); ++pos)
    tot += pos->second->n_dofs();

  return tot;      
}



GeneralSystem & EquationSystems::operator () (const std::string& name)
{
  std::map<std::string, GeneralSystem*>::iterator
    pos = _systems.find(name);
  
  if (pos == _systems.end())
    {
      std::cerr << "ERROR: system "
		<< name
		<< " not found!"
		<< std::endl;

      error();
    }

  return *pos->second;
}



const GeneralSystem & EquationSystems::operator () (const std::string& name) const
{
  std::map<std::string, GeneralSystem*>::const_iterator
    pos = _systems.find(name);
  
  if (pos == _systems.end())
    {
      std::cerr << "ERROR: system "
		<< name
		<< " not found!"
		<< std::endl;

      error();
    }

  return *pos->second;
}




const std::string & EquationSystems::name (const unsigned int num) const
{
  assert (num < n_systems());

  std::map<std::string, GeneralSystem*>::const_iterator
    pos = _systems.begin();

  // New code
#if (__GNUC__ == 2)
  std::advance (pos, static_cast<int>(num));
#else
  std::advance (pos, num);
#endif

  // Old code
//  for (unsigned int i=0; i<num; i++)
//    ++pos;

  return pos->first;
}




GeneralSystem & EquationSystems::operator () (const unsigned int num)
{
  assert (num < n_systems());

  std::map<std::string, GeneralSystem*>::iterator
    pos = _systems.begin();
  
  // New code
#if (__GNUC__ == 2)
  std::advance (pos, static_cast<int>(num));
#else
  std::advance (pos, num);
#endif

  // Old code
//  for (unsigned int i=0; i<num; i++)
//    ++pos;

  return *pos->second;
}



const GeneralSystem & EquationSystems::operator ()  (const unsigned int num) const
{
  assert (num < n_systems());

  std::map<std::string, GeneralSystem*>::const_iterator
    pos = _systems.begin();
  
    // New code
#if (__GNUC__ == 2)
  std::advance (pos, static_cast<int>(num));
#else
  std::advance (pos, num);
#endif
  
  // Old code
  //  for (unsigned int i=0; i<num; i++)
//    ++pos;
  
  return *pos->second;
}




bool EquationSystems::flag (const std::string& fl) const
{
  return (_flags.count(fl) != 0);
}



void EquationSystems::set_flag (const std::string& fl)
{
#ifdef DEBUG
  /*
  // Make sure the parameter isn't already set
  if (flags.count(fl))
    std::cerr << "WARNING: flag "
      	      << "\""
	      << fl
      	      << "\""
	      << " is already set!"
	      << std::endl;
  */
#endif

  _flags.insert (fl);
}



void EquationSystems::unset_flag (const std::string& fl)
{
  // Look for the flag in the database
  if (!_flags.count(fl))
    {
      std::cerr << "ERROR: flag " << fl
		<< " was not set!"
		<< std::endl;
      error(); 
    }

  // Remove the flag
  _flags.erase (fl);  
}



Real EquationSystems::parameter (const std::string& id) const
{
  // Look for the id in the database
  std::map<std::string, Real>::const_iterator
    pos = _parameters.find(id);
  
  if (pos == _parameters.end())
    {
      std::cerr << "ERROR: parameter " << id
		<< " was not set!"
		<< std::endl;
      error();
    }
  
  // Return the parameter value if found
  return pos->second;
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
  
  // Make sure the parameter was found
  if (pos == _parameters.end())
    {
      std::cerr << "ERROR: parameter " << id
		<< " was not set!"
		<< std::endl;
      error();
    }
  
  // Erase the entry
  _parameters.erase(pos);
}



void EquationSystems::build_variable_names (std::vector<std::string>& var_names)
{
  assert (n_systems());
  
  var_names.resize (n_vars());

  unsigned int var_num=0;
  
  for (unsigned int sys=0; sys<n_systems(); sys++)
    for (unsigned int vn=0; vn < (*this)(sys).n_vars(); vn++)
      var_names[var_num++] = (*this)(sys).variable_name(vn);	   
}



void EquationSystems::build_solution_vector (std::vector<Complex>& soln)
{
  assert (n_systems());

  const unsigned int dim = _mesh.mesh_dimension();
  const unsigned int nn  = _mesh.n_nodes();
  const unsigned int nv  = n_vars();

  if (_mesh.processor_id() == 0)
    soln.resize(nn*nv);

  std::vector<Complex> sys_soln; 
  
  unsigned int var_num=0;

  for (unsigned int sys=0; sys<n_systems(); sys++)
    {
      const GeneralSystem& system  = (*this)(sys);	      
      const unsigned int nv_sys    = system.n_vars();
      
      system.update_global_solution (sys_soln);

      if (_mesh.processor_id() == 0)
	{
	  std::vector<Complex>      elem_soln;   // The finite element solution
	  std::vector<Complex>      nodal_soln;  // The finite elemnt solution interpolated to the nodes
	  std::vector<unsigned int> dof_indices; // The DOF indices for the finite element 
	      
	  for (unsigned int var=0; var<nv_sys; var++)
	    {
	      const FEType& fe_type    = system.variable_type(var);
	      
	      for (unsigned int e=0; e<_mesh.n_elem(); e++)
		if (_mesh.elem(e)->active())
		  {
		    const Elem* elem = _mesh.elem(e);
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



std::string EquationSystems::get_info () const
{
  std::ostringstream out;
  
  out << " EquationSystems:" << std::endl
      << "  n_systems()=" << n_systems() << std::endl;
  
  for (std::map<std::string, GeneralSystem*>::const_iterator it=_systems.begin();
       it != _systems.end(); ++it)
    {
      const std::string& sys_name    = it->first;
      const GeneralSystem&  system   = *it->second;
      
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
  
  
  if (!_flags.empty())
    {  
      out << "  Flags:" << std::endl;
      
      for (std::set<std::string>::const_iterator flag = _flags.begin();
	   flag != _flags.end(); ++flag)
	out << "   "
	    << "\""
	    << *flag
	    << "\""
	    << std::endl;
    }
  
  
  if (!_parameters.empty())
    {  
      out << "  Parameters:" << std::endl;
      
      for (std::map<std::string, Real>::const_iterator
	     param = _parameters.begin(); param != _parameters.end();
	   ++param)
	out << "   "
	    << "\""
	    << param->first
	    << "\""
	    << "="
	    << param->second
	    << std::endl;
    }
  
  return out.str();
}

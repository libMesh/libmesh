// $Id: equation_systems.C,v 1.3 2003-01-20 17:06:15 jwpeterson Exp $

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
#include "petsc_vector.h"
#include "petsc_matrix.h"
#include "petsc_interface.h"
#include "system_data.h"
#include "equation_systems.h"

// Forward Declarations




// ------------------------------------------------------------
// EquationSystem class implementation
EquationSystems::EquationSystems (const Mesh& m, const bool up) :
  mesh(m),
  use_petsc(up)
{
  // Default parameters
  set_parameter("linear solver tolerance")          = 1.e-12;
  set_parameter("linear solver maximum iterations") = 5000;
};



EquationSystems::~EquationSystems ()
{
};



void EquationSystems::clear ()
{
  systems.clear ();

  flags.clear ();

  parameters.clear ();
};



void EquationSystems::init ()
{
  assert (!systems.empty());


  for (std::map<std::string, SystemData>::iterator
	 sys = systems.begin(); sys != systems.end();
       ++sys)
    {
//        /**
//         * Pass our parameters on to the system
//         */
//        for (std::map<std::string, real>::const_iterator
//  	     param = parameters.begin();
//  	   param != parameters.end(); ++param)
//  	sys->second.set_parameter(param->first) =
//  	  param->second;

//        /**
//         * Pass our flags on to the system
//         */
//        for (std::set<std::string>::const_iterator
//  	     fl = flags.begin(); fl != flags.end();
//  	   ++fl)
//  	sys->second.set_flag(*fl);

      /**
       * Initialize the system.
       */
      sys->second.init();
    };
};



void EquationSystems::add_system (const std::string& name)
{
  if (!systems.count(name))
    {  
      SystemData sd(*this, name);
      
      std::pair<std::string, SystemData>
	kv(name, sd);
      
      systems.insert (kv);
    }
  else
    {
      std::cerr << "WARNING: There was already a system"
		<< " named " << name
		<< std::endl;
    }
};



void EquationSystems::delete_system (const std::string& name)
{
  if (!systems.count(name))
    {
      std::cerr << "ERROR: no system named "
		<< name  << std::endl;

      error();
    }

  systems.erase (name);
};



unsigned int EquationSystems::n_vars () const
{
  if (systems.empty())
    return 0;

  unsigned int tot=0;
  
  for (std::map<std::string, SystemData>::const_iterator
	 pos = systems.begin(); pos != systems.end(); ++pos)
    tot += pos->second.n_vars();

  return tot;      
};



unsigned int EquationSystems::n_dofs () const
{
  if (systems.empty())
    return 0;

  unsigned int tot=0;
  
  for (std::map<std::string, SystemData>::const_iterator
	 pos = systems.begin(); pos != systems.end(); ++pos)
    tot += pos->second.n_dofs();

  return tot;      
};



SystemData & EquationSystems::operator () (const std::string& name)
{
  std::map<std::string, SystemData>::iterator
    pos = systems.find(name);
  
  if (pos == systems.end())
    {
      std::cerr << "ERROR: system "
		<< name
		<< " not found!"
		<< std::endl;

      error();
    }

  return pos->second;
};



const SystemData & EquationSystems::operator () (const std::string& name) const
{
  std::map<std::string, SystemData>::const_iterator
    pos = systems.find(name);
  
  if (pos == systems.end())
    {
      std::cerr << "ERROR: system "
		<< name
		<< " not found!"
		<< std::endl;

      error();
    }

  return pos->second;
};




const std::string & EquationSystems::name (const unsigned int num) const
{
  assert (num < n_systems());

  std::map<std::string, SystemData>::const_iterator
    pos = systems.begin();
  
  for (unsigned int i=0; i<num; i++)
    ++pos;

  return pos->first;
};




SystemData & EquationSystems::operator () (const unsigned int num)
{
  assert (num < n_systems());

  std::map<std::string, SystemData>::iterator
    pos = systems.begin();
  
  for (unsigned int i=0; i<num; i++)
    ++pos;

  return pos->second;
};



const SystemData & EquationSystems::operator ()  (const unsigned int num) const
{
  assert (num < n_systems());

  std::map<std::string, SystemData>::const_iterator
    pos = systems.begin();
  
  for (unsigned int i=0; i<num; i++)
    ++pos;

  return pos->second;
};




bool EquationSystems::flag (const std::string& fl) const
{
  return (flags.count(fl) != 0);
};



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

  flags.insert (fl);
};



void EquationSystems::unset_flag (const std::string& fl)
{
  // Look for the flag in the database
  if (!flags.count(fl))
    {
      std::cerr << "ERROR: flag " << fl
		<< " was not set!"
		<< std::endl;
      error(); 
    }

  // Remove the flag
  flags.erase (fl);  
};



real EquationSystems::parameter (const std::string& id) const
{
  // Look for the id in the database
  std::map<std::string, real>::const_iterator
    pos = parameters.find(id);
  
  if (pos == parameters.end())
    {
      std::cerr << "ERROR: parameter " << id
		<< " was not set!"
		<< std::endl;
      error();
    }
  
  // Return the parameter value if found
  return pos->second;
};



real & EquationSystems::set_parameter (const std::string& id)
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
  return parameters[id];
};



void EquationSystems::unset_parameter (const std::string& id)
{
  // Look for the id in the database
  std::map<std::string, real>::iterator
    pos = parameters.find(id);
  
  // Make sure the parameter was found
  if (pos == parameters.end())
    {
      std::cerr << "ERROR: parameter " << id
		<< " was not set!"
		<< std::endl;
      error();
    }
  
  // Erase the entry
  parameters.erase(pos);
};



void EquationSystems::build_variable_names (std::vector<std::string>& var_names)
{
  assert (n_systems());
  
  var_names.resize (n_vars());

  unsigned int var_num=0;
  
  for (unsigned int sys=0; sys<n_systems(); sys++)
    for (unsigned int vn=0; vn < (*this)(sys).n_vars(); vn++)
      var_names[var_num++] = (*this)(sys).variable_name(vn);	   
};



void EquationSystems::build_solution_vector (std::vector<number>& soln)
{
  assert (n_systems());

  const unsigned int dim = mesh.mesh_dimension();
  const unsigned int nn  = mesh.n_nodes();
  const unsigned int nv  = n_vars();

  if (mesh.processor_id() == 0)
    soln.resize(nn*nv);

  std::vector<number> sys_soln; 
  
  unsigned int var_num=0;

  for (unsigned int sys=0; sys<n_systems(); sys++)
    {
      const unsigned int nv_sys = (*this)(sys).n_vars();
      
      (*this)(sys).update_global_solution (sys_soln);

      if (mesh.processor_id() == 0)
	{
	  std::vector<number>         elem_soln; // The finite element solution
	  std::vector<number>         nodal_soln;  // The finite elemnt solution interpolated to the nodes
	  std::vector<unsigned int> dof_indices; // The DOF indices for the finite element 
	      
	  for (unsigned int var=0; var<nv_sys; var++)
	    {
	      const FEType fe_type = (*this)(sys).variable_type(var);
	      
	      for (unsigned int e=0; e<mesh.n_elem(); e++)
		if (mesh.elem(e)->active())
		  {
		    const Elem* elem = mesh.elem(e);
		    (*this)(sys).dof_map.dof_indices (e, dof_indices, var);
		    
		    elem_soln.resize(dof_indices.size());

		    for (unsigned int i=0; i<dof_indices.size(); i++)
		      elem_soln[i] = sys_soln[dof_indices[i]];
		    
		    FEInterface::nodal_soln (dim, fe_type, elem,
					     elem_soln, nodal_soln);

		    assert (nodal_soln.size() == elem->n_nodes());
		    
		    for (unsigned int n=0; n<elem->n_nodes(); n++)
		      soln[nv*(elem->node(n)) + (var + var_num)] =
			nodal_soln[n];
		  };
	    };	 
	};

      var_num += nv_sys;
    };
};



std::string EquationSystems::get_info () const
{
  std::ostringstream out;
  
  out << " EquationSystems:" << std::endl
      << "  n_systems()=" << n_systems() << std::endl;
  
  for (std::map<std::string, SystemData>::const_iterator it=systems.begin();
       it != systems.end(); ++it)
    {
      const std::string& sys_name = it->first;
      const SystemData&  system   = it->second;
      
      out << "   System \"" << sys_name << "\"" << std::endl
	  << "    Variables=";
      for (unsigned int vn=0; vn<system.n_vars(); vn++)
	out << "\"" << system.variable_name(vn) << "\" ";
     
      out << std::endl;

      out << "    Finite Element Types=";
      for (unsigned int vn=0; vn<system.n_vars(); vn++)
      {
#ifndef ENABLE_INFINITE_ELEMENTS
	out << "\"" << system.dof_map.component_type(vn).family << "\" ";
#else
	out << "(" << system.dof_map.component_type(vn).family << ",";
	out << system.dof_map.component_type(vn).base_family << ") ";
#endif
      };      

      out << std::endl;
      
      out << "    Approximation Orders=";
      for (unsigned int vn=0; vn<system.n_vars(); vn++)
	out << "\"" << system.dof_map.component_order(vn) << "\" ";
      
      out << std::endl;
      
      out << "    n_dofs()="             << system.n_dofs()             << std::endl;
      out << "    n_local_dofs()="       << system.n_local_dofs()       << std::endl;
#ifdef ENABLE_AMR
      out << "    n_constrained_dofs()=" << system.n_constrained_dofs() << std::endl;
#endif
    };
  
  
  if (!flags.empty())
    {  
      out << "  Flags:" << std::endl;
      
      for (std::set<std::string>::const_iterator flag = flags.begin();
	   flag != flags.end(); ++flag)
	out << "   "
	    << "\""
	    << *flag
	    << "\""
	    << std::endl;
    };
  
  
  if (!parameters.empty())
    {  
      out << "  Parameters:" << std::endl;
      
      for (std::map<std::string, real>::const_iterator
	     param = parameters.begin(); param != parameters.end();
	   ++param)
	out << "   "
	    << "\""
	    << param->first
	    << "\""
	    << "="
	    << param->second
	    << std::endl;
    };
  
  return out.str();
};

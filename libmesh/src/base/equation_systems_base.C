// $Id: equation_systems_base.C,v 1.1 2003-03-11 04:35:18 ddreyer Exp $

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
#include "mesh.h"
#include "equation_systems_base.h"

// Forward Declarations




// ------------------------------------------------------------
// EquationSystemsBase class implementation
EquationSystemsBase::EquationSystemsBase (Mesh& m,
					  const SolverPackage sp) :
  _mesh(m),
  _solver_package(sp)
{
  // Default parameters
  this->set_parameter("linear solver tolerance")          = 1.e-12;
  this->set_parameter("linear solver maximum iterations") = 5000;
}



EquationSystemsBase::~EquationSystemsBase ()
{
  /*
   * Do _not_ call \p this->clear(), since
   * the Destructor of EquationSystems<T_sys>
   * already invokes our clear().
   */
}



void EquationSystemsBase::clear ()
{
  _flags.clear ();

  _parameters.clear ();
}



void EquationSystemsBase::init ()
{
  assert (_mesh.is_prepared());
}




bool EquationSystemsBase::flag (const std::string& fl) const
{
  return (_flags.count(fl) != 0);
}



void EquationSystemsBase::set_flag (const std::string& fl)
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



void EquationSystemsBase::unset_flag (const std::string& fl)
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



Real EquationSystemsBase::parameter (const std::string& id) const
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



Real & EquationSystemsBase::set_parameter (const std::string& id)
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



void EquationSystemsBase::unset_parameter (const std::string& id)
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







std::string EquationSystemsBase::get_info () const
{
  std::ostringstream out;
  
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


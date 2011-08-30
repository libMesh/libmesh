// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2008 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner
  
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

#ifndef __variable_h__
#define __variable_h__

// C++ includes
#include <set>
#include <string>

// Local Includes
#include "fe_type.h"
#include "id_types.h"

namespace libMesh {

/**
 * This class defines the notion of a variable in the system.
 * A variable is one of potentially several unknowns in the 
 * problem at hand.  A variable is described by a unique 
 * name, a finite element approximation family, and 
 * (optionally) a list of subdomains to which the 
 * variable is restricted.
 */  
class Variable
{
public:
  
  /**
   * Constructor.  Omits the subdomain mapping, hence this
   * constructor creates a variable which is active on 
   * all subdomains.
   */
  Variable (const std::string &var_name,
	    const unsigned int var_number,
	    const FEType &var_type) :
    _name(var_name),
    _number(var_number),
    _type(var_type),
    _active_subdomains()
  {}
  
  /**
   * Constructor.  Takes a set which contains the subdomain
   * indices for which this variable is active.
   */ 
  Variable (const std::string &var_name,
	    const unsigned int var_number,
	    const FEType &var_type,
	    const std::set<subdomain_id_type> &var_active_subdomains) :
    _name(var_name),
    _number(var_number),
    _type(var_type),
    _active_subdomains(var_active_subdomains)
  {}
  
  /**
   * Arbitrary, user-specified name of the variable.
   */
  const std::string & name() const 
  { return _name; }

  /**
   * The rank of this variable in the system.
   */
  unsigned int number() const 
  { return _number; }

  /**
   * The \p FEType for this variable.
   */
  const FEType & type() const 
  { return _type; }

  /**
   * \p returns \p true if this variable is active on subdomain \p sid,
   * \p false otherwise.  Note that we interperet the special case of an 
   * empty \p _active_subdomains container as active everywhere, i.e. 
   * for all subdomains.
   */
  bool active_on_subdomain (const subdomain_id_type sid) const
  { return (_active_subdomains.empty() || _active_subdomains.count(sid));  }

  /**
   * \p returns \p true if this variable is active on all subdomains
   * because it has no specified activity map.  This can be used
   * to perform more efficient computations in some places.
   */
  bool implicitly_active () const
  { return _active_subdomains.empty(); }

  /**
   * Returns set of subdomain ids this variable lives on
   */
  const std::set<subdomain_id_type> & active_subdomains() const
  { return _active_subdomains; }

private:
  std::string             _name; 
  unsigned int            _number;
  FEType                  _type;
  std::set<subdomain_id_type> _active_subdomains;
};

} // namespace libMesh

#endif // #define __variable_h__

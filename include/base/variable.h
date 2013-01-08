// The libMesh Finite Element Library.
// Copyright (C) 2002-2012 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#ifndef LIBMESH_VARIABLE_H
#define LIBMESH_VARIABLE_H

// Local Includes
#include "libmesh/libmesh_common.h"
#include "libmesh/fe_type.h"
#include "libmesh/id_types.h"

// C++ includes
#include <set>
#include <string>
#include <vector>

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
	    const unsigned int first_scalar_number,
	    const FEType &var_type) :
    _name(var_name),
    _number(var_number),
    _first_scalar_number(first_scalar_number),
    _type(var_type),
    _active_subdomains()
  {}

  /**
   * Constructor.  Takes a set which contains the subdomain
   * indices for which this variable is active.
   */
  Variable (const std::string &var_name,
	    const unsigned int var_number,
	    const unsigned int first_scalar_number,
	    const FEType &var_type,
	    const std::set<subdomain_id_type> &var_active_subdomains) :
    _name(var_name),
    _number(var_number),
    _first_scalar_number(first_scalar_number),
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
   * The index of the first scalar component of this variable in the
   * system.
   */
  unsigned int first_scalar_number() const
  { return _first_scalar_number; }

  /**
   * The \p FEType for this variable.
   */
  const FEType & type() const
  { return _type; }

  /**
   * The number of components of this variable.
   */
  unsigned int n_components() const
  { return type().family == SCALAR ? _type.order : 1; }

  /**
   * \p returns \p true if this variable is active on subdomain \p sid,
   * \p false otherwise.  Note that we interperet the special case of an
   * empty \p _active_subdomains container as active everywhere, i.e.
   * for all subdomains.
   */
  bool active_on_subdomain (subdomain_id_type sid) const
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

protected:
  std::string             _name;
  unsigned int            _number;
  unsigned int            _first_scalar_number;
  FEType                  _type;
  std::set<subdomain_id_type> _active_subdomains;
};



/**
 * This class defines a logically grouped set of variables in
 * the system.  \p VariableGroup is appropriate for representing
 * several unknowns in the problem that are all approximated
 * with the same finite element approximation family and
 * (optionally) a list of subdomains to which the
 * variables are restricted.
 */
class VariableGroup : public Variable
{
public:
  /**
   * Constructor.  Omits the subdomain mapping, hence this
   * constructor creates a variable which is active on
   * all subdomains.
   */
  VariableGroup (const std::vector<std::string> &var_names,
		 const unsigned int var_number,
		 const unsigned int first_scalar_number,
		 const FEType &var_type) :
    Variable ("var_group",
	      var_number,
	      first_scalar_number,
	      var_type),
    _names(var_names)
  {}

   
  /**
   * Constructor.  Takes a set which contains the subdomain
   * indices for which this variable is active.
   */
  VariableGroup (const std::vector<std::string> &var_names,
		 const unsigned int var_number,
		 const unsigned int first_scalar_number,
		 const FEType &var_type,
		 const std::set<subdomain_id_type> &var_active_subdomains) :
     
    Variable ("var_group",
	      var_number,
	      first_scalar_number,
	      var_type,
	      var_active_subdomains),
    _names(var_names)
  {}

  /**
   * The number of variables in this \p VariableGroup
   */
  unsigned int n_variables () const
  { return _names.size(); }

  /**
   * Construct a \p Variable object for an individual member
   * of our group.
   */
  Variable variable (unsigned int v) const
  {
    libmesh_assert_less (v, this->n_variables());
    return Variable (this->name(v),
		     this->number(v),
		     this->first_scalar_number(v),
		     this->type(),
		     this->active_subdomains());
  }

  /**
   * Support vg(v) - returns a \p Variable for v.
   */
  Variable operator() (unsigned int v) const
  { return this->variable(v); }
  
  /**
   * Arbitrary, user-specified name of the variable.
   */
  const std::string & name(unsigned int v) const
  {
    libmesh_assert_less (v, this->n_variables());
    return _names[v];
  }

  /**
   * The rank of this variable in the system.
   */
  unsigned int number(unsigned int v) const
  {
    libmesh_assert_less (v, this->n_variables());
    return _number + v;
  }

  /**
   * The index of the first scalar component of this variable in the
   * system.
   */
  unsigned int first_scalar_number(unsigned int v) const
  {
    libmesh_assert_less (v, this->n_variables());
    if ((this->type().family == SCALAR) &&
	(this->n_variables() > 1)) libmesh_error(); // [BSK] I am not yet sure what this means!
    return _first_scalar_number;
  }

  /**
   * Appends a variable to the group.  Really only can be used by \p System in
   * a very limited window of opportunity - after the user specifies variables 
   * but before the system is initialized. 
   */
  void append (const std::string &var_name)
  { _names.push_back (var_name); }

protected:
  std::vector<std::string> _names;
};
  
} // namespace libMesh

#endif // LIBMESH_VARIABLE_H

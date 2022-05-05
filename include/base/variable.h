// The libMesh Finite Element Library.
// Copyright (C) 2002-2022 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

namespace libMesh
{

// Forward Declaration
class System;

/**
 * This class defines the notion of a variable in the system.
 * A variable is one of potentially several unknowns in the
 * problem at hand.  A variable is described by a unique
 * name, a finite element approximation family, and
 * (optionally) a list of subdomains to which the
 * variable is restricted.
 *
 * \author Roy Stogner
 * \date 2010
 * \brief A variable which is solved for in a System of equations.
 */
class Variable
{
public:

  /**
   * Constructor.  Omits the subdomain mapping, hence this
   * constructor creates a variable which is active on
   * all subdomains.
   */
  Variable (System * sys,
            std::string var_name,
            const unsigned int var_number,
            const unsigned int first_scalar_num,
            const FEType & var_type) :
    _sys(sys),
    _name(std::move(var_name)),
    _active_subdomains(),
    _number(var_number),
    _first_scalar_number(first_scalar_num),
    _type(var_type)
  {}

  /**
   * Constructor.  Takes a set which contains the subdomain
   * indices for which this variable is active.
   */
  Variable (System * sys,
            std::string var_name,
            const unsigned int var_number,
            const unsigned int first_scalar_num,
            const FEType & var_type,
            const std::set<subdomain_id_type> & var_active_subdomains) :
    _sys(sys),
    _name(std::move(var_name)),
    _active_subdomains(var_active_subdomains),
    _number(var_number),
    _first_scalar_number(first_scalar_num),
    _type(var_type)
  {}

  /**
   * Standard constructors.
   */
  Variable (const Variable &) = default;
  Variable & operator= (const Variable &) = default;
  Variable (Variable &&) = default;
  Variable & operator= (Variable &&) = default;

  /**
   * \returns true iff the \p other Variable has the same
   * characteristics and system numbering as this one.
   */
  bool operator== ( const Variable & other) const {
    return (_sys == other._sys) &&
           (_name == other._name) &&
           (_active_subdomains == other._active_subdomains) &&
           (_first_scalar_number == other._first_scalar_number) &&
           (_type == other._type);
  }

  /**
   * \returns A pointer to the System this Variable is part of.
   */
  System * system() const
  {
    return _sys;
  }

  /**
   * \returns The user-specified name of the variable.
   */
  const std::string & name() const
  { return _name; }

  /**
   * \returns The rank of this variable in the system.
   */
  unsigned int number() const
  { return _number; }

  /**
   * \returns The index of the first scalar component of this variable in the
   * system.
   */
  unsigned int first_scalar_number() const
  { return _first_scalar_number; }

  /**
   * \returns The \p FEType for this variable.
   */
  const FEType & type() const
  { return _type; }

  /**
   * \returns The number of components of this variable.
   */
  unsigned int n_components() const
  { return type().family == SCALAR ? _type.order.get_order() : 1; }

  /**
   * \returns \p true if this variable is active on subdomain \p sid,
   * \p false otherwise.
   *
   * \note We interpret the special case of an empty \p
   * _active_subdomains container as active everywhere, i.e. for all
   * subdomains.
   */
  bool active_on_subdomain (subdomain_id_type sid) const
  { return (_active_subdomains.empty() || _active_subdomains.count(sid));  }

  /**
   * \returns \p true if this variable is active on all subdomains
   * because it has no specified activity map.  This can be used
   * to perform more efficient computations in some places.
   */
  bool implicitly_active () const
  { return _active_subdomains.empty(); }

  /**
   * \returns The set of subdomain ids this variable lives on.
   */
  const std::set<subdomain_id_type> & active_subdomains() const
  { return _active_subdomains; }

protected:
  System *                _sys;
  std::string             _name;
  std::set<subdomain_id_type> _active_subdomains;
  unsigned int            _number;
  unsigned int            _first_scalar_number;
  FEType                  _type;
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
  VariableGroup (System * sys,
                 std::vector<std::string> var_names,
                 const unsigned int var_number,
                 const unsigned int first_scalar_num,
                 const FEType & var_type) :
    Variable (sys,
              "var_group",
              var_number,
              first_scalar_num,
              var_type),
    _names(std::move(var_names))
  {}


  /**
   * Constructor.  Takes a set which contains the subdomain
   * indices for which this variable is active.
   */
  VariableGroup (System * sys,
                 std::vector<std::string> var_names,
                 const unsigned int var_number,
                 const unsigned int first_scalar_num,
                 const FEType & var_type,
                 const std::set<subdomain_id_type> & var_active_subdomains) :

    Variable (sys,
              "var_group",
              var_number,
              first_scalar_num,
              var_type,
              var_active_subdomains),
    _names(std::move(var_names))
  {}

  /**
   * Standard constructors.
   */
  VariableGroup (const VariableGroup &) = default;
  VariableGroup & operator= (const VariableGroup &) = default;
  VariableGroup (VariableGroup &&) = default;
  VariableGroup & operator= (VariableGroup &&) = default;

  /**
   * \returns true iff the \p other VariableGroup has exactly the same
   * Variable members as this one.
   */
  bool operator== ( const VariableGroup & other) const {
    return (this->Variable::operator==(other)) &&
           (_names == other._names);
  }

  /**
   * \returns The number of variables in this \p VariableGroup
   */
  unsigned int n_variables () const
  { return cast_int<unsigned int>(_names.size()); }

  /**
   * \returns A \p Variable object constructed for an individual member
   * of our group.
   */
  Variable variable (unsigned int v) const
  {
    libmesh_assert_less (v, this->n_variables());
    return Variable (this->system(),
                     this->name(v),
                     this->number(v),
                     this->first_scalar_number(v),
                     this->type(),
                     this->active_subdomains());
  }

  /**
   * Support vg(v).
   *
   * \returns A \p Variable for v.
   */
  Variable operator() (unsigned int v) const
  { return this->variable(v); }

  /**
   * \returns The user-specified name of the variable.
   */
  const std::string & name(unsigned int v) const
  {
    libmesh_assert_less (v, this->n_variables());
    return _names[v];
  }

  /**
   * \returns The rank of this variable in the system.
   */
  unsigned int number(unsigned int v) const
  {
    libmesh_assert_less (v, this->n_variables());
    return _number + v;
  }

  // Don't let number(uint) hide number()
  using Variable::number;

  /**
   * \returns The index of the first scalar component of this variable in the
   * system.
   */
  unsigned int first_scalar_number(unsigned int v) const
  {
    libmesh_assert_less (v, this->n_variables());
    return _first_scalar_number+v;
  }

  /**
   * Appends a variable to the group.  Really only can be used by \p System in
   * a very limited window of opportunity - after the user specifies variables
   * but before the system is initialized.
   */
  void append (std::string var_name)
  { _names.push_back (std::move(var_name)); }

protected:
  std::vector<std::string> _names;
};

} // namespace libMesh

#endif // LIBMESH_VARIABLE_H

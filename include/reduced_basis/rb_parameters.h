// rbOOmit: An implementation of the Certified Reduced Basis method.
// Copyright (C) 2009, 2010 David J. Knezevic

// This file is part of rbOOmit.

// rbOOmit is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// rbOOmit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

#ifndef LIBMESH_RB_PARAMETERS_H
#define LIBMESH_RB_PARAMETERS_H

// libMesh includes
#include "libmesh/libmesh_common.h"

// C++ includes
#include <string>
#include <map>
#include <set>
#include <vector>

namespace libMesh
{

/**
 * This class is part of the rbOOmit framework.
 *
 * This class defines a set of parameters index by strings.
 *
 * \author David J. Knezevic
 * \date 2012
 */
class RBParameters
{
public:

  /**
   * The special functions can be defaulted for this class, as it
   * does not manage any memory itself.
   */
  RBParameters () = default;
  RBParameters (RBParameters &&) = default;
  RBParameters (const RBParameters &) = default;
  RBParameters & operator= (const RBParameters &) = default;
  RBParameters & operator= (RBParameters &&) = default;
  ~RBParameters() = default;

  /**
   * Constructor. Set parameters based on the std::map \p parameter_map.
   *
   * This constructor will still be supported once we switch over to
   * the vector-based storage for RBParameters objects. It will just set
   * the 0th entry of the vector corresponding to each parameter name.
   */
  RBParameters(const std::map<std::string, Real> & parameter_map);

  /**
   * Define a constant iterator for this class. This custom iterator
   * design is copied from the chunked_mapvector class.
   */
  class const_iterator
  {
  public:
    // Typedefs needed for interoperating with other STL algorithms and containers.
    typedef std::forward_iterator_tag iterator_category;
    typedef std::pair<std::string, Real> value_type;
    typedef std::ptrdiff_t difference_type;
    typedef value_type* pointer;
    typedef value_type& reference;

    // Underlying iterator type, must match the container type of _parameters.
    typedef std::map<std::string, std::vector<Real>>::const_iterator MapIter;

    // Constructor
    const_iterator(const MapIter & it,
                   const std::size_t vec_index)
      : _it(it),
        _vec_index(vec_index)
    {}

    // Copy ctor
    const_iterator(const const_iterator & i) = default;

    // Prefix increment operator "++it"
    const_iterator & operator++()
    {
      // First increment the vector index
      ++_vec_index;

      // If _vec_index goes past the current vector, start at beginning of next one
      if (_vec_index >= _it->second.size())
        {
          _vec_index = 0;
          ++_it;
        }

      return *this;
    }

    // Postfix increment operator "it++". This is actually less
    // efficient than doing the prefix increment, so if nothing
    // requires it, let's skip defining it.
    // const_iterator operator++(int)
    // {
    //   const_iterator i = *this;
    //   ++(*this);
    //   return i;
    // }

    // Dereference operator - returns a const reference to the value
    // indexed by it and _vec_index. Note: this is needed for backwards
    // compatibility but it is not the most efficient thing since we
    // need to construct a std::pair every time we dereference a
    // const_iterator.
    const value_type &
    operator*() const
    {
      _emulator = std::make_pair(_it->first, _it->second[_vec_index]);
      return _emulator;
    }

    // Equivalence comparison operator.
    bool operator==(const const_iterator & other) const
    {
      return (_it == other._it && _vec_index == other._vec_index);
    }

    // Non-equvialence comparison operator
    bool operator!=(const const_iterator & other) const
    {
      return !(other == *this);
    }

  private:
    // Give RBParameters access to our private members. At the moment
    // this is not needed since the RBParameters class does not really
    // need to know anything about the const_iterator once it has been
    // created.
    // friend class RBParameters;

    // Iterator into real container
    MapIter _it;

    // Accompanying current vector index into it->second
    std::size_t _vec_index;

    // Returned by the operator* function. Emulates dereferencing a
    // map<string, Real> iterator.
    mutable value_type _emulator;
  };

  /**
   * Clear this object.
   */
  void clear();

  /**
   * \returns true if there is a parameter named "param_name" present
   * in this class, false otherwise.
   */
  bool has_value(const std::string & param_name) const;

  /**
   * \returns true if there is an extra parameter named "param_name" present
   * in this class, false otherwise.
   */
  bool has_extra_value(const std::string & param_name) const;

  /**
   * Get the value of the specified parameter, throwing an error if it
   * does not exist.
   */
  Real get_value(const std::string & param_name) const;

  /**
   * Get the value of the specified parameter, returning the provided
   * default value if it does not exist.
   */
  Real get_value(const std::string & param_name, const Real & default_val) const;

  /**
   * Get the value of the specified parameter at the specified step,
   * throwing an error if it does not exist.
   */
  Real get_step_value(const std::string & param_name, std::size_t index) const;

  /**
   * Get the value of the specified parameter at the specified step,
   * returning the provided default value if either the parameter is
   * not defined or the step is invalid.
   */
  Real get_step_value(const std::string & param_name, std::size_t index, const Real & default_val) const;

  /**
   * Set the value of the specified parameter. If param_name
   * doesn't already exist, it is added to the RBParameters object.
   * For backwards compatibility, calling this function sets up
   * "param_name" to be a single-entry vector with "value" as the
   * only entry.
   */
  void set_value(const std::string & param_name, Real value);

  /**
   * Set the value of the specified parameter at the specified vector
   * index.  Note: each parameter is now allowed to be vector-valued,
   * it is up to the user to organize what the vector indices refer to
   * (e.g. load or time steps).
   */
  void set_value(const std::string & param_name, std::size_t index, Real value);

  /**
   * Similar to set_value(name, index, value) but instead of specifying a particular
   * index, just appends one more. Calling push_back_value() many times is more efficient
   * than calling set_value(name, index, value) many times because it takes advantage
   * of the std::vector's size-doubling t reduce allocations.
   */
  void push_back_value(const std::string & param_name, Real value);

  /**
   * Get the value of the specified extra parameter, throwing an error
   * if it does not exist.
   */
  Real get_extra_value(const std::string & param_name) const;

  /**
   * Get the value of the specified extra parameter, returning the
   * provided default value if it does not exist.
   */
  Real get_extra_value(const std::string & param_name, const Real & default_val) const;

  /**
   * Set the value of the specified extra parameter. If param_name
   * doesn't already exist, it is added to the extra parameters.
   */
  void set_extra_value(const std::string & param_name, Real value);

  /**
   * Set the value of the specified extra parameter at the specified vector
   * index.  Note: each parameter is now allowed to be vector-valued,
   * it is up to the user to organize what the vector indices refer to
   * (e.g. load or time steps).
   */
  void set_extra_value(const std::string & param_name, std::size_t index, Real value);

  /**
   * Get the number of parameters that have been added.
   */
  unsigned int n_parameters() const;

  /**
   * Returns the number of steps stored for all parameters. For
   * simplicity, we require all parameters to store the same number of
   * steps ("step" here may refer to time step or load step) and in
   * debug mode we actually verify that is the case.
   */
  unsigned int n_steps() const;

  /**
   * Fill \p param_names with the names of the parameters.
   *
   * \deprecated to avoid making it too easy to create copies that in
   * most circumstances aren't needed.  If you really need a list of
   * the parameter names, the best approach is to iterate over this
   * object using the begin()/end() APIs and build a std::set that
   * way.
   */
  void get_parameter_names(std::set<std::string> & param_names) const;

  /**
   * Fill \p param_names with the names of the extra parameters.
   *
   * \deprecated to avoid making it too easy to create copies that in
   * most circumstances aren't needed.  If you really need a list of
   * the parameter names, the best approach is to iterate over this
   * object using the extra_begin()/extra_end() APIs and build a std::set that
   * way.
   */
  void get_extra_parameter_names(std::set<std::string> & param_names) const;

  /**
   * Erase \p param_name  from _parameters. If \p param_name is not present
   * in _parameters, then do nothing.
   */
  void erase_parameter(const std::string & param_name);

  /**
   * Erase \p param_name  from _extra_parameters. If \p param_name is not present
   * in _extra_parameters, then do nothing.
   */
  void erase_extra_parameter(const std::string & param_name);

  /**
   * Get const_iterator access to the parameters stored in this RBParameters object.
   */
  const_iterator begin() const;
  const_iterator end() const;

  /**
   * Get const_iterator access to the extra parameters stored in this RBParameters object.
   */
  const_iterator extra_begin() const;
  const_iterator extra_end() const;

  /**
   * Two RBParameters are equal if they have the same _parameters map.
   */
  bool operator== (const RBParameters & rhs) const;

  /**
   * \returns !(*this == rhs).
   */
  bool operator!= (const RBParameters & rhs) const;

  /**
   * Append "rhs" to "*this".  Both RBParameters objects must have the
   * same n_steps(), otherwise an error is thrown. If some of the
   * parameter names overlap, then the values from rhs overwrite
   * *this. Both parameters and "extra" parameters are appended.
   */
  RBParameters & operator+= (const RBParameters & rhs);

  /**
   * Get a string that specifies the contents of this RBParameters object.
   * \p precision specifies the number of digits of precision we use
   * in scientific notation in the string.
   */
  std::string get_string(unsigned int precision=6) const;

  /**
   * Print the parameters.
   */
  void print() const;

private:

  /**
   * Helper function for the 3-parameter versions of set_value() and
   * set_extra_value().
   */
  void set_value_helper(std::map<std::string, std::vector<Real>> & map,
                        const std::string & param_name,
                        std::size_t index,
                        Real value);

  /**
   * The map that stores the actual parameter vectors. Each vector is
   * indexed by a name.
   */
  std::map<std::string, std::vector<Real>> _parameters;

  /**
   * The map that stores extra parameter vectors not used for RB
   * training. Each vector is indexed by a name.
   */
  std::map<std::string, std::vector<Real>> _extra_parameters;
};

} // namespace libMesh


#endif // LIBMESH_RB_PARAMETERS_H

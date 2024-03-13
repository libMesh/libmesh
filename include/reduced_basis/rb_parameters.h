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
 * Typedef for an individual RB parameter. Each parameter is now stored
 * as a vector of values (different from the vector of samples).
 */
using RBParameter = std::vector<Real>;

/**
 * This class is part of the rbOOmit framework.
 *
 * This class defines a set of parameters indexed by strings.
 * Multiple samples can be defined, where a sample is a set of values for
 * each parameter. The parameters themselves can be multi-valued, e.g.
 * for storing a matrix-type parameter.
 *
 * \author David J. Knezevic
 * \date 2012
 */
class RBParameters
{
public:

  /**
   * Constructor. Initializes the _n_samples parameter to 1 for
   * backwards compatibility, but the set_n_samples() function can
   * always be called later to update this value.
   */
  RBParameters();

  /**
   * The special functions can be defaulted for this class, as it
   * does not manage any memory itself.
   */
  RBParameters (RBParameters &&) = default;
  RBParameters (const RBParameters &) = default;
  RBParameters & operator= (const RBParameters &) = default;
  RBParameters & operator= (RBParameters &&) = default;
  ~RBParameters() = default;

  /**
   * Constructor. Set parameters based on the std::map \p parameter_map.
   *
   * It sets the values as the 0th entry of the sample-vector
   * corresponding to each parameter name.
   */
  RBParameters(const std::map<std::string, Real> & parameter_map);

  /**
   * Return const_iterators to the internal parameter map, as a convenient
   * way to access the parameter names and values.
   * For example: for (const auto & [param_name, sample_vec] : my_parameters)
   */
  std::map<std::string,std::vector<RBParameter>>::const_iterator begin() const;
  std::map<std::string,std::vector<RBParameter>>::const_iterator end() const;
  std::map<std::string,std::vector<RBParameter>>::const_iterator extra_begin() const;
  std::map<std::string,std::vector<RBParameter>>::const_iterator extra_end() const;

  /**
   * Define a constant iterator for iterating over the map of parameters.
   * This will iterate over every individual value in the map, meaning all
   * three levels (param name, sample vector, value vector).
   * This custom iterator design is copied from the chunked_mapvector class.
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
    typedef std::map<std::string, std::vector<RBParameter>>::const_iterator MapIter;

    // Constructor
    const_iterator(const MapIter & it,
                   const std::size_t sample_vec_index,
                   const std::size_t value_vec_index)
      : _it(it),
        _sample_vec_index(sample_vec_index),
        _value_vec_index(value_vec_index)
    {}

    // Copy ctor
    const_iterator(const const_iterator & i) = default;

    // Prefix increment operator "++it"
    const_iterator & operator++()
    {
      // First increment the value-vector index.
      // If _value_vec_index goes beyond the current value-vector, reset to the next one.
      ++_value_vec_index;
      if (_value_vec_index >= _it->second[_sample_vec_index].size())
        {
          _value_vec_index = 0;
          // Now increment the sample-vector index, and do the same check.
          // If we go beyond the current sample-vector, reset to the next one.
          ++_sample_vec_index;
          if (_sample_vec_index >= _it->second.size())
            {
              _sample_vec_index = 0;
              ++_it;
            }
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
    const value_type & operator*() const
    {
      _emulator = std::make_pair(_it->first, _it->second[_sample_vec_index][_value_vec_index]);
      return _emulator;
    }

    // Equivalence comparison operator.
    bool operator==(const const_iterator & other) const
    {
      return (_it == other._it && _sample_vec_index == other._sample_vec_index);
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

    // Accompanying current sample-vector index into it->second
    std::size_t _sample_vec_index;

    // Accompanying current value-vector index into it->second[_sample_vec_index]
    std::size_t _value_vec_index;

    // Returned by the operator* function. Emulates dereferencing a
    // map<string, Real> iterator.
    mutable value_type _emulator;
  };  // end const_iterator

  /**
   * Get const_iterator access to the parameters stored in this RBParameters object.
   * This gives serialized access to all the individual Real values in the
   * nested vector<vector<Real>>.
   * Use this in a for loop with the following syntax, for example:
   * for (const auto &[key,val] : as_range(rb_parameters.begin_serialized(), rb_parameters.end_serialized())
   */
  const_iterator begin_serialized() const;
  const_iterator end_serialized() const;

  /**
   * Get const_iterator access to the extra parameters stored in this RBParameters object.
   * This gives serialized access to all the individual Real values in the
   * nested vector<vector<Real>>.
   * Use this in a for loop with the following syntax, for example:
   * for (const auto &[key,val] : as_range(rb_parameters.begin_serialized_extra(), rb_parameters.end_serialized_extra())
   */
  const_iterator begin_serialized_extra() const;
  const_iterator end_serialized_extra() const;

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
   * Get the value of the specified parameter, throw an error if it does not exist.
   * Here we assume that there is only one sample, throw an error otherwise.
   * The Real-returning version also throws an error if the parameter exists
   * but contains multiple values.
   */
  Real get_value(const std::string & param_name) const;
  const RBParameter & get_vector_value(const std::string & param_name) const;

  /**
   * Get the value of the specified parameter, returning the provided
   * default value if it does not exist.
   * If the value does exist, we assume that there is only one sample,
   * and throw an error otherwise.
   * The Real-returning version also throws an error if the parameter exists
   * but contains multiple values.
   */
  Real get_value(const std::string & param_name, const Real & default_val) const;
  const RBParameter & get_vector_value(const std::string & param_name, const RBParameter & default_val) const;

  /**
   * Get the value of the specified parameter at the specified sample,
   * throwing an error if it does not exist.
   * The Real-returning version throws an error if the parameter exists
   * but contains multiple values.
   */
  Real get_sample_value(const std::string & param_name, std::size_t sample_idx) const;
  const RBParameter & get_sample_vector_value(const std::string & param_name, std::size_t sample_idx) const;

  /**
   * Get the value of the specified parameter at the specified sample,
   * returning the provided default value if either the parameter is
   * not defined or the sample is invalid.
   * The Real-returning version throws an error if the parameter exists
   * but contains multiple values.
   */
  Real get_sample_value(const std::string & param_name, std::size_t sample_idx, const Real & default_val) const;
  const RBParameter & get_sample_vector_value(const std::string & param_name, std::size_t sample_idx, const RBParameter & default_val) const;

  /**
   * Set the value of the specified parameter. If param_name
   * doesn't already exist, it is added to the RBParameters object.
   * For backwards compatibility, calling this function sets up
   * "param_name" to be a single-entry vector with "value" as the
   * only entry.
   */
  void set_value(const std::string & param_name, Real value);
  void set_value(const std::string & param_name, const RBParameter & value);

  /**
   * Set the value of the specified parameter at the specified sample
   * index. The sample index can refer to, e.g., load or time steps.
   */
  void set_value(const std::string & param_name, std::size_t index, Real value);
  void set_value(const std::string & param_name, std::size_t index, const RBParameter & value);

  /**
   * Similar to set_value(name, index, value) but instead of specifying a particular
   * index, just appends one more. Calling push_back_value() many times is more efficient
   * than calling set_value(name, index, value) many times because it takes advantage
   * of the std::vector's size-doubling t reduce allocations.
   */
  void push_back_value(const std::string & param_name, Real value);
  void push_back_value(const std::string & param_name, const RBParameter & value);

  /**
   * Same as push_back_value(), but for "extra" parameters.
   */
  void push_back_extra_value(const std::string & param_name, Real value);
  void push_back_extra_value(const std::string & param_name, const RBParameter & value);

  /**
   * Get the value of the specified extra parameter, throwing an error
   * if it does not exist.
   */
  Real get_extra_value(const std::string & param_name) const;
  const RBParameter & get_extra_vector_value(const std::string & param_name) const;

  /**
   * Get the value of the specified extra parameter, returning the
   * provided default value if it does not exist.
   */
  Real get_extra_value(const std::string & param_name, const Real & default_val) const;
  const RBParameter & get_extra_vector_value(const std::string & param_name, const RBParameter & default_val) const;

  /**
   * Get the value of the specified "extra" parameter at the specified sample index,
   * throwing an error if it does not exist.
   */
  Real get_extra_sample_value(const std::string & param_name, std::size_t sample_idx) const;
  const RBParameter & get_extra_sample_vector_value(const std::string & param_name, std::size_t sample_idx) const;

  /**
   * Get the value of the specified extra parameter at the specified sample index,
   * returning the provided default value if either the parameter is
   * not defined or the sample index is invalid.
   */
  Real get_extra_sample_value(const std::string & param_name, std::size_t sample_idx, const Real & default_val) const;
  const RBParameter & get_extra_sample_vector_value(const std::string & param_name, std::size_t sample_idx, const RBParameter & default_val) const;

  /**
   * Set the value of the specified extra parameter. If param_name
   * doesn't already exist, it is added to the extra parameters.
   */
  void set_extra_value(const std::string & param_name, Real value);
  void set_extra_value(const std::string & param_name, const RBParameter & value);

  /**
   * Set the value of the specified extra parameter at the specified sample
   * index. The sample index can refer to, e.g., load or time steps.
   */
  void set_extra_value(const std::string & param_name, std::size_t index, Real value);
  void set_extra_value(const std::string & param_name, std::size_t index, const RBParameter & value);

  /**
   * Get the number of parameters that have been added.
   */
  unsigned int n_parameters() const;

  /**
   * Set the number of samples this RBParameters object is intended to
   * represent, in the case that there are no actual parameters stored
   * on it. Note: this value will only be used in the no-parameters
   * case; if there are actual parameters specified in this class, the
   * number set via this API is ignored. All parameters stored within
   * the RBParameters object must have n_samples() samples.
   */
  void set_n_samples(unsigned int n_samples);

  /**
   * Returns the number of samples stored for all parameters. For
   * simplicity, we require all parameters to store the same number of
   * "samples" ("sample" here may refer to, e.g., time step or load step) and in
   * debug mode we actually verify that is the case.
   */
  unsigned int n_samples() const;

  /**
   * \return a set with the names of the parameters.
   *
   * Note that instead of creating a new set of strings here, it's
   * better to iterate over the RBParameters object directly, using
   * the iterators from the begin()/end() functions. The .first
   * will provide the parameter name, .second can be ignored as needed.
   */
  std::set<std::string> get_parameter_names() const;

  /**
   * \return a set with the names of the extra parameters.
   *
   * Note that instead of creating a new set of strings here, it's
   * better to iterate over the RBParameters extra object directly, using
   * the iterators from the extra_begin()/extra_end() functions. The .first
   * will provide the parameter name, .second can be ignored as needed.
   */
  std::set<std::string> get_extra_parameter_names() const;

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
   * Two RBParameters are equal if they have the same _parameters map.
   */
  bool operator== (const RBParameters & rhs) const;

  /**
   * \returns !(*this == rhs).
   */
  bool operator!= (const RBParameters & rhs) const;

  /**
   * Append "rhs" to "*this".  Both RBParameters objects must have the
   * same n_samples(), otherwise an error is thrown. If some of the
   * parameter names overlap, then the values from rhs overwrite
   * *this. Both parameters and "extra" parameters are appended.
   */
  RBParameters & operator+= (const RBParameters & rhs);

  /**
   * Get a string that specifies the contents of this RBParameters object.
   * \p precision specifies the number of digits of precision we use
   * in scientific notation in the string.
   * \p max_values is the max number of values to print out if the parameter
   * is vector-valued. Set to negative value to print all.
   */
  std::string get_string(unsigned precision=6, int max_values=5) const;

  /**
   * Print the parameters.
   */
  void print(unsigned precision=6, int max_values=5) const;

private:

  /**
   * Helper function for the 3-parameter versions of set_value() and
   * set_extra_value().
   */
  void set_value_helper(std::map<std::string, std::vector<RBParameter>> & map,
                        const std::string & param_name,
                        const std::size_t index,
                        RBParameter value);

  /**
   * The number of samples represented by this RBParameters object, in
   * the case where there are no parameters actually stored on it. If
   * there are parameters stored on this RBParameters object, then the
   * n_samples() API returns that number of samples instead.
   */
  unsigned int _n_samples;

  /**
   * Actual parameter values (in std::vector<RBParameter> form) across a vector of samples.
   * Each vector is indexed by a name.
   * Note that the number of samples in the outer vector should be the same
   * across all parameters, however, this is not necessary for the inner
   * "value-vector".
   */
  std::map<std::string, std::vector<RBParameter>> _parameters;

  /**
   * Extra parameter vectors not used for RB training.
   * Each vector is indexed by a name.
   */
  std::map<std::string, std::vector<RBParameter>> _extra_parameters;
};

} // namespace libMesh


#endif // LIBMESH_RB_PARAMETERS_H

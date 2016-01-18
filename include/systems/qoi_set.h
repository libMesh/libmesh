// The libMesh Finite Element Library.
// Copyright (C) 2002-2016 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_QOI_SET_H
#define LIBMESH_QOI_SET_H


// Local Includes -----------------------------------
#include "libmesh/libmesh_common.h"

// C++ Includes   -----------------------------------
#include <vector>

namespace libMesh
{

// Forward Declarations -----------------------------
class System;

/**
 * Data structure for specifying which Quantities of Interest
 * should be calculated in an adjoint or a parameter sensitivity
 * calculation.
 */
class QoISet
{
public:
  class iterator
  {
  public:
    iterator(unsigned int i, const std::vector<bool> & v) : _i(i), _vecbool(v)
    {
      while (_i < _vecbool.size() && !_vecbool[_i])
        _i++;
    }

    unsigned int operator*() const { return _i; }

    iterator & operator++()
    {
      do {
        _i++;
      } while (_i < _vecbool.size() && !_vecbool[_i]);
      return *this;
    }

    iterator operator++(int) {
      iterator it = *this;
      ++(*this);
      return it;
    }

    bool operator==(const iterator & other) const {
      libmesh_assert_equal_to (&_vecbool, &other._vecbool);
      return _i == other._i;
    }

    bool operator!=(const iterator & other) const {
      libmesh_assert_equal_to (&_vecbool, &other._vecbool);
      return _i != other._i;
    }

  private:

    unsigned int _i;

    const std::vector<bool> & _vecbool;
  };

  /**
   * Empty constructor: "calculate all QoIs in the System"
   *
   * No further changes to this special QoISet should be made;
   * it doesn't even know how many QoIs your system has, it
   * just knows to instruct a function to use all of them.
   */
  QoISet() : _indices(), _weights() {}

  /**
   * Default constructor: "calculate all QoIs in the System",
   * "give every QoI weight 1.0"
   */
  explicit
  QoISet(const System & sys);

  /**
   * Constructor-from-vector-of-bool: "calculate the QoIs for which
   * \p indices[q] is true"
   */
  explicit
  QoISet(const std::vector<bool> & indices) :
    _indices(indices), _weights() {}

  /**
   * Constructor-from-vector: "calculate the listed QoIs", "give every
   * QoI weight 1.0"
   */
  explicit
  QoISet(const std::vector<unsigned int> & indices);

  /**
   * Resets to "calculate all QoIs, give every QoI weight 1.0"
   */
  void clear() { _indices.clear(); _weights.clear(); }

  /**
   * Returns the number of QoIs that would be computed for the
   * System \p sys
   */
  unsigned int size(const System & sys) const;

  /**
   * Add this indices to the set to be calculated
   */
  void add_indices(const std::vector<unsigned int> & indices);

  /**
   * Add this index to the set to be calculated
   */
  void add_index(unsigned int);

  /**
   * Remove these indices from the set to be calculated
   */
  void remove_indices(const std::vector<unsigned int> & indices);

  /**
   * Remove this index from the set to be calculated
   */
  void remove_index(unsigned int);

  /**
   * Set the weight for this index
   */
  void set_weight(unsigned int, Real);

  /**
   * Get the weight for this index (default 1.0)
   */
  Real weight(unsigned int) const;

  /**
   * Return whether or not this index is in the set to be calculated
   */
  bool has_index(unsigned int) const;

  /**
   * Return an iterator pointing to the first index in the set
   */
  iterator begin() const { return iterator(0, _indices); }

private:
  /**
   * Interpret _indices.empty() to mean "calculate all indices"
   */
  std::vector<bool> _indices;

  /**
   * Interpret _weights.size() <= i to mean "weight i = 1.0"
   */
  std::vector<Real> _weights;
};



// ------------------------------------------------------------
// QoISet inline methods



inline
QoISet::QoISet(const std::vector<unsigned int> & indices) :
  _indices(), _weights()
{
  this->add_indices(indices);
}



inline
void QoISet::add_index(unsigned int i)
{
  if (i >= _indices.size())
    _indices.resize(i+1, true);
  _indices[i] = true;
}



inline
void QoISet::remove_index(unsigned int i)
{
  if (i >= _indices.size())
    _indices.resize(i+1, true);
  _indices[i] = false;
}



inline
bool QoISet::has_index(unsigned int i) const
{
  return (_indices.size() <= i || _indices[i]);
}



inline
void QoISet::set_weight(unsigned int i, Real w)
{
  if (_weights.size() <= i)
    _weights.resize(i+1, 1.0);

  _weights[i] = w;
}



inline
Real QoISet::weight(unsigned int i) const
{
  if (_weights.size() <= i)
    return 1.0;
  return _weights[i];
}

} // namespace libMesh

#endif // LIBMESH_QOI_SET_H

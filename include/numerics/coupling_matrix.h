// The libMesh Finite Element Library.
// Copyright (C) 2002-2014 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_COUPLING_MATRIX_H
#define LIBMESH_COUPLING_MATRIX_H

// Local Includes
#include "libmesh/libmesh_common.h"

// C++ includes
#include <algorithm>
#include <limits>
#include <utility> // std::pair
#include <vector>

namespace libMesh
{

// Forward declarations
class ConstCouplingAccessor;
class CouplingAccessor;

/**
 * This class defines a coupling matrix.  A coupling
 * matrix is simply a matrix of ones and zeros describing
 * how different components in a system couple with each
 * other.  A coupling matrix is necessarily square but not
 * necessarily symmetric.
 */
class CouplingMatrix
{
public:

  /**
   * Constructor.
   */
  explicit
  CouplingMatrix (const unsigned int n=0);

  /**
   * @returns the (i,j) entry of the matrix.
   */
  bool operator() (const unsigned int i,
                   const unsigned int j) const;

  /**
   * @returns the (i,j) entry of the matrix as
   * a smart-reference.
   */
  CouplingAccessor operator() (const unsigned int i,
                               const unsigned int j);

  /**
   * @returns the size of the matrix, i.e. N for an
   * NxN matrix.
   */
  unsigned int size() const;

  /**
   * Resizes the matrix and initializes
   * all entries to be 0.
   */
  void resize(const unsigned int n);

  /**
   * Clears the matrix.
   */
  void clear();

  /**
   * @returns true if the matrix is empty.
   */
  bool empty() const;

private:

  friend class ConstCouplingAccessor;
  friend class CouplingAccessor;

  /**
   * Coupling matrices are typically either full or very sparse, and
   * all values are only zero or one.
   *
   * We store non-zeros as ranges: the first entry of each range pair
   * is the location of the first non-zero, and the second is the
   * location of the last subsequent non-zero (*not* the next
   * subsequent zero; we drop empty ranges).
   * 
   * We store locations (i,j) as long integers i*_size+j
   */
  typedef std::pair<std::size_t, std::size_t> range_type;
  typedef std::vector<range_type> rc_type;
  rc_type _ranges;

  /**
   * The size of the matrix.
   */
  unsigned int _size;
};



/**
 * This accessor class allows simple access to CouplingMatrix values.
 */
class ConstCouplingAccessor
{
public:
  ConstCouplingAccessor(std::size_t loc_in,
                        const CouplingMatrix& mat_in) :
    _location(loc_in), _mat(mat_in)
  {
    libmesh_assert_less(_location, _mat.size() * _mat.size());
  }

  operator bool() const {
    const std::size_t max_size = std::numeric_limits<std::size_t>::max();

    // Find the range that might contain i,j
    // lower_bound isn't *quite* what we want
    CouplingMatrix::rc_type::const_iterator lb = std::upper_bound
      (_mat._ranges.begin(), _mat._ranges.end(),
       std::make_pair(_location, max_size));
    if (lb!=_mat._ranges.begin())
      --lb;
    else
      lb=_mat._ranges.end();

    // If no range might contain i,j then it's 0
    if (lb == _mat._ranges.end())
      return false;

    const std::size_t firstloc = lb->first;
    const std::size_t lastloc  = lb->second;
    libmesh_assert_less_equal(firstloc, lastloc);
    libmesh_assert_less_equal(firstloc, _location);

#ifdef DEBUG
    CouplingMatrix::rc_type::const_iterator next = lb;
    next++;
    if (next != _mat._ranges.end())
      {
        // Ranges should be sorted and should not touch
        libmesh_assert_greater(next->first, lastloc+1);
      }
#endif

    return (lastloc >= _location);
  }

protected:

  std::size_t _location;
  const CouplingMatrix& _mat;
};


/**
 * This accessor class allows simple setting of CouplingMatrix values.
 */
class CouplingAccessor : public ConstCouplingAccessor
{
public:
  CouplingAccessor(std::size_t loc_in,
                   CouplingMatrix& mat_in) :
    ConstCouplingAccessor(loc_in, mat_in), _my_mat(mat_in) {}

  template <typename T>
  CouplingAccessor& operator = (T new_value)
  {
    // For backwards compatibility we take integer arguments,
    // but coupling matrix entries are really all zero or one.
    const bool as_bool = new_value;
    libmesh_assert_equal_to(new_value, as_bool);

    *this = as_bool;
    return *this;
  }

  CouplingAccessor& operator = (bool new_value)
  {
    const std::size_t max_size = std::numeric_limits<std::size_t>::max();

    // Find the range that might contain i,j
    // lower_bound isn't *quite* what we want
    CouplingMatrix::rc_type::iterator lb =
      std::upper_bound (_my_mat._ranges.begin(), _my_mat._ranges.end(),
                        std::make_pair(_location, max_size));
    if (lb!=_my_mat._ranges.begin())
      --lb;
    else
      lb=_my_mat._ranges.end();

    // If no range might contain i,j then we might need to make a new
    // one.
    if (lb == _my_mat._ranges.end())
      {
        if (new_value == true)
          _my_mat._ranges.insert(_my_mat._ranges.begin(),
                                 std::make_pair(_location, _location));
      }
    else
      {
        const std::size_t firstloc = lb->first;
        const std::size_t lastloc  = lb->second;
        libmesh_assert_less_equal(firstloc, lastloc);
        libmesh_assert_less_equal(firstloc, _location);

#ifdef DEBUG
        CouplingMatrix::rc_type::const_iterator next = lb;
        next++;
        if (next != _my_mat._ranges.end())
          {
            // Ranges should be sorted and should not touch
            libmesh_assert_greater(next->first, lastloc+1);
          }
#endif

        // If we're in this range, we might need to shorten or remove
        // or split it
        if (new_value == false)
          {
            if (_location == firstloc)
              {
                if (_location == lastloc)
                  {
                    _my_mat._ranges.erase(lb);
                  }
                else
                  {
                    libmesh_assert_less (lb->first, lastloc);
                    lb->first++;
                  }
              }
            else if (_location == lastloc)
              {
                libmesh_assert_less (firstloc, lb->second);

                lb->second--;
              }
            else if (_location < lastloc)
              {
                libmesh_assert_less_equal(_location+1, lastloc);

                lb->first = _location+1;

                libmesh_assert_less_equal(firstloc, _location-1);

                _my_mat._ranges.insert
                  (lb, std::make_pair(firstloc, _location-1));
              }
          }

        // If we're not in this range, we might need to extend it or
        // join it with its neighbor or add a new one.
        else // new_value == true
          {
            CouplingMatrix::rc_type::iterator next = lb;
            next++;
            const std::size_t nextloc =
              (next == _my_mat._ranges.end()) ?
              std::numeric_limits<std::size_t>::max() :
              next->first;

            // Ranges should be sorted and should not touch
            libmesh_assert_greater(nextloc, lastloc+1);

            if (_location > lastloc)
              {
                if (_location == lastloc + 1)
                  {
                    if (_location == nextloc - 1)
                      {
                        next->first = firstloc;
                        _my_mat._ranges.erase(lb);
                      }
                    else
                      lb->second++;
                  }
                else
                  {
                    if (_location == nextloc - 1)
                      next->first--;
                    else
                      _my_mat._ranges.insert
                        (next, std::make_pair(_location, _location));
                  }
              }
          }
      }

    return *this;
  }

private:

  CouplingMatrix& _my_mat;
};

//--------------------------------------------------
// CouplingMatrix inline methods
inline
CouplingMatrix::CouplingMatrix (const unsigned int n) :
  _ranges(), _size(n)
{
  this->resize(n);
}



inline
bool CouplingMatrix::operator() (const unsigned int i,
                                 const unsigned int j) const
{
  libmesh_assert_less (i, _size);
  libmesh_assert_less (j, _size);

  const std::size_t location = std::size_t(i)*_size + j;

  return bool(ConstCouplingAccessor(location, *this));
}




inline
CouplingAccessor CouplingMatrix::operator() (const unsigned int i,
                                             const unsigned int j)
{
  const std::size_t location = std::size_t(i)*_size + j;

  return CouplingAccessor(location, *this);
}



inline
unsigned int CouplingMatrix::size() const
{
  return _size;
}



inline
void CouplingMatrix::resize(const unsigned int n)
{
  _size = n;

  _ranges.clear();
}



inline
void CouplingMatrix::clear()
{
  this->resize(0);
}



inline
bool CouplingMatrix::empty() const
{
  return (_size == 0);
}


} // namespace libMesh


#endif // LIBMESH_COUPLING_MATRIX_H

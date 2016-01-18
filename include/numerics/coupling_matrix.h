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

class ConstCouplingRow;

class ConstCouplingRowConstIterator;

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
  friend class ConstCouplingRow;
  friend class ConstCouplingRowConstIterator;

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
                        const CouplingMatrix & mat_in) :
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

    const std::size_t lastloc  = lb->second;

#ifdef DEBUG
    const std::size_t firstloc = lb->first;
    libmesh_assert_less_equal(firstloc, lastloc);
    libmesh_assert_less_equal(firstloc, _location);

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
  const CouplingMatrix & _mat;
};


/**
 * This accessor class allows simple setting of CouplingMatrix values.
 */
class CouplingAccessor : public ConstCouplingAccessor
{
public:
  CouplingAccessor(std::size_t loc_in,
                   CouplingMatrix & mat_in) :
    ConstCouplingAccessor(loc_in, mat_in), _my_mat(mat_in) {}

  template <typename T>
  CouplingAccessor & operator = (T new_value)
  {
    // For backwards compatibility we take integer arguments,
    // but coupling matrix entries are really all zero or one.
    const bool as_bool = new_value;
    libmesh_assert_equal_to(new_value, as_bool);

    *this = as_bool;
    return *this;
  }

  CouplingAccessor & operator = (bool new_value)
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

  CouplingMatrix & _my_mat;
};



/**
 * This proxy class acts like a container of indices from a single
 * coupling row
 */
class ConstCouplingRow
{
public:
  ConstCouplingRow(unsigned int row_in,
                   const CouplingMatrix & mat_in) :
    _row_i(row_in), _mat(mat_in)
  {
    libmesh_assert_less(_row_i, _mat.size());

    // Location for i,N
    _begin_location = _row_i*_mat.size()+_mat.size()-1;

    const std::size_t max_size = std::numeric_limits<std::size_t>::max();

    // Find the range that might contain i,N
    // lower_bound isn't *quite* what we want
    _begin_it = std::upper_bound
      (_mat._ranges.begin(), _mat._ranges.end(),
       std::make_pair(_begin_location, max_size));
    if (_begin_it !=_mat._ranges.begin())
      --_begin_it;
    else
      _begin_it=_mat._ranges.end();

    // If that range doesn't exist then we're an empty row
    if (_begin_it == _mat._ranges.end())
      _begin_location = max_size;
    else
      {
        const std::size_t lastloc  = _begin_it->second;
#ifdef DEBUG
        const std::size_t firstloc = _begin_it->first;
        libmesh_assert_less_equal(firstloc, lastloc);
#endif

        // If that range ends before i,0 then we're an empty row
        std::size_t zero_location = _row_i*_mat.size();
        if (zero_location > lastloc)
          {
            _begin_location = max_size;
            _begin_it = _mat._ranges.end();
          }
        else
          // We have *some* entry(s) in this row, we just need to find
          // the earliest
          {
            while (_begin_it != _mat._ranges.begin())
              {
                CouplingMatrix::rc_type::const_iterator prev =
                  _begin_it;
                --prev;

                if (prev->second < zero_location)
                  break;

                _begin_it = prev;
              }
            if (_begin_it->first < zero_location)
              _begin_location = zero_location;
            else
              _begin_location = _begin_it->first;
          }
      }
  }

  /*
   * A forward iterator type for looping over indices in this row
   */
  typedef ConstCouplingRowConstIterator const_iterator;

  /*
   * An iterator to the first index in this row, or to end() for an
   * empty row
   */
  const_iterator begin() const;

  /*
   * An iterator representing past-the-end of this row
   */
  const_iterator end() const;

  bool operator== (const ConstCouplingRow & other) const
  {
    // Thinking that rows from different matrix objects are equal is
    // not even wrong
    libmesh_assert(&_mat == &other._mat);

    return ((_begin_location == other._begin_location) &&
            (_begin_it == other._begin_it));
  }

  bool operator!= (const ConstCouplingRow & other) const
  {
    return !(*this == other);
  }
protected:

  friend class ConstCouplingRowConstIterator;

  unsigned int _row_i;
  const CouplingMatrix & _mat;

  // The location (i*size+j) corresponding to the first entry in this
  // row, or numeric_limits<size_t>::max() for an empty row.
  std::size_t _begin_location;

  // Iterator to the range containing the first row element, or
  // _row._mat._ranges.end() if no CouplingMatrix values are true for
  // this row
  CouplingMatrix::rc_type::const_iterator _begin_it;
};



class ConstCouplingRowConstIterator
{
public:
  ConstCouplingRowConstIterator (const ConstCouplingRow & row_in,
                                 std::size_t loc_in,
                                 CouplingMatrix::rc_type::const_iterator it_in) :
    _location(loc_in),
    _row(row_in),
    _it(it_in)
  {
#ifndef NDEBUG
    if (_it != _row._mat._ranges.end())
      {
        libmesh_assert_less_equal(_it->first, _location);
        libmesh_assert_less_equal(_location, _it->second);
      }
    else
      {
        libmesh_assert_equal_to
          (_location, std::numeric_limits<size_t>::max());
      }
#endif
  }

  unsigned int operator* ()
  {
    libmesh_assert_not_equal_to
      (_location, std::numeric_limits<std::size_t>::max());
    return _location % _row._mat.size();
  }

  ConstCouplingRowConstIterator & operator++ ()
  {
    libmesh_assert_not_equal_to
      (_location, std::numeric_limits<std::size_t>::max());

    if (_location == _it->second)
      {
        ++_it;

        // Are we past the end of the matrix?
        if (_it == _row._mat._ranges.end())
          _location = std::numeric_limits<std::size_t>::max();
        else
          {
            _location = _it->first;
            // Are we past the end of the row?
            if (_location >= _row._mat.size()*(_row._row_i+1))
              {
                _location = std::numeric_limits<std::size_t>::max();
                _it = _row._mat._ranges.end();
              }
          }
      }
    else
      ++_location;

    return *this;
  }

  bool operator== (const ConstCouplingRowConstIterator & other) const
  {
    // Thinking that iterators from different row objects are equal
    // is not even wrong
    libmesh_assert(_row == other._row);

    return ((_location == other._location) &&
            (_it == other._it));
  }

  bool operator!= (const ConstCouplingRowConstIterator & other) const
  {
    return !(*this == other);
  }

private:
  // The location (i*size+j) corresponding to this iterator, or
  // numeric_limits<size_t>::max() to signify end()
  std::size_t _location;
  const ConstCouplingRow & _row;
  // The range containing this iterator location, or
  // _row._mat._ranges.end() if no range contains this location
  CouplingMatrix::rc_type::const_iterator _it;
};



//--------------------------------------------------
// ConstCouplingRow inline methods
inline
ConstCouplingRow::const_iterator ConstCouplingRow::begin() const {
  return const_iterator (*this, _begin_location, _begin_it);
}

inline
ConstCouplingRow::const_iterator ConstCouplingRow::end() const {
  return const_iterator
    (*this, std::numeric_limits<std::size_t>::max(),
     _mat._ranges.end());
}



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

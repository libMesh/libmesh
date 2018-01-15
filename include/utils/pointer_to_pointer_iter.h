// The libMesh Finite Element Library.
// Copyright (C) 2002-2018 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_PPITER_H
#define LIBMESH_PPITER_H


namespace libMesh
{

/**
 * The \p PointerToPointerIter templated class is intended to wrap
 * pointer-to-pointer iterators in an interface which works more like
 * a standard iterator, by returning a value rather than a pointer.
 *
 * \author  Roy H. Stogner
 */

template <typename T>
class PointerToPointerIter
{
public:
  PointerToPointerIter (T * const * it) : _it(it) {}

  T & operator* () const { return **_it; }

  const PointerToPointerIter & operator++ ()
  {
    ++_it;
    return *this;
  }

  PointerToPointerIter operator++ (int)
  {
    PointerToPointerIter returnval(*this);
    ++_it;
    return returnval;
  }

  bool operator== (const PointerToPointerIter & j) const
  {
    return ( _it == j._it );
  }

  bool operator!= (const PointerToPointerIter & j) const
  {
    return !(*this == j);
  }

private:
  T * const * _it;
};

} // namespace libMesh

#endif // LIBMESH_PPITER_H

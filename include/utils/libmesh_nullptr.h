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
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA

#ifndef LIBMESH_NULLPTR_H
#define LIBMESH_NULLPTR_H

#include "libmesh/libmesh_config.h"

#if defined(LIBMESH_HAVE_CXX11_NULLPTR)

// Use the built-in nullptr keyword
#define libmesh_nullptr nullptr

#elif defined(LIBMESH_HAVE_CXX11_NULLPTR_WORKAROUND)

/**
 * A C++03-compatible nullptr type implemented as a library solution
 * rather than a language keyword.  Originally due to Scott Meyers,
 * Effective C++, 2nd Edition, Item 25.  The main benefits of using
 * this type over NULL are
 * 1.) Expressiveness:
 * 2.) Disambiguating overloads, e.g.
 *     void f(int);
 *     void f(double *);
 *     f(NULL); // ambiguous
 *     f(nullptr); // unambiguous - calls the version of f taking a pointer
 *
 * We are specifically *not* attempting to redefine the 'nullptr'
 * keyword itself -- that can lead to issues when C++11 is disabled,
 * yet the compiler still has 'nullptr' defined.  Therefore, client
 * code that wishes to be more expressive/correct should use
 * 'libmesh_nullptr' for maximum backwards compatibility.
 *
 * Note: the leading const applies to the following nullptr
 * instance.
 *
 * Code originally coped from:
 * https://en.wikibooks.org/wiki/More_C%2B%2B_Idioms/nullptr#Solution_and_Sample_Code
 */
const class libmesh_nullptr_t
{
public:
  // convertible to any type of null non-member pointer...
  template<class T>
  inline operator T * () const { return 0; }

  // or any type of null member pointer...
  template<class C, class T>
  inline operator T C::*() const { return 0; }

private:
  // Can't take address of nullptr
  void operator & () const;
} libmesh_nullptr = {};

#else

// No nullptr and workaround doesn't compile.  This should not be
// common, but in this case, We have to fall back on the C++03
// definition of NULL.
#define libmesh_nullptr NULL

#endif // LIBMESH_HAVE_CXX11_NULLPTR

#endif // LIBMESH_NULLPTR_H

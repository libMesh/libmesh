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

#ifndef LIBMESH_ELEM_HASH_H
#define LIBMESH_ELEM_HASH_H

#include "elem.h"
#include LIBMESH_INCLUDE_UNORDERED_MAP
#include LIBMESH_INCLUDE_UNORDERED_SET

// This header defines some typedefs that are useful for working with
// "unordered" containers of Elem * that use Elem::key() as a hash
// function.
namespace libMesh
{

// The ElemHashUtils struct defines functions used for the "Hash" and
// "Pred" template arguments of the various "unordered" containers,
// e.g.
// template <class Key,                         // unordered_multiset::key_type/value_type
//           class Hash = hash<Key>,            // unordered_multiset::hasher
//           class Pred = equal_to<Key>,        // unordered_multiset::key_equal
//           class Alloc = allocator<Key>       // unordered_multiset::allocator_type
//           > class unordered_multiset;
struct ElemHashUtils
{
public:
  // The "Hash" template argument.  A custom hash functor that can be
  // used with the "unordered" container types.  Simply returns
  // elem->key() as the hash.
  inline
  std::size_t operator()(const Elem * elem) const
  {
    return cast_int<std::size_t>(elem->key());
  }

  // The "Pred" template argument.  A binary predicate that takes two
  // arguments of the same type as the elements and returns a bool.
  // We need to specify this in order to use the unordered_multiset,
  // otherwise it just uses std::equal_to to compare two pointers...
  inline
  bool operator()(const Elem * lhs, const Elem * rhs) const
  {
    return lhs->key() == rhs->key();
  }
};

// A convenient type for working with unordered_multiset<Elem *>
typedef LIBMESH_BEST_UNORDERED_MULTISET<Elem *, ElemHashUtils, ElemHashUtils> unordered_multiset_elem;

}

#endif

// The libMesh Finite Element Library.
// Copyright (C) 2002-2019 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#ifndef LIBMESH_MESH_ITERATORS_IMPL_H
#define LIBMESH_MESH_ITERATORS_IMPL_H

// C++ includes

// Local includes
#include "libmesh/distributed_mesh.h"
#include "libmesh/elem.h"
#include "libmesh/replicated_mesh.h"

#include "libmesh/ignore_warnings.h" // Ignore warnings about variadic macros

namespace libMesh
{

// This file contains the implementation of all the different iterator
// functions for the mesh class.

// This macro generates four iterator accessor function definitions
// (const/non-const and begin/end) for both Replicated and DistributedMesh
// given the Predicate PRED, which may be passed an arbitrary number
// of arguments.
#define INSTANTIATE_ELEM_ACCESSORS(RealType, FUNC_PREFIX, PRED, FUNC_ARG, ...) \
  template <>                                                            \
  typename ReplicatedMeshTempl<RealType>::element_iterator              \
  ReplicatedMeshTempl<RealType>::FUNC_PREFIX##_begin (FUNC_ARG)                        \
  {                                                                     \
    return element_iterator(_elements.begin(), _elements.end(), Predicates::PRED<elem_iterator_imp>(__VA_ARGS__)); \
  }                                                                     \
  template <>                                                           \
  typename ReplicatedMeshTempl<RealType>::const_element_iterator                                \
  ReplicatedMeshTempl<RealType>::FUNC_PREFIX##_begin (FUNC_ARG) const   \
  {                                                                     \
    return const_element_iterator(_elements.begin(), _elements.end(), Predicates::PRED<const_elem_iterator_imp>(__VA_ARGS__)); \
  }                                                                     \
  template <>                                                           \
  typename ReplicatedMeshTempl<RealType>::element_iterator                                      \
  ReplicatedMeshTempl<RealType>::FUNC_PREFIX##_end (FUNC_ARG)           \
  {                                                                     \
    return element_iterator(_elements.end(), _elements.end(), Predicates::PRED<elem_iterator_imp>(__VA_ARGS__)); \
  }                                                                     \
  template <>                                                           \
  typename ReplicatedMeshTempl<RealType>::const_element_iterator                                \
  ReplicatedMeshTempl<RealType>::FUNC_PREFIX##_end (FUNC_ARG) const     \
  {                                                                     \
    return const_element_iterator(_elements.end(), _elements.end(), Predicates::PRED<const_elem_iterator_imp>(__VA_ARGS__)); \
  }                                                                     \
  template <>                                                           \
  typename DistributedMeshTempl<RealType>::element_iterator                                     \
  DistributedMeshTempl<RealType>::FUNC_PREFIX##_begin (FUNC_ARG)                       \
  {                                                                     \
    return element_iterator(_elements.begin(), _elements.end(), Predicates::PRED<elem_iterator_imp>(__VA_ARGS__)); \
  }                                                                     \
  template <>                                                           \
  typename DistributedMeshTempl<RealType>::const_element_iterator                               \
  DistributedMeshTempl<RealType>::FUNC_PREFIX##_begin (FUNC_ARG) const                 \
  {                                                                     \
    return const_element_iterator(_elements.begin(), _elements.end(), Predicates::PRED<const_elem_iterator_imp>(__VA_ARGS__)); \
  }                                                                     \
  template <>                                                           \
  typename DistributedMeshTempl<RealType>::element_iterator                                     \
  DistributedMeshTempl<RealType>::FUNC_PREFIX##_end (FUNC_ARG)                         \
  {                                                                     \
    return element_iterator(_elements.end(), _elements.end(), Predicates::PRED<elem_iterator_imp>(__VA_ARGS__)); \
  }                                                                     \
  template <>                                                           \
  typename DistributedMeshTempl<RealType>::const_element_iterator                               \
  DistributedMeshTempl<RealType>::FUNC_PREFIX##_end (FUNC_ARG) const                   \
  {                                                                     \
    return const_element_iterator(_elements.end(), _elements.end(), Predicates::PRED<const_elem_iterator_imp>(__VA_ARGS__)); \
  }



// This macro is similar to the one above except that it generates
// node iterator accessor functions.
#define INSTANTIATE_NODE_ACCESSORS(RealType, FUNC_PREFIX, PRED, FUNC_ARG, ...) \
  template <>                                                            \
  typename ReplicatedMeshTempl<RealType>::node_iterator                                         \
  ReplicatedMeshTempl<RealType>::FUNC_PREFIX##_begin (FUNC_ARG)                        \
  {                                                                     \
    return node_iterator(_nodes.begin(), _nodes.end(), Predicates::PRED<node_iterator_imp>(__VA_ARGS__)); \
  }                                                                     \
  template <>                                                           \
  typename ReplicatedMeshTempl<RealType>::const_node_iterator           \
  ReplicatedMeshTempl<RealType>::FUNC_PREFIX##_begin (FUNC_ARG) const   \
  {                                                                     \
    return const_node_iterator(_nodes.begin(), _nodes.end(), Predicates::PRED<const_node_iterator_imp>(__VA_ARGS__)); \
  }                                                                     \
  template <>                                                           \
  typename ReplicatedMeshTempl<RealType>::node_iterator                 \
  ReplicatedMeshTempl<RealType>::FUNC_PREFIX##_end (FUNC_ARG)           \
  {                                                                     \
    return node_iterator(_nodes.end(), _nodes.end(), Predicates::PRED<node_iterator_imp>(__VA_ARGS__)); \
  }                                                                     \
  template <>                                                           \
  typename ReplicatedMeshTempl<RealType>::const_node_iterator           \
  ReplicatedMeshTempl<RealType>::FUNC_PREFIX##_end (FUNC_ARG) const     \
  {                                                                     \
    return const_node_iterator(_nodes.end(), _nodes.end(), Predicates::PRED<const_node_iterator_imp>(__VA_ARGS__)); \
  }                                                                     \
  template <>                                                           \
  typename DistributedMeshTempl<RealType>::node_iterator                                        \
  DistributedMeshTempl<RealType>::FUNC_PREFIX##_begin (FUNC_ARG)                       \
  {                                                                     \
    return node_iterator(_nodes.begin(), _nodes.end(), Predicates::PRED<node_iterator_imp>(__VA_ARGS__)); \
  }                                                                     \
  template <>                                                           \
  typename DistributedMeshTempl<RealType>::const_node_iterator                                  \
  DistributedMeshTempl<RealType>::FUNC_PREFIX##_begin (FUNC_ARG) const                 \
  {                                                                     \
    return const_node_iterator(_nodes.begin(), _nodes.end(), Predicates::PRED<const_node_iterator_imp>(__VA_ARGS__)); \
  }                                                                     \
  template <>                                                           \
  typename DistributedMeshTempl<RealType>::node_iterator                                        \
  DistributedMeshTempl<RealType>::FUNC_PREFIX##_end (FUNC_ARG)                         \
  {                                                                     \
    return node_iterator(_nodes.end(), _nodes.end(), Predicates::PRED<node_iterator_imp>(__VA_ARGS__)); \
  }                                                                     \
  template <>                                                           \
  typename DistributedMeshTempl<RealType>::const_node_iterator                                  \
  DistributedMeshTempl<RealType>::FUNC_PREFIX##_end (FUNC_ARG) const                   \
  {                                                                     \
    return const_node_iterator(_nodes.end(), _nodes.end(), Predicates::PRED<const_node_iterator_imp>(__VA_ARGS__)); \
  }

// Use an empty preprocessor token to silence older compilers that
// still warn about empty macro arguments.
#define EMPTY

// Use a second macro layer to allow us to pass commas into a macro
// argument
#define LIBMESH_COMMA ,

} // namespace libMesh

#endif // LIBMESH_MESH_ITERATORS_IMPL_H

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



// C++ includes

// Local includes
#include "libmesh/serial_mesh.h"
#include "libmesh/parallel_mesh.h"
#include "libmesh/elem.h"

namespace libMesh
{

// This file contains the implementation of all the different iterator
// functions for the mesh class.

// This macro generates four iterator accessor function definitions
// (const/non-const and begin/end) for both Serial and ParallelMesh
// given the Predicate PRED, which may be passed an arbitrary number
// of arguments.
#define INSTANTIATE_ELEM_ACCESSORS(FUNC_PREFIX, PRED, FUNC_ARG, ...)    \
  SerialMesh::element_iterator                                          \
  SerialMesh::FUNC_PREFIX##_begin (FUNC_ARG)                            \
  {                                                                     \
    return element_iterator(_elements.begin(), _elements.end(), Predicates::PRED<elem_iterator_imp>(__VA_ARGS__)); \
  }                                                                     \
  SerialMesh::const_element_iterator                                    \
  SerialMesh::FUNC_PREFIX##_begin (FUNC_ARG) const                      \
  {                                                                     \
    return const_element_iterator(_elements.begin(), _elements.end(), Predicates::PRED<const_elem_iterator_imp>(__VA_ARGS__)); \
  }                                                                     \
  SerialMesh::element_iterator                                          \
  SerialMesh::FUNC_PREFIX##_end (FUNC_ARG)                              \
  {                                                                     \
    return element_iterator(_elements.end(), _elements.end(), Predicates::PRED<elem_iterator_imp>(__VA_ARGS__)); \
  }                                                                     \
  SerialMesh::const_element_iterator                                    \
  SerialMesh::FUNC_PREFIX##_end (FUNC_ARG) const                        \
  {                                                                     \
    return const_element_iterator(_elements.end(), _elements.end(), Predicates::PRED<const_elem_iterator_imp>(__VA_ARGS__)); \
  }                                                                     \
  ParallelMesh::element_iterator                                        \
  ParallelMesh::FUNC_PREFIX##_begin (FUNC_ARG)                          \
  {                                                                     \
    return element_iterator(_elements.begin(), _elements.end(), Predicates::PRED<elem_iterator_imp>(__VA_ARGS__)); \
  }                                                                     \
  ParallelMesh::const_element_iterator                                  \
  ParallelMesh::FUNC_PREFIX##_begin (FUNC_ARG) const                    \
  {                                                                     \
    return const_element_iterator(_elements.begin(), _elements.end(), Predicates::PRED<const_elem_iterator_imp>(__VA_ARGS__)); \
  }                                                                     \
  ParallelMesh::element_iterator                                        \
  ParallelMesh::FUNC_PREFIX##_end (FUNC_ARG)                            \
  {                                                                     \
    return element_iterator(_elements.end(), _elements.end(), Predicates::PRED<elem_iterator_imp>(__VA_ARGS__)); \
  }                                                                     \
  ParallelMesh::const_element_iterator                                  \
  ParallelMesh::FUNC_PREFIX##_end (FUNC_ARG) const                      \
  {                                                                     \
    return const_element_iterator(_elements.end(), _elements.end(), Predicates::PRED<const_elem_iterator_imp>(__VA_ARGS__)); \
  }



// This macro is similar to the one above except that it generates
// node iterator accessor functions.
#define INSTANTIATE_NODE_ACCESSORS(FUNC_PREFIX, PRED, FUNC_ARG, ...)    \
  SerialMesh::node_iterator                                             \
  SerialMesh::FUNC_PREFIX##_begin (FUNC_ARG)                            \
  {                                                                     \
    return node_iterator(_nodes.begin(), _nodes.end(), Predicates::PRED<node_iterator_imp>(__VA_ARGS__)); \
  }                                                                     \
  SerialMesh::const_node_iterator                                       \
  SerialMesh::FUNC_PREFIX##_begin (FUNC_ARG) const                      \
  {                                                                     \
    return const_node_iterator(_nodes.begin(), _nodes.end(), Predicates::PRED<const_node_iterator_imp>(__VA_ARGS__)); \
  }                                                                     \
  SerialMesh::node_iterator                                             \
  SerialMesh::FUNC_PREFIX##_end (FUNC_ARG)                              \
  {                                                                     \
    return node_iterator(_nodes.end(), _nodes.end(), Predicates::PRED<node_iterator_imp>(__VA_ARGS__)); \
  }                                                                     \
  SerialMesh::const_node_iterator                                       \
  SerialMesh::FUNC_PREFIX##_end (FUNC_ARG) const                        \
  {                                                                     \
    return const_node_iterator(_nodes.end(), _nodes.end(), Predicates::PRED<const_node_iterator_imp>(__VA_ARGS__)); \
  }                                                                     \
  ParallelMesh::node_iterator                                           \
  ParallelMesh::FUNC_PREFIX##_begin (FUNC_ARG)                          \
  {                                                                     \
    return node_iterator(_nodes.begin(), _nodes.end(), Predicates::PRED<node_iterator_imp>(__VA_ARGS__)); \
  }                                                                     \
  ParallelMesh::const_node_iterator                                     \
  ParallelMesh::FUNC_PREFIX##_begin (FUNC_ARG) const                    \
  {                                                                     \
    return const_node_iterator(_nodes.begin(), _nodes.end(), Predicates::PRED<const_node_iterator_imp>(__VA_ARGS__)); \
  }                                                                     \
  ParallelMesh::node_iterator                                           \
  ParallelMesh::FUNC_PREFIX##_end (FUNC_ARG)                            \
  {                                                                     \
    return node_iterator(_nodes.end(), _nodes.end(), Predicates::PRED<node_iterator_imp>(__VA_ARGS__)); \
  }                                                                     \
  ParallelMesh::const_node_iterator                                     \
  ParallelMesh::FUNC_PREFIX##_end (FUNC_ARG) const                      \
  {                                                                     \
    return const_node_iterator(_nodes.end(), _nodes.end(), Predicates::PRED<const_node_iterator_imp>(__VA_ARGS__)); \
  }


// Instantiate various element iterator accessor functions.
INSTANTIATE_ELEM_ACCESSORS(elements,                        NotNull,              /*none*/,                       /*none*/)
INSTANTIATE_ELEM_ACCESSORS(active_elements,                 Active,               /*none*/,                       /*none*/)
INSTANTIATE_ELEM_ACCESSORS(not_active_elements,             NotActive,            /*none*/,                       /*none*/)
INSTANTIATE_ELEM_ACCESSORS(ancestor_elements,               Ancestor,             /*none*/,                       /*none*/)
INSTANTIATE_ELEM_ACCESSORS(not_ancestor_elements,           NotAncestor,          /*none*/,                       /*none*/)
INSTANTIATE_ELEM_ACCESSORS(subactive_elements,              SubActive,            /*none*/,                       /*none*/)
INSTANTIATE_ELEM_ACCESSORS(not_subactive_elements,          NotSubActive,         /*none*/,                       /*none*/)
INSTANTIATE_ELEM_ACCESSORS(local_elements,                  Local,                /*none*/,                       this->processor_id())
INSTANTIATE_ELEM_ACCESSORS(semilocal_elements,              SemiLocal,            /*none*/,                       this->processor_id())
INSTANTIATE_ELEM_ACCESSORS(facelocal_elements,              FaceLocal,            /*none*/,                       this->processor_id())
INSTANTIATE_ELEM_ACCESSORS(not_local_elements,              NotLocal,             /*none*/,                       this->processor_id())
INSTANTIATE_ELEM_ACCESSORS(active_local_elements,           ActiveLocal,          /*none*/,                       this->processor_id())
INSTANTIATE_ELEM_ACCESSORS(active_not_local_elements,       ActiveNotLocal,       /*none*/,                       this->processor_id())
INSTANTIATE_ELEM_ACCESSORS(level_elements,                  Level,                unsigned int level,             level)
INSTANTIATE_ELEM_ACCESSORS(not_level_elements,              NotLevel,             unsigned int level,             level)
INSTANTIATE_ELEM_ACCESSORS(pid_elements,                    PID,                  processor_id_type proc_id,      proc_id)
INSTANTIATE_ELEM_ACCESSORS(type_elements,                   Type,                 ElemType type,                  type)
INSTANTIATE_ELEM_ACCESSORS(active_type_elements,            ActiveType,           ElemType type,                  type)
INSTANTIATE_ELEM_ACCESSORS(active_pid_elements,             ActivePID,            processor_id_type proc_id,      proc_id)
INSTANTIATE_ELEM_ACCESSORS(active_subdomain_elements,       ActiveSubdomain,      subdomain_id_type subdomain_id, subdomain_id)
INSTANTIATE_ELEM_ACCESSORS(ghost_elements,                  Ghost,                /*none*/,                       this->processor_id())
INSTANTIATE_ELEM_ACCESSORS(unpartitioned_elements,          PID,                  /*none*/,                       DofObject::invalid_processor_id)
INSTANTIATE_ELEM_ACCESSORS(local_level_elements,            LocalLevel,           unsigned int level,             this->processor_id(), level)
INSTANTIATE_ELEM_ACCESSORS(local_not_level_elements,        LocalNotLevel,        unsigned int level,             this->processor_id(), level)
INSTANTIATE_ELEM_ACCESSORS(active_local_subdomain_elements, ActiveLocalSubdomain, subdomain_id_type subdomain_id, this->processor_id(), subdomain_id)

// Instantiate various node iterator accessor functions.
INSTANTIATE_NODE_ACCESSORS(nodes,        NotNull, /*none*/,                  /*none*/)
INSTANTIATE_NODE_ACCESSORS(active_nodes, Active,  /*none*/,                  /*none*/)
INSTANTIATE_NODE_ACCESSORS(local_nodes,  Local,   /*none*/,                  this->processor_id())
INSTANTIATE_NODE_ACCESSORS(pid_nodes,    PID,     processor_id_type proc_id, proc_id)
INSTANTIATE_NODE_ACCESSORS(bnd_nodes,    BND,     /*none*/,                  this->get_boundary_info())
INSTANTIATE_NODE_ACCESSORS(bid_nodes,    BID,     boundary_id_type bndry_id, bndry_id, this->get_boundary_info())

} // namespace libMesh

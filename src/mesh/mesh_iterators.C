// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#define INSTANTIATE_ELEM_ACCESSORS(FUNC_PREFIX, PRED, FUNC_ARG, ...)    \
  ReplicatedMesh::element_iterator                                      \
  ReplicatedMesh::FUNC_PREFIX##_begin (FUNC_ARG)                        \
  {                                                                     \
    return element_iterator(_elements.begin(), _elements.end(), Predicates::PRED<elem_iterator_imp>(__VA_ARGS__)); \
  }                                                                     \
  ReplicatedMesh::const_element_iterator                                \
  ReplicatedMesh::FUNC_PREFIX##_begin (FUNC_ARG) const                  \
  {                                                                     \
    return const_element_iterator(_elements.begin(), _elements.end(), Predicates::PRED<const_elem_iterator_imp>(__VA_ARGS__)); \
  }                                                                     \
  ReplicatedMesh::element_iterator                                      \
  ReplicatedMesh::FUNC_PREFIX##_end (FUNC_ARG)                          \
  {                                                                     \
    return element_iterator(_elements.end(), _elements.end(), Predicates::PRED<elem_iterator_imp>(__VA_ARGS__)); \
  }                                                                     \
  ReplicatedMesh::const_element_iterator                                \
  ReplicatedMesh::FUNC_PREFIX##_end (FUNC_ARG) const                    \
  {                                                                     \
    return const_element_iterator(_elements.end(), _elements.end(), Predicates::PRED<const_elem_iterator_imp>(__VA_ARGS__)); \
  }                                                                     \
  DistributedMesh::element_iterator                                     \
  DistributedMesh::FUNC_PREFIX##_begin (FUNC_ARG)                       \
  {                                                                     \
    return element_iterator(_elements.begin(), _elements.end(), Predicates::PRED<elem_iterator_imp>(__VA_ARGS__)); \
  }                                                                     \
  DistributedMesh::const_element_iterator                               \
  DistributedMesh::FUNC_PREFIX##_begin (FUNC_ARG) const                 \
  {                                                                     \
    return const_element_iterator(_elements.begin(), _elements.end(), Predicates::PRED<const_elem_iterator_imp>(__VA_ARGS__)); \
  }                                                                     \
  DistributedMesh::element_iterator                                     \
  DistributedMesh::FUNC_PREFIX##_end (FUNC_ARG)                         \
  {                                                                     \
    return element_iterator(_elements.end(), _elements.end(), Predicates::PRED<elem_iterator_imp>(__VA_ARGS__)); \
  }                                                                     \
  DistributedMesh::const_element_iterator                               \
  DistributedMesh::FUNC_PREFIX##_end (FUNC_ARG) const                   \
  {                                                                     \
    return const_element_iterator(_elements.end(), _elements.end(), Predicates::PRED<const_elem_iterator_imp>(__VA_ARGS__)); \
  }



// This macro is similar to the one above except that it generates
// node iterator accessor functions.
#define INSTANTIATE_NODE_ACCESSORS(FUNC_PREFIX, PRED, FUNC_ARG, ...)    \
  ReplicatedMesh::node_iterator                                         \
  ReplicatedMesh::FUNC_PREFIX##_begin (FUNC_ARG)                        \
  {                                                                     \
    return node_iterator(_nodes.begin(), _nodes.end(), Predicates::PRED<node_iterator_imp>(__VA_ARGS__)); \
  }                                                                     \
  ReplicatedMesh::const_node_iterator                                   \
  ReplicatedMesh::FUNC_PREFIX##_begin (FUNC_ARG) const                  \
  {                                                                     \
    return const_node_iterator(_nodes.begin(), _nodes.end(), Predicates::PRED<const_node_iterator_imp>(__VA_ARGS__)); \
  }                                                                     \
  ReplicatedMesh::node_iterator                                         \
  ReplicatedMesh::FUNC_PREFIX##_end (FUNC_ARG)                          \
  {                                                                     \
    return node_iterator(_nodes.end(), _nodes.end(), Predicates::PRED<node_iterator_imp>(__VA_ARGS__)); \
  }                                                                     \
  ReplicatedMesh::const_node_iterator                                   \
  ReplicatedMesh::FUNC_PREFIX##_end (FUNC_ARG) const                    \
  {                                                                     \
    return const_node_iterator(_nodes.end(), _nodes.end(), Predicates::PRED<const_node_iterator_imp>(__VA_ARGS__)); \
  }                                                                     \
  DistributedMesh::node_iterator                                        \
  DistributedMesh::FUNC_PREFIX##_begin (FUNC_ARG)                       \
  {                                                                     \
    return node_iterator(_nodes.begin(), _nodes.end(), Predicates::PRED<node_iterator_imp>(__VA_ARGS__)); \
  }                                                                     \
  DistributedMesh::const_node_iterator                                  \
  DistributedMesh::FUNC_PREFIX##_begin (FUNC_ARG) const                 \
  {                                                                     \
    return const_node_iterator(_nodes.begin(), _nodes.end(), Predicates::PRED<const_node_iterator_imp>(__VA_ARGS__)); \
  }                                                                     \
  DistributedMesh::node_iterator                                        \
  DistributedMesh::FUNC_PREFIX##_end (FUNC_ARG)                         \
  {                                                                     \
    return node_iterator(_nodes.end(), _nodes.end(), Predicates::PRED<node_iterator_imp>(__VA_ARGS__)); \
  }                                                                     \
  DistributedMesh::const_node_iterator                                  \
  DistributedMesh::FUNC_PREFIX##_end (FUNC_ARG) const                   \
  {                                                                     \
    return const_node_iterator(_nodes.end(), _nodes.end(), Predicates::PRED<const_node_iterator_imp>(__VA_ARGS__)); \
  }

// Use an empty preprocessor token to silence older compilers that
// still warn about empty macro arguments.
#define EMPTY

// Use a second macro layer to allow us to pass commas into a macro
// argument
#define LIBMESH_COMMA ,

// Instantiate various element iterator accessor functions.
INSTANTIATE_ELEM_ACCESSORS(elements,                        NotNull,              EMPTY,                          EMPTY)
INSTANTIATE_ELEM_ACCESSORS(active_elements,                 Active,               EMPTY,                          EMPTY)
INSTANTIATE_ELEM_ACCESSORS(not_active_elements,             NotActive,            EMPTY,                          EMPTY)
INSTANTIATE_ELEM_ACCESSORS(ancestor_elements,               Ancestor,             EMPTY,                          EMPTY)
INSTANTIATE_ELEM_ACCESSORS(not_ancestor_elements,           NotAncestor,          EMPTY,                          EMPTY)
INSTANTIATE_ELEM_ACCESSORS(subactive_elements,              SubActive,            EMPTY,                          EMPTY)
INSTANTIATE_ELEM_ACCESSORS(not_subactive_elements,          NotSubActive,         EMPTY,                          EMPTY)
INSTANTIATE_ELEM_ACCESSORS(local_elements,                  Local,                EMPTY,                          this->processor_id())
INSTANTIATE_ELEM_ACCESSORS(semilocal_elements,              ActiveSemiLocal,      EMPTY,                          this->processor_id())
INSTANTIATE_ELEM_ACCESSORS(active_semilocal_elements,       ActiveSemiLocal,      EMPTY,                          this->processor_id())
INSTANTIATE_ELEM_ACCESSORS(facelocal_elements,              FaceLocal,            EMPTY,                          this->processor_id())
INSTANTIATE_ELEM_ACCESSORS(not_local_elements,              NotLocal,             EMPTY,                          this->processor_id())
INSTANTIATE_ELEM_ACCESSORS(active_local_elements,           ActiveLocal,          EMPTY,                          this->processor_id())
INSTANTIATE_ELEM_ACCESSORS(active_not_local_elements,       ActiveNotLocal,       EMPTY,                          this->processor_id())
INSTANTIATE_ELEM_ACCESSORS(level_elements,                  Level,                unsigned int level,             level)
INSTANTIATE_ELEM_ACCESSORS(not_level_elements,              NotLevel,             unsigned int level,             level)
INSTANTIATE_ELEM_ACCESSORS(pid_elements,                    PID,                  processor_id_type proc_id,      proc_id)
INSTANTIATE_ELEM_ACCESSORS(type_elements,                   Type,                 ElemType type,                  type)
INSTANTIATE_ELEM_ACCESSORS(active_type_elements,            ActiveType,           ElemType type,                  type)
INSTANTIATE_ELEM_ACCESSORS(active_pid_elements,             ActivePID,            processor_id_type proc_id,      proc_id)
INSTANTIATE_ELEM_ACCESSORS(active_subdomain_elements,       ActiveSubdomain,      subdomain_id_type subdomain_id, subdomain_id)
INSTANTIATE_ELEM_ACCESSORS(active_subdomain_set_elements,   ActiveSubdomainSet,   std::set<subdomain_id_type> ss, ss)
INSTANTIATE_ELEM_ACCESSORS(ghost_elements,                  Ghost,                EMPTY,                          this->processor_id())
INSTANTIATE_ELEM_ACCESSORS(evaluable_elements,              Evaluable,            const DofMap & dof_map LIBMESH_COMMA unsigned int var_num, this->processor_id(), dof_map, var_num)
INSTANTIATE_ELEM_ACCESSORS(unpartitioned_elements,          PID,                  EMPTY,                          DofObject::invalid_processor_id)
INSTANTIATE_ELEM_ACCESSORS(active_unpartitioned_elements,   ActivePID,            EMPTY,                          DofObject::invalid_processor_id)

#ifdef LIBMESH_ENABLE_AMR
INSTANTIATE_ELEM_ACCESSORS(flagged_elements,                Flagged,              unsigned char rflag,            rflag)
INSTANTIATE_ELEM_ACCESSORS(flagged_pid_elements,            FlaggedPID,           unsigned char rflag LIBMESH_COMMA processor_id_type pid,    rflag, pid)
#endif

INSTANTIATE_ELEM_ACCESSORS(local_level_elements,            LocalLevel,           unsigned int level,             this->processor_id(), level)
INSTANTIATE_ELEM_ACCESSORS(local_not_level_elements,        LocalNotLevel,        unsigned int level,             this->processor_id(), level)
INSTANTIATE_ELEM_ACCESSORS(active_local_subdomain_elements, ActiveLocalSubdomain, subdomain_id_type subdomain_id, this->processor_id(), subdomain_id)

// Instantiate various node iterator accessor functions.
INSTANTIATE_NODE_ACCESSORS(nodes,        NotNull, EMPTY,                               EMPTY)
INSTANTIATE_NODE_ACCESSORS(active_nodes, Active,  EMPTY,                               EMPTY)
INSTANTIATE_NODE_ACCESSORS(local_nodes,  Local,   EMPTY,                               this->processor_id())
INSTANTIATE_NODE_ACCESSORS(pid_nodes,    PID,     processor_id_type proc_id,           proc_id)
INSTANTIATE_NODE_ACCESSORS(bnd_nodes,    BND,     EMPTY,                               this->get_boundary_info())
INSTANTIATE_NODE_ACCESSORS(bid_nodes,    BID,     boundary_id_type bndry_id, bndry_id, this->get_boundary_info())
INSTANTIATE_NODE_ACCESSORS(evaluable_nodes, Evaluable, const DofMap & dof_map LIBMESH_COMMA unsigned int var_num, this->processor_id(), dof_map, var_num)

} // namespace libMesh

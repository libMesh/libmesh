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
#include "libmesh/boundary_info.h"
#include "libmesh/distributed_mesh.h"
#include "libmesh/mesh_base.h"
#include "libmesh/node.h"
#include "libmesh/parallel.h"
#include "libmesh/parallel_mesh.h"

// Helper functions in anonymous namespace

namespace
{
using namespace libMesh;

#ifdef LIBMESH_ENABLE_UNIQUE_ID
static const unsigned int header_size = 3;
#else
static const unsigned int header_size = 2;
#endif

// use "(a+b-1)/b" trick to get a/b to round up
static const unsigned int idtypes_per_Real =
  (sizeof(Real) + sizeof(largest_id_type) - 1) / sizeof(largest_id_type);

#ifndef NDEBUG
// Currently this constant is only used for debugging.
static const largest_id_type node_magic_header = 1234567890;
#endif
}


namespace libMesh
{

namespace Parallel
{

template <>
template <>
unsigned int
Packing<const Node *>::packable_size (const Node * const & node,
                                      const MeshBase * mesh)
{
  return
#ifndef NDEBUG
    1 + // add an int for the magic header when testing
#endif
    header_size + LIBMESH_DIM*idtypes_per_Real +
    node->packed_indexing_size() +
    1 + mesh->get_boundary_info().n_boundary_ids(node);
}



template <>
template <>
unsigned int
Packing<const Node *>::packed_size (const std::vector<largest_id_type>::const_iterator in)
{
  const unsigned int pre_indexing_size =
#ifndef NDEBUG
    1 + // add an int for the magic header when testing
#endif
    header_size + LIBMESH_DIM*idtypes_per_Real;

  const unsigned int indexing_size =
    DofObject::unpackable_indexing_size(in+pre_indexing_size);

  const int n_bcs = cast_int<int>
    (*(in + pre_indexing_size + indexing_size));
  libmesh_assert_greater_equal (n_bcs, 0);

  return pre_indexing_size + indexing_size + 1 + n_bcs;
}



template <>
template <>
unsigned int
Packing<const Node *>::packed_size (const std::vector<largest_id_type>::iterator in)
{
  return packed_size(std::vector<largest_id_type>::const_iterator(in));
}



template <>
template <>
unsigned int
Packing<const Node *>::packable_size (const Node * const & node,
                                      const DistributedMesh * mesh)
{
  return packable_size(node, static_cast<const MeshBase *>(mesh));
}



template <>
template <>
unsigned int
Packing<const Node *>::packable_size (const Node * const & node,
                                      const ParallelMesh * mesh)
{
  return packable_size(node, static_cast<const MeshBase *>(mesh));
}



template <>
template <>
void
Packing<const Node *>::pack (const Node * const & node,
                             std::back_insert_iterator<std::vector<largest_id_type>> data_out,
                             const MeshBase * mesh)
{
  libmesh_assert(node);

#ifndef NDEBUG
  *data_out++ = (node_magic_header);
#endif

  *data_out++ = (static_cast<largest_id_type>(node->processor_id()));
  *data_out++ = (static_cast<largest_id_type>(node->id()));

#ifdef LIBMESH_ENABLE_UNIQUE_ID
  if (node->valid_unique_id())
    *data_out++ = (static_cast<largest_id_type>(node->unique_id()));
  else
    // OK to send invalid unique id, we must not own this DOF
    *data_out++ = (static_cast<largest_id_type>(DofObject::invalid_unique_id));
#endif

  // use "(a+b-1)/b" trick to get a/b to round up
  static const unsigned int idtypes_per_Real =
    (sizeof(Real) + sizeof(largest_id_type) - 1) / sizeof(largest_id_type);

  for (unsigned int i=0; i != LIBMESH_DIM; ++i)
    {
      const largest_id_type * Real_as_idtypes =
        reinterpret_cast<const largest_id_type *>(&((*node)(i)));
      for (unsigned int j=0; j != idtypes_per_Real; ++j)
        {
          *data_out++ =(Real_as_idtypes[j]);
        }
    }

  // Add any DofObject indices
  node->pack_indexing(data_out);

  // Add any nodal boundary condition ids
  std::vector<boundary_id_type> bcs;
  mesh->get_boundary_info().boundary_ids(node, bcs);

  libmesh_assert(bcs.size() < std::numeric_limits<largest_id_type>::max());

  *data_out++ =(bcs.size());

  for (std::vector<boundary_id_type>::iterator bc_it=bcs.begin();
       bc_it != bcs.end(); ++bc_it)
    *data_out++ =(*bc_it);
}



template <>
template <>
void
Packing<const Node *>::pack (const Node * const & node,
                             std::back_insert_iterator<std::vector<largest_id_type>> data_out,
                             const DistributedMesh * mesh)
{
  pack(node, data_out, static_cast<const MeshBase*>(mesh));
}



template <>
template <>
void
Packing<const Node *>::pack (const Node * const & node,
                             std::back_insert_iterator<std::vector<largest_id_type>> data_out,
                             const ParallelMesh * mesh)
{
  pack(node, data_out, static_cast<const MeshBase*>(mesh));
}



template <>
template <>
Node *
Packing<Node *>::unpack (std::vector<largest_id_type>::const_iterator in,
                         MeshBase * mesh)
{
#ifndef NDEBUG
  const std::vector<largest_id_type>::const_iterator original_in = in;
  const largest_id_type incoming_header = *in++;
  libmesh_assert_equal_to (incoming_header, node_magic_header);
#endif

  const processor_id_type processor_id = cast_int<processor_id_type>(*in++);
  libmesh_assert(processor_id == DofObject::invalid_processor_id ||
                 processor_id < mesh->n_processors());

  const dof_id_type id = cast_int<dof_id_type>(*in++);

#ifdef LIBMESH_ENABLE_UNIQUE_ID
  const unique_id_type unique_id = cast_int<unique_id_type>(*in++);
#endif

  Node * node = mesh->query_node_ptr(id);

  if (node)
    {
      libmesh_assert_equal_to (node->processor_id(), processor_id);

      // We currently don't communicate mesh motion via packed Nodes,
      // so it should usually be safe to assume (and assert) that Node
      // locations are consistent between processors.
      //
      // There may be exceptions due to rounding in file I/O, so we'll
      // only assert equality to within a tight tolerance, and we'll
      // believe the sender's node locations over our own if we don't
      // own the node.
      for (unsigned int i=0; i != LIBMESH_DIM; ++i)
        {
          const Real & idtypes_as_Real = *(reinterpret_cast<const Real *>(&(*in)));
          libmesh_assert_less_equal ((*node)(i), idtypes_as_Real + (std::max(Real(1),idtypes_as_Real)*TOLERANCE*TOLERANCE));
          libmesh_assert_greater_equal ((*node)(i), idtypes_as_Real - (std::max(Real(1),idtypes_as_Real)*TOLERANCE*TOLERANCE));

          if (processor_id != mesh->processor_id())
            (*node)(i) = idtypes_as_Real;
          in += idtypes_per_Real;
        }

      if (!node->has_dofs())
        {
          node->unpack_indexing(in);
          libmesh_assert_equal_to (DofObject::unpackable_indexing_size(in),
                                   node->packed_indexing_size());
          in += node->packed_indexing_size();
        }
      else
        {
          // FIXME: We should add some debug mode tests to ensure that
          // the encoded indexing is consistent
          in += DofObject::unpackable_indexing_size(in);
        }
    }
  else
    {
      // If we don't already have it, we need to allocate it
      node = new Node();

      for (unsigned int i=0; i != LIBMESH_DIM; ++i)
        {
          const Real * idtypes_as_Real = reinterpret_cast<const Real *>(&(*in));
          (*node)(i) = *idtypes_as_Real;
          in += idtypes_per_Real;
        }

      node->set_id() = id;
#ifdef LIBMESH_ENABLE_UNIQUE_ID
      node->set_unique_id() = unique_id;
#endif
      node->processor_id() = processor_id;

      node->unpack_indexing(in);
      libmesh_assert_equal_to (DofObject::unpackable_indexing_size(in),
                               node->packed_indexing_size());
      in += node->packed_indexing_size();
    }

  // FIXME: We should add some debug mode tests to ensure that the
  // encoded boundary conditions are consistent

  // Add any nodal boundary condition ids
  const largest_id_type num_bcs = *in++;
  // libmesh_assert_greater_equal (num_bcs, 0);

  for (largest_id_type bc_it=0; bc_it < num_bcs; bc_it++)
    mesh->get_boundary_info().add_node
      (node, cast_int<boundary_id_type>(*in++));

#ifndef NDEBUG
  libmesh_assert (in - original_in ==
                  cast_int<int>
                  (Parallel::Packing<const Node*>::packed_size(original_in)));
#endif

  return node;
}



template <>
template <>
Node *
Packing<Node *>::unpack (std::vector<largest_id_type>::const_iterator in,
                         DistributedMesh * mesh)
{
  return unpack(in, static_cast<MeshBase*>(mesh));
}



template <>
template <>
Node *
Packing<Node *>::unpack (std::vector<largest_id_type>::const_iterator in,
                         ParallelMesh * mesh)
{
  return unpack(in, static_cast<MeshBase*>(mesh));
}

} // namespace Parallel

} // namespace libMesh

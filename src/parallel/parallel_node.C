// The libMesh Finite Element Library.
// Copyright (C) 2002-2012 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

  static const largest_id_type node_magic_header = 1234567890;
}


namespace libMesh
{

#ifdef LIBMESH_HAVE_MPI

namespace Parallel
{

template <>
unsigned int packable_size (const Node* node, const MeshBase* mesh)
{
  return
#ifndef NDEBUG
         1 + // add an int for the magic header when testing
#endif
	 header_size + LIBMESH_DIM*idtypes_per_Real +
         node->packed_indexing_size() +
         1 + mesh->boundary_info->n_boundary_ids(node);
}



template <>
unsigned int packed_size (const Node*,
			  const std::vector<largest_id_type>::const_iterator in)
{
  const unsigned int pre_indexing_size =
#ifndef NDEBUG
    1 + // add an int for the magic header when testing
#endif
    header_size + LIBMESH_DIM*idtypes_per_Real;

  const unsigned int indexing_size =
    DofObject::unpackable_indexing_size(in+pre_indexing_size);

  const int n_bcs =
    *(in + pre_indexing_size + indexing_size);
  libmesh_assert_greater_equal (n_bcs, 0);

  return pre_indexing_size + indexing_size + 1 + n_bcs;
}



template <>
unsigned int packed_size (const Node* n,
			  const std::vector<largest_id_type>::iterator in)
{
  return packed_size(n, std::vector<largest_id_type>::const_iterator(in));
}



template <>
unsigned int packable_size (const Node* node, const ParallelMesh* mesh)
{
  return packable_size(node, static_cast<const MeshBase*>(mesh));
}



template <>
void pack (const Node* node,
           std::vector<largest_id_type>& data,
           const MeshBase* mesh)
{
  libmesh_assert(node);

  // This should be redundant when used with Parallel::pack_range()
  // data.reserve (data.size() + Parallel::packable_size(node, mesh));

#ifndef NDEBUG
  data.push_back (node_magic_header);
#endif

  data.push_back (static_cast<largest_id_type>(node->processor_id()));
  data.push_back (static_cast<largest_id_type>(node->id()));

#ifdef LIBMESH_ENABLE_UNIQUE_ID
  if (node->valid_unique_id())
    data.push_back (static_cast<largest_id_type>(node->unique_id()));
  else
    // OK to send invalid unique id, we must not own this DOF
    data.push_back (static_cast<largest_id_type>(DofObject::invalid_unique_id));
#endif

  // use "(a+b-1)/b" trick to get a/b to round up
  static const unsigned int idtypes_per_Real =
    (sizeof(Real) + sizeof(largest_id_type) - 1) / sizeof(largest_id_type);

  for (unsigned int i=0; i != LIBMESH_DIM; ++i)
    {
      const largest_id_type* Real_as_idtypes =
        reinterpret_cast<const largest_id_type*>(&((*node)(i)));
      for (unsigned int j=0; j != idtypes_per_Real; ++j)
        {
          data.push_back(Real_as_idtypes[j]);
        }
    }

#ifndef NDEBUG
  const std::size_t start_indices = data.size();
#endif
  // Add any DofObject indices
  node->pack_indexing(std::back_inserter(data));

  libmesh_assert(node->packed_indexing_size() ==
		 DofObject::unpackable_indexing_size(data.begin() +
						     start_indices));

  libmesh_assert_equal_to (node->packed_indexing_size(),
		          data.size() - start_indices);

  // Add any nodal boundary condition ids
  std::vector<boundary_id_type> bcs =
    mesh->boundary_info->boundary_ids(node);

  libmesh_assert(bcs.size() < std::numeric_limits<largest_id_type>::max());

  data.push_back(bcs.size());

  for (std::size_t bc_it=0; bc_it < bcs.size(); bc_it++)
    data.push_back(bcs[bc_it]);
}



template <>
void pack (const Node* node,
           std::vector<largest_id_type>& data,
           const ParallelMesh* mesh)
{
  pack(node, data, static_cast<const MeshBase*>(mesh));
}



template <>
void unpack (std::vector<largest_id_type>::const_iterator in,
             Node** out,
             MeshBase* mesh)
{
#ifndef NDEBUG
  const std::vector<largest_id_type>::const_iterator original_in = in;
  const largest_id_type incoming_header = *in++;
  libmesh_assert_equal_to (incoming_header, node_magic_header);
#endif

  const processor_id_type processor_id = static_cast<processor_id_type>(*in++);
  libmesh_assert(processor_id == DofObject::invalid_processor_id ||
                 processor_id < mesh->n_processors());

  const dof_id_type id = static_cast<dof_id_type>(*in++);

#ifdef LIBMESH_ENABLE_UNIQUE_ID
  const unique_id_type unique_id = static_cast<unique_id_type>(*in++);
#endif

  Node *node = mesh->query_node_ptr(id);

  if (node)
    {
      libmesh_assert_equal_to (node->processor_id(), processor_id);

      // We currently don't communicate mesh motion via packed Nodes,
      // so it should be safe to assume (and assert) that Node
      // locations are consistent between processors
#ifndef NDEBUG
      for (unsigned int i=0; i != LIBMESH_DIM; ++i)
        {
          const Real* idtypes_as_Real = reinterpret_cast<const Real*>(&(*in));
          libmesh_assert_equal_to ((*node)(i), *idtypes_as_Real);
          in += idtypes_per_Real;
        }
#else
      in += LIBMESH_DIM * idtypes_per_Real;
#endif // !NDEBUG

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

      *out = node;
    }
  else
    {
      // If we don't already have it, we need to allocate it
      node = new Node();

      for (unsigned int i=0; i != LIBMESH_DIM; ++i)
        {
          const Real* idtypes_as_Real = reinterpret_cast<const Real*>(&(*in));
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

  for(largest_id_type bc_it=0; bc_it < num_bcs; bc_it++)
    mesh->boundary_info->add_node (node, *in++);

  *out = node;

#ifndef NDEBUG
  libmesh_assert (in - original_in ==
		  static_cast<int>
                    (Parallel::packed_size(node, original_in)));
#endif
}



template <>
void unpack (std::vector<largest_id_type>::const_iterator in,
             Node** out,
             ParallelMesh* mesh)
{
  unpack(in, out, static_cast<MeshBase*>(mesh));
}

} // namespace Parallel

#endif // #ifdef LIBMESH_HAVE_MPI

} // namespace libMesh

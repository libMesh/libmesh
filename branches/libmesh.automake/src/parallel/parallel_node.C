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
#include "boundary_info.h"
#include "mesh_base.h"
#include "node.h"
#include "parallel.h"
#include "parallel_mesh.h"

// Helper functions in anonymous namespace

namespace 
{
  static const unsigned int header_size = 2;

  // use "(a+b-1)/b" trick to get a/b to round up
  static const unsigned int ints_per_Real = 
    (sizeof(Real) + sizeof(int) - 1) / sizeof(int);

  static const int node_magic_header = 1234567890;
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
	 header_size + LIBMESH_DIM*ints_per_Real +
         node->packed_indexing_size() +
         1 + mesh->boundary_info->n_boundary_ids(node);
}



template <>
unsigned int packed_size (const Node*,
			  std::vector<int>::const_iterator in)
{
  const unsigned int pre_indexing_size = 
#ifndef NDEBUG
    1 + // add an int for the magic header when testing
#endif
    header_size + LIBMESH_DIM*ints_per_Real;

  const unsigned int indexing_size = 
    DofObject::unpackable_indexing_size(in+pre_indexing_size);

  const int n_bcs =
    *(in + pre_indexing_size + indexing_size);
  libmesh_assert (n_bcs >= 0);

  return pre_indexing_size + indexing_size + 1 + n_bcs;
}



template <>
unsigned int packable_size (const Node* node, const ParallelMesh* mesh)
{
  return packable_size(node, static_cast<const MeshBase*>(mesh));
}



template <>
void pack (const Node* node,
           std::vector<int>& data,
           const MeshBase* mesh)
{
  libmesh_assert (node != NULL);

  // This should be redundant when used with Parallel::pack_range()
  // data.reserve (data.size() + Parallel::packable_size(node, mesh));

#ifndef NDEBUG
  data.push_back (node_magic_header);
#endif

  data.push_back (static_cast<int>(node->processor_id()));
  data.push_back (static_cast<int>(node->id()));

  // use "(a+b-1)/b" trick to get a/b to round up
  static const unsigned int ints_per_Real = 
    (sizeof(Real) + sizeof(int) - 1) / sizeof(int);

  for (unsigned int i=0; i != LIBMESH_DIM; ++i)
    {
      const int* Real_as_ints = reinterpret_cast<const int*>(&((*node)(i)));
      for (unsigned int j=0; j != ints_per_Real; ++j)
        {
          data.push_back(Real_as_ints[j]);
        }
    }

#ifndef NDEBUG
  const int start_indices = data.size();
#endif
  // Add any DofObject indices
  node->pack_indexing(std::back_inserter(data));

  libmesh_assert(node->packed_indexing_size() == 
		 DofObject::unpackable_indexing_size(data.begin() +
						     start_indices));

  libmesh_assert(node->packed_indexing_size() ==
		 data.size() - start_indices);

  // Add any nodal boundary condition ids
  std::vector<boundary_id_type> bcs =
    mesh->boundary_info->boundary_ids(node);

  data.push_back(bcs.size());

  for(unsigned int bc_it=0; bc_it < bcs.size(); bc_it++)
    data.push_back(bcs[bc_it]);
}



template <>
void pack (const Node* node,
           std::vector<int>& data,
           const ParallelMesh* mesh)
{
  pack(node, data, static_cast<const MeshBase*>(mesh));
}



template <>
void unpack (std::vector<int>::const_iterator in,
             Node** out,
             MeshBase* mesh)
{
#ifndef NDEBUG
  const std::vector<int>::const_iterator original_in = in;
  const int incoming_header = *in++;
  libmesh_assert (incoming_header == node_magic_header);
#endif

  const unsigned int processor_id = static_cast<unsigned int>(*in++);
  libmesh_assert(processor_id == DofObject::invalid_processor_id ||
                 processor_id < libMesh::n_processors());

  const unsigned int id = static_cast<unsigned int>(*in++);

  Node *node = mesh->query_node_ptr(id);

  if (node)
    {
      libmesh_assert(node->processor_id() == processor_id);

      // We currently don't communicate mesh motion via packed Nodes,
      // so it should be safe to assume (and assert) that Node
      // locations are consistent between processors
#ifndef NDEBUG
      for (unsigned int i=0; i != LIBMESH_DIM; ++i)
        {
          const Real* ints_as_Real = reinterpret_cast<const Real*>(&(*in));
          libmesh_assert((*node)(i) == *ints_as_Real);
          in += ints_per_Real;
        }
#endif // !NDEBUG

      if (!node->has_dofs())
        {
	  node->unpack_indexing(in);
	  libmesh_assert (DofObject::unpackable_indexing_size(in) ==
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
          const Real* ints_as_Real = reinterpret_cast<const Real*>(&(*in));
          (*node)(i) = *ints_as_Real;
          in += ints_per_Real;
        }

      node->set_id() = id;
      node->processor_id() = processor_id;

      node->unpack_indexing(in);
      libmesh_assert (DofObject::unpackable_indexing_size(in) ==
		      node->packed_indexing_size());
      in += node->packed_indexing_size();
    }

  // FIXME: We should add some debug mode tests to ensure that the
  // encoded boundary conditions are consistent

  // Add any nodal boundary condition ids
  const int num_bcs = *in++;
  libmesh_assert(num_bcs >= 0);

  for(int bc_it=0; bc_it < num_bcs; bc_it++)
    mesh->boundary_info->add_node (node, *in++);

  *out = node;

#ifndef NDEBUG
  libmesh_assert (in - original_in == 
		  Parallel::packed_size(node, original_in));
#endif
}



template <>
void unpack (std::vector<int>::const_iterator in,
             Node** out,
             ParallelMesh* mesh)
{
  unpack(in, out, static_cast<MeshBase*>(mesh));
}
  
} // namespace Parallel

#endif // #ifdef LIBMESH_HAVE_MPI

} // namespace libMesh

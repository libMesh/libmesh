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
#include "libmesh/elem.h"
#include "libmesh/mesh_base.h"
#include "libmesh/parallel.h"
#include "libmesh/parallel_mesh.h"
#include "libmesh/remote_elem.h"

// Helper functions in anonymous namespace

namespace
{
#ifdef LIBMESH_ENABLE_UNIQUE_ID
  static const unsigned int header_size = 11;
#else
  static const unsigned int header_size = 10;
#endif

  static const largest_id_type elem_magic_header = 987654321;
}


namespace libMesh
{

#ifdef LIBMESH_HAVE_MPI

namespace Parallel
{

template <>
unsigned int packed_size (const Elem*,
			  std::vector<largest_id_type>::const_iterator in)
{
#ifndef NDEBUG
  const largest_id_type packed_header = *in++;
  libmesh_assert_equal_to (packed_header, elem_magic_header);
#endif

  // int 0: level
  const unsigned int level =
    static_cast<unsigned int>(*in);

  // int 4: element type
  const int typeint = *(in+4);
  libmesh_assert_greater_equal (typeint, 0);
  libmesh_assert_less (typeint, INVALID_ELEM);
  const ElemType type =
    static_cast<ElemType>(typeint);

  const unsigned int n_nodes =
    Elem::type_to_n_nodes_map[type];

  const unsigned int n_sides =
    Elem::type_to_n_sides_map[type];

  const unsigned int pre_indexing_size =
    header_size + n_nodes + n_sides;

  const unsigned int indexing_size =
    DofObject::unpackable_indexing_size(in+pre_indexing_size);

  unsigned int total_packed_bc_data = 0;
  if (level == 0)
    {
      for (unsigned int s = 0; s != n_sides; ++s)
        {
	  const int n_bcs =
	    *(in + pre_indexing_size + indexing_size +
	      total_packed_bc_data++);
	  libmesh_assert_greater_equal (n_bcs, 0);
          total_packed_bc_data += n_bcs;
        }
    }

  return
#ifndef NDEBUG
    1 + // Account for magic header
#endif
    pre_indexing_size + indexing_size + total_packed_bc_data;
}



template <>
unsigned int packed_size (const Elem* e,
			  std::vector<largest_id_type>::iterator in)
{
  return packed_size(e, std::vector<largest_id_type>::const_iterator(in));
}



template <>
unsigned int packable_size (const Elem* elem, const MeshBase* mesh)
{
  unsigned int total_packed_bcs = 0;
  if (elem->level() == 0)
    {
      total_packed_bcs += elem->n_sides();
      for (unsigned int s = 0; s != elem->n_sides(); ++s)
        total_packed_bcs += mesh->boundary_info->n_boundary_ids(elem,s);
    }

  return
#ifndef NDEBUG
         1 + // add an int for the magic header when testing
#endif
	 header_size + elem->n_nodes() +
         elem->n_neighbors() +
         elem->packed_indexing_size() + total_packed_bcs;
}



template <>
unsigned int packable_size (const Elem* elem, const ParallelMesh* mesh)
{
  return packable_size(elem, static_cast<const MeshBase*>(mesh));
}



template <>
void pack (const Elem* elem,
           std::vector<largest_id_type>& data,
           const MeshBase* mesh)
{
  libmesh_assert(elem);

  // This should be redundant when used with Parallel::pack_range()
  // data.reserve (data.size() + Parallel::packable_size(elem, mesh));

#ifndef NDEBUG
  data.push_back (elem_magic_header);
#endif

#ifdef LIBMESH_ENABLE_AMR
  data.push_back (static_cast<largest_id_type>(elem->level()));
  data.push_back (static_cast<largest_id_type>(elem->p_level()));
  data.push_back (static_cast<largest_id_type>(elem->refinement_flag()));
  data.push_back (static_cast<largest_id_type>(elem->p_refinement_flag()));
#else
  data.push_back (0);
  data.push_back (0);
  data.push_back (0);
  data.push_back (0);
#endif
  data.push_back (static_cast<largest_id_type>(elem->type()));
  data.push_back (elem->processor_id());
  data.push_back (elem->subdomain_id());
  data.push_back (elem->id());

#ifdef LIBMESH_ENABLE_UNIQUE_ID
  if (elem->valid_unique_id())
    data.push_back (static_cast<largest_id_type>(elem->unique_id()));
  else
    // OK to send invalid unique id, we must not own this DOF
    data.push_back (static_cast<largest_id_type>(DofObject::invalid_unique_id));
#endif

#ifdef LIBMESH_ENABLE_AMR
  // use parent_ID of -1 to indicate a level 0 element
  if (elem->level() == 0)
    {
      data.push_back(-1);
      data.push_back(-1);
    }
  else
    {
      data.push_back(elem->parent()->id());
      data.push_back(elem->parent()->which_child_am_i(elem));
    }
#else
  data.push_back (-1);
  data.push_back (-1);
#endif

  for (unsigned int n=0; n<elem->n_nodes(); n++)
    data.push_back (elem->node(n));

  for (unsigned int n=0; n<elem->n_neighbors(); n++)
    {
      const Elem *neigh = elem->neighbor(n);
      if (neigh)
        data.push_back (neigh->id());
      else
        data.push_back (-1);
    }

#ifndef NDEBUG
  const std::size_t start_indices = data.size();
#endif
  // Add any DofObject indices
  elem->pack_indexing(std::back_inserter(data));

  libmesh_assert(elem->packed_indexing_size() ==
		 DofObject::unpackable_indexing_size(data.begin() +
						     start_indices));

  libmesh_assert_equal_to (elem->packed_indexing_size(),
		          data.size() - start_indices);


  // If this is a coarse element,
  // Add any element side boundary condition ids
  if (elem->level() == 0)
    for (unsigned int s = 0; s != elem->n_sides(); ++s)
      {
        std::vector<boundary_id_type> bcs =
          mesh->boundary_info->boundary_ids(elem, s);

        data.push_back(bcs.size());

        for(unsigned int bc_it=0; bc_it < bcs.size(); bc_it++)
          data.push_back(bcs[bc_it]);
      }
}



template <>
void pack (const Elem* elem,
           std::vector<largest_id_type>& data,
           const ParallelMesh* mesh)
{
  pack(elem, data, static_cast<const MeshBase*>(mesh));
}



// FIXME - this needs serious work to be 64-bit compatible
template <>
void unpack(std::vector<largest_id_type>::const_iterator in,
            Elem** out,
            MeshBase* mesh)
{
#ifndef NDEBUG
  const std::vector<largest_id_type>::const_iterator original_in = in;

  const largest_id_type incoming_header = *in++;
  libmesh_assert_equal_to (incoming_header, elem_magic_header);
#endif

  // int 0: level
  const unsigned int level =
    static_cast<unsigned int>(*in++);

#ifdef LIBMESH_ENABLE_AMR
  // int 1: p level
  const unsigned int p_level =
    static_cast<unsigned int>(*in++);

  // int 2: refinement flag
  const int rflag = *in++;
  libmesh_assert_greater_equal (rflag, 0);
  libmesh_assert_less (rflag, Elem::INVALID_REFINEMENTSTATE);
  const Elem::RefinementState refinement_flag =
    static_cast<Elem::RefinementState>(rflag);

  // int 3: p refinement flag
  const int pflag = *in++;
  libmesh_assert_greater_equal (pflag, 0);
  libmesh_assert_less (pflag, Elem::INVALID_REFINEMENTSTATE);
  const Elem::RefinementState p_refinement_flag =
    static_cast<Elem::RefinementState>(pflag);
#else
  in += 3;
#endif // LIBMESH_ENABLE_AMR

  // int 4: element type
  const int typeint = *in++;
  libmesh_assert_greater_equal (typeint, 0);
  libmesh_assert_less (typeint, INVALID_ELEM);
  const ElemType type =
    static_cast<ElemType>(typeint);

  const unsigned int n_nodes =
    Elem::type_to_n_nodes_map[type];

  // int 5: processor id
  const processor_id_type processor_id =
    static_cast<processor_id_type>(*in++);
  libmesh_assert (processor_id < mesh->n_processors() ||
                  processor_id == DofObject::invalid_processor_id);

  // int 6: subdomain id
  const subdomain_id_type subdomain_id =
    static_cast<subdomain_id_type>(*in++);

  // int 7: dof object id
  const dof_id_type id =
    static_cast<dof_id_type>(*in++);
  libmesh_assert_not_equal_to (id, DofObject::invalid_id);

#ifdef LIBMESH_ENABLE_UNIQUE_ID
  // int 8: dof object unique id
  const unique_id_type unique_id =
    static_cast<unique_id_type>(*in++);
#endif

#ifdef LIBMESH_ENABLE_AMR
  // int 9: parent dof object id
  const dof_id_type parent_id =
    static_cast<dof_id_type>(*in++);
  libmesh_assert (level == 0 || parent_id != DofObject::invalid_id);
  libmesh_assert (level != 0 || parent_id == DofObject::invalid_id);

  // int 10: local child id
  const unsigned int which_child_am_i =
    static_cast<unsigned int>(*in++);
#else
  in += 2;
#endif // LIBMESH_ENABLE_AMR

  // Make sure we don't miscount above when adding the "magic" header
  // plus the real data header
  libmesh_assert_equal_to (in - original_in, header_size + 1);

  Elem *elem = mesh->query_elem(id);

  // if we already have this element, make sure its
  // properties match, and update any missing neighbor
  // links, but then go on
  if (elem)
    {
      libmesh_assert_equal_to (elem->level(), level);
      libmesh_assert_equal_to (elem->id(), id);
//#ifdef LIBMESH_ENABLE_UNIQUE_ID
      // No check for unqiue id sanity
//#endif
      libmesh_assert_equal_to (elem->processor_id(), processor_id);
      libmesh_assert_equal_to (elem->subdomain_id(), subdomain_id);
      libmesh_assert_equal_to (elem->type(), type);
      libmesh_assert_equal_to (elem->n_nodes(), n_nodes);

#ifndef NDEBUG
      // All our nodes should be correct
      for (unsigned int i=0; i != n_nodes; ++i)
        libmesh_assert(elem->node(i) ==
                       static_cast<dof_id_type>(*in++));
#else
      in += n_nodes;
#endif

#ifdef LIBMESH_ENABLE_AMR
      libmesh_assert_equal_to (elem->p_level(), p_level);
      libmesh_assert_equal_to (elem->refinement_flag(), refinement_flag);
      libmesh_assert_equal_to (elem->p_refinement_flag(), p_refinement_flag);

      libmesh_assert (!level || elem->parent() != NULL);
      libmesh_assert (!level || elem->parent()->id() == parent_id);
      libmesh_assert (!level || elem->parent()->child(which_child_am_i) == elem);
#endif

      // Our neighbor links should be "close to" correct - we may have
      // to update them, but we can check for some inconsistencies.
      for (unsigned int n=0; n != elem->n_neighbors(); ++n)
        {
          const dof_id_type neighbor_id =
            static_cast<dof_id_type>(*in++);

	  // If the sending processor sees a domain boundary here,
	  // we'd better agree.
          if (neighbor_id == DofObject::invalid_id)
            {
              libmesh_assert (!(elem->neighbor(n)));
              continue;
            }

	  // If the sending processor has a remote_elem neighbor here,
	  // then all we know is that we'd better *not* have a domain
	  // boundary.
          if (neighbor_id == remote_elem->id())
            {
              libmesh_assert(elem->neighbor(n));
              continue;
            }

          Elem *neigh = mesh->query_elem(neighbor_id);

          // The sending processor sees a neighbor here, so if we
          // don't have that neighboring element, then we'd better
          // have a remote_elem signifying that fact.
          if (!neigh)
            {
              libmesh_assert_equal_to (elem->neighbor(n), remote_elem);
              continue;
            }

          // The sending processor has a neighbor here, and we have
          // that element, but that does *NOT* mean we're already
	  // linking to it.  Perhaps we initially received both elem
	  // and neigh from processors on which their mutual link was
	  // remote?
          libmesh_assert(elem->neighbor(n) == neigh ||
			 elem->neighbor(n) == remote_elem);

	  // If the link was originally remote, we should update it,
	  // and make sure the appropriate parts of its family link
	  // back to us.
	  if (elem->neighbor(n) == remote_elem)
            {
              elem->set_neighbor(n, neigh);

              elem->make_links_to_me_local(n);
	    }
	}

      // FIXME: We should add some debug mode tests to ensure that the
      // encoded indexing and boundary conditions are consistent.
    }
  else
    {
      // We don't already have the element, so we need to create it.

      // Find the parent if necessary
      Elem *parent = NULL;
#ifdef LIBMESH_ENABLE_AMR
      // Find a child element's parent
      if (level > 0)
        {
	  // Note that we must be very careful to construct the send
	  // connectivity so that parents are encountered before
	  // children.  If we get here and can't find the parent that
	  // is a fatal error.
          parent = mesh->elem(parent_id);
        }
      // Or assert that the sending processor sees no parent
      else
        libmesh_assert_equal_to (parent_id, static_cast<dof_id_type>(-1));
#else
      // No non-level-0 elements without AMR
      libmesh_assert_equal_to (level, 0);
#endif

      elem = Elem::build(type,parent).release();
      libmesh_assert (elem);

#ifdef LIBMESH_ENABLE_AMR
      if (level != 0)
        {
          // Since this is a newly created element, the parent must
          // have previously thought of this child as a remote element.
          libmesh_assert_equal_to (parent->child(which_child_am_i), remote_elem);

          parent->add_child(elem, which_child_am_i);
        }

      // Assign the refinement flags and levels
      elem->set_p_level(p_level);
      elem->set_refinement_flag(refinement_flag);
      elem->set_p_refinement_flag(p_refinement_flag);
      libmesh_assert_equal_to (elem->level(), level);

      // If this element definitely should have children, assign
      // remote_elem to all of them for now, for consistency.  Later
      // unpacked elements may overwrite that.
      if (!elem->active())
        for (unsigned int c=0; c != elem->n_children(); ++c)
          elem->add_child(const_cast<RemoteElem*>(remote_elem), c);

#endif // LIBMESH_ENABLE_AMR

      // Assign the IDs
      elem->subdomain_id()  = subdomain_id;
      elem->processor_id()  = processor_id;
      elem->set_id()        = id;
#ifdef LIBMESH_ENABLE_UNIQUE_ID
      elem->set_unique_id() = unique_id;
#endif

      // Assign the connectivity
      libmesh_assert_equal_to (elem->n_nodes(), n_nodes);

      for (unsigned int n=0; n != n_nodes; n++)
        elem->set_node(n) =
          mesh->node_ptr
	    (static_cast<dof_id_type>(*in++));

      for (unsigned int n=0; n<elem->n_neighbors(); n++)
        {
          const dof_id_type neighbor_id =
            static_cast<dof_id_type>(*in++);

          if (neighbor_id == DofObject::invalid_id)
	    continue;

          // We may be unpacking an element that was a ghost element on the
          // sender, in which case the element's neighbors may not all be
          // known by the packed element.  We'll have to set such
          // neighbors to remote_elem ourselves and wait for a later
          // packed element to give us better information.
          if (neighbor_id == remote_elem->id())
            {
              elem->set_neighbor(n, const_cast<RemoteElem*>(remote_elem));
	      continue;
	    }

          // If we don't have the neighbor element, then it's a
          // remote_elem until we get it.
          Elem *neigh = mesh->query_elem(neighbor_id);
          if (!neigh)
            {
              elem->set_neighbor(n, const_cast<RemoteElem*>(remote_elem));
	      continue;
	    }

          // If we have the neighbor element, then link to it, and
          // make sure the appropriate parts of its family link back
          // to us.
          elem->set_neighbor(n, neigh);

          elem->make_links_to_me_local(n);
        }

      elem->unpack_indexing(in);
    }

  in += elem->packed_indexing_size();

  // If this is a coarse element,
  // add any element side boundary condition ids
  if (level == 0)
    for (unsigned int s = 0; s != elem->n_sides(); ++s)
      {
        const int num_bcs = *in++;
        libmesh_assert_greater_equal (num_bcs, 0);

        for(int bc_it=0; bc_it < num_bcs; bc_it++)
          mesh->boundary_info->add_side (elem, s, *in++);
      }

  // Return the new element
  *out = elem;
}



template <>
void unpack(std::vector<largest_id_type>::const_iterator in,
            Elem** out,
            ParallelMesh* mesh)
{
  unpack(in, out, static_cast<MeshBase*>(mesh));
}

} // namespace Parallel

#endif // LIBMESH_HAVE_MPI

} // namespace libMesh

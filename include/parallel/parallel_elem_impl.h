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

#ifndef LIBMESH_PARALLEL_ELEM_IMPL_H
#define LIBMESH_PARALLEL_ELEM_IMPL_H

// C++ includes

// Local includes
#include "libmesh/boundary_info.h"
#include "libmesh/distributed_mesh.h"
#include "libmesh/elem.h"
#include "libmesh/mesh_base.h"
#include "libmesh/parallel_elem.h"
#include "libmesh/parallel_mesh.h"
#include "libmesh/remote_elem.h"

// Helper functions in anonymous namespace

namespace
{
using namespace libMesh;

#ifdef LIBMESH_ENABLE_UNIQUE_ID
static const unsigned int header_size = 12;
#else
static const unsigned int header_size = 11;
#endif

#ifndef NDEBUG
// Currently this constant is only used for debugging.
static const largest_id_type elem_magic_header = 987654321;
#endif
}


namespace libMesh
{

namespace Parallel
{

namespace
{
template <typename, typename> struct PackedSizeShim;

template <typename RealType>
struct PackedSizeShim<const ElemTempl<RealType> *, typename std::vector<largest_id_type>::const_iterator>
{
  static unsigned int packed_size(typename std::vector<largest_id_type>::const_iterator in)
    {
#ifndef NDEBUG
      const largest_id_type packed_header = *in++;
      libmesh_assert_equal_to (packed_header, elem_magic_header);
#endif

      // int 0: level
      const unsigned int level =
        cast_int<unsigned int>(*in);

      // int 4: element type
      const int typeint = cast_int<int>(*(in+4));
      libmesh_assert_greater_equal (typeint, 0);
      libmesh_assert_less (typeint, INVALID_ELEM);
      const ElemType type =
        cast_int<ElemType>(typeint);

      const unsigned int n_nodes =
        Elem::type_to_n_nodes_map[type];

      const unsigned int n_sides =
        Elem::type_to_n_sides_map[type];

      const unsigned int n_edges =
        Elem::type_to_n_edges_map[type];

      const unsigned int pre_indexing_size =
        header_size + n_nodes + n_sides*2;

      const unsigned int indexing_size =
        DofObject::unpackable_indexing_size(in+pre_indexing_size);

      unsigned int total_packed_bc_data = 0;
      if (level == 0)
      {
        for (unsigned int s = 0; s != n_sides; ++s)
        {
          const int n_bcs = cast_int<int>
            (*(in + pre_indexing_size + indexing_size +
               total_packed_bc_data++));
          libmesh_assert_greater_equal (n_bcs, 0);
          total_packed_bc_data += n_bcs;
        }

        for (unsigned int e = 0; e != n_edges; ++e)
        {
          const int n_bcs = cast_int<int>
            (*(in + pre_indexing_size + indexing_size +
               total_packed_bc_data++));
          libmesh_assert_greater_equal (n_bcs, 0);
          total_packed_bc_data += n_bcs;
        }

        for (unsigned short sf=0; sf != 2; ++sf)
        {
          const int n_bcs = cast_int<int>
            (*(in + pre_indexing_size + indexing_size +
               total_packed_bc_data++));
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
};

template <typename RealType>
struct PackedSizeShim<const ElemTempl<RealType> *, typename std::vector<largest_id_type>::iterator>
{
  static unsigned int packed_size(typename std::vector<largest_id_type>::iterator in)
    {
      return PackedSizeShim<const ElemTempl<RealType> *, typename std::vector<largest_id_type>::const_iterator>::packed_size(in);
    }
};

template <typename, typename> struct PackableSizeShim;

template <typename RealType>
struct PackableSizeShim<const ElemTempl<RealType> *, MeshBaseTempl<RealType>>
{
  typedef ElemTempl<RealType> Elem;
  typedef MeshBaseTempl<RealType> MeshBase;

  static unsigned int packable_size(const Elem * const & elem,
                                    const MeshBase * mesh)
    {
      unsigned int total_packed_bcs = 0;
      const unsigned short n_sides = elem->n_sides();

      if (elem->level() == 0)
      {
        total_packed_bcs += n_sides;
        for (unsigned short s = 0; s != n_sides; ++s)
          total_packed_bcs +=
            mesh->get_boundary_info().n_boundary_ids(elem,s);

        const unsigned short n_edges = elem->n_edges();
        total_packed_bcs += n_edges;
        for (unsigned short e = 0; e != n_edges; ++e)
          total_packed_bcs +=
            mesh->get_boundary_info().n_edge_boundary_ids(elem,e);

        total_packed_bcs += 2; // shellfaces
        for (unsigned short sf=0; sf != 2; ++sf)
          total_packed_bcs +=
            mesh->get_boundary_info().n_shellface_boundary_ids(elem,sf);
      }

      return
#ifndef NDEBUG
        1 + // add an int for the magic header when testing
#endif
        header_size + elem->n_nodes() + n_sides*2 +
        elem->packed_indexing_size() + total_packed_bcs;
    }
};

template <typename RealType>
struct PackableSizeShim<const ElemTempl<RealType> *, DistributedMeshTempl<RealType>>
{
  typedef ElemTempl<RealType> Elem;
  typedef MeshBaseTempl<RealType> MeshBase;
  typedef DistributedMeshTempl<RealType> DistributedMesh;

  static unsigned int packable_size (const Elem * const & elem,
                                     const DistributedMeshTempl<RealType> * mesh)
    {
      return PackableSizeShim<const Elem *, MeshBase>::packable_size(
        elem,
        static_cast<const MeshBase *>(mesh));
    }
};

template <typename RealType>
struct PackableSizeShim<const ElemTempl<RealType> *, ParallelMeshTempl<RealType>>
{
  typedef ElemTempl<RealType> Elem;
  typedef MeshBaseTempl<RealType> MeshBase;
  typedef ParallelMeshTempl<RealType> ParallelMesh;

  static unsigned int packable_size (const Elem * const & elem,
                                     const ParallelMeshTempl<RealType> * mesh)
    {
      return PackableSizeShim<const Elem *, MeshBase>::packable_size(
        elem,
        static_cast<const MeshBase *>(mesh));
    }
};

template <typename, typename, typename>
struct PackShim;

template <typename RealType>
struct PackShim<const ElemTempl<RealType> *,
                std::back_insert_iterator<std::vector<largest_id_type>>,
                MeshBaseTempl<RealType>>
{
  typedef ElemTempl<RealType> Elem;
  typedef NodeTempl<RealType> Node;
  typedef MeshBaseTempl<RealType> MeshBase;

  static void pack (const Elem * const & elem,
                    std::back_insert_iterator<std::vector<largest_id_type>> data_out,
                    const MeshBase * mesh)
    {
      libmesh_assert(elem);

#ifndef NDEBUG
      *data_out++ = elem_magic_header;
#endif

#ifdef LIBMESH_ENABLE_AMR
      *data_out++ = (static_cast<largest_id_type>(elem->level()));
      *data_out++ = (static_cast<largest_id_type>(elem->p_level()));

      // Encode both the refinement flag and whether the element has
      // children together.  This coding is unambiguous because our
      // refinement state encoding starts at 0 and ends at
      // INVALID_REFINEMENTSTATE
      largest_id_type refinement_info =
        static_cast<largest_id_type>(elem->refinement_flag());
      if (elem->has_children())
        refinement_info +=
          static_cast<largest_id_type>(Elem::INVALID_REFINEMENTSTATE) + 1;
      *data_out++ = (refinement_info);

      *data_out++ = (static_cast<largest_id_type>(elem->p_refinement_flag()));
#else
      *data_out++ = (0);
      *data_out++ = (0);
      *data_out++ = (0);
      *data_out++ = (0);
#endif
      *data_out++ = (static_cast<largest_id_type>(elem->type()));
      *data_out++ = (elem->processor_id());
      *data_out++ = (elem->subdomain_id());
      *data_out++ = (elem->id());

#ifdef LIBMESH_ENABLE_UNIQUE_ID
      if (elem->valid_unique_id())
        *data_out++ = (static_cast<largest_id_type>(elem->unique_id()));
      else
        // OK to send invalid unique id, we must not own this DOF
        *data_out++ = (static_cast<largest_id_type>(DofObject::invalid_unique_id));
#endif

#ifdef LIBMESH_ENABLE_AMR
      // use parent_ID of invalid_id to indicate a level 0 element
      if (elem->level() == 0)
      {
        *data_out++ =(DofObject::invalid_id);
        *data_out++ =(DofObject::invalid_id);
      }
      else
      {
        *data_out++ =(elem->parent()->id());
        *data_out++ =(elem->parent()->which_child_am_i(elem));
      }
#else
      *data_out++ = (DofObject::invalid_id);
      *data_out++ = (DofObject::invalid_id);
#endif

      if ((elem->dim() < LIBMESH_DIM) &&
          elem->interior_parent())
        *data_out++ =(elem->interior_parent()->id());
      else
        *data_out++ =(DofObject::invalid_id);

      for (const Node & node : elem->node_ref_range())
        *data_out++ = node.id();

      // Add the id of and the side for any return link from each neighbor
      for (auto neigh : elem->neighbor_ptr_range())
      {
        if (neigh)
        {
          *data_out++ = (neigh->id());
          if (neigh == RemoteElem::get_instance())
            *data_out++ = (DofObject::invalid_id);
          else
            *data_out++ = neigh->which_neighbor_am_i(elem);
        }
        else
        {
          *data_out++ = (DofObject::invalid_id);
          *data_out++ = (DofObject::invalid_id);
        }
      }

      // Add any DofObject indices
      elem->pack_indexing(data_out);

      // If this is a coarse element,
      // Add any element side boundary condition ids
      if (elem->level() == 0)
      {
        std::vector<boundary_id_type> bcs;
        for (auto s : elem->side_index_range())
        {
          mesh->get_boundary_info().boundary_ids(elem, s, bcs);

          *data_out++ =(bcs.size());

          for (const auto & bid : bcs)
            *data_out++ = bid;
        }

        for (auto e : elem->edge_index_range())
        {
          mesh->get_boundary_info().edge_boundary_ids(elem, e, bcs);

          *data_out++ =(bcs.size());

          for (const auto & bid : bcs)
            *data_out++ = bid;
        }

        for (unsigned short sf=0; sf != 2; ++sf)
        {
          mesh->get_boundary_info().shellface_boundary_ids(elem, sf, bcs);

          *data_out++ =(bcs.size());

          for (const auto & bid : bcs)
            *data_out++ = bid;
        }
      }
    }
};

template <typename RealType>
struct PackShim<const ElemTempl<RealType> *,
                std::back_insert_iterator<std::vector<largest_id_type>>,
                DistributedMeshTempl<RealType>>
{
  typedef ElemTempl<RealType> Elem;
  typedef MeshBaseTempl<RealType> MeshBase;
  typedef DistributedMeshTempl<RealType> DistributedMesh;

  static void pack(const Elem * const & elem,
                   std::back_insert_iterator<std::vector<largest_id_type>> data_out,
                   const DistributedMesh * mesh)
    {
      PackShim<const Elem *,
               std::back_insert_iterator<std::vector<largest_id_type>>,
               MeshBase>::pack(elem, data_out, static_cast<const MeshBase*>(mesh));
    }
};

template <typename RealType>
struct PackShim<const ElemTempl<RealType> *,
                std::back_insert_iterator<std::vector<largest_id_type>>,
                ParallelMeshTempl<RealType>>
{
  typedef ElemTempl<RealType> Elem;
  typedef MeshBaseTempl<RealType> MeshBase;
  typedef ParallelMeshTempl<RealType> ParallelMesh;

  static void pack(const Elem * const & elem,
                   std::back_insert_iterator<std::vector<largest_id_type>> data_out,
                   const ParallelMesh * mesh)
    {
      PackShim<const Elem *,
               std::back_insert_iterator<std::vector<largest_id_type>>,
               MeshBase>::pack(elem, data_out, static_cast<const MeshBase*>(mesh));
    }
};

template <typename, typename, typename>
struct UnpackShim;

template <typename RealType>
struct UnpackShim<ElemTempl<RealType> *,
                  std::vector<largest_id_type>::const_iterator,
                  MeshBaseTempl<RealType>>
{
  typedef ElemTempl<RealType> Elem;
  typedef MeshBaseTempl<RealType> MeshBase;

  static Elem * unpack (std::vector<largest_id_type>::const_iterator in,
                        MeshBase * mesh)
    {
#ifndef NDEBUG
      const std::vector<largest_id_type>::const_iterator original_in = in;

      const largest_id_type incoming_header = *in++;
      libmesh_assert_equal_to (incoming_header, elem_magic_header);
#endif

      // int 0: level
      const unsigned int level =
        cast_int<unsigned int>(*in++);

#ifdef LIBMESH_ENABLE_AMR
      // int 1: p level
      const unsigned int p_level =
        cast_int<unsigned int>(*in++);

      // int 2: refinement flag and encoded has_children
      const int rflag = cast_int<int>(*in++);
      const int invalid_rflag =
        cast_int<int>(Elem::INVALID_REFINEMENTSTATE);
      libmesh_assert_greater_equal (rflag, 0);

      libmesh_assert_less (rflag, invalid_rflag*2+1);

      const bool has_children = (rflag > invalid_rflag);

      const typename Elem::RefinementState refinement_flag = has_children ?
        cast_int<typename Elem::RefinementState>(rflag - invalid_rflag - 1) :
        cast_int<typename Elem::RefinementState>(rflag);

      // int 3: p refinement flag
      const int pflag = cast_int<int>(*in++);
      libmesh_assert_greater_equal (pflag, 0);
      libmesh_assert_less (pflag, Elem::INVALID_REFINEMENTSTATE);
      const typename Elem::RefinementState p_refinement_flag =
        cast_int<typename Elem::RefinementState>(pflag);
#else
      in += 3;
#endif // LIBMESH_ENABLE_AMR

      // int 4: element type
      const int typeint = cast_int<int>(*in++);
      libmesh_assert_greater_equal (typeint, 0);
      libmesh_assert_less (typeint, INVALID_ELEM);
      const ElemType type =
        cast_int<ElemType>(typeint);

      const unsigned int n_nodes =
        Elem::type_to_n_nodes_map[type];

      // int 5: processor id
      const processor_id_type processor_id =
        cast_int<processor_id_type>(*in++);
      libmesh_assert (processor_id < mesh->n_processors() ||
                      processor_id == DofObject::invalid_processor_id);

      // int 6: subdomain id
      const subdomain_id_type subdomain_id =
        cast_int<subdomain_id_type>(*in++);

      // int 7: dof object id
      const dof_id_type id =
        cast_int<dof_id_type>(*in++);
      libmesh_assert_not_equal_to (id, DofObject::invalid_id);

#ifdef LIBMESH_ENABLE_UNIQUE_ID
      // int 8: dof object unique id
      const unique_id_type unique_id =
        cast_int<unique_id_type>(*in++);
#endif

#ifdef LIBMESH_ENABLE_AMR
      // int 9: parent dof object id.
      // Note: If level==0, then (*in) == invalid_id.  In
      // this case, the equality check in cast_int<unsigned>(*in) will
      // never succeed.  Therefore, we should only attempt the more
      // rigorous cast verification in cases where level != 0.
      const dof_id_type parent_id =
        (level == 0)
        ? static_cast<dof_id_type>(*in++)
        : cast_int<dof_id_type>(*in++);
      libmesh_assert (level == 0 || parent_id != DofObject::invalid_id);
      libmesh_assert (level != 0 || parent_id == DofObject::invalid_id);

      // int 10: local child id
      // Note: If level==0, then which_child_am_i is not valid, so don't
      // do the more rigorous cast verification.
      const unsigned int which_child_am_i =
        (level == 0)
        ? static_cast<unsigned int>(*in++)
        : cast_int<unsigned int>(*in++);
#else
      in += 2;
#endif // LIBMESH_ENABLE_AMR

      const dof_id_type interior_parent_id =
        static_cast<dof_id_type>(*in++);

      // Make sure we don't miscount above when adding the "magic" header
      // plus the real data header
      libmesh_assert_equal_to (in - original_in, header_size + 1);

      Elem * elem = mesh->query_elem_ptr(id);

      // if we already have this element, make sure its
      // properties match, and update any missing neighbor
      // links, but then go on
      if (elem)
      {
        libmesh_assert_equal_to (elem->level(), level);
        libmesh_assert_equal_to (elem->id(), id);
        //#ifdef LIBMESH_ENABLE_UNIQUE_ID
        // No check for unique id sanity
        //#endif
        libmesh_assert_equal_to (elem->processor_id(), processor_id);
        libmesh_assert_equal_to (elem->subdomain_id(), subdomain_id);
        libmesh_assert_equal_to (elem->type(), type);
        libmesh_assert_equal_to (elem->n_nodes(), n_nodes);

#ifndef NDEBUG
        // All our nodes should be correct
        for (unsigned int i=0; i != n_nodes; ++i)
          libmesh_assert(elem->node_id(i) ==
                         cast_int<dof_id_type>(*in++));
#else
        in += n_nodes;
#endif

#ifdef LIBMESH_ENABLE_AMR
        libmesh_assert_equal_to (elem->refinement_flag(), refinement_flag);
        libmesh_assert_equal_to (elem->has_children(), has_children);

#ifdef DEBUG
        if (elem->active())
        {
          libmesh_assert_equal_to (elem->p_level(), p_level);
          libmesh_assert_equal_to (elem->p_refinement_flag(), p_refinement_flag);
        }
#endif

        libmesh_assert (!level || elem->parent() != nullptr);
        libmesh_assert (!level || elem->parent()->id() == parent_id);
        libmesh_assert (!level || elem->parent()->child_ptr(which_child_am_i) == elem);
#endif
        // Our interior_parent link should be "close to" correct - we
        // may have to update it, but we can check for some
        // inconsistencies.
        {
          // If the sending processor sees no interior_parent here, we'd
          // better agree.
          if (interior_parent_id == DofObject::invalid_id)
          {
            if (elem->dim() < LIBMESH_DIM)
              libmesh_assert (!(elem->interior_parent()));
          }

          // If the sending processor has a remote_elem interior_parent,
          // then all we know is that we'd better have *some*
          // interior_parent
          else if (interior_parent_id == RemoteElem::get_instance()->id())
          {
            libmesh_assert(elem->interior_parent());
          }
          else
          {
            Elem * ip = mesh->query_elem_ptr(interior_parent_id);

            // The sending processor sees an interior parent here, so
            // if we don't have that interior element, then we'd
            // better have a remote_elem signifying that fact.
            if (!ip)
              libmesh_assert_equal_to (elem->interior_parent(), RemoteElem::get_instance());
            else
            {
              // The sending processor has an interior_parent here,
              // and we have that element, but that does *NOT* mean
              // we're already linking to it.  Perhaps we initially
              // received elem from a processor on which the
              // interior_parent link was remote?
              libmesh_assert(elem->interior_parent() == ip ||
                             elem->interior_parent() == RemoteElem::get_instance());

              // If the link was originally remote, update it
              if (elem->interior_parent() == RemoteElem::get_instance())
              {
                elem->set_interior_parent(ip);
              }
            }
          }
        }

        // Our neighbor links should be "close to" correct - we may have
        // to update a remote_elem link, and we can check for possible
        // inconsistencies along the way.
        //
        // For subactive elements, we don't bother keeping neighbor
        // links in good shape, so there's nothing we need to set or can
        // safely assert here.
        if (!elem->subactive())
          for (auto n : elem->side_index_range())
          {
            const dof_id_type neighbor_id =
              cast_int<dof_id_type>(*in++);

            const dof_id_type neighbor_side =
              cast_int<dof_id_type>(*in++);

            // If the sending processor sees a domain boundary here,
            // we'd better agree.
            if (neighbor_id == DofObject::invalid_id)
            {
              libmesh_assert (!(elem->neighbor_ptr(n)));
              continue;
            }

            // If the sending processor has a remote_elem neighbor here,
            // then all we know is that we'd better *not* have a domain
            // boundary.
            if (neighbor_id == RemoteElem::get_instance()->id())
            {
              libmesh_assert(elem->neighbor_ptr(n));
              continue;
            }

            Elem * neigh = mesh->query_elem_ptr(neighbor_id);

            // The sending processor sees a neighbor here, so if we
            // don't have that neighboring element, then we'd better
            // have a remote_elem signifying that fact.
            if (!neigh)
            {
              libmesh_assert_equal_to (elem->neighbor_ptr(n), RemoteElem::get_instance());
              continue;
            }

            // The sending processor has a neighbor here, and we have
            // that element, but that does *NOT* mean we're already
            // linking to it.  Perhaps we initially received both elem
            // and neigh from processors on which their mutual link was
            // remote?
            libmesh_assert(elem->neighbor_ptr(n) == neigh ||
                           elem->neighbor_ptr(n) == RemoteElem::get_instance());

            // If the link was originally remote, we should update it,
            // and make sure the appropriate parts of its family link
            // back to us.
            if (elem->neighbor_ptr(n) == RemoteElem::get_instance())
            {
              elem->set_neighbor(n, neigh);

              elem->make_links_to_me_local(n, neighbor_side);
            }
            else
              libmesh_assert(neigh->level() < elem->level() ||
                             neigh->neighbor_ptr(neighbor_side) == elem);
          }

        // Our p level and refinement flags should be "close to" correct
        // if we're not an active element - we might have a p level
        // increased or decreased by changes in remote_elem children.
        //
        // But if we have remote_elem children, then we shouldn't be
        // doing a projection on this inactive element on this
        // processor, so we won't need correct p settings.  Couldn't
        // hurt to update, though.
#ifdef LIBMESH_ENABLE_AMR
        if (elem->processor_id() != mesh->processor_id())
        {
          elem->hack_p_level(p_level);
          elem->set_p_refinement_flag(p_refinement_flag);
        }
#endif // LIBMESH_ENABLE_AMR

        // FIXME: We should add some debug mode tests to ensure that the
        // encoded indexing and boundary conditions are consistent.
      }
      else
      {
        // We don't already have the element, so we need to create it.

        // Find the parent if necessary
        Elem * parent = nullptr;
#ifdef LIBMESH_ENABLE_AMR
        // Find a child element's parent
        if (level > 0)
        {
          // Note that we must be very careful to construct the send
          // connectivity so that parents are encountered before
          // children.  If we get here and can't find the parent that
          // is a fatal error.
          parent = mesh->elem_ptr(parent_id);
        }
        // Or assert that the sending processor sees no parent
        else
          libmesh_assert_equal_to (parent_id, DofObject::invalid_id);
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
          libmesh_assert_equal_to (parent->child_ptr(which_child_am_i), RemoteElem::get_instance());

          parent->add_child(elem, which_child_am_i);
        }

        // Assign the refinement flags and levels
        elem->set_p_level(p_level);
        elem->set_refinement_flag(refinement_flag);
        elem->set_p_refinement_flag(p_refinement_flag);
        libmesh_assert_equal_to (elem->level(), level);

        // If this element should have children, assign remote_elem to
        // all of them for now, for consistency.  Later unpacked
        // elements may overwrite that.
        if (has_children)
        {
          const unsigned int nc = elem->n_children();
          for (unsigned int c=0; c != nc; ++c)
            elem->add_child(const_cast<RemoteElem *>(RemoteElem::get_instance()), c);
        }

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
            (cast_int<dof_id_type>(*in++));

        // Set interior_parent if found
        {
          // We may be unpacking an element that was a ghost element on the
          // sender, in which case the element's interior_parent may not be
          // known by the packed element.  We'll have to set such
          // interior_parents to remote_elem ourselves and wait for a
          // later packed element to give us better information.
          if (interior_parent_id == RemoteElem::get_instance()->id())
          {
            elem->set_interior_parent
              (const_cast<RemoteElem *>(RemoteElem::get_instance()));
          }
          else if (interior_parent_id != DofObject::invalid_id)
          {
            // If we don't have the interior parent element, then it's
            // a remote_elem until we get it.
            Elem * ip = mesh->query_elem_ptr(interior_parent_id);
            if (!ip )
              elem->set_interior_parent
                (const_cast<RemoteElem *>(RemoteElem::get_instance()));
            else
              elem->set_interior_parent(ip);
          }
        }

        for (auto n : elem->side_index_range())
        {
          const dof_id_type neighbor_id =
            cast_int<dof_id_type>(*in++);

          const dof_id_type neighbor_side =
            cast_int<dof_id_type>(*in++);

          if (neighbor_id == DofObject::invalid_id)
            continue;

          // We may be unpacking an element that was a ghost element on the
          // sender, in which case the element's neighbors may not all be
          // known by the packed element.  We'll have to set such
          // neighbors to remote_elem ourselves and wait for a later
          // packed element to give us better information.
          if (neighbor_id == RemoteElem::get_instance()->id())
          {
            elem->set_neighbor(n, const_cast<RemoteElem *>(RemoteElem::get_instance()));
            continue;
          }

          // If we don't have the neighbor element, then it's a
          // remote_elem until we get it.
          Elem * neigh = mesh->query_elem_ptr(neighbor_id);
          if (!neigh)
          {
            elem->set_neighbor(n, const_cast<RemoteElem *>(RemoteElem::get_instance()));
            continue;
          }

          // If we have the neighbor element, then link to it, and
          // make sure any appropriate parts of its family link back
          // to us.
          elem->set_neighbor(n, neigh);

          elem->make_links_to_me_local(n, neighbor_side);
        }

        elem->unpack_indexing(in);
      }

      in += elem->packed_indexing_size();

      // If this is a coarse element,
      // add any element side or edge boundary condition ids
      if (level == 0)
      {
        for (auto s : elem->side_index_range())
        {
          const boundary_id_type num_bcs =
            cast_int<boundary_id_type>(*in++);

          for (boundary_id_type bc_it=0; bc_it < num_bcs; bc_it++)
            mesh->get_boundary_info().add_side
              (elem, s, cast_int<boundary_id_type>(*in++));
        }

        for (auto e : elem->edge_index_range())
        {
          const boundary_id_type num_bcs =
            cast_int<boundary_id_type>(*in++);

          for (boundary_id_type bc_it=0; bc_it < num_bcs; bc_it++)
            mesh->get_boundary_info().add_edge
              (elem, e, cast_int<boundary_id_type>(*in++));
        }

        for (unsigned short sf=0; sf != 2; ++sf)
        {
          const boundary_id_type num_bcs =
            cast_int<boundary_id_type>(*in++);

          for (boundary_id_type bc_it=0; bc_it < num_bcs; bc_it++)
            mesh->get_boundary_info().add_shellface
              (elem, sf, cast_int<boundary_id_type>(*in++));
        }
      }

      // Return the new element
      return elem;
    }
};

template <typename RealType>
struct UnpackShim<ElemTempl<RealType> *,
                  std::vector<largest_id_type>::const_iterator,
                  DistributedMeshTempl<RealType>>
{
  typedef ElemTempl<RealType> Elem;
  typedef DistributedMeshTempl<RealType> DistributedMesh;
  typedef MeshBaseTempl<RealType> MeshBase;

  static Elem * unpack (std::vector<largest_id_type>::const_iterator in,
                        DistributedMesh * mesh)
    {
      return UnpackShim<Elem *,
                        std::vector<largest_id_type>::const_iterator,
                        MeshBase>::unpack(in, static_cast<MeshBase*>(mesh));
    }
};

template <typename RealType>
struct UnpackShim<ElemTempl<RealType> *,
                  std::vector<largest_id_type>::const_iterator,
                  ParallelMeshTempl<RealType>>
{
  typedef ElemTempl<RealType> Elem;
  typedef ParallelMeshTempl<RealType> ParallelMesh;
  typedef MeshBaseTempl<RealType> MeshBase;

  static Elem * unpack (std::vector<largest_id_type>::const_iterator in,
                        ParallelMesh * mesh)
    {
      return UnpackShim<Elem *,
                        std::vector<largest_id_type>::const_iterator,
                        MeshBase>::unpack(in, static_cast<MeshBase*>(mesh));
    }
};

}

template <typename RealType>
template <typename BufferIter>
unsigned int
Packing<const ElemTempl<RealType> *>::packed_size (BufferIter in)
{
  return PackedSizeShim<const ElemTempl<RealType> *, BufferIter>::packed_size(in);
}

template <typename RealType>
template <typename Context>
unsigned int
Packing<const ElemTempl<RealType> *>::packable_size (const Elem * const & elem,
                                                     const Context * ctx)
{
  return PackableSizeShim<const ElemTempl<RealType> *, Context>::packable_size(elem, ctx);
}


template <typename RealType>
template <typename OutputIter, typename Context>
void
Packing<const ElemTempl<RealType> *>::pack (const Elem * const & elem,
                                            OutputIter data_out,
                                            const Context * context)
{
  PackShim<const Elem *, OutputIter, Context>::pack(elem, data_out, context);
}


// FIXME - this needs serious work to be 64-bit compatible
template <typename RealType>
template <typename BufferIter, typename Context>
ElemTempl<RealType> *
Packing<ElemTempl<RealType> *>::unpack (BufferIter in,
                                        Context * ctx)
{
  return UnpackShim<Elem *, BufferIter, Context>::unpack(in, ctx);
}

} // namespace Parallel

} // namespace libMesh

#endif // LIBMESH_PARALLEL_ELEM_IMPL_H

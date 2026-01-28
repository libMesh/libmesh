// The libMesh Finite Element Library.
// Copyright (C) 2002-2025 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



// Local includes
#include "libmesh/boundary_info.h"
#include "libmesh/ghosting_functor.h"
#include "libmesh/ghost_point_neighbors.h"
#include "libmesh/unstructured_mesh.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/elem.h"
#include "libmesh/mesh_tools.h" // For n_levels
#include "libmesh/parallel.h"
#include "libmesh/remote_elem.h"
#include "libmesh/namebased_io.h"
#include "libmesh/partitioner.h"
#include "libmesh/enum_order.h"
#include "libmesh/mesh_communication.h"
#include "libmesh/enum_to_string.h"
#include "libmesh/mesh_serializer.h"
#include "libmesh/utility.h"

#ifdef LIBMESH_HAVE_NANOFLANN
#include "libmesh/nanoflann.hpp"
#endif

// C++ includes
#include <algorithm> // std::all_of
#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>
#include <unordered_map>

// for disjoint neighbors
#include "libmesh/periodic_boundaries.h"
#include "libmesh/periodic_boundary.h"

namespace {

using namespace libMesh;

// Helper functions for all_second_order, all_complete_order

std::map<std::vector<dof_id_type>, Node *>::iterator
map_hi_order_node(unsigned int hon,
                  const Elem & hi_elem,
                  std::map<std::vector<dof_id_type>, Node *> & adj_vertices_to_ho_nodes)
{
  /*
   * form a vector that will hold the node id's of
   * the vertices that are adjacent to the nth
   * higher-order node.
   */
  const unsigned int n_adjacent_vertices =
    hi_elem.n_second_order_adjacent_vertices(hon);

  std::vector<dof_id_type> adjacent_vertices_ids(n_adjacent_vertices);

  for (unsigned int v=0; v<n_adjacent_vertices; v++)
    adjacent_vertices_ids[v] =
      hi_elem.node_id( hi_elem.second_order_adjacent_vertex(hon,v) );

  /*
   * \p adjacent_vertices_ids is now in order of the current
   * side.  sort it, so that comparisons  with the
   * \p adjacent_vertices_ids created through other elements'
   * sides can match
   */
  std::sort(adjacent_vertices_ids.begin(),
            adjacent_vertices_ids.end());

  // Does this set of vertices already have a mid-node added?  If not
  // we'll want to add it.
  return adj_vertices_to_ho_nodes.try_emplace(adjacent_vertices_ids, nullptr).first;
}

void transfer_elem(Elem & lo_elem,
                   std::unique_ptr<Elem> hi_elem,
#ifdef LIBMESH_ENABLE_UNIQUE_ID
                   unique_id_type max_unique_id,
                   unique_id_type max_new_nodes_per_elem,
#endif
                   UnstructuredMesh & mesh,
                   std::map<std::vector<dof_id_type>, Node *> & adj_vertices_to_ho_nodes,
                   std::unordered_map<Elem *, std::vector<Elem *>> & exterior_children_of)
{
  libmesh_assert_equal_to (lo_elem.n_vertices(), hi_elem->n_vertices());

  const processor_id_type my_pid = mesh.processor_id();
  const processor_id_type lo_pid = lo_elem.processor_id();

  /*
   * Now handle the additional higher-order nodes.  This
   * is simply handled through a map that remembers
   * the already-added nodes.  This map maps the global
   * ids of the vertices (that uniquely define this
   * higher-order node) to the new node.
   * Notation: hon = high-order node
   */
  const unsigned int hon_begin = lo_elem.n_nodes();
  const unsigned int hon_end   = hi_elem->n_nodes();

  for (unsigned int hon=hon_begin; hon<hon_end; hon++)
    {
      auto pos = map_hi_order_node(hon, *hi_elem, adj_vertices_to_ho_nodes);

      // no, not added yet
      if (!pos->second)
        {
          const auto & adjacent_vertices_ids = pos->first;

          /*
           * for this set of vertices, there is no
           * second_order node yet.  Add it.
           *
           * compute the location of the new node as
           * the average over the adjacent vertices.
           */
          Point new_location = 0;
          for (dof_id_type vertex_id : adjacent_vertices_ids)
            new_location += mesh.point(vertex_id);

          new_location /= static_cast<Real>(adjacent_vertices_ids.size());

          /* Add the new point to the mesh.
           *
           * If we are on a serialized mesh, then we're doing this
           * all in sync, and the node processor_id will be
           * consistent between processors.
           *
           * If we are on a distributed mesh, we can fix
           * inconsistent processor ids later, but only if every
           * processor gives new nodes a *locally* consistent
           * processor id, so we'll give the new node the
           * processor id of an adjacent element for now and then
           * we'll update that later if appropriate.
           */
          Node * hi_node = mesh.add_point
            (new_location, DofObject::invalid_id, lo_pid);

          /* Come up with a unique unique_id for a potentially new
           * node.  On a distributed mesh we don't yet know what
           * processor_id will definitely own it, so we can't let
           * the pid determine the unique_id.  But we're not
           * adding unpartitioned nodes in sync, so we can't let
           * the mesh autodetermine a unique_id for a new
           * unpartitioned node either.  So we have to pick unique
           * unique_id values manually.
           *
           * We don't have to pick the *same* unique_id value as
           * will be picked on other processors, though; we'll
           * sync up each node later.  We just need to make sure
           * we don't duplicate any unique_id that might be chosen
           * by the same process elsewhere.
           */
#ifdef LIBMESH_ENABLE_UNIQUE_ID
          unique_id_type new_unique_id = max_unique_id +
            max_new_nodes_per_elem * lo_elem.id() +
            hon - hon_begin;

          hi_node->set_unique_id(new_unique_id);
#endif

          /*
           * insert the new node with its defining vertex
           * set into the map, and relocate pos to this
           * new entry, so that the hi_elem can use
           * \p pos for inserting the node
           */
          pos->second = hi_node;

          hi_elem->set_node(hon, hi_node);
        }
      // yes, already added.
      else
        {
          Node * hi_node = pos->second;
          libmesh_assert(hi_node);
          libmesh_assert_equal_to(mesh.node_ptr(hi_node->id()), hi_node);

          hi_elem->set_node(hon, hi_node);

          // We need to ensure that the processor who should own a
          // node *knows* they own the node.  And because
          // Node::choose_processor_id() may depend on Node id,
          // which may not yet be authoritative, we still have to
          // use a dumb-but-id-independent partitioning heuristic.
          processor_id_type chosen_pid =
            std::min (hi_node->processor_id(), lo_pid);

          // Plus, if we just discovered that we own this node,
          // then on a distributed mesh we need to make sure to
          // give it a valid id, not just a placeholder id!
          if (!mesh.is_replicated() &&
              hi_node->processor_id() != my_pid &&
              chosen_pid == my_pid)
            mesh.own_node(*hi_node);

          hi_node->processor_id() = chosen_pid;
        }
    }

  /*
   * find_neighbors relies on remote_elem neighbor links being
   * properly maintained.  Our own code here relies on ordinary
   * neighbor links being properly maintained, so let's just keep
   * everything up to date.
   */
  for (auto s : lo_elem.side_index_range())
    {
      Elem * neigh = lo_elem.neighbor_ptr(s);
      if (!neigh)
        continue;

      if (neigh != remote_elem)
        {
          // We don't support AMR even outside our own range yet.
          libmesh_assert_equal_to (neigh->level(), 0);

          const unsigned int ns = neigh->which_neighbor_am_i(&lo_elem);
          libmesh_assert_not_equal_to(ns, libMesh::invalid_uint);

          neigh->set_neighbor(ns, hi_elem.get());
        }

      hi_elem->set_neighbor(s, neigh);
    }

  /**
   * If the old element has an interior_parent(), transfer it to the
   * new element ... and if the interior_parent itself might be
   * getting upgraded, make sure we later consider the new element to
   * be its exterior child, not the old element.
   */
  Elem * interior_p = lo_elem.interior_parent();
  if (interior_p)
    hi_elem->set_interior_parent(interior_p);

  if (auto parent_exterior_it = exterior_children_of.find(interior_p);
      parent_exterior_it != exterior_children_of.end())
    {
      auto & exteriors = parent_exterior_it->second;
      for (std::size_t i : index_range(exteriors))
        if (exteriors[i] == &lo_elem)
          {
            exteriors[i] = hi_elem.get();
            break;
          }
    }

  /**
   * If we had interior_parent() links to the old element, transfer
   * them to the new element.
   */
  if (auto exterior_it = exterior_children_of.find(&lo_elem);
      exterior_it != exterior_children_of.end())
    {
      for (Elem * exterior_elem : exterior_it->second)
        {
          libmesh_assert(exterior_elem->interior_parent() == &lo_elem);
          exterior_elem->set_interior_parent(hi_elem.get());
        }
    }

  /**
   * If the old element had any boundary conditions they
   * should be transferred to the second-order element.  The old
   * boundary conditions will be removed from the BoundaryInfo
   * data structure by insert_elem.
   *
   * Also, prepare_for_use() will reconstruct most of our neighbor
   * links, but if we have any remote_elem links in a distributed
   * mesh, they need to be preserved.  We do that in the same loop
   * here.
   */
  mesh.get_boundary_info().copy_boundary_ids
    (mesh.get_boundary_info(), &lo_elem, hi_elem.get());

  /*
   * The new second-order element is ready.
   * Inserting it into the mesh will replace and delete
   * the first-order element.
   */
  hi_elem->set_id(lo_elem.id());
#ifdef LIBMESH_ENABLE_UNIQUE_ID
  hi_elem->set_unique_id(lo_elem.unique_id());
#endif

  const unsigned int nei = lo_elem.n_extra_integers();
  hi_elem->add_extra_integers(nei);
  for (unsigned int i=0; i != nei; ++i)
    hi_elem->set_extra_integer(i, lo_elem.get_extra_integer(i));

  hi_elem->inherit_data_from(lo_elem);

  mesh.insert_elem(std::move(hi_elem));
}


template <typename ElemTypeConverter>
void
all_increased_order_range (UnstructuredMesh & mesh,
                           const SimpleRange<MeshBase::element_iterator> & range,
                           const unsigned int max_new_nodes_per_elem,
                           const ElemTypeConverter & elem_type_converter)
{
  // This function must be run on all processors at once
  timpi_parallel_only(mesh.comm());

  /*
   * The maximum number of new higher-order nodes we might be adding,
   * for use when picking unique unique_id values later. This variable
   * is not used unless unique ids are enabled, so libmesh_ignore() it
   * to avoid warnings in that case.
   */
  libmesh_ignore(max_new_nodes_per_elem);

  /*
   * The mesh should at least be consistent enough for us to add new
   * nodes consistently.
   */
  libmesh_assert(mesh.comm().verify(mesh.n_elem()));
  libmesh_assert(mesh.comm().verify(mesh.max_elem_id()));

  /*
   * If the mesh is empty then we have nothing to do
   */
  if (!mesh.n_elem())
    return;

  // If every element in the range _on every proc_ is already of the
  // requested higher order then we have nothing to do. However, if
  // any proc has some lower-order elements in the range, then _all_
  // processors need to continue this function because it is
  // parallel_only().
  //
  // Note: std::all_of() returns true for an empty range, which can
  // happen for example in the DistributedMesh case when there are
  // more processors than elements. In the case of an empty range we
  // therefore set already_second_order to true on that proc.
  auto is_higher_order = [&elem_type_converter](const Elem * elem) {
    ElemType old_type = elem->type();
    ElemType new_type = elem_type_converter(old_type);
    return old_type == new_type;
  };

  bool already_higher_order =
    std::all_of(range.begin(), range.end(), is_higher_order);

  // Check with other processors and possibly return early
  mesh.comm().min(already_higher_order);
  if (already_higher_order)
    return;

  /*
   * this map helps in identifying higher order
   * nodes.  Namely, a higher-order node:
   * - edge node
   * - face node
   * - bubble node
   * is uniquely defined through a set of adjacent
   * vertices.  This set of adjacent vertices is
   * used to identify already added higher-order
   * nodes.  We are safe to use node id's since we
   * make sure that these are correctly numbered.
   *
   * We lazily use an ordered map here to avoid having to implement a
   * good hash for vector<dof_id_type>
   */
  std::map<std::vector<dof_id_type>, Node *> adj_vertices_to_ho_nodes;

  /*
   * This map helps us reset any interior_parent() values from the
   * lower order element to its higher order replacement.  Unlike with
   * neighbor pointers, we don't have backlinks here, so we have to
   * iterate over the mesh to track forward links.
   */
  std::unordered_map<Elem *, std::vector<Elem *>> exterior_children_of;

  /*
   * max_new_nodes_per_elem is the maximum number of new higher order
   * nodes we might be adding, for use when picking unique unique_id
   * values later. This variable is not used unless unique ids are
   * enabled.
   */
#ifdef LIBMESH_ENABLE_UNIQUE_ID
  unique_id_type max_unique_id = mesh.parallel_max_unique_id();
#endif

  /**
   * On distributed meshes we currently only support unpartitioned
   * meshes (where we'll add every node in sync) or
   * completely-partitioned meshes (where we'll sync nodes later);
   * let's keep track to make sure we're not in any in-between state.
   */
  dof_id_type n_unpartitioned_elem = 0,
              n_partitioned_elem = 0;

  /**
   * Loop over the elements in the given range.  If any are
   * already at higher than first-order, track their higher-order
   * nodes in case we need them for neighboring elements later.
   *
   * In this way we can use this method to "fix up" a mesh which has
   * otherwise inconsistent neighbor pairs of lower and higher order
   * geometric elements.
   *
   * If any elements are not at the desired order yet, we need to
   * check their neighbors and even their edge neighbors for higher
   * order; we may need to share elements with a neighbor not in the
   * range.
   */
  auto track_if_necessary = [&adj_vertices_to_ho_nodes,
                             &exterior_children_of,
                             &elem_type_converter](Elem * elem) {
    if (elem && elem != remote_elem)
      {
        if (elem->default_order() != FIRST)
          for (unsigned int hon : make_range(elem->n_vertices(), elem->n_nodes()))
            {
              auto pos = map_hi_order_node(hon, *elem, adj_vertices_to_ho_nodes);
              pos->second = elem->node_ptr(hon);
            }

        const ElemType old_type = elem->type();
        const ElemType new_type = elem_type_converter(old_type);
        if (old_type != new_type)
          exterior_children_of.emplace(elem, std::vector<Elem *>());
      }
  };

  // If we're in the common case then just track everything; otherwise
  // find point neighbors to track
  if (range.begin() == mesh.elements_begin() &&
      range.end() == mesh.elements_end())
    {
      for (auto & elem : range)
        track_if_necessary(elem);
    }
  else
    {
      GhostingFunctor::map_type point_neighbor_elements;

      GhostPointNeighbors point_neighbor_finder(mesh);
      point_neighbor_finder(range.begin(), range.end(),
                            mesh.n_processors(),
                            point_neighbor_elements);

      for (auto & [elem, coupling_map] : point_neighbor_elements)
        {
          libmesh_ignore(coupling_map);
          track_if_necessary(const_cast<Elem *>(elem));
        }
    }

  /**
   * Loop over all mesh elements to look for interior_parent links we
   * need to upgrade later.
   */
  for (auto & elem : mesh.element_ptr_range())
    if (auto exterior_map_it = exterior_children_of.find(elem->interior_parent());
        exterior_map_it != exterior_children_of.end())
      exterior_map_it->second.push_back(elem);

  /**
   * Loop over the low-ordered elements in the _elements vector.
   * First make sure they _are_ indeed low-order, and then replace
   * them with an equivalent second-order element.  Don't
   * forget to delete the low-order element, or else it will leak!
   */
  for (auto & lo_elem : range)
    {
      // Now we can skip the elements in the range that are already
      // higher-order.
      const ElemType old_type = lo_elem->type();
      const ElemType new_type = elem_type_converter(old_type);

      if (old_type == new_type)
        continue;

      // this does _not_ work for refined elements
      libmesh_assert_equal_to (lo_elem->level(), 0);

      if (lo_elem->processor_id() == DofObject::invalid_processor_id)
        ++n_unpartitioned_elem;
      else
        ++n_partitioned_elem;

      /*
       * Build the higher-order equivalent; add to
       * the new_elements list.
       */
      auto ho_elem = Elem::build (new_type);

      libmesh_assert_equal_to (lo_elem->n_vertices(), ho_elem->n_vertices());

      /*
       * By definition the initial nodes of the lower and higher order
       * element are identically numbered.  Transfer these.
       */
      for (unsigned int v=0, lnn=lo_elem->n_nodes(); v < lnn; v++)
        ho_elem->set_node(v, lo_elem->node_ptr(v));

      transfer_elem(*lo_elem, std::move(ho_elem),
#ifdef LIBMESH_ENABLE_UNIQUE_ID
                    max_unique_id, max_new_nodes_per_elem,
#endif
                    mesh, adj_vertices_to_ho_nodes,
                    exterior_children_of);
    } // end for (auto & lo_elem : range)

  // we can clear the map at this point.
  adj_vertices_to_ho_nodes.clear();

#ifdef LIBMESH_ENABLE_UNIQUE_ID
  const unique_id_type new_max_unique_id = max_unique_id +
    max_new_nodes_per_elem * mesh.n_elem();
  mesh.set_next_unique_id(new_max_unique_id);
#endif

  // On a DistributedMesh our ghost node processor ids may be bad,
  // the ids of nodes touching remote elements may be inconsistent,
  // unique_ids of newly added non-local nodes remain unset, and our
  // partitioning of new nodes may not be well balanced.
  //
  // make_nodes_parallel_consistent() will fix all this.
  if (!mesh.is_replicated())
    {
      dof_id_type max_unpartitioned_elem = n_unpartitioned_elem;
      mesh.comm().max(max_unpartitioned_elem);
      if (max_unpartitioned_elem)
        {
          // We'd better be effectively serialized here.  In theory we
          // could support more complicated cases but in practice we
          // only support "completely partitioned" and/or "serialized"
          if (!mesh.comm().verify(n_unpartitioned_elem) ||
              !mesh.comm().verify(n_partitioned_elem) ||
              !mesh.is_serial())
            libmesh_not_implemented();
        }
      else
        {
          MeshCommunication().make_nodes_parallel_consistent (mesh);
        }
    }

  // renumber nodes, repartition nodes, etc.  We may no longer need a
  // find_neighbors() here since we're keeping neighbor links intact
  // ourselves, *except* that if we're not already prepared we may
  // have user code that was expecting this call to prepare neighbors.
  const bool old_find_neighbors = mesh.allow_find_neighbors();
  if (mesh.is_prepared())
    mesh.allow_find_neighbors(false);
  mesh.prepare_for_use();
  mesh.allow_find_neighbors(old_find_neighbors);
}


} // anonymous namespace


namespace libMesh
{

// This class adapts a vector of Nodes (represented by a pair of a Point and a dof_id_type)
// for use in a nanoflann KD-Tree

class VectorOfNodesAdaptor
{
private:
  const std::vector<std::pair<Point, dof_id_type>> _nodes;

public:
  VectorOfNodesAdaptor(const std::vector<std::pair<Point, dof_id_type>> & nodes) :
    _nodes(nodes)
  {}

  /**
   * Must return the number of data points
   */
  inline size_t kdtree_get_point_count() const { return _nodes.size(); }

  /**
   * \returns The dim'th component of the idx'th point in the class:
   * Since this is inlined and the "dim" argument is typically an immediate value, the
   *  "if's" are actually solved at compile time.
   */
  inline Real kdtree_get_pt(const size_t idx, int dim) const
    {
      libmesh_assert_less (idx, _nodes.size());
      libmesh_assert_less (dim, 3);

      const Point & p(_nodes[idx].first);

      if (dim==0) return p(0);
      if (dim==1) return p(1);
      return p(2);
    }

  /*
   * Optional bounding-box computation
   */
  template <class BBOX>
  bool kdtree_get_bbox(BBOX & /* bb */) const { return false; }
};


// ------------------------------------------------------------
// UnstructuredMesh class member functions
UnstructuredMesh::UnstructuredMesh (const Parallel::Communicator & comm_in,
                                    unsigned char d) :
  MeshBase (comm_in,d)
{
  libmesh_assert (libMesh::initialized());
}



UnstructuredMesh::UnstructuredMesh (const MeshBase & other_mesh) :
  MeshBase (other_mesh)
{
  libmesh_assert (libMesh::initialized());
}



void UnstructuredMesh::copy_nodes_and_elements(const MeshBase & other_mesh,
                                               const bool skip_find_neighbors,
                                               dof_id_type element_id_offset,
                                               dof_id_type node_id_offset,
                                               unique_id_type
#ifdef LIBMESH_ENABLE_UNIQUE_ID
                                                 unique_id_offset
#endif
                                               ,
                                               std::unordered_map<subdomain_id_type, subdomain_id_type> *
                                                 id_remapping,
                                               const bool skip_preparation)
{
  LOG_SCOPE("copy_nodes_and_elements()", "UnstructuredMesh");

  // If we're asked to skip all preparation, we should be skipping
  // find_neighbors specifically.
  libmesh_assert(!skip_preparation || skip_find_neighbors);

  std::pair<std::vector<unsigned int>, std::vector<unsigned int>>
    extra_int_maps = this->merge_extra_integer_names(other_mesh);

  const unsigned int n_old_node_ints = extra_int_maps.second.size(),
                     n_new_node_ints = _node_integer_names.size(),
                     n_old_elem_ints = extra_int_maps.first.size(),
                     n_new_elem_ints = _elem_integer_names.size();

  // If we are partitioned into fewer parts than the incoming mesh has
  // processors to handle, then we need to "wrap" the other Mesh's
  // processor ids to fit within our range. This can happen, for
  // example, while stitching meshes with small numbers of elements in
  // parallel...
  bool wrap_proc_ids = (this->n_processors() <
                        other_mesh.n_partitions());

  // We're assuming the other mesh has proper element number ordering,
  // so that we add parents before their children, and that the other
  // mesh is consistently partitioned.  We're not assuming that node
  // proc ids are topologically consistent, so we don't just
  // libmesh_assert_valid_procids.
#ifdef DEBUG
  MeshTools::libmesh_assert_valid_amr_elem_ids(other_mesh);
  MeshTools::libmesh_assert_parallel_consistent_procids<Node>(other_mesh);
#endif

  //Copy in Nodes
  {
    //Preallocate Memory if necessary
    this->reserve_nodes(other_mesh.n_nodes());

    for (const auto & oldn : other_mesh.node_ptr_range())
      {
        processor_id_type added_pid = cast_int<processor_id_type>
          (wrap_proc_ids ? oldn->processor_id() % this->n_processors() : oldn->processor_id());

        // Add new nodes in old node Point locations
        Node * newn =
          this->add_point(*oldn,
                          oldn->id() + node_id_offset,
                          added_pid);

        newn->add_extra_integers(n_new_node_ints);
        for (unsigned int i = 0; i != n_old_node_ints; ++i)
          newn->set_extra_integer(extra_int_maps.second[i],
                                  oldn->get_extra_integer(i));

#ifdef LIBMESH_ENABLE_UNIQUE_ID
        newn->set_unique_id(oldn->unique_id() + unique_id_offset);
#endif
      }
  }

  //Copy in Elements
  {
    //Preallocate Memory if necessary
    this->reserve_elem(other_mesh.n_elem());

    // Declare a map linking old and new elements, needed to copy the neighbor lists
    typedef std::unordered_map<const Elem *, Elem *> map_type;
    map_type old_elems_to_new_elems, ip_map;

    // Loop over the elements
    for (const auto & old : other_mesh.element_ptr_range())
      {
        // Build a new element
        Elem * newparent = old->parent() ?
          this->elem_ptr(old->parent()->id() + element_id_offset) :
          nullptr;
        auto el = old->disconnected_clone();
        el->set_parent(newparent);

        subdomain_id_type sbd_id = old->subdomain_id();
        if (id_remapping)
          {
            auto remapping_it = id_remapping->find(sbd_id);
            if (remapping_it != id_remapping->end())
              sbd_id = remapping_it->second;
          }
        el->subdomain_id() = sbd_id;

        // Hold off on trying to set the interior parent because we may actually
        // add lower dimensional elements before their interior parents
        if (old->interior_parent())
          ip_map[old] = el.get();

#ifdef LIBMESH_ENABLE_AMR
        if (old->has_children())
          for (unsigned int c = 0, nc = old->n_children(); c != nc; ++c)
            if (old->child_ptr(c) == remote_elem)
              el->add_child(const_cast<RemoteElem *>(remote_elem), c);

        //Create the parent's child pointers if necessary
        if (newparent)
          {
            unsigned int oldc = old->parent()->which_child_am_i(old);
            newparent->add_child(el.get(), oldc);
          }

        // Copy the refinement flags
        el->set_refinement_flag(old->refinement_flag());

        // Use hack_p_level since we may not have sibling elements
        // added yet
        el->hack_p_level(old->p_level());

        el->set_p_refinement_flag(old->p_refinement_flag());
#endif // #ifdef LIBMESH_ENABLE_AMR

        //Assign all the nodes
        for (auto i : el->node_index_range())
          el->set_node(i,
            this->node_ptr(old->node_id(i) + node_id_offset));

        // And start it off with the same processor id (mod _n_parts).
        el->processor_id() = cast_int<processor_id_type>
          (wrap_proc_ids ? old->processor_id() % this->n_processors() : old->processor_id());

        // Give it the same element and unique ids
        el->set_id(old->id() + element_id_offset);

        el->add_extra_integers(n_new_elem_ints);
        for (unsigned int i = 0; i != n_old_elem_ints; ++i)
          el->set_extra_integer(extra_int_maps.first[i],
                                old->get_extra_integer(i));

#ifdef LIBMESH_ENABLE_UNIQUE_ID
        el->set_unique_id(old->unique_id() + unique_id_offset);
#endif

        //Hold onto it
        if (!skip_find_neighbors)
          {
            for (auto s : old->side_index_range())
              if (old->neighbor_ptr(s) == remote_elem)
                el->set_neighbor(s, const_cast<RemoteElem *>(remote_elem));
            this->add_elem(std::move(el));
          }
        else
          {
            Elem * new_el = this->add_elem(std::move(el));
            old_elems_to_new_elems[old] = new_el;
          }
      }

    // If the other_mesh had some interior parents, we may need to
    // copy those pointers (if they're to elements in a third mesh),
    // or create new equivalent pointers (if they're to elements we
    // just copied), or scream and die (if the other mesh had interior
    // parents from a third mesh but we already have interior parents
    // that aren't to that same third mesh.
    if (!ip_map.empty())
      {
        bool existing_interior_parents = false;
        for (const auto & elem : this->element_ptr_range())
          if (elem->interior_parent())
            {
              existing_interior_parents = true;
              break;
            }

        MeshBase * other_interior_mesh =
          const_cast<MeshBase *>(&other_mesh.interior_mesh());

        // If we don't already have interior parents, then we can just
        // use whatever interior_mesh we need for the incoming
        // elements.
        if (!existing_interior_parents)
          {
            if (other_interior_mesh == &other_mesh)
              this->set_interior_mesh(*this);
            else
              this->set_interior_mesh(*other_interior_mesh);
          }

        if (other_interior_mesh == &other_mesh &&
            _interior_mesh == this)
          for (auto & elem_pair : ip_map)
            elem_pair.second->set_interior_parent(
              this->elem_ptr(elem_pair.first->interior_parent()->id() + element_id_offset));
        else if (other_interior_mesh == _interior_mesh)
          for (auto & elem_pair : ip_map)
            {
              Elem * ip = const_cast<Elem *>(elem_pair.first->interior_parent());
              libmesh_assert(ip == remote_elem ||
                             ip == other_interior_mesh->elem_ptr(ip->id()));
              elem_pair.second->set_interior_parent(ip);
            }
        else
          libmesh_error_msg("Cannot copy boundary elements between meshes with different interior meshes");
      }

    // Loop (again) over the elements to fill in the neighbors
    if (skip_find_neighbors)
      {
        old_elems_to_new_elems[remote_elem] = const_cast<RemoteElem*>(remote_elem);

        for (const auto & old_elem : other_mesh.element_ptr_range())
          {
            Elem * new_elem = old_elems_to_new_elems[old_elem];
            for (auto s : old_elem->side_index_range())
              {
                const Elem * old_neighbor = old_elem->neighbor_ptr(s);
                Elem * new_neighbor = old_elems_to_new_elems[old_neighbor];
                new_elem->set_neighbor(s, new_neighbor);
              }
          }
      }
  }

#ifdef LIBMESH_ENABLE_UNIQUE_ID
  // We set the unique ids of nodes after adding them to the mesh such that our value of
  // _next_unique_id may be wrong. So we amend that here
  this->set_next_unique_id(other_mesh.parallel_max_unique_id() + unique_id_offset + 1);
#endif

  // Finally, partially prepare the new Mesh for use, if that isn't
  // being skipped.
  // Even the default behavior here is for backwards compatibility,
  // and we don't want to prepare everything.

  if (!skip_preparation)
    {
      // Keep the same numbering and partitioning and distribution
      // status for now, but save our original policies to restore
      // later.
      const bool allowed_renumbering = this->allow_renumbering();
      const bool allowed_find_neighbors = this->allow_find_neighbors();
      const bool allowed_elem_removal = this->allow_remote_element_removal();
      this->allow_renumbering(false);
      this->allow_remote_element_removal(false);
      this->allow_find_neighbors(!skip_find_neighbors);

      // We should generally be able to skip *all* partitioning here
      // because we're only adding one already-consistent mesh to
      // another.
      const bool skipped_partitioning = this->skip_partitioning();
      this->skip_partitioning(true);

      const bool was_prepared = this->is_prepared();
      this->prepare_for_use();

      //But in the long term, don't change our policies.
      this->allow_find_neighbors(allowed_find_neighbors);
      this->allow_renumbering(allowed_renumbering);
      this->allow_remote_element_removal(allowed_elem_removal);
      this->skip_partitioning(skipped_partitioning);

      // That prepare_for_use() call marked us as prepared, but we
      // specifically avoided some important preparation, so we might not
      // actually be prepared now.
      if (skip_find_neighbors ||
          !was_prepared || !other_mesh.is_prepared())
        this->unset_is_prepared();
    }

  // In general we've just invalidated just about everything, and we'd
  // like to unset_is_prepared(), but specific use cases might know a
  // priori that they're still partitioned well, or that they've
  // copied in a disjoint mesh component and don't need new neighbor
  // pointers, or that they're not adding anything that would change
  // cached subdomain/element/boundary sets, etc., so we'll rely on
  // users of the "advanced" skip_preparation option to also set what
  // preparation they still need.

  // else
    // this->unset_is_prepared();
}



UnstructuredMesh::~UnstructuredMesh ()
{
  //  this->clear ();  // Nothing to clear at this level

  libmesh_exceptionless_assert (!libMesh::closed());
}





void UnstructuredMesh::find_neighbors (const bool reset_remote_elements,
                                       const bool reset_current_list,
                                       const bool assert_valid,
                                       const bool check_non_remote)
{
  // We might actually want to run this on an empty mesh
  // (e.g. the boundary mesh for a nonexistent bcid!)
  // libmesh_assert_not_equal_to (this->n_nodes(), 0);
  // libmesh_assert_not_equal_to (this->n_elem(), 0);

  // This function must be run on all processors at once
  parallel_object_only();

  LOG_SCOPE("find_neighbors()", "Mesh");

  //TODO:[BSK] This should be removed later?!
  if (reset_current_list)
    for (const auto & e : this->element_ptr_range())
      for (auto s : e->side_index_range())
        if (e->neighbor_ptr(s) != remote_elem || reset_remote_elements)
          e->set_neighbor(s, nullptr);

  // Find neighboring elements by first finding elements
  // with identical side keys and then check to see if they
  // are neighbors
  {
    // data structures -- Use the hash_multimap if available
    typedef dof_id_type                     key_type;
    typedef std::pair<Elem *, unsigned char> val_type;
    typedef std::unordered_multimap<key_type, val_type> map_type;

    // A map from side keys to corresponding elements & side numbers
    map_type side_to_elem_map;

    // Pull objects out of the loop to reduce heap operations
    std::unique_ptr<Elem> my_side, their_side;

    for (const auto & element : this->element_ptr_range())
      {
        for (auto ms : element->side_index_range())
          {
          next_side:
            // If we haven't yet found a neighbor on this side, try.
            // Even if we think our neighbor is remote, that
            // information may be out of date.
            //
            // If we're only checking remote neighbors, after a
            // redistribution, then we'll skip the non-remote ones
            if ((element->neighbor_ptr(ms) == nullptr && check_non_remote) ||
                element->neighbor_ptr(ms) == remote_elem)
              {
                // Get the key for the side of this element.  Use the
                // low_order_key so we can find neighbors in
                // mixed-order meshes if necessary.
                const dof_id_type key = element->low_order_key(ms);

                // Look for elements that have an identical side key
                auto bounds = side_to_elem_map.equal_range(key);

                // May be multiple keys, check all the possible
                // elements which _might_ be neighbors.
                if (bounds.first != bounds.second)
                  {
                    // Get the side for this element
                    element->side_ptr(my_side, ms);

                    // Look at all the entries with an equivalent key
                    while (bounds.first != bounds.second)
                      {
                        // Get the potential element
                        Elem * neighbor = bounds.first->second.first;

                        // Get the side for the neighboring element
                        const unsigned int ns = bounds.first->second.second;
                        neighbor->side_ptr(their_side, ns);
                        //libmesh_assert(my_side.get());
                        //libmesh_assert(their_side.get());

                        // If found a match with my side
                        //
                        // In 1D, since parents and children have an
                        // equal side (i.e. a node) we need to check
                        // for matching level() to avoid setting our
                        // neighbor pointer to any of our neighbor's
                        // descendants.
                        if ((*my_side == *their_side) &&
                            (element->level() == neighbor->level()))
                          {
                            // So share a side.  Is this a mixed pair
                            // of subactive and active/ancestor
                            // elements?
                            // If not, then we're neighbors.
                            // If so, then the subactive's neighbor is

                            if (element->subactive() ==
                                neighbor->subactive())
                              {
                                // an element is only subactive if it has
                                // been coarsened but not deleted
                                element->set_neighbor (ms,neighbor);
                                neighbor->set_neighbor(ns,element);
                              }
                            else if (element->subactive())
                              {
                                element->set_neighbor(ms,neighbor);
                              }
                            else if (neighbor->subactive())
                              {
                                neighbor->set_neighbor(ns,element);
                              }
                            side_to_elem_map.erase (bounds.first);

                            // get out of this nested crap
                            goto next_side;
                          }

                        ++bounds.first;
                      }
                  }

                // didn't find a match...
                // Build the map entry for this element
                side_to_elem_map.emplace
                  (key, std::make_pair(element, cast_int<unsigned char>(ms)));
              }
          }
      }
  }

#ifdef LIBMESH_ENABLE_PERIODIC
  // Get the disjoint neighbor boundary pairs object (from periodic BCs)
  auto * db = this->get_disjoint_neighbor_boundary_pairs();

  if (db)
    {
      // Obtain a point locator
      std::unique_ptr<PointLocatorBase> point_locator = this->sub_point_locator();

      for (const auto & element : this->element_ptr_range())
        {
          for (auto ms : element->side_index_range())
            {
              // Skip if this side already has a valid neighbor (including remote neighbors)
              if (element->neighbor_ptr(ms) != nullptr &&
                  element->neighbor_ptr(ms) != remote_elem)
                continue;

              for (const auto & [id, boundary_ptr] : *db)
                {
                  if (!this->get_boundary_info().has_boundary_id(element, ms, id))
                    continue;

                  unsigned int neigh_side;
                  const Elem * neigh =
                    db->neighbor(id, *point_locator, element, ms, &neigh_side);

                  if (neigh && neigh != remote_elem && neigh != element)
                    {
                      auto neigh_changeable = this->elem_ptr(neigh->id());
                      element->set_neighbor(ms, neigh_changeable);
                      neigh_changeable->set_neighbor(neigh_side, element);
                    }
                }
            }
        }
    }
#endif // LIBMESH_ENABLE_PERIODIC

#ifdef LIBMESH_ENABLE_AMR

  /**
   * Here we look at all of the child elements which
   * don't already have valid neighbors.
   *
   * If a child element has a nullptr neighbor it is
   * either because it is on the boundary or because
   * its neighbor is at a different level.  In the
   * latter case we must get the neighbor from the
   * parent.
   *
   * If a child element has a remote_elem neighbor
   * on a boundary it shares with its parent, that
   * info may have become out-dated through coarsening
   * of the neighbor's parent.  In this case, if the
   * parent's neighbor is active then the child should
   * share it.
   *
   * Furthermore, that neighbor better be active,
   * otherwise we missed a child somewhere.
   *
   *
   * We also need to look through children ordered by increasing
   * refinement level in order to add new interior_parent() links in
   * boundary elements which have just been generated by refinement,
   * and fix links in boundary elements whose previous
   * interior_parent() has just been coarsened away.
   */
  const unsigned int n_levels = MeshTools::n_levels(*this);
  for (unsigned int level = 1; level < n_levels; ++level)
    {
      for (auto & current_elem : as_range(level_elements_begin(level),
                                          level_elements_end(level)))
        {
          libmesh_assert(current_elem);
          Elem * parent = current_elem->parent();
          libmesh_assert(parent);
          const unsigned int my_child_num = parent->which_child_am_i(current_elem);

          for (auto s : current_elem->side_index_range())
            {
              if (current_elem->neighbor_ptr(s) == nullptr ||
                  (current_elem->neighbor_ptr(s) == remote_elem &&
                   parent->is_child_on_side(my_child_num, s)))
                {
                  Elem * neigh = parent->neighbor_ptr(s);

                  // If neigh was refined and had non-subactive children
                  // made remote earlier, then our current elem should
                  // actually have one of those remote children as a
                  // neighbor
                  if (neigh &&
                      (neigh->ancestor() ||
                       // If neigh has subactive children which should have
                       // matched as neighbors of the current element but
                       // did not, then those likewise must be remote
                       // children.
                       (current_elem->subactive() && neigh->has_children() &&
                        (neigh->level()+1) == current_elem->level())))
                    {
#ifdef DEBUG
                      // Let's make sure that "had children made remote"
                      // situation is actually the case
                      libmesh_assert(neigh->has_children());
                      bool neigh_has_remote_children = false;
                      for (auto & child : neigh->child_ref_range())
                        if (&child == remote_elem)
                          neigh_has_remote_children = true;
                      libmesh_assert(neigh_has_remote_children);

                      // And let's double-check that we don't have
                      // a remote_elem neighboring an active local element
                      if (current_elem->active())
                        libmesh_assert_not_equal_to (current_elem->processor_id(),
                                                     this->processor_id());
#endif // DEBUG
                      neigh = const_cast<RemoteElem *>(remote_elem);
                    }
                  // If neigh and current_elem are more than one level
                  // apart, figuring out whether we have a remote
                  // neighbor here becomes much harder.
                  else if (neigh && (current_elem->subactive() &&
                                     neigh->has_children()))
                    {
                      // Find the deepest descendant of neigh which
                      // we could consider for a neighbor.  If we run
                      // out of neigh children, then that's our
                      // neighbor.  If we find a potential neighbor
                      // with remote_children and we don't find any
                      // potential neighbors among its non-remote
                      // children, then our neighbor must be remote.
                      while (neigh != remote_elem &&
                             neigh->has_children())
                        {
                          bool found_neigh = false;
                          for (unsigned int c = 0, nc = neigh->n_children();
                               !found_neigh && c != nc; ++c)
                            {
                              Elem * child = neigh->child_ptr(c);
                              if (child == remote_elem)
                                continue;
                              for (auto ncn : child->neighbor_ptr_range())
                                {
                                  if (ncn != remote_elem &&
                                      ncn->is_ancestor_of(current_elem))
                                    {
                                      neigh = ncn;
                                      found_neigh = true;
                                      break;
                                    }
                                }
                            }
                          if (!found_neigh)
                            neigh = const_cast<RemoteElem *>(remote_elem);
                        }
                    }
                  current_elem->set_neighbor(s, neigh);
#ifdef DEBUG
                  if (neigh != nullptr && neigh != remote_elem)
                    // We ignore subactive elements here because
                    // we don't care about neighbors of subactive element.
                    if ((!neigh->active()) && (!current_elem->subactive()))
                      {
                        libMesh::err << "On processor " << this->processor_id()
                                     << std::endl;
                        libMesh::err << "Bad element ID = " << current_elem->id()
                                     << ", Side " << s << ", Bad neighbor ID = " << neigh->id() << std::endl;
                        libMesh::err << "Bad element proc_ID = " << current_elem->processor_id()
                                     << ", Bad neighbor proc_ID = " << neigh->processor_id() << std::endl;
                        libMesh::err << "Bad element size = " << current_elem->hmin()
                                     << ", Bad neighbor size = " << neigh->hmin() << std::endl;
                        libMesh::err << "Bad element center = " << current_elem->vertex_average()
                                     << ", Bad neighbor center = " << neigh->vertex_average() << std::endl;
                        libMesh::err << "ERROR: "
                                     << (current_elem->active()?"Active":"Ancestor")
                                     << " Element at level "
                                     << current_elem->level() << std::endl;
                        libMesh::err << "with "
                                     << (parent->active()?"active":
                                         (parent->subactive()?"subactive":"ancestor"))
                                     << " parent share "
                                     << (neigh->subactive()?"subactive":"ancestor")
                                     << " neighbor at level " << neigh->level()
                                     << std::endl;
                        NameBasedIO(*this).write ("bad_mesh.gmv");
                        libmesh_error_msg("Problematic mesh written to bad_mesh.gmv.");
                      }
#endif // DEBUG
                }
            }

          // We can skip to the next element if we're full-dimension
          // and therefore don't have any interior parents
          if (current_elem->dim() >= LIBMESH_DIM)
            continue;

          // We have no interior parents unless we can find one later
          current_elem->set_interior_parent(nullptr);

          Elem * pip = parent->interior_parent();

          if (!pip)
            continue;

          // If there's no interior_parent children, whether due to a
          // remote element or a non-conformity, then there's no
          // children to search.
          if (pip == remote_elem || pip->active())
            {
              current_elem->set_interior_parent(pip);
              continue;
            }

          // For node comparisons we'll need a sensible tolerance
          Real node_tolerance = current_elem->hmin() * TOLERANCE;

          // Otherwise our interior_parent should be a child of our
          // parent's interior_parent.
          for (auto & child : pip->child_ref_range())
            {
              // If we have a remote_elem, that might be our
              // interior_parent.  We'll set it provisionally now and
              // keep trying to find something better.
              if (&child == remote_elem)
                {
                  current_elem->set_interior_parent
                    (const_cast<RemoteElem *>(remote_elem));
                  continue;
                }

              bool child_contains_our_nodes = true;
              for (auto & n : current_elem->node_ref_range())
                {
                  bool child_contains_this_node = false;
                  for (auto & cn : child.node_ref_range())
                    if (cn.absolute_fuzzy_equals
                        (n, node_tolerance))
                      {
                        child_contains_this_node = true;
                        break;
                      }
                  if (!child_contains_this_node)
                    {
                      child_contains_our_nodes = false;
                      break;
                    }
                }
              if (child_contains_our_nodes)
                {
                  current_elem->set_interior_parent(&child);
                  break;
                }
            }

          // We should have found *some* interior_parent at this
          // point, whether semilocal or remote.
          libmesh_assert(current_elem->interior_parent());
        }
    }
#endif // AMR

#ifdef DEBUG
  if (assert_valid)
    {
      MeshTools::libmesh_assert_valid_neighbors(*this,
                                                !reset_remote_elements);
      MeshTools::libmesh_assert_valid_amr_interior_parents(*this);
    }
#else
  libmesh_ignore(assert_valid);
#endif

  this->_preparation.has_neighbor_ptrs = true;
}



void UnstructuredMesh::read (const std::string & name,
                             void *,
                             bool skip_renumber_nodes_and_elements,
                             bool skip_find_neighbors)
{
  // Set the skip_renumber_nodes_and_elements flag on all processors
  // if necessary.
  // This ensures that renumber_nodes_and_elements is *not* called
  // during prepare_for_use() for certain types of mesh files.
  // This is required in cases where there is an associated solution
  // file which expects a certain ordering of the nodes.
  if (Utility::ends_with(name, ".gmv"))
    this->allow_renumbering(false);

  NameBasedIO(*this).read(name);

  if (skip_renumber_nodes_and_elements)
    {
      // Use MeshBase::allow_renumbering() yourself instead.
      libmesh_deprecated();
      this->allow_renumbering(false);
    }

  // Done reading the mesh.  Now prepare it for use.
  const bool old_allow_find_neighbors = this->allow_find_neighbors();
  this->allow_find_neighbors(!skip_find_neighbors);
  this->prepare_for_use();
  this->allow_find_neighbors(old_allow_find_neighbors);
}



void UnstructuredMesh::write (const std::string & name) const
{
  LOG_SCOPE("write()", "Mesh");

  NameBasedIO(*this).write(name);
}



void UnstructuredMesh::write (const std::string & name,
                              const std::vector<Number> & v,
                              const std::vector<std::string> & vn) const
{
  LOG_SCOPE("write()", "Mesh");

  NameBasedIO(*this).write_nodal_data(name, v, vn);
}





void UnstructuredMesh::create_pid_mesh(UnstructuredMesh & pid_mesh,
                                       const processor_id_type pid) const
{

  // Issue a warning if the number the number of processors
  // currently available is less that that requested for
  // partitioning.  This is not necessarily an error since
  // you may run on one processor and still partition the
  // mesh into several partitions.
#ifdef DEBUG
  if (this->n_processors() < pid)
    {
      libMesh::out << "WARNING:  You are creating a "
                   << "mesh for a processor id (="
                   << pid
                   << ") greater than "
                   << "the number of processors available for "
                   << "the calculation. (="
                   << this->n_processors()
                   << ")."
                   << std::endl;
    }
#endif

  this->create_submesh (pid_mesh,
                        this->active_pid_elements_begin(pid),
                        this->active_pid_elements_end(pid));
}







void UnstructuredMesh::create_submesh (UnstructuredMesh & new_mesh,
                                       const const_element_iterator & it,
                                       const const_element_iterator & it_end) const
{
  // Just in case the subdomain_mesh already has some information
  // in it, get rid of it.
  new_mesh.clear();

  // If we're not serial, our submesh isn't either.
  // There are no remote elements to delete on an empty mesh, but
  // calling the method to do so marks the mesh as parallel.
  if (!this->is_serial())
    new_mesh.delete_remote_elements();

  // Fail if (*this == new_mesh), we cannot create a submesh inside ourself!
  // This may happen if the user accidentally passes the original mesh into
  // this function!  We will check this by making sure we did not just
  // clear ourself.
  libmesh_assert_not_equal_to (this->n_nodes(), 0);
  libmesh_assert_not_equal_to (this->n_elem(), 0);

  // Container to catch boundary IDs handed back by BoundaryInfo
  std::vector<boundary_id_type> bc_ids;

  // Put any extra integers on the new mesh too
  new_mesh.merge_extra_integer_names(*this);
  const unsigned int n_node_ints = _node_integer_names.size();

  for (const auto & old_elem : as_range(it, it_end))
    {
      // Add an equivalent element type to the new_mesh.
      // disconnected_clone() copies ids, extra element integers, etc.
      auto uelem = old_elem->disconnected_clone();
      Elem * new_elem = new_mesh.add_elem(std::move(uelem));
      libmesh_assert(new_elem);

      // Loop over the nodes on this element.
      for (auto n : old_elem->node_index_range())
        {
          const dof_id_type this_node_id = old_elem->node_id(n);

          // Add this node to the new mesh if it's not there already
          if (!new_mesh.query_node_ptr(this_node_id))
            {
              Node * newn =
                new_mesh.add_point (old_elem->point(n),
                                    this_node_id,
                                    old_elem->node_ptr(n)->processor_id());

              newn->add_extra_integers(n_node_ints);
              for (unsigned int i = 0; i != n_node_ints; ++i)
                newn->set_extra_integer(i, old_elem->node_ptr(n)->get_extra_integer(i));

#ifdef LIBMESH_ENABLE_UNIQUE_ID
              newn->set_unique_id(old_elem->node_ptr(n)->unique_id());
#endif
            }

          // Define this element's connectivity on the new mesh
          new_elem->set_node(n, new_mesh.node_ptr(this_node_id));
        }

      // Maybe add boundary conditions for this element
      for (auto s : old_elem->side_index_range())
        {
          this->get_boundary_info().boundary_ids(old_elem, s, bc_ids);
          new_mesh.get_boundary_info().add_side (new_elem, s, bc_ids);
        }
    } // end loop over elements

  // Prepare the new_mesh for use
  new_mesh.prepare_for_use();
}



#ifdef LIBMESH_ENABLE_AMR
bool UnstructuredMesh::contract ()
{
  LOG_SCOPE ("contract()", "Mesh");

  // Flag indicating if this call actually changes the mesh
  bool mesh_changed = false;

#ifdef DEBUG
  for (const auto & elem : this->element_ptr_range())
    libmesh_assert(elem->active() || elem->subactive() || elem->ancestor());
#endif

  // Loop over the elements.
  for (auto & elem : this->element_ptr_range())
    {
      // Delete all the subactive ones
      if (elem->subactive())
        {
          // No level-0 element should be subactive.
          // Note that we CAN'T test elem->level(), as that
          // touches elem->parent()->dim(), and elem->parent()
          // might have already been deleted!
          libmesh_assert(elem->parent());

          // Delete the element
          // This just sets a pointer to nullptr, and doesn't
          // invalidate any iterators
          this->delete_elem(elem);

          // the mesh has certainly changed
          mesh_changed = true;
        }
      else
        {
          // Compress all the active ones
          if (elem->active())
            elem->contract();
          else
            libmesh_assert (elem->ancestor());
        }
    }

  // Strip any newly-created nullptr voids out of the element array
  this->renumber_nodes_and_elements();

  // FIXME: Need to understand why deleting subactive children
  // invalidates the point locator.  For now we will clear it explicitly
  this->clear_point_locator();

  // Allow our GhostingFunctor objects to reinit if necessary.
  for (auto & gf : as_range(this->ghosting_functors_begin(),
                            this->ghosting_functors_end()))
    {
      libmesh_assert(gf);
      gf->mesh_reinit();
    }

  return mesh_changed;
}
#endif // #ifdef LIBMESH_ENABLE_AMR



void UnstructuredMesh::all_first_order ()
{
  LOG_SCOPE("all_first_order()", "Mesh");

  /**
   * Prepare to identify (and then delete) a bunch of no-longer-used nodes.
   */
  std::vector<bool> node_touched_by_me(this->max_node_id(), false);

  // Loop over the high-ordered elements.
  // First make sure they _are_ indeed high-order, and then replace
  // them with an equivalent first-order element.
  for (auto & so_elem : element_ptr_range())
    {
      libmesh_assert(so_elem);

      /*
       * build the first-order equivalent, add to
       * the new_elements list.
       */
      auto lo_elem = Elem::build
        (Elem::first_order_equivalent_type
         (so_elem->type()), so_elem->parent());

      const unsigned short n_sides = so_elem->n_sides();

      for (unsigned short s=0; s != n_sides; ++s)
        if (so_elem->neighbor_ptr(s) == remote_elem)
          lo_elem->set_neighbor(s, const_cast<RemoteElem *>(remote_elem));

#ifdef LIBMESH_ENABLE_AMR
      /*
       * Reset the parent links of any child elements
       */
      if (so_elem->has_children())
        for (unsigned int c = 0, nc = so_elem->n_children(); c != nc; ++c)
          {
            Elem * child = so_elem->child_ptr(c);
            if (child != remote_elem)
              child->set_parent(lo_elem.get());
            lo_elem->add_child(child, c);
          }

      /*
       * Reset the child link of any parent element
       */
      if (so_elem->parent())
        {
          unsigned int c =
            so_elem->parent()->which_child_am_i(so_elem);
          lo_elem->parent()->replace_child(lo_elem.get(), c);
        }

      /*
       * Copy as much data to the new element as makes sense
       */
      lo_elem->set_p_level(so_elem->p_level());
      lo_elem->set_refinement_flag(so_elem->refinement_flag());
      lo_elem->set_p_refinement_flag(so_elem->p_refinement_flag());
#endif

      libmesh_assert_equal_to (lo_elem->n_vertices(), so_elem->n_vertices());

      /*
       * By definition the vertices of the linear and
       * second order element are identically numbered.
       * transfer these.
       */
      for (unsigned int v=0, snv=so_elem->n_vertices(); v < snv; v++)
        {
          lo_elem->set_node(v, so_elem->node_ptr(v));
          node_touched_by_me[lo_elem->node_id(v)] = true;
        }

      /*
       * find_neighbors relies on remote_elem neighbor links being
       * properly maintained.
       */
      for (unsigned short s=0; s != n_sides; s++)
        {
          if (so_elem->neighbor_ptr(s) == remote_elem)
            lo_elem->set_neighbor(s, const_cast<RemoteElem*>(remote_elem));
        }

      /**
       * If the second order element had any boundary conditions they
       * should be transferred to the first-order element.  The old
       * boundary conditions will be removed from the BoundaryInfo
       * data structure by insert_elem.
       */
      this->get_boundary_info().copy_boundary_ids
        (this->get_boundary_info(), so_elem, lo_elem.get());

      /*
       * The new first-order element is ready.
       * Inserting it into the mesh will replace and delete
       * the second-order element.
       */
      lo_elem->set_id(so_elem->id());
#ifdef LIBMESH_ENABLE_UNIQUE_ID
      lo_elem->set_unique_id(so_elem->unique_id());
#endif

      const unsigned int nei = so_elem->n_extra_integers();
      lo_elem->add_extra_integers(nei);
      for (unsigned int i=0; i != nei; ++i)
        lo_elem->set_extra_integer(i, so_elem->get_extra_integer(i));

      lo_elem->inherit_data_from(*so_elem);

      this->insert_elem(std::move(lo_elem));
    }

  // Deleting nodes does not invalidate iterators, so this is safe.
  for (const auto & node : this->node_ptr_range())
    if (!node_touched_by_me[node->id()])
      this->delete_node(node);

  // If crazy people applied boundary info to non-vertices and then
  // deleted those non-vertices, we should make sure their boundary id
  // caches are correct.
  this->get_boundary_info().regenerate_id_sets();

  // On hanging nodes that used to also be second order nodes, we
  // might now have an invalid nodal processor_id()
  Partitioner::set_node_processor_ids(*this);

  // delete or renumber nodes if desired
  this->prepare_for_use();
}



void
UnstructuredMesh::all_second_order_range (const SimpleRange<element_iterator> & range,
                                          const bool full_ordered)
{
  LOG_SCOPE("all_second_order_range()", "Mesh");

  /*
   * The maximum number of new second order nodes we might be adding,
   * for use when picking unique unique_id values later. This variable
   * is not used unless unique ids are enabled.
   */
  unsigned int max_new_nodes_per_elem;

  /*
   * For speed-up of the \p add_point() method, we
   * can reserve memory.  Guess the number of additional
   * nodes based on the element spatial dimensions and the
   * total number of nodes in the mesh as an upper bound.
   */
  switch (this->mesh_dimension())
    {
    case 1:
      /*
       * in 1D, there can only be order-increase from Edge2
       * to Edge3.  Something like 1/2 of n_nodes() have
       * to be added
       */
      max_new_nodes_per_elem = 3 - 2;
      this->reserve_nodes(static_cast<unsigned int>
                          (1.5*static_cast<double>(this->n_nodes())));
      break;

    case 2:
      /*
       * in 2D, either refine from Tri3 to Tri6 (double the nodes)
       * or from Quad4 to Quad8 (again, double) or Quad9 (2.25 that much)
       */
      max_new_nodes_per_elem = 9 - 4;
      this->reserve_nodes(static_cast<unsigned int>
                          (2*static_cast<double>(this->n_nodes())));
      break;


    case 3:
      /*
       * in 3D, either refine from Tet4 to Tet10 (factor = 2.5) up to
       * Hex8 to Hex27 (something  > 3).  Since in 3D there _are_ already
       * quite some nodes, and since we do not want to overburden the memory by
       * a too conservative guess, use the lower bound
       */
      max_new_nodes_per_elem = 27 - 8;
      this->reserve_nodes(static_cast<unsigned int>
                          (2.5*static_cast<double>(this->n_nodes())));
      break;

    default:
      // Hm?
      libmesh_error_msg("Unknown mesh dimension " << this->mesh_dimension());
    }

  // All the real work is done in the helper function
  all_increased_order_range(*this, range, max_new_nodes_per_elem,
    [full_ordered](ElemType t) {
      return Elem::second_order_equivalent_type(t, full_ordered);
    });
}



void UnstructuredMesh::all_complete_order_range(const SimpleRange<element_iterator> & range)
{
  LOG_SCOPE("all_complete_order()", "Mesh");

  /*
   * The maximum number of new higher-order nodes we might be adding,
   * for use when picking unique unique_id values later. This variable
   * is not used unless unique ids are enabled.
   */
  unsigned int max_new_nodes_per_elem;

  /*
   * for speed-up of the \p add_point() method, we
   * can reserve memory.  Guess the number of additional
   * nodes based on the element spatial dimensions and the
   * total number of nodes in the mesh as an upper bound.
   */
  switch (this->mesh_dimension())
    {
    case 1:
      /*
       * in 1D, there can only be order-increase from Edge2
       * to Edge3.  Something like 1/2 of n_nodes() have
       * to be added
       */
      max_new_nodes_per_elem = 3 - 2;
      this->reserve_nodes(static_cast<unsigned int>
                          (1.5*static_cast<double>(this->n_nodes())));
      break;

    case 2:
      /*
       * in 2D, we typically refine from Tri6 to Tri7 (1.1667 times
       * the nodes) but might refine from Quad4 to Quad9
       * (2.25 times the nodes)
       */
      max_new_nodes_per_elem = 9 - 4;
      this->reserve_nodes(static_cast<unsigned int>
                          (2*static_cast<double>(this->n_nodes())));
      break;


    case 3:
      /*
       * in 3D, we typically refine from Tet10 to Tet14 (factor = 1.4)
       * but may go Hex8 to Hex27 (something  > 3).  Since in 3D there
       * _are_ already quite some nodes, and since we do not want to
       * overburden the memory by a too conservative guess, use a
       * moderate bound
       */
      max_new_nodes_per_elem = 27 - 8;
      this->reserve_nodes(static_cast<unsigned int>
                          (2.5*static_cast<double>(this->n_nodes())));
      break;

    default:
      // Hm?
      libmesh_error_msg("Unknown mesh dimension " << this->mesh_dimension());
    }

  // All the real work is done in the helper function
  all_increased_order_range(*this, range, max_new_nodes_per_elem,
    [](ElemType t) {
      return Elem::complete_order_equivalent_type(t);
    });
}


std::size_t
UnstructuredMesh::stitch_meshes (const MeshBase & other_mesh,
                                 boundary_id_type this_mesh_boundary_id,
                                 boundary_id_type other_mesh_boundary_id,
                                 Real tol,
                                 bool clear_stitched_boundary_ids,
                                 bool verbose,
                                 bool use_binary_search,
                                 bool enforce_all_nodes_match_on_boundaries,
                                 bool merge_boundary_nodes_all_or_nothing,
                                 bool remap_subdomain_ids,
                                 bool prepare_after_stitching)
{
  LOG_SCOPE("stitch_meshes()", "UnstructuredMesh");
  return stitching_helper(&other_mesh,
                          this_mesh_boundary_id,
                          other_mesh_boundary_id,
                          tol,
                          clear_stitched_boundary_ids,
                          verbose,
                          use_binary_search,
                          enforce_all_nodes_match_on_boundaries,
                          true,
                          merge_boundary_nodes_all_or_nothing,
                          remap_subdomain_ids,
                          prepare_after_stitching);
}


std::size_t
UnstructuredMesh::stitch_surfaces (boundary_id_type boundary_id_1,
                                   boundary_id_type boundary_id_2,
                                   Real tol,
                                   bool clear_stitched_boundary_ids,
                                   bool verbose,
                                   bool use_binary_search,
                                   bool enforce_all_nodes_match_on_boundaries,
                                   bool merge_boundary_nodes_all_or_nothing,
                                   bool prepare_after_stitching)

{
  return stitching_helper(nullptr,
                          boundary_id_1,
                          boundary_id_2,
                          tol,
                          clear_stitched_boundary_ids,
                          verbose,
                          use_binary_search,
                          enforce_all_nodes_match_on_boundaries,
                          /* skip_find_neighbors = */ true,
                          merge_boundary_nodes_all_or_nothing,
                          /* remap_subdomain_ids = */ false,
                          prepare_after_stitching);
}


std::size_t
UnstructuredMesh::stitching_helper (const MeshBase * other_mesh,
                                    boundary_id_type this_mesh_boundary_id,
                                    boundary_id_type other_mesh_boundary_id,
                                    Real tol,
                                    bool clear_stitched_boundary_ids,
                                    bool verbose,
                                    bool use_binary_search,
                                    bool enforce_all_nodes_match_on_boundaries,
                                    bool skip_find_neighbors,
                                    bool merge_boundary_nodes_all_or_nothing,
                                    bool remap_subdomain_ids,
                                    bool prepare_after_stitching)
{
#ifdef DEBUG
  // We rely on neighbor links here
  MeshTools::libmesh_assert_valid_neighbors(*this);
#endif

  bool is_valid_disjoint_pair_to_stitch = false;

#ifdef LIBMESH_ENABLE_PERIODIC
  auto * this_db  = this->get_disjoint_neighbor_boundary_pairs();
  auto * other_db = (other_mesh ? other_mesh->get_disjoint_neighbor_boundary_pairs() : nullptr);
  const bool have_disc_bdys =
    (this_db && !this_db->empty()) || (other_db && !other_db->empty());

  if (have_disc_bdys)
    {
      const boundary_id_type a = this_mesh_boundary_id;
      const boundary_id_type b = other_mesh_boundary_id;

      auto get_pb = [](const PeriodicBoundaries * db, boundary_id_type id)
        {
          return db ? db->boundary(id) : nullptr;
        };

      // this mesh
      const auto * pb_this_a = get_pb(this_db, a);
      const auto * pb_this_b = get_pb(this_db, b);
      const bool in_this =
        (pb_this_a && pb_this_a->pairedboundary == b) ||
        (pb_this_b && pb_this_b->pairedboundary == a);

      // other mesh
      const auto * pb_other_b = get_pb(other_db, b);
      const auto * pb_other_a = get_pb(other_db, a);
      const bool in_other =
        (pb_other_b && pb_other_b->pairedboundary == a) ||
        (pb_other_a && pb_other_a->pairedboundary == b);

      // Conflict conditions:
      // Case 1: On "this" mesh, a or b exist but are not paired,
      //         while the other mesh pairs them.
      if (!in_this && (pb_this_a || pb_this_b) && in_other)
        libmesh_error_msg("Disjoint neighbor boundary pairing mismatch: on 'this' mesh, "
                          "boundary (" << a << " or " << b
                          << ") exists but is not paired; on 'other' mesh the pair is present.");

      // Case 2: On "other" mesh, a or b exist but are not paired,
      //         while this mesh pairs them.
      if (!in_other && (pb_other_a || pb_other_b) && in_this)
        libmesh_error_msg("Disjoint neighbor boundary pairing mismatch: on 'other' mesh, "
                          "boundary (" << a << " or " << b
                          << ") exists but is not paired; on 'this' mesh the pair is present.");

      // Legal conditions: either side has a correct pairing
      if (in_this || in_other)
        is_valid_disjoint_pair_to_stitch = true;
    }
#endif // LIBMESH_ENABLE_PERIODIC

  // We can't even afford any unset neighbor links here.
  if (!this->is_prepared())
    this->find_neighbors();

  // FIXME: make distributed mesh support efficient.
  // Yes, we currently suck.
  MeshSerializer serialize(*this);

  // *Badly*.
  std::unique_ptr<MeshSerializer> serialize_other;
  if (other_mesh)
    serialize_other = std::make_unique<MeshSerializer>
      (*const_cast<MeshBase *>(other_mesh));

  std::map<dof_id_type, dof_id_type> node_to_node_map, other_to_this_node_map; // The second is the inverse map of the first
  std::map<dof_id_type, std::vector<dof_id_type>> node_to_elems_map;

  typedef dof_id_type                     key_type;
  typedef std::pair<const Elem *, unsigned char> val_type;
  typedef std::pair<key_type, val_type>   key_val_pair;
  typedef std::unordered_multimap<key_type, val_type> map_type;
  // Mapping between all side keys in this mesh and elements+side numbers relevant to the boundary in this mesh as well.
  map_type side_to_elem_map;

  // If there is only one mesh (i.e. other_mesh == nullptr), then loop over this mesh twice
  if (!other_mesh)
    {
      other_mesh = this;
    }

  if ((this_mesh_boundary_id  != BoundaryInfo::invalid_id) &&
      (other_mesh_boundary_id != BoundaryInfo::invalid_id))
    {
      LOG_SCOPE("stitch_meshes node merging", "UnstructuredMesh");

      // While finding nodes on the boundary, also find the minimum edge length
      // of all faces on both boundaries.  This will later be used in relative
      // distance checks when stitching nodes.
      Real h_min = std::numeric_limits<Real>::max();
      bool h_min_updated = false;

      // Loop below fills in these sets for the two meshes.
      std::set<dof_id_type> this_boundary_node_ids, other_boundary_node_ids;

      // Pull objects out of the loop to reduce heap operations
      std::unique_ptr<const Elem> side;

      {
        // Make temporary fixed-size arrays for loop
        boundary_id_type id_array[2]         = {this_mesh_boundary_id, other_mesh_boundary_id};
        std::set<dof_id_type> * set_array[2] = {&this_boundary_node_ids, &other_boundary_node_ids};
        const MeshBase * mesh_array[2] = {this, other_mesh};

        for (unsigned i=0; i<2; ++i)
          {
            // First we deal with node boundary IDs.  We only enter
            // this loop if we have at least one nodeset. Note that we
            // do not attempt to make an h_min determination here.
            // The h_min determination is done while looping over the
            // Elems and checking their sides and edges for boundary
            // information, below.
            if (mesh_array[i]->get_boundary_info().n_nodeset_conds() > 0)
              {
                // build_node_list() returns a vector of (node-id, bc-id) tuples
                for (const auto & t : mesh_array[i]->get_boundary_info().build_node_list())
                  {
                    boundary_id_type node_bc_id = std::get<1>(t);
                    if (node_bc_id == id_array[i])
                      {
                        dof_id_type this_node_id = std::get<0>(t);
                        set_array[i]->insert( this_node_id );
                      }
                  }
              }

            // Container to catch boundary IDs passed back from BoundaryInfo.
            std::vector<boundary_id_type> bc_ids;

            // Pointers to boundary NodeElems encountered while looping over the entire Mesh
            // and checking side and edge boundary ids. The Nodes associated with NodeElems
            // may be in a boundary nodeset, but not connected to any other Elems. In this
            // case, we also consider the "minimum node separation distance" amongst all
            // NodeElems when determining the relevant h_min value for this mesh.
            std::vector<const Elem *> boundary_node_elems;

            for (auto & el : mesh_array[i]->element_ptr_range())
              {
                // Now check whether elem has a face on the specified boundary
                for (auto side_id : el->side_index_range())
                  {
                    bool should_stitch_this_side =
                      (el->neighbor_ptr(side_id) == nullptr) ||
                      (is_valid_disjoint_pair_to_stitch &&
                      mesh_array[i]->get_boundary_info().has_boundary_id(el, side_id, id_array[i]));

                    if (should_stitch_this_side)
                      {
                        // Get *all* boundary IDs on this side, not just the first one!
                        mesh_array[i]->get_boundary_info().boundary_ids (el, side_id, bc_ids);

                        if (std::find(bc_ids.begin(), bc_ids.end(), id_array[i]) != bc_ids.end())
                          {
                            el->build_side_ptr(side, side_id);
                            for (auto & n : side->node_ref_range())
                              set_array[i]->insert(n.id());

                            h_min = std::min(h_min, side->hmin());
                            h_min_updated = true;

                            // This side is on the boundary, add its information to side_to_elem
                            if (skip_find_neighbors && (i==0))
                              {
                                key_type key = el->low_order_key(side_id);
                                val_type val;
                                val.first = el;
                                val.second = cast_int<unsigned char>(side_id);

                                key_val_pair kvp;
                                kvp.first = key;
                                kvp.second = val;
                                side_to_elem_map.insert (kvp);
                              }
                          }

                        // Also, check the edges on this side. We don't have to worry about
                        // updating neighbor info in this case since elements don't store
                        // neighbor info on edges.
                        for (auto edge_id : el->edge_index_range())
                          {
                            if (el->is_edge_on_side(edge_id, side_id))
                              {
                                // Get *all* boundary IDs on this edge, not just the first one!
                                mesh_array[i]->get_boundary_info().edge_boundary_ids (el, edge_id, bc_ids);

                                if (std::find(bc_ids.begin(), bc_ids.end(), id_array[i]) != bc_ids.end())
                                  {
                                    std::unique_ptr<const Elem> edge (el->build_edge_ptr(edge_id));
                                    for (auto & n : edge->node_ref_range())
                                      set_array[i]->insert( n.id() );

                                    h_min = std::min(h_min, edge->hmin());
                                    h_min_updated = true;
                                  }
                              }
                          } // end for (edge_id)
                      } // end if (should_stitch_this_side)
                  } // end for (side_id)

                // Alternatively, is this a boundary NodeElem? If so,
                // add it to a list of NodeElems that will later be
                // used to set h_min based on the minimum node
                // separation distance between all pairs of boundary
                // NodeElems.
                if (el->type() == NODEELEM)
                  {
                    mesh_array[i]->get_boundary_info().boundary_ids(el->node_ptr(0), bc_ids);
                    if (std::find(bc_ids.begin(), bc_ids.end(), id_array[i]) != bc_ids.end())
                      {
                        boundary_node_elems.push_back(el);

                        // Debugging:
                        // libMesh::out << "Elem " << el->id() << " is a NodeElem on boundary " << id_array[i] << std::endl;
                      }
                  } // end if (el->type() == NODEELEM)
              } // end for (el)

            // Compute the minimum node separation distance amongst
            // all boundary NodeElem pairs.
            {
              const auto N = boundary_node_elems.size();
              for (auto node_elem_i : make_range(N))
                for (auto node_elem_j : make_range(node_elem_i+1, N))
                  {
                    Real node_sep =
                      (boundary_node_elems[node_elem_i]->point(0) - boundary_node_elems[node_elem_j]->point(0)).norm();

                    // We only want to consider non-coincident
                    // boundary NodeElem pairs when determining the
                    // minimum node separation distance.
                    if (node_sep > 0.)
                      {
                        h_min = std::min(h_min, node_sep);
                        h_min_updated = true;
                      }
                  } // end for (node_elem_j)
            } // end minimum NodeElem separation scope
          } // end for (i)
      } // end scope

      if (verbose)
        {
          libMesh::out << "In UnstructuredMesh::stitch_meshes:\n"
                       << "This mesh has "  << this_boundary_node_ids.size()
                       << " nodes on boundary `"
                       << this->get_boundary_info().get_sideset_name(this_mesh_boundary_id)
                       << "' (" << this_mesh_boundary_id  << ").\n"
                       << "Other mesh has " << other_boundary_node_ids.size()
                       << " nodes on boundary `"
                       << this->get_boundary_info().get_sideset_name(other_mesh_boundary_id)
                       << "' (" << other_mesh_boundary_id  << ").\n";

          if (h_min_updated)
            {
              libMesh::out << "Minimum edge length on both surfaces is " << h_min << ".\n";
            }
          else
            {
              libMesh::out << "No minimum edge length determined on specified surfaces." << std::endl;
            }
        }

      // At this point, if h_min==0 it means that there were at least two coincident
      // nodes on the surfaces being stitched, and we don't currently support that case.
      // (It might be possible to support, but getting it exactly right would be tricky
      // and probably not worth the extra complications to the "normal" case.)
      libmesh_error_msg_if(h_min < std::numeric_limits<Real>::epsilon(),
                           "Coincident nodes detected on source and/or target "
                           "surface, stitching meshes is not possible.");

      // We require nanoflann for the "binary search" (really kd-tree)
      // option to work. If it's not available, turn that option off,
      // warn the user, and fall back on the N^2 search algorithm.
      if (use_binary_search)
        {
#ifndef LIBMESH_HAVE_NANOFLANN
          use_binary_search = false;
          libmesh_warning("The use_binary_search option in the "
                          "UnstructuredMesh stitching algorithms requires nanoflann "
                          "support. Falling back on N^2 search algorithm.");
#endif
        }

      if (!this_boundary_node_ids.empty())
      {
        if (use_binary_search)
        {
#ifdef LIBMESH_HAVE_NANOFLANN
          typedef nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<Real, VectorOfNodesAdaptor>,
            VectorOfNodesAdaptor, 3, std::size_t> kd_tree_t;

          // Create the dataset needed to build the kd tree with nanoflann
          std::vector<std::pair<Point, dof_id_type>> this_mesh_nodes(this_boundary_node_ids.size());

          for (auto [it, ctr] = std::make_tuple(this_boundary_node_ids.begin(), 0u);
               it != this_boundary_node_ids.end(); ++it, ++ctr)
          {
            this_mesh_nodes[ctr].first = this->point(*it);
            this_mesh_nodes[ctr].second = *it;
          }

          VectorOfNodesAdaptor vec_nodes_adaptor(this_mesh_nodes);

          kd_tree_t this_kd_tree(3, vec_nodes_adaptor, 10);
          this_kd_tree.buildIndex();

          // Storage for nearest neighbor in the loop below
          std::size_t ret_index;
          Real ret_dist_sqr;

          // Loop over other mesh. For each node, find its nearest neighbor in this mesh, and fill in the maps.
          for (const auto & node_id : other_boundary_node_ids)
          {
            const auto & p = other_mesh->point(node_id);
            const Real query_pt[] = {p(0), p(1), p(2)};
            this_kd_tree.knnSearch(&query_pt[0], 1, &ret_index, &ret_dist_sqr);

            // TODO: here we should use the user's specified tolerance
            // and the previously determined value of h_min in the
            // distance comparison, not just TOLERANCE^2.
            if (ret_dist_sqr < TOLERANCE*TOLERANCE)
            {
              node_to_node_map[this_mesh_nodes[ret_index].second] = node_id;
              other_to_this_node_map[node_id] = this_mesh_nodes[ret_index].second;
            }
          }

          // If the two maps don't have the same size, it means one
          // node in this mesh is the nearest neighbor of several
          // nodes in other mesh. Since the stitching is ambiguous in
          // this case, we throw an error.
          libmesh_error_msg_if(node_to_node_map.size() != other_to_this_node_map.size(),
                               "Error: Found multiple matching nodes in stitch_meshes");
#endif
        }
        else // !use_binary_search
        {
          // In the unlikely event that two meshes composed entirely of
          // NodeElems are being stitched together, we will not have
          // selected a valid h_min value yet, and the distance
          // comparison below will be true for essentially any two
          // nodes. In this case we simply fall back on an absolute
          // distance check.
          if (!h_min_updated)
            {
              libmesh_warning("No valid h_min value was found, falling back on "
                              "absolute distance check in the N^2 search algorithm.");
              h_min = 1.;
            }

          // Otherwise, use a simple N^2 search to find the closest matching points. This can be helpful
          // in the case that we have tolerance issues which cause mismatch between the two surfaces
          // that are being stitched.
          for (const auto & this_node_id : this_boundary_node_ids)
          {
            Node & this_node = this->node_ref(this_node_id);

            bool found_matching_nodes = false;

            for (const auto & other_node_id : other_boundary_node_ids)
            {
              const Node & other_node = other_mesh->node_ref(other_node_id);

              Real node_distance = (this_node - other_node).norm();

              if (node_distance < tol*h_min)
              {
                // Make sure we didn't already find a matching node!
                libmesh_error_msg_if(found_matching_nodes,
                                     "Error: Found multiple matching nodes in stitch_meshes");

                node_to_node_map[this_node_id] = other_node_id;
                other_to_this_node_map[other_node_id] = this_node_id;

                found_matching_nodes = true;
              }
            }
          }
        }
      }

      // Build up the node_to_elems_map, using only one loop over other_mesh
      for (auto & el : other_mesh->element_ptr_range())
        {
          // For each node on the element, find the corresponding node
          // on "this" Mesh, 'this_node_id', if it exists, and push
          // the current element ID back onto node_to_elems_map[this_node_id].
          // For that we will use the reverse mapping we created at
          // the same time as the forward mapping.
          for (auto & n : el->node_ref_range())
            if (const auto it = other_to_this_node_map.find(/*other_node_id=*/n.id());
                it != other_to_this_node_map.end())
              node_to_elems_map[/*this_node_id=*/it->second].push_back( el->id() );
        }

      if (verbose)
        {
          libMesh::out << "In UnstructuredMesh::stitch_meshes:\n"
                       << "Found " << node_to_node_map.size()
                       << " matching nodes.\n"
                       << std::endl;
        }

      if (enforce_all_nodes_match_on_boundaries)
        {
          std::size_t n_matching_nodes = node_to_node_map.size();
          std::size_t this_mesh_n_nodes = this_boundary_node_ids.size();
          std::size_t other_mesh_n_nodes = other_boundary_node_ids.size();
          libmesh_error_msg_if((n_matching_nodes != this_mesh_n_nodes) || (n_matching_nodes != other_mesh_n_nodes),
                               "Error: We expected the number of nodes to match.");
        }

      if (merge_boundary_nodes_all_or_nothing)
        {
          std::size_t n_matching_nodes = node_to_node_map.size();
          std::size_t this_mesh_n_nodes = this_boundary_node_ids.size();
          std::size_t other_mesh_n_nodes = other_boundary_node_ids.size();
          if ((n_matching_nodes != this_mesh_n_nodes) || (n_matching_nodes != other_mesh_n_nodes))
            {
              if (verbose)
                {
                  libMesh::out << "Skipping node merging in "
                                  "UnstructuredMesh::stitch_meshes because not "
                                  "all boundary nodes were matched."
                               << std::endl;
                }
              node_to_node_map.clear();
              other_to_this_node_map.clear();
              node_to_elems_map.clear();
            }
        }
    }
  else
    {
      if (verbose)
        {
          libMesh::out << "Skip node merging in UnstructuredMesh::stitch_meshes:" << std::endl;
        }
    }

  dof_id_type node_delta = this->max_node_id();
  dof_id_type elem_delta = this->max_elem_id();

  unique_id_type unique_delta =
#ifdef LIBMESH_ENABLE_UNIQUE_ID
    this->parallel_max_unique_id();
#else
    0;
#endif

  // If other_mesh != nullptr, then we have to do a bunch of work
  // in order to copy it to this mesh
  if (this!=other_mesh)
    {
      LOG_SCOPE("stitch_meshes copying", "UnstructuredMesh");

#ifdef LIBMESH_ENABLE_PERIODIC
        // Copy disjoint neighbor boundary pairs (PeriodicBoundary objects)
        // from `other_mesh` to `this` mesh
        if (other_db && !other_db->empty())
          {
            for (const auto & [bdy_id, pb_ptr] : *other_db)
              {
                const auto & pb = *pb_ptr;
                const boundary_id_type a = pb.myboundary;
                const boundary_id_type b = pb.pairedboundary;

                if (this_db)
                  {
                    // Skip if identical pair already exists
                    if (const auto * existing_pb = this_db->boundary(a))
                      if ((existing_pb->myboundary == a && existing_pb->pairedboundary == b) ||
                          (existing_pb->myboundary == b && existing_pb->pairedboundary == a))
                        continue;

                    // If both boundary ids exist on this mesh but aren't paired here, refuse to create a new pair
                    const auto & bdy_ids = this->get_boundary_info().get_boundary_ids();
                    const bool a_exists = bdy_ids.count(a);
                    const bool b_exists = bdy_ids.count(b);
                    // If a and b already exist on `this`, we should be screaming and dying
                    // unless they already have PeriodicBoundary objects connecting them too
                    if (a_exists && b_exists && !this_db->boundary(a))
                      libmesh_error_msg("Conflict: boundaries " << a << " and " << b
                                        << " already exist on this mesh but are not paired.");
                  }

                this->add_disjoint_neighbor_boundary_pairs(a, b, pb.get_corresponding_pos(Point(0.0,0.0,0.0)));
              }
          }
#endif // LIBMESH_ENABLE_PERIODIC


      // Increment the node_to_node_map and node_to_elems_map
      // to account for id offsets
      for (auto & pr : node_to_node_map)
        pr.second += node_delta;

      for (auto & pr : node_to_elems_map)
        for (auto & entry : pr.second)
          entry += elem_delta;

      // We run into problems when the libMesh subdomain standard (the
      // id defines the subdomain; the name was an afterthought) and
      // the MOOSE standard (the name defines the subdomain; the id
      // might be autogenerated) clash.
      //
      // Subdomain ids with the same name in both meshes are surely
      // meant to represent the same subdomain.  We can just merge
      // them.
      //
      // Subdomain ids which don't have a name in either mesh are
      // almost surely meant to represent the same subdomain.  We'll
      // just merge them.
      //
      // Subdomain ids with different names in different meshes, or
      // names with different ids in different meshes, are trickier.
      // For backwards compatibility we default to the old "just copy
      // all the subdomain ids over" behavior, but if requested we'll
      // remap any ids that appear to be clear conflicts, and we'll
      // scream and die if we see any ids that are ambiguous due to
      // being named in one mesh but not the other.
      std::unordered_map<subdomain_id_type, subdomain_id_type> id_remapping;
      if (remap_subdomain_ids)
        {
          const auto & this_map = this->get_subdomain_name_map();
          const auto & other_map = other_mesh->get_subdomain_name_map();
          std::unordered_map<std::string, subdomain_id_type> other_map_reversed;
          for (auto & [sid, sname] : other_map)
            other_map_reversed.emplace(sname, sid);

          std::unordered_map<std::string, subdomain_id_type> this_map_reversed;
          for (auto & [sid, sname] : this_map)
            this_map_reversed.emplace(sname, sid);

          // We don't require either mesh to be prepared, but that
          // means we need to check for subdomains manually.
          auto get_subdomains = [](const MeshBase & mesh) {
            std::set<subdomain_id_type> all_subdomains;
            for (auto & el : mesh.element_ptr_range())
              all_subdomains.insert(el->subdomain_id());
            return all_subdomains;
          };

          const auto this_subdomains = get_subdomains(*this);
          const auto other_subdomains = get_subdomains(*other_mesh);

          for (auto & [sid, sname] : this_map)
            {
              // The same name with the same id means we're fine.  The
              // same name with another id means we remap their id to
              // ours
              if (const auto other_reverse_it = other_map_reversed.find(sname);
                  other_reverse_it != other_map_reversed.end() && other_reverse_it->second != sid)
                id_remapping[other_reverse_it->second] = sid;

              // The same id with a different name, we'll get to
              // later.  The same id without any name means we don't
              // know what the user wants.
              if (other_subdomains.count(sid) && !other_map.count(sid))
                libmesh_error_msg("Can't safely stitch with a mesh sharing subdomain id "
                                  << sid << " but not subdomain name " << sname);
            }

          subdomain_id_type next_free_id = 0;
          // We might try to stitch empty meshes ...
          if (!this_subdomains.empty())
            next_free_id = *this_subdomains.rbegin() + 1;
          if (!other_subdomains.empty())
            next_free_id =
              std::max(next_free_id,
                       cast_int<subdomain_id_type>
                         (*other_subdomains.rbegin() + 1));

          for (auto & [sid, sname] : other_map)
            {
              // At this point we've figured out any remapping
              // necessary for an sname that we share.  And we don't
              // need to remap any sid we don't share.
              if (!this_map_reversed.count(sname))
                {
                  // But if we don't have this sname and we do have this
                  // sid then we can't just merge into that.
                  if (this_subdomains.count(sid))
                    {
                      // If we have this sid with no name, we don't
                      // know what the user wants.
                      if (!this_map.count(sid))
                        libmesh_error_msg("Can't safely stitch with a mesh sharing subdomain id "
                                          << sid << " but under subdomain name " << sname);

                      // We have this sid under a different name, so
                      // we just need to give the other elements a new
                      // id.

                      // Users might have done crazy things with id
                      // choice so let's make sure they didn't get too
                      // crazy.
                      libmesh_error_msg_if ((!this_subdomains.empty() &&
                                             next_free_id < *this_subdomains.rbegin()) ||
                                            (!other_subdomains.empty() &&
                                             next_free_id < *other_subdomains.rbegin()),
                                            "Subdomain id overflow");

                      id_remapping[sid] = next_free_id++;
                      this->subdomain_name(next_free_id) = sname;
                    }
                  // If we don't have this subdomain id, well, we're
                  // about to, so we should have its name too.
                  else
                    this->subdomain_name(sid) = sname;
                }
            }
        }

      // Copy mesh data. If we skip the call to find_neighbors(), the lists
      // of neighbors will be copied verbatim from the other mesh
      this->copy_nodes_and_elements(*other_mesh, skip_find_neighbors,
                                    elem_delta, node_delta,
                                    unique_delta, &id_remapping);

      // Copy BoundaryInfo from other_mesh too.  We do this via the
      // list APIs rather than element-by-element for speed.
      BoundaryInfo & boundary = this->get_boundary_info();
      const BoundaryInfo & other_boundary = other_mesh->get_boundary_info();

      for (const auto & t : other_boundary.build_node_list())
        boundary.add_node(std::get<0>(t) + node_delta,
                          std::get<1>(t));

      for (const auto & t : other_boundary.build_side_list())
        boundary.add_side(std::get<0>(t) + elem_delta,
                          std::get<1>(t),
                          std::get<2>(t));

      for (const auto & t : other_boundary.build_edge_list())
        boundary.add_edge(std::get<0>(t) + elem_delta,
                          std::get<1>(t),
                          std::get<2>(t));

      for (const auto & t : other_boundary.build_shellface_list())
        boundary.add_shellface(std::get<0>(t) + elem_delta,
                               std::get<1>(t),
                               std::get<2>(t));

      const auto & other_ns_id_to_name = other_boundary.get_nodeset_name_map();
      auto & ns_id_to_name = boundary.set_nodeset_name_map();
      ns_id_to_name.insert(other_ns_id_to_name.begin(), other_ns_id_to_name.end());

      const auto & other_ss_id_to_name = other_boundary.get_sideset_name_map();
      auto & ss_id_to_name = boundary.set_sideset_name_map();
      ss_id_to_name.insert(other_ss_id_to_name.begin(), other_ss_id_to_name.end());

      const auto & other_es_id_to_name = other_boundary.get_edgeset_name_map();
      auto & es_id_to_name = boundary.set_edgeset_name_map();
      es_id_to_name.insert(other_es_id_to_name.begin(), other_es_id_to_name.end());

      // Merge other_mesh's elemset information with ours. Throw an
      // error if this and other_mesh have overlapping elemset codes
      // that refer to different elemset ids.
      std::vector<dof_id_type> this_elemset_codes = this->get_elemset_codes();
      MeshBase::elemset_type this_id_set_to_fill, other_id_set_to_fill;
      for (const auto & elemset_code : other_mesh->get_elemset_codes())
        {
          // Get the elemset ids for this elemset_code on other_mesh
          other_mesh->get_elemsets(elemset_code, other_id_set_to_fill);

          // Check that this elemset code does not already exist
          // in this mesh, or if it does, that it has the same elemset
          // ids associated with it.
          //
          // Note: get_elemset_codes() is guaranteed to return a
          // sorted vector, so we can binary search in it.
          auto it = Utility::binary_find(this_elemset_codes.begin(),
                                         this_elemset_codes.end(),
                                         elemset_code);

          if (it != this_elemset_codes.end())
            {
              // This mesh has the same elemset code. Does it refer to
              // the same elemset ids?
              this->get_elemsets(elemset_code, this_id_set_to_fill);

              // Throw an error if they don't match, otherwise we
              // don't need to do anything
              libmesh_error_msg_if(other_id_set_to_fill != this_id_set_to_fill,
                                   "Attempted to stitch together meshes with conflicting elemset codes.");
            }
          else
            {
              // Add other_mesh's elemset code to this mesh
              this->add_elemset_code(elemset_code, other_id_set_to_fill);
            }
        }

    } // end if (other_mesh)

  // Finally, we need to "merge" the overlapping nodes
  // We do this by iterating over node_to_elems_map and updating
  // the elements so that they "point" to the nodes that came
  // from this mesh, rather than from other_mesh.
  // Then we iterate over node_to_node_map and delete the
  // duplicate nodes that came from other_mesh.

  {
    LOG_SCOPE("stitch_meshes node updates", "UnstructuredMesh");

    // Container to catch boundary IDs passed back from BoundaryInfo.
    std::vector<boundary_id_type> bc_ids;

    for (const auto & [target_node_id, elem_vec] : node_to_elems_map)
      {
        dof_id_type other_node_id = node_to_node_map[target_node_id];
        Node & target_node = this->node_ref(target_node_id);

        std::size_t n_elems = elem_vec.size();
        for (std::size_t i=0; i<n_elems; i++)
          {
            dof_id_type elem_id = elem_vec[i];
            Elem * el = this->elem_ptr(elem_id);

            // find the local node index that we want to update
            unsigned int local_node_index = el->local_node(other_node_id);
            libmesh_assert_not_equal_to(local_node_index, libMesh::invalid_uint);

            // We also need to copy over the nodeset info here,
            // because the node will get deleted below
            this->get_boundary_info().boundary_ids(el->node_ptr(local_node_index), bc_ids);
            el->set_node(local_node_index, &target_node);
            this->get_boundary_info().add_node(&target_node, bc_ids);
          }
      }
  }

  {
    LOG_SCOPE("stitch_meshes node deletion", "UnstructuredMesh");
    for (const auto & [other_node_id, this_node_id] : node_to_node_map)
      {
        // In the case that this==other_mesh, the two nodes might be the same (e.g. if
        // we're stitching a "sliver"), hence we need to skip node deletion in that case.
        if ((this == other_mesh) && (this_node_id == other_node_id))
          continue;

        this->delete_node( this->node_ptr(this_node_id) );
      }
  }

  // If find_neighbors() wasn't called in prepare_for_use(), we need to
  // manually loop once more over all elements adjacent to the stitched boundary
  // and fix their lists of neighbors.
  // This is done according to the following steps:
  //   1. Loop over all copied elements adjacent to the boundary using node_to_elems_map (trying to avoid duplicates)
  //   2. Look at all their sides with a nullptr neighbor and update them using side_to_elem_map if necessary
  //   3. Update the corresponding side in side_to_elem_map as well
  if (skip_find_neighbors)
    {
      LOG_SCOPE("stitch_meshes neighbor fixes", "UnstructuredMesh");

      // Pull objects out of the loop to reduce heap operations
      std::unique_ptr<const Elem> my_side, their_side;

      std::set<dof_id_type> fixed_elems;
      for (const auto & pr : node_to_elems_map)
        {
          std::size_t n_elems = pr.second.size();
          for (std::size_t i=0; i<n_elems; i++)
            {
              dof_id_type elem_id = pr.second[i];
              if (!fixed_elems.count(elem_id))
                {
                  Elem * el = this->elem_ptr(elem_id);
                  fixed_elems.insert(elem_id);
                  for (auto s : el->side_index_range())
                    {
                      bool has_real_neighbor = (el->neighbor_ptr(s) != nullptr);
                      bool has_disdjoint_neighbor = is_valid_disjoint_pair_to_stitch &&
                      (this->get_boundary_info().has_boundary_id(el, s, this_mesh_boundary_id)
                      || this->get_boundary_info().has_boundary_id(el, s, other_mesh_boundary_id));

                      if (!has_real_neighbor || has_disdjoint_neighbor)
                        {
                          key_type key = el->low_order_key(s);
                          auto bounds = side_to_elem_map.equal_range(key);

                          if (bounds.first != bounds.second)
                            {
                              // Get the side for this element
                              el->side_ptr(my_side, s);

                              // Look at all the entries with an equivalent key
                              while (bounds.first != bounds.second)
                                {
                                  // Get the potential element
                                  Elem * neighbor = const_cast<Elem *>(bounds.first->second.first);

                                  // Get the side for the neighboring element
                                  const unsigned int ns = bounds.first->second.second;
                                  neighbor->side_ptr(their_side, ns);
                                  //libmesh_assert(my_side.get());
                                  //libmesh_assert(their_side.get());

                                  // If found a match with my side
                                  //
                                  // We need special tests here for 1D:
                                  // since parents and children have an equal
                                  // side (i.e. a node), we need to check
                                  // ns != ms, and we also check level() to
                                  // avoid setting our neighbor pointer to
                                  // any of our neighbor's descendants
                                  if ((*my_side == *their_side) &&
                                      (el->level() == neighbor->level()) &&
                                      ((el->dim() != 1) || (ns != s)))
                                    {
                                      // So share a side.  Is this a mixed pair
                                      // of subactive and active/ancestor
                                      // elements?
                                      // If not, then we're neighbors.
                                      // If so, then the subactive's neighbor is

                                      if (el->subactive() ==
                                          neighbor->subactive())
                                        {
                                          // an element is only subactive if it has
                                          // been coarsened but not deleted
                                          el->set_neighbor (s,neighbor);
                                          neighbor->set_neighbor(ns,el);
                                        }
                                      else if (el->subactive())
                                        {
                                          el->set_neighbor(s,neighbor);
                                        }
                                      else if (neighbor->subactive())
                                        {
                                          neighbor->set_neighbor(ns,el);
                                        }
                                      // It's OK to invalidate the
                                      // bounds.first iterator here,
                                      // as we are immediately going
                                      // to break out of this while
                                      // loop. bounds.first will
                                      // therefore not be used for
                                      // anything else.
                                      side_to_elem_map.erase (bounds.first);
                                      break;
                                    }

                                  ++bounds.first;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

#ifdef LIBMESH_ENABLE_PERIODIC
  // Remove only the disjoint pair that was actually stitched.
  // Safe because `is_valid_disjoint_pair_to_stitch` is true
  // only if this exact (a,b) pair exists in the registry.
  // Other disjoint pairs remain untouched.
  if (is_valid_disjoint_pair_to_stitch)
    this->remove_disjoint_boundary_pair(this_mesh_boundary_id, other_mesh_boundary_id);
#endif

  if (prepare_after_stitching)
    {
      // We set our new neighbor pointers already
      const bool old_allow_find_neighbors = this->allow_find_neighbors();
      this->allow_find_neighbors(!skip_find_neighbors);

      // We haven't newly remoted any elements
      const bool old_allow_remote_element_removal = this->allow_remote_element_removal();
      this->allow_remote_element_removal(false);

      this->prepare_for_use();

      this->allow_find_neighbors(old_allow_find_neighbors);
      this->allow_remote_element_removal(old_allow_remote_element_removal);
    }

  // After the stitching, we may want to clear boundary IDs from element
  // faces that are now internal to the mesh
  if (clear_stitched_boundary_ids)
    {
      LOG_SCOPE("stitch_meshes clear bcids", "UnstructuredMesh");

      this->get_boundary_info().clear_stitched_boundary_side_ids(
          this_mesh_boundary_id, other_mesh_boundary_id, /*clear_nodeset_data=*/true);
    }

  // Return the number of nodes which were merged.
  return node_to_node_map.size();
}


} // namespace libMesh

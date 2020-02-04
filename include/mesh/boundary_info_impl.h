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

#ifndef LIBMESH_BOUNDARY_INFO_IMPL_H
#define LIBMESH_BOUNDARY_INFO_IMPL_H

// C++ includes
#include <iterator>  // std::distance

// Local includes
#include "libmesh/libmesh_config.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/boundary_info.h"
#include "libmesh/distributed_mesh.h"
#include "libmesh/elem.h"
#include "libmesh/mesh_communication.h"
#include "libmesh/mesh_serializer.h"
#include "libmesh/parallel.h"
#include "libmesh/partitioner.h"
#include "libmesh/remote_elem.h"
#include "libmesh/unstructured_mesh.h"

namespace
{

// Templated helper function for removing a subset of keys from a
// multimap that further satisfy a given predicate on the
// corresponding values.
template <class Key, class T, class Pred>
void erase_if(std::multimap<Key,T> & map, Key k, Pred pred)
{
  auto rng = map.equal_range(k);
  auto it = rng.first;
  while (it != rng.second)
    {
      if (pred(it->second))
        it = map.erase(it);
      else
        ++it;
    }
}

// Similar to the helper function above but doesn't take a key,
// instead it applies the predicate to every value in the map.
template <class Key, class T, class Pred>
void erase_if(std::multimap<Key,T> & map, Pred pred)
{
  auto it = map.begin();
  while (it != map.end())
    {
      if (pred(it->second))
        it = map.erase(it);
      else
        ++it;
    }
}

}

namespace libMesh
{

//------------------------------------------------------
// BoundaryInfo functions
template <typename RealType>
BoundaryInfoTempl<RealType>::BoundaryInfoTempl(MeshBase & m) :
  ParallelObject(m.comm()),
  _mesh (m)
{
}

template <typename RealType>
BoundaryInfoTempl<RealType> & BoundaryInfoTempl<RealType>::operator=(const BoundaryInfo & other_boundary_info)
{
  // Overwrite any preexisting boundary info
  this->clear();

  /**
   * We're going to attempt to pull _new_ pointers out of the mesh
   * assigned to this boundary info.
   *
   * This will only work if the mesh assigned to this BoundaryInfo is
   * the same mesh object as other_boundary_info _or_ was constructed
   * in exactly the same way (or constructed as a copy, or a refined
   * copy without renumbering, etc.).
   */

  // Copy node boundary info
  for (const auto & pr : other_boundary_info._boundary_node_id)
    _boundary_node_id.insert(std::make_pair(_mesh.node_ptr(pr.first->id()),
                                            pr.second));

  // Copy edge boundary info
  for (const auto & pr : other_boundary_info._boundary_edge_id)
    _boundary_edge_id.insert(std::make_pair(_mesh.elem_ptr(pr.first->id()),
                                            pr.second));

  // Copy shellface boundary info
  for (const auto & pr : other_boundary_info._boundary_shellface_id)
    _boundary_shellface_id.insert(std::make_pair(_mesh.elem_ptr(pr.first->id()),
                                                 pr.second));

  // Copy side boundary info
  for (const auto & pr : other_boundary_info._boundary_side_id)
    _boundary_side_id.insert(std::make_pair(_mesh.elem_ptr(pr.first->id()),
                                            pr.second));

  _boundary_ids = other_boundary_info._boundary_ids;
  _side_boundary_ids = other_boundary_info._side_boundary_ids;
  _node_boundary_ids = other_boundary_info._node_boundary_ids;
  _edge_boundary_ids = other_boundary_info._edge_boundary_ids;
  _shellface_boundary_ids = other_boundary_info._shellface_boundary_ids;

  return *this;
}


template <typename RealType>
BoundaryInfoTempl<RealType>::~BoundaryInfoTempl()
{
  this->clear();
}



template <typename RealType>
void BoundaryInfoTempl<RealType>::clear()
{
  _boundary_node_id.clear();
  _boundary_side_id.clear();
  _boundary_edge_id.clear();
  _boundary_shellface_id.clear();
  _boundary_ids.clear();
  _side_boundary_ids.clear();
  _node_boundary_ids.clear();
  _edge_boundary_ids.clear();
  _shellface_boundary_ids.clear();
}



template <typename RealType>
void BoundaryInfoTempl<RealType>::regenerate_id_sets()
{
  // Clear the old caches
  _boundary_ids.clear();
  _side_boundary_ids.clear();
  _node_boundary_ids.clear();
  _edge_boundary_ids.clear();
  _shellface_boundary_ids.clear();

  // Loop over id maps to regenerate each set.
  for (const auto & pr : _boundary_node_id)
    {
      const boundary_id_type id = pr.second;
      _boundary_ids.insert(id);
      _node_boundary_ids.insert(id);
    }

  for (const auto & pr : _boundary_edge_id)
    {
      const boundary_id_type id = pr.second.second;
      _boundary_ids.insert(id);
      _edge_boundary_ids.insert(id);
    }

  for (const auto & pr : _boundary_side_id)
    {
      const boundary_id_type id = pr.second.second;
      _boundary_ids.insert(id);
      _side_boundary_ids.insert(id);
    }

  for (const auto & pr : _boundary_shellface_id)
    {
      const boundary_id_type id = pr.second.second;
      _boundary_ids.insert(id);
      _shellface_boundary_ids.insert(id);
    }
}



template <typename RealType>
void BoundaryInfoTempl<RealType>::sync (UnstructuredMesh & boundary_mesh)
{
  std::set<boundary_id_type> request_boundary_ids(_boundary_ids);
  request_boundary_ids.insert(invalid_id);
  if (!_mesh.is_serial())
    this->comm().set_union(request_boundary_ids);

  this->sync(request_boundary_ids,
             boundary_mesh);
}



template <typename RealType>
void BoundaryInfoTempl<RealType>::sync (const std::set<boundary_id_type> & requested_boundary_ids,
                         UnstructuredMesh & boundary_mesh)
{
  // Call the 3 argument version of this function with a dummy value for the third set.
  std::set<subdomain_id_type> subdomains_relative_to;
  subdomains_relative_to.insert(Elem::invalid_subdomain_id);

  this->sync(requested_boundary_ids,
             boundary_mesh,
             subdomains_relative_to);
}



template <typename RealType>
void BoundaryInfoTempl<RealType>::sync (const std::set<boundary_id_type> & requested_boundary_ids,
                         UnstructuredMesh & boundary_mesh,
                         const std::set<subdomain_id_type> & subdomains_relative_to)
{
  LOG_SCOPE("sync()", "BoundaryInfo");

  boundary_mesh.clear();

  /**
   * Deleting 0 elements seems weird, but it's better encapsulating
   * than exposing a set_is_serial(false) capability that might be
   * easily misused.
   */
  if (!_mesh.is_serial())
    boundary_mesh.delete_remote_elements();

  /**
   * If the boundary_mesh is still serial, that means we *can't*
   * parallelize it, so to make sure we can construct it in full on
   * every processor we'll serialize the interior mesh.  Use a
   * temporary serializer here.
   */
  MeshSerializer serializer
    (const_cast<MeshBase &>(_mesh), boundary_mesh.is_serial());

  /**
   * Re-create the boundary mesh.
   */

  boundary_mesh.set_n_partitions() = _mesh.n_partitions();

  std::map<dof_id_type, dof_id_type> node_id_map;
  std::map<std::pair<dof_id_type, unsigned char>, dof_id_type> side_id_map;

  this->_find_id_maps(requested_boundary_ids, 0, &node_id_map, 0, &side_id_map, subdomains_relative_to);

  // Let's add all the boundary nodes we found to the boundary mesh
  for (const auto & node : _mesh.node_ptr_range())
    {
      dof_id_type node_id = node->id();
      if (node_id_map.count(node_id))
        {
          boundary_mesh.add_point(*node, node_id_map[node_id], node->processor_id());

          // Copy over all the node's boundary IDs to boundary_mesh
          std::vector<boundary_id_type> node_boundary_ids;
          this->boundary_ids(node, node_boundary_ids);
          for (const auto & node_bid : node_boundary_ids)
            boundary_mesh.get_boundary_info().add_node(node_id_map[node_id], node_bid);
        }
    }

  // Let's add the elements
  this->add_elements (requested_boundary_ids, boundary_mesh, subdomains_relative_to);

  // The new elements are currently using the interior mesh's nodes;
  // we want them to use the boundary mesh's nodes instead.

  // This side's Node pointers still point to the nodes of the original mesh.
  // We need to re-point them to the boundary mesh's nodes!  Since we copied *ALL* of
  // the original mesh's nodes over, we should be guaranteed to have the same ordering.
  for (auto & new_elem : boundary_mesh.element_ptr_range())
    {
      for (auto nn : new_elem->node_index_range())
        {
          // Get the correct node pointer, based on the id()
          Node * new_node =
            boundary_mesh.node_ptr(node_id_map[new_elem->node_id(nn)]);

          // sanity check: be sure that the new Node exists and its
          // global id really matches
          libmesh_assert (new_node);
          libmesh_assert_equal_to (new_node->id(),
                                   node_id_map[new_elem->node_id(nn)]);

          // Assign the new node pointer
          new_elem->set_node(nn) = new_node;
        }
    }

  // Don't repartition this mesh; we want it to stay in sync with the
  // interior partitioning.
  boundary_mesh.partitioner().reset(nullptr);

  // Make boundary_mesh nodes and elements contiguous
  boundary_mesh.prepare_for_use(/*skip_renumber =*/ false);

  // and finally distribute element partitioning to the nodes
  Partitioner::set_node_processor_ids(boundary_mesh);
}


template <typename RealType>
void BoundaryInfoTempl<RealType>::get_side_and_node_maps (UnstructuredMesh & boundary_mesh,
                                           std::map<dof_id_type, dof_id_type> & node_id_map,
                                           std::map<dof_id_type, unsigned char> & side_id_map,
                                           Real tolerance)
{
  LOG_SCOPE("get_side_and_node_maps()", "BoundaryInfo");

  node_id_map.clear();
  side_id_map.clear();

  // Pull objects out of the loop to reduce heap operations
  std::unique_ptr<const Elem> interior_parent_side;

  for (const auto & boundary_elem : boundary_mesh.active_element_ptr_range())
    {
      const Elem * interior_parent = boundary_elem->interior_parent();

      // Find out which side of interior_parent boundary_elem corresponds to.
      // Use centroid comparison as a way to check.
      unsigned char interior_parent_side_index = 0;
      bool found_matching_sides = false;
      for (auto side : interior_parent->side_index_range())
        {
          interior_parent->build_side_ptr(interior_parent_side, side);
          Real centroid_distance = (boundary_elem->centroid() - interior_parent_side->centroid()).norm();

          if (centroid_distance < (tolerance * boundary_elem->hmin()))
            {
              interior_parent_side_index = cast_int<unsigned char>(side);
              found_matching_sides = true;
              break;
            }
        }

      if (!found_matching_sides)
        libmesh_error_msg("No matching side found within the specified tolerance");

      side_id_map[boundary_elem->id()] = interior_parent_side_index;

      for (auto local_node_index : boundary_elem->node_index_range())
        {
          dof_id_type boundary_node_id = boundary_elem->node_id(local_node_index);
          dof_id_type interior_node_id = interior_parent_side->node_id(local_node_index);

          node_id_map[interior_node_id] = boundary_node_id;
        }
    }
}



template <typename RealType>
void BoundaryInfoTempl<RealType>::add_elements(const std::set<boundary_id_type> & requested_boundary_ids,
                                UnstructuredMesh & boundary_mesh)
{
  // Call the 3 argument version of this function with a dummy value for the third arg.
  std::set<subdomain_id_type> subdomains_relative_to;
  subdomains_relative_to.insert(Elem::invalid_subdomain_id);

  this->add_elements(requested_boundary_ids,
                     boundary_mesh,
                     subdomains_relative_to);
}



template <typename RealType>
void BoundaryInfoTempl<RealType>::add_elements(const std::set<boundary_id_type> & requested_boundary_ids,
                                UnstructuredMesh & boundary_mesh,
                                const std::set<subdomain_id_type> & subdomains_relative_to)
{
  LOG_SCOPE("add_elements()", "BoundaryInfo");

  // We're not prepared to mix serial and distributed meshes in this
  // method, so make sure they match from the start.
  libmesh_assert_equal_to(_mesh.is_serial(),
                          boundary_mesh.is_serial());

  std::map<std::pair<dof_id_type, unsigned char>, dof_id_type> side_id_map;
  this->_find_id_maps(requested_boundary_ids,
                      0,
                      nullptr,
                      boundary_mesh.max_elem_id(),
                      &side_id_map,
                      subdomains_relative_to);

  // We have to add sides *outside* any element loop, because if
  // boundary_mesh and _mesh are the same then those additions can
  // invalidate our element iterators.  So we just use the element
  // loop to make a list of sides to add.
  typedef std::vector<std::pair<dof_id_type, unsigned char>>
    side_container;
  side_container sides_to_add;

  for (const auto & elem : _mesh.element_ptr_range())
    {
      // If the subdomains_relative_to container has the
      // invalid_subdomain_id, we fall back on the "old" behavior of
      // adding sides regardless of this Elem's subdomain. Otherwise,
      // if the subdomains_relative_to container doesn't contain the
      // current Elem's subdomain_id(), we won't add any sides from
      // it.
      if (!subdomains_relative_to.count(Elem::invalid_subdomain_id) &&
          !subdomains_relative_to.count(elem->subdomain_id()))
        continue;

      // Get the top-level parent for this element
      const Elem * top_parent = elem->top_parent();

      // Find all the boundary side ids for this Elem.
      auto bounds = _boundary_side_id.equal_range(top_parent);

      for (auto s : elem->side_index_range())
        {
          bool add_this_side = false;
          boundary_id_type this_bcid = invalid_id;

          for (const auto & pr : as_range(bounds))
            {
              this_bcid = pr.second.second;

              // if this side is flagged with a boundary condition
              // and the user wants this id
              if ((pr.second.first == s) &&
                  (requested_boundary_ids.count(this_bcid)))
                {
                  add_this_side = true;
                  break;
                }
            }

          // We may still want to add this side if the user called
          // sync() with no requested_boundary_ids. This corresponds
          // to the "old" style of calling sync() in which the entire
          // boundary was copied to the BoundaryMesh, and handles the
          // case where elements on the geometric boundary are not in
          // any sidesets.
          if (bounds.first == bounds.second            &&
              requested_boundary_ids.count(invalid_id) &&
              elem->neighbor_ptr(s) == nullptr)
            add_this_side = true;

          if (add_this_side)
            sides_to_add.push_back(std::make_pair(elem->id(), s));
        }
    }

#ifdef LIBMESH_ENABLE_UNIQUE_ID
  unique_id_type old_max_unique_id = boundary_mesh.parallel_max_unique_id();
#endif

  for (const auto & pr : sides_to_add)
    {
      const dof_id_type elem_id = pr.first;
      const unsigned char s = pr.second;
      Elem * elem = _mesh.elem_ptr(elem_id);

      // Build the side - do not use a "proxy" element here:
      // This will be going into the boundary_mesh and needs to
      // stand on its own.
      std::unique_ptr<Elem> side (elem->build_side_ptr(s, false));

      side->processor_id() = elem->processor_id();

      const std::pair<dof_id_type, unsigned char> side_pair(elem_id, s);

      libmesh_assert(side_id_map.count(side_pair));

      const dof_id_type new_side_id = side_id_map[side_pair];

      side->set_id(new_side_id);

#ifdef LIBMESH_ENABLE_UNIQUE_ID
      side->set_unique_id() = old_max_unique_id + new_side_id;
#endif

      // Add the side
      Elem * new_elem = boundary_mesh.add_elem(side.release());

#ifdef LIBMESH_ENABLE_AMR
      // Set parent links
      if (elem->parent())
        {
          const std::pair<dof_id_type, unsigned char> parent_side_pair(elem->parent()->id(), s);

          libmesh_assert(side_id_map.count(parent_side_pair));

          Elem * side_parent = boundary_mesh.elem_ptr(side_id_map[parent_side_pair]);

          libmesh_assert(side_parent);

          new_elem->set_parent(side_parent);

          side_parent->set_refinement_flag(Elem::INACTIVE);

          // Figuring out which child we are of our parent
          // is a trick.  Due to libMesh child numbering
          // conventions, if we are an element on a vertex,
          // then we share that vertex with our parent, with
          // the same local index.
          bool found_child = false;
          for (auto v : IntRange<unsigned int>(0, new_elem->n_vertices()))
            if (new_elem->node_ptr(v) == side_parent->node_ptr(v))
              {
                side_parent->add_child(new_elem, v);
                found_child = true;
              }

          // If we don't share any vertex with our parent,
          // then we're the fourth child (index 3) of a
          // triangle.
          if (!found_child)
            {
              libmesh_assert_equal_to (new_elem->n_vertices(), 3);
              side_parent->add_child(new_elem, 3);
            }
        }
#endif

      new_elem->set_interior_parent (elem);

      // On non-local elements on DistributedMesh we might have
      // RemoteElem neighbor links to construct
      if (!_mesh.is_serial() &&
          (elem->processor_id() != this->processor_id()))
        {
          const unsigned short n_nodes = elem->n_nodes();

          const unsigned short bdy_n_sides = new_elem->n_sides();
          const unsigned short bdy_n_nodes = new_elem->n_nodes();

          // Check every interior side for a RemoteElem
          for (auto interior_side : elem->side_index_range())
            {
              // Might this interior side have a RemoteElem that
              // needs a corresponding Remote on a boundary side?
              if (elem->neighbor_ptr(interior_side) != RemoteElem::get_instance())
                continue;

              // Which boundary side?
              for (unsigned short boundary_side = 0;
                   boundary_side != bdy_n_sides; ++boundary_side)
                {
                  // Look for matching node points.  This is safe in
                  // *this* context.
                  bool found_all_nodes = true;
                  for (unsigned short boundary_node = 0;
                       boundary_node != bdy_n_nodes; ++boundary_node)
                    {
                      if (!new_elem->is_node_on_side(boundary_node,
                                                     boundary_side))
                        continue;

                      bool found_this_node = false;
                      for (unsigned short interior_node = 0;
                           interior_node != n_nodes; ++interior_node)
                        {
                          if (!elem->is_node_on_side(interior_node,
                                                     interior_side))
                            continue;

                          if (new_elem->point(boundary_node) ==
                              elem->point(interior_node))
                            {
                              found_this_node = true;
                              break;
                            }
                        }
                      if (!found_this_node)
                        {
                          found_all_nodes = false;
                          break;
                        }
                    }

                  if (found_all_nodes)
                    {
                      new_elem->set_neighbor
                        (boundary_side,
                         const_cast<RemoteElem *>(RemoteElem::get_instance()));
                      break;
                    }
                }
            }
        }
    }

  // We haven't been bothering to keep unique ids consistent on ghost
  // elements
  if (!boundary_mesh.is_serial())
    MeshCommunication().make_node_unique_ids_parallel_consistent(boundary_mesh);

  // Make sure we didn't add ids inconsistently
#ifdef DEBUG
# ifdef LIBMESH_HAVE_RTTI
  DistributedMesh * parmesh = dynamic_cast<DistributedMesh *>(&boundary_mesh);
  if (parmesh)
    parmesh->libmesh_assert_valid_parallel_ids();
# endif
#endif
}



template <typename RealType>
void BoundaryInfoTempl<RealType>::add_node(const dof_id_type node_id,
                            const boundary_id_type id)
{
  const Node * node_ptr = _mesh.query_node_ptr(node_id);

  // The user could easily ask for an invalid node id, so let's throw
  // an easy-to-understand error message when this happens.
  if (!node_ptr)
    libmesh_error_msg("BoundaryInfo::add_node(): Could not retrieve pointer for node " << node_id << ", no boundary id was added.");

  this->add_node (node_ptr, id);
}



template <typename RealType>
void BoundaryInfoTempl<RealType>::add_node(const Node * node,
                            const boundary_id_type id)
{
  if (id == invalid_id)
    libmesh_error_msg("ERROR: You may not set a boundary ID of "   \
                      << invalid_id                                \
                      << "\n That is reserved for internal use.");

  // Don't add the same ID twice
  for (const auto & pr : as_range(_boundary_node_id.equal_range(node)))
    if (pr.second == id)
      return;

  _boundary_node_id.insert(std::make_pair(node, id));
  _boundary_ids.insert(id);
  _node_boundary_ids.insert(id); // Also add this ID to the set of node boundary IDs
}



template <typename RealType>
void BoundaryInfoTempl<RealType>::add_node(const Node * node,
                            const std::vector<boundary_id_type> & ids)
{
  if (ids.empty())
    return;

  libmesh_assert(node);

  // Don't add the same ID twice
  auto bounds = _boundary_node_id.equal_range(node);

  // The entries in the ids vector may be non-unique.  If we expected
  // *lots* of ids, it might be fastest to construct a std::set from
  // the entries, but for a small number of entries, which is more
  // typical, it is probably faster to copy the vector and do sort+unique.
  // http://stackoverflow.com/questions/1041620/whats-the-most-efficient-way-to-erase-duplicates-and-sort-a-vector
  std::vector<boundary_id_type> unique_ids(ids.begin(), ids.end());
  std::sort(unique_ids.begin(), unique_ids.end());
  std::vector<boundary_id_type>::iterator new_end =
    std::unique(unique_ids.begin(), unique_ids.end());

  for (auto & id : as_range(unique_ids.begin(), new_end))
    {
      if (id == invalid_id)
        libmesh_error_msg("ERROR: You may not set a boundary ID of "    \
                          << invalid_id                                 \
                          << "\n That is reserved for internal use.");

      bool already_inserted = false;
      for (const auto & pr : as_range(bounds))
        if (pr.second == id)
          {
            already_inserted = true;
            break;
          }
      if (already_inserted)
        continue;

      _boundary_node_id.insert(std::make_pair(node,id));
      _boundary_ids.insert(id);
      _node_boundary_ids.insert(id); // Also add this ID to the set of node boundary IDs
    }
}



template <typename RealType>
void BoundaryInfoTempl<RealType>::clear_boundary_node_ids()
{
  _boundary_node_id.clear();
}

template <typename RealType>
void BoundaryInfoTempl<RealType>::add_edge(const dof_id_type e,
                            const unsigned short int edge,
                            const boundary_id_type id)
{
  this->add_edge (_mesh.elem_ptr(e), edge, id);
}



template <typename RealType>
void BoundaryInfoTempl<RealType>::add_edge(const Elem * elem,
                            const unsigned short int edge,
                            const boundary_id_type id)
{
  libmesh_assert(elem);

  // Only add BCs for level-0 elements.
  libmesh_assert_equal_to (elem->level(), 0);

  if (id == invalid_id)
    libmesh_error_msg("ERROR: You may not set a boundary ID of "        \
                      << invalid_id                                     \
                      << "\n That is reserved for internal use.");

  // Don't add the same ID twice
  for (const auto & pr : as_range(_boundary_edge_id.equal_range(elem)))
    if (pr.second.first == edge &&
        pr.second.second == id)
      return;

  _boundary_edge_id.insert(std::make_pair(elem, std::make_pair(edge, id)));
  _boundary_ids.insert(id);
  _edge_boundary_ids.insert(id); // Also add this ID to the set of edge boundary IDs
}



template <typename RealType>
void BoundaryInfoTempl<RealType>::add_edge(const Elem * elem,
                            const unsigned short int edge,
                            const std::vector<boundary_id_type> & ids)
{
  if (ids.empty())
    return;

  libmesh_assert(elem);

  // Only add BCs for level-0 elements.
  libmesh_assert_equal_to (elem->level(), 0);

  // Don't add the same ID twice
  auto bounds = _boundary_edge_id.equal_range(elem);

  // The entries in the ids vector may be non-unique.  If we expected
  // *lots* of ids, it might be fastest to construct a std::set from
  // the entries, but for a small number of entries, which is more
  // typical, it is probably faster to copy the vector and do sort+unique.
  // http://stackoverflow.com/questions/1041620/whats-the-most-efficient-way-to-erase-duplicates-and-sort-a-vector
  std::vector<boundary_id_type> unique_ids(ids.begin(), ids.end());
  std::sort(unique_ids.begin(), unique_ids.end());
  std::vector<boundary_id_type>::iterator new_end =
    std::unique(unique_ids.begin(), unique_ids.end());

  for (auto & id : as_range(unique_ids.begin(), new_end))
    {
      if (id == invalid_id)
        libmesh_error_msg("ERROR: You may not set a boundary ID of "   \
                          << invalid_id                                \
                          << "\n That is reserved for internal use.");

      bool already_inserted = false;
      for (const auto & pr : as_range(bounds))
        if (pr.second.first == edge &&
            pr.second.second == id)
          {
            already_inserted = true;
            break;
          }
      if (already_inserted)
        continue;

      _boundary_edge_id.insert(std::make_pair(elem, std::make_pair(edge, id)));
      _boundary_ids.insert(id);
      _edge_boundary_ids.insert(id); // Also add this ID to the set of edge boundary IDs
    }
}



template <typename RealType>
void BoundaryInfoTempl<RealType>::add_shellface(const dof_id_type e,
                                 const unsigned short int shellface,
                                 const boundary_id_type id)
{
  this->add_shellface (_mesh.elem_ptr(e), shellface, id);
}



template <typename RealType>
void BoundaryInfoTempl<RealType>::add_shellface(const Elem * elem,
                                 const unsigned short int shellface,
                                 const boundary_id_type id)
{
  libmesh_assert(elem);

  // Only add BCs for level-0 elements.
  libmesh_assert_equal_to (elem->level(), 0);

  // Shells only have 2 faces
  libmesh_assert_less(shellface, 2);

  if (id == invalid_id)
    libmesh_error_msg("ERROR: You may not set a boundary ID of "        \
                      << invalid_id                                     \
                      << "\n That is reserved for internal use.");

  // Don't add the same ID twice
  for (const auto & pr : as_range(_boundary_shellface_id.equal_range(elem)))
    if (pr.second.first == shellface &&
        pr.second.second == id)
      return;

  _boundary_shellface_id.insert(std::make_pair(elem, std::make_pair(shellface, id)));
  _boundary_ids.insert(id);
  _shellface_boundary_ids.insert(id); // Also add this ID to the set of shellface boundary IDs
}



template <typename RealType>
void BoundaryInfoTempl<RealType>::add_shellface(const Elem * elem,
                                 const unsigned short int shellface,
                                 const std::vector<boundary_id_type> & ids)
{
  if (ids.empty())
    return;

  libmesh_assert(elem);

  // Only add BCs for level-0 elements.
  libmesh_assert_equal_to (elem->level(), 0);

  // Shells only have 2 faces
  libmesh_assert_less(shellface, 2);

  // Don't add the same ID twice
  auto bounds = _boundary_shellface_id.equal_range(elem);

  // The entries in the ids vector may be non-unique.  If we expected
  // *lots* of ids, it might be fastest to construct a std::set from
  // the entries, but for a small number of entries, which is more
  // typical, it is probably faster to copy the vector and do sort+unique.
  // http://stackoverflow.com/questions/1041620/whats-the-most-efficient-way-to-erase-duplicates-and-sort-a-vector
  std::vector<boundary_id_type> unique_ids(ids.begin(), ids.end());
  std::sort(unique_ids.begin(), unique_ids.end());
  std::vector<boundary_id_type>::iterator new_end =
    std::unique(unique_ids.begin(), unique_ids.end());

  for (auto & id : as_range(unique_ids.begin(), new_end))
    {
      if (id == invalid_id)
        libmesh_error_msg("ERROR: You may not set a boundary ID of "   \
                          << invalid_id                                \
                          << "\n That is reserved for internal use.");

      bool already_inserted = false;
      for (const auto & pr : as_range(bounds))
        if (pr.second.first == shellface &&
            pr.second.second == id)
          {
            already_inserted = true;
            break;
          }
      if (already_inserted)
        continue;

      _boundary_shellface_id.insert(std::make_pair(elem, std::make_pair(shellface, id)));
      _boundary_ids.insert(id);
      _shellface_boundary_ids.insert(id); // Also add this ID to the set of shellface boundary IDs
    }
}


template <typename RealType>
void BoundaryInfoTempl<RealType>::add_side(const dof_id_type e,
                            const unsigned short int side,
                            const boundary_id_type id)
{
  this->add_side (_mesh.elem_ptr(e), side, id);
}



template <typename RealType>
void BoundaryInfoTempl<RealType>::add_side(const Elem * elem,
                            const unsigned short int side,
                            const boundary_id_type id)
{
  libmesh_assert(elem);

  // Only add BCs for level-0 elements.
  libmesh_assert_equal_to (elem->level(), 0);

  if (id == invalid_id)
    libmesh_error_msg("ERROR: You may not set a boundary ID of "        \
                      << invalid_id                                     \
                      << "\n That is reserved for internal use.");

  // Don't add the same ID twice
  for (const auto & pr : as_range(_boundary_side_id.equal_range(elem)))
    if (pr.second.first == side &&
        pr.second.second == id)
      return;

  _boundary_side_id.insert(std::make_pair(elem, std::make_pair(side, id)));
  _boundary_ids.insert(id);
  _side_boundary_ids.insert(id); // Also add this ID to the set of side boundary IDs
}



template <typename RealType>
void BoundaryInfoTempl<RealType>::add_side(const Elem * elem,
                            const unsigned short int side,
                            const std::vector<boundary_id_type> & ids)
{
  if (ids.empty())
    return;

  libmesh_assert(elem);

  // Only add BCs for level-0 elements.
  libmesh_assert_equal_to (elem->level(), 0);

  // Don't add the same ID twice
  auto bounds = _boundary_side_id.equal_range(elem);

  // The entries in the ids vector may be non-unique.  If we expected
  // *lots* of ids, it might be fastest to construct a std::set from
  // the entries, but for a small number of entries, which is more
  // typical, it is probably faster to copy the vector and do sort+unique.
  // http://stackoverflow.com/questions/1041620/whats-the-most-efficient-way-to-erase-duplicates-and-sort-a-vector
  std::vector<boundary_id_type> unique_ids(ids.begin(), ids.end());
  std::sort(unique_ids.begin(), unique_ids.end());
  std::vector<boundary_id_type>::iterator new_end =
    std::unique(unique_ids.begin(), unique_ids.end());

  for (auto & id : as_range(unique_ids.begin(), new_end))
    {
      if (id == invalid_id)
        libmesh_error_msg("ERROR: You may not set a boundary ID of "    \
                          << invalid_id                                 \
                          << "\n That is reserved for internal use.");

      bool already_inserted = false;
      for (const auto & pr : as_range(bounds))
        if (pr.second.first == side && pr.second.second == id)
          {
            already_inserted = true;
            break;
          }
      if (already_inserted)
        continue;

      _boundary_side_id.insert(std::make_pair(elem, std::make_pair(side, id)));
      _boundary_ids.insert(id);
      _side_boundary_ids.insert(id); // Also add this ID to the set of side boundary IDs
    }
}



template <typename RealType>
bool BoundaryInfoTempl<RealType>::has_boundary_id(const Node * const node,
                                   const boundary_id_type id) const
{
  for (const auto & pr : as_range(_boundary_node_id.equal_range(node)))
    if (pr.second == id)
      return true;

  return false;
}



#ifdef LIBMESH_ENABLE_DEPRECATED
template <typename RealType>
std::vector<boundary_id_type> BoundaryInfoTempl<RealType>::boundary_ids(const Node * node) const
{
  libmesh_deprecated();

  std::vector<boundary_id_type> ids;
  this->boundary_ids(node, ids);
  return ids;
}
#endif



template <typename RealType>
void BoundaryInfoTempl<RealType>::boundary_ids (const Node * node,
                                 std::vector<boundary_id_type> & vec_to_fill) const
{
  // Clear out any previous contents
  vec_to_fill.clear();

  for (const auto & pr : as_range(_boundary_node_id.equal_range(node)))
    vec_to_fill.push_back(pr.second);
}



template <typename RealType>
unsigned int BoundaryInfoTempl<RealType>::n_boundary_ids(const Node * node) const
{
  auto pos = _boundary_node_id.equal_range(node);
  return cast_int<unsigned int>(std::distance(pos.first, pos.second));
}



#ifdef LIBMESH_ENABLE_DEPRECATED
template <typename RealType>
std::vector<boundary_id_type> BoundaryInfoTempl<RealType>::edge_boundary_ids (const Elem * const elem,
                                                               const unsigned short int edge) const
{
  libmesh_deprecated();

  std::vector<boundary_id_type> ids;
  this->edge_boundary_ids(elem, edge, ids);
  return ids;
}
#endif



template <typename RealType>
void BoundaryInfoTempl<RealType>::edge_boundary_ids (const Elem * const elem,
                                      const unsigned short int edge,
                                      std::vector<boundary_id_type> & vec_to_fill) const
{
  libmesh_assert(elem);

  // Clear out any previous contents
  vec_to_fill.clear();

  // Only level-0 elements store BCs.  If this is not a level-0
  // element get its level-0 parent and infer the BCs.
  const Elem * searched_elem = elem;
#ifdef LIBMESH_ENABLE_AMR
  if (elem->level() != 0)
    {
      // Find all the sides that contain edge. If one of those is a boundary
      // side, then this must be a boundary edge. In that case, we just use the
      // top-level parent.
      bool found_boundary_edge = false;
      for (auto side : elem->side_index_range())
        {
          if (elem->is_edge_on_side(edge,side))
            {
              if (elem->neighbor_ptr(side) == nullptr)
                {
                  searched_elem = elem->top_parent ();
                  found_boundary_edge = true;
                  break;
                }
            }
        }

      if (!found_boundary_edge)
        {
          // Child element is not on external edge, but it may have internal
          // "boundary" IDs.  We will walk up the tree, at each level checking that
          // the current child is actually on the same edge of the parent that is
          // currently being searched for (i.e. that was passed in as "edge").
          while (searched_elem->parent() != nullptr)
            {
              const Elem * parent = searched_elem->parent();
              if (parent->is_child_on_edge(parent->which_child_am_i(searched_elem), edge) == false)
                return;
              searched_elem = parent;
            }
        }
    }
#endif

  // Check each element in the range to see if its edge matches the requested edge.
  for (const auto & pr : as_range(_boundary_edge_id.equal_range(searched_elem)))
    if (pr.second.first == edge)
      vec_to_fill.push_back(pr.second.second);
}



template <typename RealType>
unsigned int BoundaryInfoTempl<RealType>::n_edge_boundary_ids (const Elem * const elem,
                                                const unsigned short int edge) const
{
  std::vector<boundary_id_type> ids;
  this->edge_boundary_ids(elem, edge, ids);
  return cast_int<unsigned int>(ids.size());
}



#ifdef LIBMESH_ENABLE_DEPRECATED
template <typename RealType>
std::vector<boundary_id_type> BoundaryInfoTempl<RealType>::raw_edge_boundary_ids (const Elem * const elem,
                                                                   const unsigned short int edge) const
{
  libmesh_deprecated();

  std::vector<boundary_id_type> ids;
  this->raw_edge_boundary_ids(elem, edge, ids);
  return ids;
}
#endif



template <typename RealType>
void BoundaryInfoTempl<RealType>::raw_edge_boundary_ids (const Elem * const elem,
                                          const unsigned short int edge,
                                          std::vector<boundary_id_type> & vec_to_fill) const
{
  libmesh_assert(elem);

  // Clear out any previous contents
  vec_to_fill.clear();

  // Only level-0 elements store BCs.
  if (elem->parent())
    return;

  // Check each element in the range to see if its edge matches the requested edge.
  for (const auto & pr : as_range(_boundary_edge_id.equal_range(elem)))
    if (pr.second.first == edge)
      vec_to_fill.push_back(pr.second.second);
}



template <typename RealType>
void BoundaryInfoTempl<RealType>::shellface_boundary_ids (const Elem * const elem,
                                           const unsigned short int shellface,
                                           std::vector<boundary_id_type> & vec_to_fill) const
{
  libmesh_assert(elem);

  // Shells only have 2 faces
  libmesh_assert_less(shellface, 2);

  // Clear out any previous contents
  vec_to_fill.clear();

  // Only level-0 elements store BCs.  If this is not a level-0
  // element get its level-0 parent and infer the BCs.
  const Elem * searched_elem = elem;
#ifdef LIBMESH_ENABLE_AMR
  if (elem->level() != 0)
    {
      while (searched_elem->parent() != nullptr)
        {
          const Elem * parent = searched_elem->parent();
          searched_elem = parent;
        }
    }
#endif

  // Check each element in the range to see if its shellface matches the requested shellface.
  for (const auto & pr : as_range(_boundary_shellface_id.equal_range(searched_elem)))
    if (pr.second.first == shellface)
      vec_to_fill.push_back(pr.second.second);
}



template <typename RealType>
unsigned int BoundaryInfoTempl<RealType>::n_shellface_boundary_ids (const Elem * const elem,
                                                     const unsigned short int shellface) const
{
  std::vector<boundary_id_type> ids;
  this->shellface_boundary_ids(elem, shellface, ids);
  return cast_int<unsigned int>(ids.size());
}



template <typename RealType>
void BoundaryInfoTempl<RealType>::raw_shellface_boundary_ids (const Elem * const elem,
                                               const unsigned short int shellface,
                                               std::vector<boundary_id_type> & vec_to_fill) const
{
  libmesh_assert(elem);

  // Shells only have 2 faces
  libmesh_assert_less(shellface, 2);

  // Clear out any previous contents
  vec_to_fill.clear();

  // Only level-0 elements store BCs.
  if (elem->parent())
    return;

  // Check each element in the range to see if its shellface matches the requested shellface.
  for (const auto & pr : as_range(_boundary_shellface_id.equal_range(elem)))
    if (pr.second.first == shellface)
      vec_to_fill.push_back(pr.second.second);
}


#ifdef LIBMESH_ENABLE_DEPRECATED
template <typename RealType>
boundary_id_type BoundaryInfoTempl<RealType>::boundary_id(const Elem * const elem,
                                           const unsigned short int side) const
{
  libmesh_deprecated();

  std::vector<boundary_id_type> ids;
  this->boundary_ids(elem, side, ids);

  // If the set is empty, return invalid_id
  if (ids.empty())
    return invalid_id;

  // Otherwise, just return the first id we came across for this
  // element on this side.
  return *(ids.begin());
}
#endif



template <typename RealType>
bool BoundaryInfoTempl<RealType>::has_boundary_id(const Elem * const elem,
                                   const unsigned short int side,
                                   const boundary_id_type id) const
{
  std::vector<boundary_id_type> ids;
  this->boundary_ids(elem, side, ids);
  return (std::find(ids.begin(), ids.end(), id) != ids.end());
}



#ifdef LIBMESH_ENABLE_DEPRECATED
template <typename RealType>
std::vector<boundary_id_type> BoundaryInfoTempl<RealType>::boundary_ids (const Elem * const elem,
                                                          const unsigned short int side) const
{
  libmesh_deprecated();

  std::vector<boundary_id_type> ids;
  this->boundary_ids(elem, side, ids);
  return ids;
}
#endif



template <typename RealType>
void BoundaryInfoTempl<RealType>::boundary_ids (const Elem * const elem,
                                 const unsigned short int side,
                                 std::vector<boundary_id_type> & vec_to_fill) const
{
  libmesh_assert(elem);

  // Clear out any previous contents
  vec_to_fill.clear();

  // Only level-0 elements store BCs.  If this is not a level-0
  // element get its level-0 parent and infer the BCs.
  const Elem * searched_elem = elem;
  if (elem->level() != 0)
    {
      if (elem->neighbor_ptr(side) == nullptr)
        searched_elem = elem->top_parent ();
#ifdef LIBMESH_ENABLE_AMR
      else
        while (searched_elem->parent() != nullptr)
          {
            const Elem * parent = searched_elem->parent();
            if (parent->is_child_on_side(parent->which_child_am_i(searched_elem), side) == false)
              return;
            searched_elem = parent;
          }
#endif
    }

  // Check each element in the range to see if its side matches the requested side.
  for (const auto & pr : as_range(_boundary_side_id.equal_range(searched_elem)))
    if (pr.second.first == side)
      vec_to_fill.push_back(pr.second.second);
}




template <typename RealType>
unsigned int BoundaryInfoTempl<RealType>::n_boundary_ids (const Elem * const elem,
                                           const unsigned short int side) const
{
  std::vector<boundary_id_type> ids;
  this->boundary_ids(elem, side, ids);
  return cast_int<unsigned int>(ids.size());
}



#ifdef LIBMESH_ENABLE_DEPRECATED
template <typename RealType>
std::vector<boundary_id_type> BoundaryInfoTempl<RealType>::raw_boundary_ids (const Elem * const elem,
                                                              const unsigned short int side) const
{
  libmesh_deprecated();

  std::vector<boundary_id_type> ids;
  this->raw_boundary_ids(elem, side, ids);
  return ids;
}
#endif



template <typename RealType>
void BoundaryInfoTempl<RealType>::raw_boundary_ids (const Elem * const elem,
                                     const unsigned short int side,
                                     std::vector<boundary_id_type> & vec_to_fill) const
{
  libmesh_assert(elem);

  // Clear out any previous contents
  vec_to_fill.clear();

  // Only level-0 elements store BCs.
  if (elem->parent())
    return;

  // Check each element in the range to see if its side matches the requested side.
  for (const auto & pr : as_range(_boundary_side_id.equal_range(elem)))
    if (pr.second.first == side)
      vec_to_fill.push_back(pr.second.second);
}



template <typename RealType>
void BoundaryInfoTempl<RealType>::copy_boundary_ids (const BoundaryInfo & old_boundary_info,
                                      const Elem * const old_elem,
                                      const Elem * const new_elem)
{
  libmesh_assert_equal_to (old_elem->n_sides(), new_elem->n_sides());
  libmesh_assert_equal_to (old_elem->n_edges(), new_elem->n_edges());

  std::vector<boundary_id_type> bndry_ids;

  for (auto s : old_elem->side_index_range())
    {
      old_boundary_info.raw_boundary_ids (old_elem, s, bndry_ids);
      this->add_side (new_elem, s, bndry_ids);
    }

  for (auto e : old_elem->edge_index_range())
    {
      old_boundary_info.raw_edge_boundary_ids (old_elem, e, bndry_ids);
      this->add_edge (new_elem, e, bndry_ids);
    }

  for (unsigned short sf=0; sf != 2; sf++)
    {
      old_boundary_info.raw_shellface_boundary_ids (old_elem, sf, bndry_ids);
      this->add_shellface (new_elem, sf, bndry_ids);
    }
}



template <typename RealType>
void BoundaryInfoTempl<RealType>::remove (const Node * node)
{
  libmesh_assert(node);

  // Erase everything associated with node
  _boundary_node_id.erase (node);
}



template <typename RealType>
void BoundaryInfoTempl<RealType>::remove_node (const Node * node,
                                const boundary_id_type id)
{
  libmesh_assert(node);

  // Erase (node, id) entry from map.
  erase_if(_boundary_node_id, node,
           [id](typename decltype(_boundary_node_id)::mapped_type & val)
           {return val == id;});
}



template <typename RealType>
void BoundaryInfoTempl<RealType>::remove (const Elem * elem)
{
  libmesh_assert(elem);

  // Erase everything associated with elem
  _boundary_edge_id.erase (elem);
  _boundary_side_id.erase (elem);
  _boundary_shellface_id.erase (elem);
}



template <typename RealType>
void BoundaryInfoTempl<RealType>::remove_edge (const Elem * elem,
                                const unsigned short int edge)
{
  libmesh_assert(elem);

  // Only level 0 elements are stored in BoundaryInfo.
  libmesh_assert_equal_to (elem->level(), 0);

  // Erase (elem, edge, *) entries from map.
  erase_if(_boundary_edge_id, elem,
           [edge](typename decltype(_boundary_edge_id)::mapped_type & pr)
           {return pr.first == edge;});
}



template <typename RealType>
void BoundaryInfoTempl<RealType>::remove_edge (const Elem * elem,
                                const unsigned short int edge,
                                const boundary_id_type id)
{
  libmesh_assert(elem);

  // Only level 0 elements are stored in BoundaryInfo.
  libmesh_assert_equal_to (elem->level(), 0);

  // Erase (elem, edge, id) entries from map.
  erase_if(_boundary_edge_id, elem,
           [edge, id](typename decltype(_boundary_edge_id)::mapped_type & pr)
           {return pr.first == edge && pr.second == id;});
}


template <typename RealType>
void BoundaryInfoTempl<RealType>::remove_shellface (const Elem * elem,
                                     const unsigned short int shellface)
{
  libmesh_assert(elem);

  // Only level 0 elements are stored in BoundaryInfo.
  libmesh_assert_equal_to (elem->level(), 0);

  // Shells only have 2 faces
  libmesh_assert_less(shellface, 2);

  // Erase (elem, shellface, *) entries from map.
  erase_if(_boundary_shellface_id, elem,
           [shellface](typename decltype(_boundary_shellface_id)::mapped_type & pr)
           {return pr.first == shellface;});
}



template <typename RealType>
void BoundaryInfoTempl<RealType>::remove_shellface (const Elem * elem,
                                     const unsigned short int shellface,
                                     const boundary_id_type id)
{
  libmesh_assert(elem);

  // Only level 0 elements are stored in BoundaryInfo.
  libmesh_assert_equal_to (elem->level(), 0);

  // Shells only have 2 faces
  libmesh_assert_less(shellface, 2);

  // Erase (elem, shellface, id) entries from map.
  erase_if(_boundary_shellface_id, elem,
           [shellface, id](typename decltype(_boundary_shellface_id)::mapped_type & pr)
           {return pr.first == shellface && pr.second == id;});
}

template <typename RealType>
void BoundaryInfoTempl<RealType>::remove_side (const Elem * elem,
                                const unsigned short int side)
{
  libmesh_assert(elem);

  // Only level 0 elements are stored in BoundaryInfo.
  libmesh_assert_equal_to (elem->level(), 0);

  // Erase (elem, side, *) entries from map.
  erase_if(_boundary_side_id, elem,
           [side](typename decltype(_boundary_side_id)::mapped_type & pr)
           {return pr.first == side;});
}



template <typename RealType>
void BoundaryInfoTempl<RealType>::remove_side (const Elem * elem,
                                const unsigned short int side,
                                const boundary_id_type id)
{
  libmesh_assert(elem);

  // Erase (elem, side, id) entries from map.
  erase_if(_boundary_side_id, elem,
           [side, id](typename decltype(_boundary_side_id)::mapped_type & pr)
           {return pr.first == side && pr.second == id;});
}



template <typename RealType>
void BoundaryInfoTempl<RealType>::remove_id (boundary_id_type id)
{
  // Erase id from ids containers
  _boundary_ids.erase(id);
  _side_boundary_ids.erase(id);
  _edge_boundary_ids.erase(id);
  _shellface_boundary_ids.erase(id);
  _node_boundary_ids.erase(id);
  _ss_id_to_name.erase(id);
  _ns_id_to_name.erase(id);
  _es_id_to_name.erase(id);

  // Erase (*, id) entries from map.
  erase_if(_boundary_node_id,
           [id](typename decltype(_boundary_node_id)::mapped_type & val)
           {return val == id;});

  // Erase (*, *, id) entries from map.
  erase_if(_boundary_edge_id,
           [id](typename decltype(_boundary_edge_id)::mapped_type & pr)
           {return pr.second == id;});

  // Erase (*, *, id) entries from map.
  erase_if(_boundary_shellface_id,
           [id](typename decltype(_boundary_shellface_id)::mapped_type & pr)
           {return pr.second == id;});

  // Erase (*, *, id) entries from map.
  erase_if(_boundary_side_id,
           [id](typename decltype(_boundary_side_id)::mapped_type & pr)
           {return pr.second == id;});
}



template <typename RealType>
unsigned int BoundaryInfoTempl<RealType>::side_with_boundary_id(const Elem * const elem,
                                                 const boundary_id_type boundary_id_in) const
{
  const Elem * searched_elem = elem;
  if (elem->level() != 0)
    searched_elem = elem->top_parent();

  // elem may have zero or multiple occurrences
  for (const auto & pr : as_range(_boundary_side_id.equal_range(searched_elem)))
    {
      // if this is true we found the requested boundary_id
      // of the element and want to return the side
      if (pr.second.second == boundary_id_in)
        {
          unsigned int side = pr.second.first;

          // If we're on this external boundary then we share this
          // external boundary id
          if (elem->neighbor_ptr(side) == nullptr)
            return side;

          // If we're on an internal boundary then we need to be sure
          // it's the same internal boundary as our top_parent
          const Elem * p = elem;

#ifdef LIBMESH_ENABLE_AMR

          while (p != nullptr)
            {
              const Elem * parent = p->parent();
              if (!parent->is_child_on_side(parent->which_child_am_i(p), side))
                break;
              p = parent;
            }
#endif
          // We're on that side of our top_parent; return it
          if (!p)
            return side;
        }
    }

  // if we get here, we found elem in the data structure but not
  // the requested boundary id, so return the default value
  return libMesh::invalid_uint;
}

template <typename RealType>
void
BoundaryInfoTempl<RealType>::build_node_boundary_ids(std::vector<boundary_id_type> & b_ids) const
{
  b_ids.clear();

  for (const auto & pr : _boundary_node_id)
    {
      boundary_id_type id = pr.second;

      if (std::find(b_ids.begin(),b_ids.end(),id) == b_ids.end())
        b_ids.push_back(id);
    }
}

template <typename RealType>
void
BoundaryInfoTempl<RealType>::build_side_boundary_ids(std::vector<boundary_id_type> & b_ids) const
{
  b_ids.clear();

  for (const auto & pr : _boundary_side_id)
    {
      boundary_id_type id = pr.second.second;

      if (std::find(b_ids.begin(),b_ids.end(),id) == b_ids.end())
        b_ids.push_back(id);
    }
}

template <typename RealType>
void
BoundaryInfoTempl<RealType>::build_shellface_boundary_ids(std::vector<boundary_id_type> & b_ids) const
{
  b_ids.clear();

  for (const auto & pr :_boundary_shellface_id)
    {
      boundary_id_type id = pr.second.second;

      if (std::find(b_ids.begin(),b_ids.end(),id) == b_ids.end())
        b_ids.push_back(id);
    }
}

template <typename RealType>
std::size_t BoundaryInfoTempl<RealType>::n_boundary_conds () const
{
  // in serial we know the number of bcs from the
  // size of the container
  if (_mesh.is_serial())
    return _boundary_side_id.size();

  // in parallel we need to sum the number of local bcs
  parallel_object_only();

  std::size_t nbcs=0;

  for (const auto & pr : _boundary_side_id)
    if (pr.first->processor_id() == this->processor_id())
      nbcs++;

  this->comm().sum (nbcs);

  return nbcs;
}

template <typename RealType>
std::size_t BoundaryInfoTempl<RealType>::n_edge_conds () const
{
  // in serial we know the number of nodesets from the
  // size of the container
  if (_mesh.is_serial())
    return _boundary_edge_id.size();

  // in parallel we need to sum the number of local nodesets
  parallel_object_only();

  std::size_t n_edge_bcs=0;

  for (const auto & pr : _boundary_edge_id)
    if (pr.first->processor_id() == this->processor_id())
      n_edge_bcs++;

  this->comm().sum (n_edge_bcs);

  return n_edge_bcs;
}


template <typename RealType>
std::size_t BoundaryInfoTempl<RealType>::n_shellface_conds () const
{
  // in serial we know the number of nodesets from the
  // size of the container
  if (_mesh.is_serial())
    return _boundary_shellface_id.size();

  // in parallel we need to sum the number of local nodesets
  parallel_object_only();

  std::size_t n_shellface_bcs=0;

  for (const auto & pr : _boundary_shellface_id)
    if (pr.first->processor_id() == this->processor_id())
      n_shellface_bcs++;

  this->comm().sum (n_shellface_bcs);

  return n_shellface_bcs;
}


template <typename RealType>
std::size_t BoundaryInfoTempl<RealType>::n_nodeset_conds () const
{
  // in serial we know the number of nodesets from the
  // size of the container
  if (_mesh.is_serial())
    return _boundary_node_id.size();

  // in parallel we need to sum the number of local nodesets
  parallel_object_only();

  std::size_t n_nodesets=0;

  for (const auto & pr : _boundary_node_id)
    if (pr.first->processor_id() == this->processor_id())
      n_nodesets++;

  this->comm().sum (n_nodesets);

  return n_nodesets;
}



#ifdef LIBMESH_ENABLE_DEPRECATED
template <typename RealType>
void BoundaryInfoTempl<RealType>::build_node_list (std::vector<dof_id_type> & nl,
                                    std::vector<boundary_id_type> & il) const
{
  libmesh_deprecated();

  // Call the non-deprecated version of this function.
  auto bc_tuples = this->build_node_list();

  // Clear the input vectors, just in case they were used for
  // something else recently...
  nl.clear();
  il.clear();

  // Reserve the size, then use push_back
  nl.reserve (bc_tuples.size());
  il.reserve (bc_tuples.size());

  for (const auto & t : bc_tuples)
    {
      nl.push_back(std::get<0>(t));
      il.push_back(std::get<1>(t));
    }
}
#endif


template <typename RealType>
std::vector<std::tuple<dof_id_type, boundary_id_type>>
BoundaryInfoTempl<RealType>::build_node_list() const
{
  std::vector<std::tuple<dof_id_type, boundary_id_type>> bc_tuples;
  bc_tuples.reserve(_boundary_node_id.size());

  for (const auto & pr : _boundary_node_id)
    bc_tuples.emplace_back(pr.first->id(), pr.second);

  // This list is currently in memory address (arbitrary) order, so
  // sort to make it consistent on all procs.
  std::sort(bc_tuples.begin(), bc_tuples.end());

  return bc_tuples;
}


template <typename RealType>
void
BoundaryInfoTempl<RealType>::build_node_list_from_side_list()
{
  // If we're on a distributed mesh, even the owner of a node is not
  // guaranteed to be able to properly assign its new boundary id(s)!
  // Nodal neighbors are not always ghosted, and a nodal neighbor
  // might have a boundary side.
  const bool mesh_is_serial = _mesh.is_serial();

  typedef std::set<std::pair<dof_id_type, boundary_id_type>> set_type;

  const processor_id_type n_proc     = this->n_processors();
  const processor_id_type my_proc_id = this->processor_id();
  std::vector<set_type> nodes_to_push(n_proc);

  // Pull objects out of the loop to reduce heap operations
  std::unique_ptr<const Elem> side;

  // Loop over the side list
  for (const auto & pr : _boundary_side_id)
    {
      // Don't add remote sides
      if (pr.first->is_remote())
        continue;

      // Need to loop over the sides of any possible children
      std::vector<const Elem *> family;
#ifdef LIBMESH_ENABLE_AMR
      pr.first->active_family_tree_by_side (family, pr.second.first);
#else
      family.push_back(pr.first);
#endif

      for (const auto & cur_elem : family)
        {
          cur_elem->build_side_ptr(side, pr.second.first);

          // Add each node node on the side with the side's boundary id
          for (auto i : side->node_index_range())
            {
              const boundary_id_type bcid = pr.second.second;
              this->add_node(side->node_ptr(i), bcid);
              if (!mesh_is_serial)
                {
                  const processor_id_type proc_id =
                    side->node_ptr(i)->processor_id();
                  if (proc_id != my_proc_id)
                    nodes_to_push[proc_id].insert
                      (std::make_pair(side->node_id(i), bcid));
                }
            }
        }
    }

  // If we're on a serial mesh then we're done.
  if (mesh_is_serial)
    return;

  // Otherwise we need to push ghost node bcids to their owners, then
  // pull ghost node bcids from their owners.
  Parallel::MessageTag
    node_pushes_tag = this->comm().get_unique_tag(),
    node_pulls_tag = this->comm().get_unique_tag(),
    node_responses_tag = this->comm().get_unique_tag();

  std::vector<Parallel::Request> node_push_requests(n_proc-1);

  for (processor_id_type p = 0; p != n_proc; ++p)
    {
      if (p == my_proc_id)
        continue;

      Parallel::Request &request =
        node_push_requests[p - (p > my_proc_id)];

      this->comm().send
        (p, nodes_to_push[p], request, node_pushes_tag);
    }

  for (processor_id_type p = 1; p != n_proc; ++p)
    {
      set_type received_nodes;

      this->comm().receive
        (Parallel::any_source, received_nodes, node_pushes_tag);

      for (const auto & pr : received_nodes)
        this->add_node(_mesh.node_ptr(pr.first), pr.second);
    }

  // At this point we should know all the BCs for our own nodes; now
  // we need BCs for ghost nodes.
  //
  // FIXME - parallel_ghost_sync.h doesn't work here because it
  // assumes a fixed size datum on each node.
  std::vector<std::vector<dof_id_type>> node_ids_requested(n_proc);

  // Determine what nodes we need to request
  for (const auto & node : _mesh.node_ptr_range())
    {
      const processor_id_type pid = node->processor_id();
      if (pid != my_proc_id)
        node_ids_requested[pid].push_back(node->id());
    }

  typedef std::vector<std::pair<dof_id_type, boundary_id_type>> vec_type;

  std::vector<Parallel::Request>
    node_pull_requests(n_proc-1),
    node_response_requests(n_proc-1);

  // Make all requests
  for (processor_id_type p = 0; p != n_proc; ++p)
    {
      if (p == my_proc_id)
        continue;

      Parallel::Request &request =
        node_pull_requests[p - (p > my_proc_id)];

      this->comm().send
        (p, node_ids_requested[p], request, node_pulls_tag);
    }

  // Process all incoming requests
  std::vector<vec_type> responses(n_proc-1);

  for (processor_id_type p = 1; p != n_proc; ++p)
    {
      std::vector<dof_id_type> requested_nodes;

      Parallel::Status
        status(this->comm().probe (Parallel::any_source, node_pulls_tag));
      const processor_id_type
        source_pid = cast_int<processor_id_type>(status.source());

      this->comm().receive
        (source_pid, requested_nodes, node_pulls_tag);

      Parallel::Request &request =
        node_response_requests[p-1];

      std::vector<boundary_id_type> bcids;

      for (const auto & id : requested_nodes)
        {
          this->boundary_ids(_mesh.node_ptr(id), bcids);

          for (const auto & b : bcids)
            responses[p-1].push_back(std::make_pair(id, b));
        }

      this->comm().send
        (source_pid, responses[p-1], request, node_responses_tag);
    }

  // Process all incoming responses
  for (processor_id_type p = 1; p != n_proc; ++p)
    {
      Parallel::Status
        status(this->comm().probe (Parallel::any_source, node_responses_tag));
      const processor_id_type
        source_pid = cast_int<processor_id_type>(status.source());

      vec_type response;

      this->comm().receive
        (source_pid, response, node_responses_tag);

      for (const auto & pr : response)
        this->add_node(_mesh.node_ptr(pr.first), pr.second);
    }

  Parallel::wait (node_push_requests);
  Parallel::wait (node_pull_requests);
  Parallel::wait (node_response_requests);
}




template <typename RealType>
void BoundaryInfoTempl<RealType>::build_side_list_from_node_list()
{
  // Check for early return
  if (_boundary_node_id.empty())
    {
      libMesh::out << "No boundary node IDs have been added: cannot build side list!" << std::endl;
      return;
    }

  // Pull objects out of the loop to reduce heap operations
  std::unique_ptr<const Elem> side_elem;

  for (const auto & elem : _mesh.active_element_ptr_range())
    for (auto side : elem->side_index_range())
      {
        elem->build_side_ptr(side_elem, side);

        // map from nodeset_id to count for that ID
        std::map<boundary_id_type, unsigned> nodesets_node_count;

        // For each nodeset that this node is a member of, increment the associated
        // nodeset ID count
        for (const auto & node : side_elem->node_ref_range())
          for (const auto & pr : as_range(_boundary_node_id.equal_range(&node)))
            nodesets_node_count[pr.second]++;

        // Now check to see what nodeset_counts have the correct
        // number of nodes in them.  For any that do, add this side to
        // the sideset, making sure the sideset inherits the
        // nodeset's name, if there is one.
        for (const auto & pr : nodesets_node_count)
          if (pr.second == side_elem->n_nodes())
            {
              add_side(elem, side, pr.first);

              // Let the sideset inherit any non-empty name from the nodeset
              std::string & nset_name = nodeset_name(pr.first);

              if (nset_name != "")
                sideset_name(pr.first) = nset_name;
            }
      } // end for side
}




#ifdef LIBMESH_ENABLE_DEPRECATED
template <typename RealType>
void BoundaryInfoTempl<RealType>::build_side_list (std::vector<dof_id_type> & el,
                                    std::vector<unsigned short int> & sl,
                                    std::vector<boundary_id_type> & il) const
{
  libmesh_deprecated();

  // Call the non-deprecated version of this function.
  auto bc_tuples = this->build_side_list();

  // Clear the input vectors, just in case they were used for
  // something else recently...
  el.clear();
  sl.clear();
  il.clear();

  // Reserve the size, then use push_back
  el.reserve (bc_tuples.size());
  sl.reserve (bc_tuples.size());
  il.reserve (bc_tuples.size());

  for (const auto & t : bc_tuples)
    {
      el.push_back(std::get<0>(t));
      sl.push_back(std::get<1>(t));
      il.push_back(std::get<2>(t));
    }
}
#endif


template <typename RealType>
std::vector<std::tuple<dof_id_type, unsigned short int, boundary_id_type>>
BoundaryInfoTempl<RealType>::build_side_list() const
{
  std::vector<std::tuple<dof_id_type, unsigned short int, boundary_id_type>> bc_triples;
  bc_triples.reserve(_boundary_side_id.size());

  for (const auto & pr : _boundary_side_id)
    bc_triples.emplace_back(pr.first->id(), pr.second.first, pr.second.second);

  // bc_triples is currently in whatever order the Elem pointers in
  // the _boundary_side_id multimap are in, and in particular might be
  // in different orders on different processors. To avoid this
  // inconsistency, we'll sort using the default operator< for tuples.
  std::sort(bc_triples.begin(), bc_triples.end());

  return bc_triples;
}



#ifdef LIBMESH_ENABLE_DEPRECATED
template <typename RealType>
void BoundaryInfoTempl<RealType>::build_active_side_list (std::vector<dof_id_type> & el,
                                           std::vector<unsigned short int> & sl,
                                           std::vector<boundary_id_type> & il) const
{
  libmesh_deprecated();

  // Call the non-deprecated version of this function.
  auto bc_tuples = this->build_active_side_list();

  // Clear the input vectors, just in case they were used for
  // something else recently...
  el.clear();
  sl.clear();
  il.clear();

  // Reserve the size, then use push_back
  el.reserve (bc_tuples.size());
  sl.reserve (bc_tuples.size());
  il.reserve (bc_tuples.size());

  for (const auto & t : bc_tuples)
    {
      el.push_back(std::get<0>(t));
      sl.push_back(std::get<1>(t));
      il.push_back(std::get<2>(t));
    }
}
#endif


template <typename RealType>
std::vector<std::tuple<dof_id_type, unsigned short int, boundary_id_type>>
BoundaryInfoTempl<RealType>::build_active_side_list () const
{
  std::vector<std::tuple<dof_id_type, unsigned short int, boundary_id_type>> bc_triples;
  bc_triples.reserve(_boundary_side_id.size());

  for (const auto & pr : _boundary_side_id)
    {
      // Don't add remote sides
      if (pr.first->is_remote())
        continue;

      // Loop over the sides of possible children
      std::vector<const Elem *> family;
#ifdef LIBMESH_ENABLE_AMR
      pr.first->active_family_tree_by_side(family, pr.second.first);
#else
      family.push_back(pr.first);
#endif

      // Populate the list items
      for (const auto & elem : family)
        bc_triples.emplace_back(elem->id(), pr.second.first, pr.second.second);
    }

  // This list is currently in memory address (arbitrary) order, so
  // sort to make it consistent on all procs.
  std::sort(bc_triples.begin(), bc_triples.end());

  return bc_triples;
}


#ifdef LIBMESH_ENABLE_DEPRECATED
template <typename RealType>
void BoundaryInfoTempl<RealType>::build_edge_list (std::vector<dof_id_type> & el,
                                    std::vector<unsigned short int> & sl,
                                    std::vector<boundary_id_type> & il) const
{
  libmesh_deprecated();

  // Call the non-deprecated version of this function.
  auto bc_tuples = this->build_edge_list();

  // Clear the input vectors, just in case they were used for
  // something else recently...
  el.clear();
  sl.clear();
  il.clear();

  // Reserve the size, then use push_back
  el.reserve (bc_tuples.size());
  sl.reserve (bc_tuples.size());
  il.reserve (bc_tuples.size());

  for (const auto & t : bc_tuples)
    {
      el.push_back(std::get<0>(t));
      sl.push_back(std::get<1>(t));
      il.push_back(std::get<2>(t));
    }
}
#endif


template <typename RealType>
std::vector<std::tuple<dof_id_type, unsigned short int, boundary_id_type>>
BoundaryInfoTempl<RealType>::build_edge_list() const
{
  std::vector<std::tuple<dof_id_type, unsigned short int, boundary_id_type>> bc_triples;
  bc_triples.reserve(_boundary_edge_id.size());

  for (const auto & pr : _boundary_edge_id)
    bc_triples.emplace_back(pr.first->id(), pr.second.first, pr.second.second);

  // This list is currently in memory address (arbitrary) order, so
  // sort to make it consistent on all procs.
  std::sort(bc_triples.begin(), bc_triples.end());

  return bc_triples;
}


#ifdef LIBMESH_ENABLE_DEPRECATED
template <typename RealType>
void BoundaryInfoTempl<RealType>::build_shellface_list (std::vector<dof_id_type> & el,
                                         std::vector<unsigned short int> & sl,
                                         std::vector<boundary_id_type> & il) const
{
  libmesh_deprecated();

  // Call the non-deprecated version of this function.
  auto bc_tuples = this->build_shellface_list();

  // Clear the input vectors, just in case they were used for
  // something else recently...
  el.clear();
  sl.clear();
  il.clear();

  // Reserve the size, then use push_back
  el.reserve (bc_tuples.size());
  sl.reserve (bc_tuples.size());
  il.reserve (bc_tuples.size());

  for (const auto & t : bc_tuples)
    {
      el.push_back(std::get<0>(t));
      sl.push_back(std::get<1>(t));
      il.push_back(std::get<2>(t));
    }
}
#endif


template <typename RealType>
std::vector<std::tuple<dof_id_type, unsigned short int, boundary_id_type>>
BoundaryInfoTempl<RealType>::build_shellface_list() const
{
  std::vector<std::tuple<dof_id_type, unsigned short int, boundary_id_type>> bc_triples;
  bc_triples.reserve(_boundary_shellface_id.size());

  for (const auto & pr : _boundary_shellface_id)
    bc_triples.emplace_back(pr.first->id(), pr.second.first, pr.second.second);

  // This list is currently in memory address (arbitrary) order, so
  // sort to make it consistent on all procs.
  std::sort(bc_triples.begin(), bc_triples.end());

  return bc_triples;
}


template <typename RealType>
void BoundaryInfoTempl<RealType>::print_info(std::ostream & out_stream) const
{
  // Print out the nodal BCs
  if (!_boundary_node_id.empty())
    {
      out_stream << "Nodal Boundary conditions:" << std::endl
                 << "--------------------------" << std::endl
                 << "  (Node No., ID)               " << std::endl;

      for (const auto & pr : _boundary_node_id)
        out_stream << "  (" << pr.first->id()
                   << ", "  << pr.second
                   << ")"  << std::endl;
    }

  // Print out the element edge BCs
  if (!_boundary_edge_id.empty())
    {
      out_stream << std::endl
                 << "Edge Boundary conditions:" << std::endl
                 << "-------------------------" << std::endl
                 << "  (Elem No., Edge No., ID)      " << std::endl;

      for (const auto & pr : _boundary_edge_id)
        out_stream << "  (" << pr.first->id()
                   << ", "  << pr.second.first
                   << ", "  << pr.second.second
                   << ")"   << std::endl;
    }

  // Print out the element shellface BCs
  if (!_boundary_shellface_id.empty())
    {
      out_stream << std::endl
                 << "Shell-face Boundary conditions:" << std::endl
                 << "-------------------------" << std::endl
                 << "  (Elem No., Shell-face No., ID)      " << std::endl;

      for (const auto & pr : _boundary_shellface_id)
        out_stream << "  (" << pr.first->id()
                   << ", "  << pr.second.first
                   << ", "  << pr.second.second
                   << ")"   << std::endl;
    }

  // Print out the element side BCs
  if (!_boundary_side_id.empty())
    {
      out_stream << std::endl
                 << "Side Boundary conditions:" << std::endl
                 << "-------------------------" << std::endl
                 << "  (Elem No., Side No., ID)      " << std::endl;

      for (const auto & pr : _boundary_side_id)
        out_stream << "  (" << pr.first->id()
                   << ", "  << pr.second.first
                   << ", "  << pr.second.second
                   << ")"   << std::endl;
    }
}



template <typename RealType>
void BoundaryInfoTempl<RealType>::print_summary(std::ostream & out_stream) const
{
  // Print out the nodal BCs
  if (!_boundary_node_id.empty())
    {
      out_stream << "Nodal Boundary conditions:" << std::endl
                 << "--------------------------" << std::endl
                 << "  (ID, number of nodes)   " << std::endl;

      std::map<boundary_id_type, std::size_t> ID_counts;

      for (const auto & pr : _boundary_node_id)
        ID_counts[pr.second]++;

      for (const auto & pr : ID_counts)
        out_stream << "  (" << pr.first
                   << ", "  << pr.second
                   << ")"  << std::endl;
    }

  // Print out the element edge BCs
  if (!_boundary_edge_id.empty())
    {
      out_stream << std::endl
                 << "Edge Boundary conditions:" << std::endl
                 << "-------------------------" << std::endl
                 << "  (ID, number of edges)   " << std::endl;

      std::map<boundary_id_type, std::size_t> ID_counts;

      for (const auto & pr : _boundary_edge_id)
        ID_counts[pr.second.second]++;

      for (const auto & pr : ID_counts)
        out_stream << "  (" << pr.first
                   << ", "  << pr.second
                   << ")"  << std::endl;
    }


  // Print out the element edge BCs
  if (!_boundary_shellface_id.empty())
    {
      out_stream << std::endl
                 << "Shell-face Boundary conditions:" << std::endl
                 << "-------------------------" << std::endl
                 << "  (ID, number of shellfaces)   " << std::endl;

      std::map<boundary_id_type, std::size_t> ID_counts;

      for (const auto & pr : _boundary_shellface_id)
        ID_counts[pr.second.second]++;

      for (const auto & pr : ID_counts)
        out_stream << "  (" << pr.first
                   << ", "  << pr.second
                   << ")"  << std::endl;
    }

  // Print out the element side BCs
  if (!_boundary_side_id.empty())
    {
      out_stream << std::endl
                 << "Side Boundary conditions:" << std::endl
                 << "-------------------------" << std::endl
                 << "  (ID, number of sides)   " << std::endl;

      std::map<boundary_id_type, std::size_t> ID_counts;

      for (const auto & pr : _boundary_side_id)
        ID_counts[pr.second.second]++;

      for (const auto & pr : ID_counts)
        out_stream << "  (" << pr.first
                   << ", "  << pr.second
                   << ")"  << std::endl;
    }
}


template <typename RealType>
const std::string & BoundaryInfoTempl<RealType>::get_sideset_name(boundary_id_type id) const
{
  static const std::string empty_string;
  std::map<boundary_id_type, std::string>::const_iterator it =
    _ss_id_to_name.find(id);
  if (it == _ss_id_to_name.end())
    return empty_string;
  else
    return it->second;
}


template <typename RealType>
std::string & BoundaryInfoTempl<RealType>::sideset_name(boundary_id_type id)
{
  return _ss_id_to_name[id];
}

template <typename RealType>
const std::string & BoundaryInfoTempl<RealType>::get_nodeset_name(boundary_id_type id) const
{
  static const std::string empty_string;
  std::map<boundary_id_type, std::string>::const_iterator it =
    _ns_id_to_name.find(id);
  if (it == _ns_id_to_name.end())
    return empty_string;
  else
    return it->second;
}

template <typename RealType>
std::string & BoundaryInfoTempl<RealType>::nodeset_name(boundary_id_type id)
{
  return _ns_id_to_name[id];
}

template <typename RealType>
const std::string & BoundaryInfoTempl<RealType>::get_edgeset_name(boundary_id_type id) const
{
  static const std::string empty_string;
  std::map<boundary_id_type, std::string>::const_iterator it =
    _es_id_to_name.find(id);
  if (it == _es_id_to_name.end())
    return empty_string;
  else
    return it->second;
}


template <typename RealType>
std::string & BoundaryInfoTempl<RealType>::edgeset_name(boundary_id_type id)
{
  return _es_id_to_name[id];
}

template <typename RealType>
boundary_id_type BoundaryInfoTempl<RealType>::get_id_by_name(const std::string & name) const
{
  // Search sidesets
  for (const auto & pr : _ss_id_to_name)
    if (pr.second == name)
      return pr.first;

  // Search nodesets
  for (const auto & pr : _ns_id_to_name)
    if (pr.second == name)
      return pr.first;

  // Search edgesets
  for (const auto & pr : _es_id_to_name)
    if (pr.second == name)
      return pr.first;

  // If we made it here without returning, we don't have a sideset,
  // nodeset, or edgeset by the requested name, so return invalid_id
  return invalid_id;
}



template <typename RealType>
void BoundaryInfoTempl<RealType>::_find_id_maps(const std::set<boundary_id_type> & requested_boundary_ids,
                                 dof_id_type first_free_node_id,
                                 std::map<dof_id_type, dof_id_type> * node_id_map,
                                 dof_id_type first_free_elem_id,
                                 std::map<std::pair<dof_id_type, unsigned char>, dof_id_type> * side_id_map,
                                 const std::set<subdomain_id_type> & subdomains_relative_to)
{
  // We'll do the same modulus trick that DistributedMesh uses to avoid
  // id conflicts between different processors
  dof_id_type
    next_node_id = first_free_node_id + this->processor_id(),
    next_elem_id = first_free_elem_id + this->processor_id();

  // Pull objects out of the loop to reduce heap operations
  std::unique_ptr<const Elem> side;

  // We'll pass through the mesh once first to build
  // the maps and count boundary nodes and elements.
  // To find local boundary nodes, we have to examine all elements
  // here rather than just local elements, because it's possible to
  // have a local boundary node that's not on a local boundary
  // element, e.g. at the tip of a triangle.

  // We'll loop through two different ranges here: first all elements,
  // looking for local nodes, and second through unpartitioned
  // elements, looking for all remaining nodes.
  const typename MeshBase::const_element_iterator end_el = _mesh.elements_end();
  bool hit_end_el = false;
  const typename MeshBase::const_element_iterator end_unpartitioned_el =
    _mesh.pid_elements_end(DofObject::invalid_processor_id);

  for (typename MeshBase::const_element_iterator el = _mesh.elements_begin();
       !hit_end_el || (el != end_unpartitioned_el); ++el)
    {
      if ((el == end_el) && !hit_end_el)
        {
          // Note that we're done with local nodes and just looking
          // for remaining unpartitioned nodes
          hit_end_el = true;

          // Join up the local results from other processors
          if (side_id_map)
            this->comm().set_union(*side_id_map);
          if (node_id_map)
            this->comm().set_union(*node_id_map);

          // Finally we'll pass through any unpartitioned elements to add them
          // to the maps and counts.
          next_node_id = first_free_node_id + this->n_processors();
          next_elem_id = first_free_elem_id + this->n_processors();

          el = _mesh.pid_elements_begin(DofObject::invalid_processor_id);
          if (el == end_unpartitioned_el)
            break;
        }

      const Elem * elem = *el;

      // If the subdomains_relative_to container has the
      // invalid_subdomain_id, we fall back on the "old" behavior of
      // adding sides regardless of this Elem's subdomain. Otherwise,
      // if the subdomains_relative_to container doesn't contain the
      // current Elem's subdomain_id(), we won't add any sides from
      // it.
      if (!subdomains_relative_to.count(Elem::invalid_subdomain_id) &&
          !subdomains_relative_to.count(elem->subdomain_id()))
        continue;

      // Get the top-level parent for this element. This is used for
      // searching for boundary sides on this element.
      const Elem * top_parent = elem->top_parent();

      // Find all the boundary side ids for this Elem.
      auto bounds = _boundary_side_id.equal_range(top_parent);

      for (auto s : elem->side_index_range())
        {
          bool add_this_side = false;
          boundary_id_type this_bcid = invalid_id;

          for (const auto & pr : as_range(bounds))
            {
              this_bcid = pr.second.second;

              // if this side is flagged with a boundary condition
              // and the user wants this id
              if ((pr.second.first == s) &&
                  (requested_boundary_ids.count(this_bcid)))
                {
                  add_this_side = true;
                  break;
                }
            }

          // We may still want to add this side if the user called
          // sync() with no requested_boundary_ids. This corresponds
          // to the "old" style of calling sync() in which the entire
          // boundary was copied to the BoundaryMesh, and handles the
          // case where elements on the geometric boundary are not in
          // any sidesets.
          if (bounds.first == bounds.second            &&
              requested_boundary_ids.count(invalid_id) &&
              elem->neighbor_ptr(s) == nullptr)
            add_this_side = true;

          if (add_this_side)
            {
              // We only assign ids for our own and for
              // unpartitioned elements
              if (side_id_map &&
                  ((elem->processor_id() == this->processor_id()) ||
                   (elem->processor_id() ==
                    DofObject::invalid_processor_id)))
                {
                  std::pair<dof_id_type, unsigned char> side_pair(elem->id(), s);
                  libmesh_assert (!side_id_map->count(side_pair));
                  (*side_id_map)[side_pair] = next_elem_id;
                  next_elem_id += this->n_processors() + 1;
                }

              elem->build_side_ptr(side, s);
              for (auto n : side->node_index_range())
                {
                  const Node & node = side->node_ref(n);

                  // In parallel we don't know enough to number
                  // others' nodes ourselves.
                  if (!hit_end_el &&
                      (node.processor_id() != this->processor_id()))
                    continue;

                  dof_id_type node_id = node.id();
                  if (node_id_map && !node_id_map->count(node_id))
                    {
                      (*node_id_map)[node_id] = next_node_id;
                      next_node_id += this->n_processors() + 1;
                    }
                }
            }
        }
    }

  // FIXME: ought to renumber side/node_id_map image to be contiguous
  // to save memory, also ought to reserve memory
}


} // namespace libMesh

#endif // LIBMESH_BOUNDARY_INFO_IMPL_H

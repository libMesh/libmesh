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
#include <iterator>  // std::distance

// Local includes
#include "libmesh/libmesh_config.h"

#include "libmesh/boundary_info.h"
#include "libmesh/distributed_mesh.h"
#include "libmesh/elem.h"
#include "libmesh/mesh_communication.h"
#include "libmesh/mesh_serializer.h"
#include "libmesh/parallel.h"
#include "libmesh/partitioner.h"
#include "libmesh/remote_elem.h"
#include "libmesh/unstructured_mesh.h"

namespace libMesh
{



//------------------------------------------------------
// BoundaryInfo static member initializations
const boundary_id_type BoundaryInfo::invalid_id = -123;



//------------------------------------------------------
// BoundaryInfo functions
BoundaryInfo::BoundaryInfo(MeshBase & m) :
  ParallelObject(m.comm()),
  _mesh (m)
{
}

BoundaryInfo & BoundaryInfo::operator=(const BoundaryInfo & other_boundary_info)
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
  {
    boundary_node_iter it = other_boundary_info._boundary_node_id.begin();
    const boundary_node_iter end = other_boundary_info._boundary_node_id.end();

    for (; it != end; ++it)
      {
        const Node * other_node = it->first;
        _boundary_node_id.insert(std::make_pair(_mesh.node_ptr(other_node->id()),
                                                it->second));
      }
  }

  // Copy edge boundary info
  {
    boundary_edge_iter it = other_boundary_info._boundary_edge_id.begin();
    const boundary_edge_iter end = other_boundary_info._boundary_edge_id.end();

    for (; it != end; ++it)
      {
        const Elem * other_elem = it->first;
        _boundary_edge_id.insert(std::make_pair(_mesh.elem_ptr(other_elem->id()),
                                                it->second));
      }
  }

  // Copy shellface boundary info
  {
    boundary_shellface_iter it = other_boundary_info._boundary_shellface_id.begin();
    const boundary_shellface_iter end = other_boundary_info._boundary_shellface_id.end();

    for (; it != end; ++it)
      {
        const Elem * other_elem = it->first;
        _boundary_shellface_id.insert(std::make_pair(_mesh.elem_ptr(other_elem->id()),
                                                     it->second));
      }
  }

  // Copy side boundary info
  {
    boundary_side_iter it = other_boundary_info._boundary_side_id.begin();
    const boundary_side_iter end = other_boundary_info._boundary_side_id.end();

    for (; it != end; ++it)
      {
        const Elem * other_elem = it->first;
        _boundary_side_id.insert(std::make_pair(_mesh.elem_ptr(other_elem->id()),
                                                it->second));
      }
  }

  _boundary_ids = other_boundary_info._boundary_ids;
  _side_boundary_ids = other_boundary_info._side_boundary_ids;
  _node_boundary_ids = other_boundary_info._node_boundary_ids;
  _edge_boundary_ids = other_boundary_info._edge_boundary_ids;
  _shellface_boundary_ids = other_boundary_info._shellface_boundary_ids;

  return *this;
}


BoundaryInfo::~BoundaryInfo()
{
  this->clear();
}



void BoundaryInfo::clear()
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



void BoundaryInfo::regenerate_id_sets()
{
  // Clear the old caches
  _boundary_ids.clear();
  _side_boundary_ids.clear();
  _node_boundary_ids.clear();
  _edge_boundary_ids.clear();
  _shellface_boundary_ids.clear();

  // Loop over id maps to regenerate each set.
  for (boundary_node_iter it = _boundary_node_id.begin(),
         end = _boundary_node_id.end();
       it != end; ++it)
    {
      const boundary_id_type id = it->second;
      _boundary_ids.insert(id);
      _node_boundary_ids.insert(id);
    }

  for (boundary_edge_iter it = _boundary_edge_id.begin(),
         end = _boundary_edge_id.end();
       it != end; ++it)
    {
      const boundary_id_type id = it->second.second;
      _boundary_ids.insert(id);
      _edge_boundary_ids.insert(id);
    }

  for (boundary_side_iter it = _boundary_side_id.begin(),
         end = _boundary_side_id.end();
       it != end; ++it)
    {
      const boundary_id_type id = it->second.second;
      _boundary_ids.insert(id);
      _side_boundary_ids.insert(id);
    }

  for (boundary_shellface_iter it = _boundary_shellface_id.begin(),
         end = _boundary_shellface_id.end();
       it != end; ++it)
    {
      const boundary_id_type id = it->second.second;
      _boundary_ids.insert(id);
      _shellface_boundary_ids.insert(id);
    }
}



void BoundaryInfo::sync (UnstructuredMesh & boundary_mesh)
{
  std::set<boundary_id_type> request_boundary_ids(_boundary_ids);
  request_boundary_ids.insert(invalid_id);
  if (!_mesh.is_serial())
    this->comm().set_union(request_boundary_ids);

  this->sync(request_boundary_ids,
             boundary_mesh);
}


void BoundaryInfo::sync (const std::set<boundary_id_type> & requested_boundary_ids,
                         UnstructuredMesh & boundary_mesh)
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

  this->_find_id_maps(requested_boundary_ids, 0, &node_id_map, 0, &side_id_map);

  // Let's add all the boundary nodes we found to the boundary mesh

  MeshBase::const_node_iterator n_end  = _mesh.nodes_end();

  for(MeshBase::const_node_iterator n_it = _mesh.nodes_begin();
      n_it != n_end; ++n_it)
    {
      const Node * node = *n_it;
      dof_id_type node_id = node->id();
      if (node_id_map.count(node_id))
        {
          boundary_mesh.add_point(*node, node_id_map[node_id], node->processor_id());

          // Copy over all the node's boundary IDs to boundary_mesh
          std::vector<boundary_id_type> node_boundary_ids;
          this->boundary_ids(node, node_boundary_ids);
          for (std::size_t index=0; index<node_boundary_ids.size(); index++)
            {
              boundary_mesh.boundary_info->add_node(node_id_map[node_id],
                                                    node_boundary_ids[index]);
            }
        }
    }

  // Let's add the elements
  this->add_elements (requested_boundary_ids, boundary_mesh);

  // The new elements are currently using the interior mesh's nodes;
  // we want them to use the boundary mesh's nodes instead.

  // This side's Node pointers still point to the nodes of the original mesh.
  // We need to re-point them to the boundary mesh's nodes!  Since we copied *ALL* of
  // the original mesh's nodes over, we should be guaranteed to have the same ordering.

  const MeshBase::element_iterator end_bdy_el =
    boundary_mesh.elements_end();

  for (MeshBase::element_iterator el = boundary_mesh.elements_begin();
       el != end_bdy_el; ++el)
    {
      Elem * new_elem = *el;

      for (unsigned int nn=0; nn<new_elem->n_nodes(); ++nn)
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
  boundary_mesh.partitioner().reset(libmesh_nullptr);

  // Make boundary_mesh nodes and elements contiguous
  boundary_mesh.prepare_for_use(/*skip_renumber =*/ false);

  // and finally distribute element partitioning to the nodes
  Partitioner::set_node_processor_ids(boundary_mesh);
}


void BoundaryInfo::get_side_and_node_maps (UnstructuredMesh & boundary_mesh,
                                           std::map<dof_id_type, dof_id_type> & node_id_map,
                                           std::map<dof_id_type, unsigned char> & side_id_map,
                                           Real tolerance)
{
  LOG_SCOPE("get_side_and_node_maps()", "BoundaryInfo");

  node_id_map.clear();
  side_id_map.clear();

  MeshBase::const_element_iterator el =
    boundary_mesh.active_elements_begin();
  const MeshBase::const_element_iterator end_el =
    boundary_mesh.active_elements_end();

  for (; el != end_el; ++el)
    {
      const Elem * boundary_elem = *el;
      const Elem * interior_parent = boundary_elem->interior_parent();

      // Find out which side of interior_parent boundary_elem correponds to.
      // Use centroid comparison as a way to check.
      unsigned char interior_parent_side_index = 0;
      bool found_matching_sides = false;
      for(unsigned char side=0; side<interior_parent->n_sides(); side++)
        {
          UniquePtr<const Elem> interior_parent_side = interior_parent->build_side_ptr(side);
          Real centroid_distance = (boundary_elem->centroid() - interior_parent_side->centroid()).norm();

          if( centroid_distance < (tolerance * boundary_elem->hmin()) )
            {
              interior_parent_side_index = side;
              found_matching_sides = true;
              break;
            }
        }

      if(!found_matching_sides)
        {
          libmesh_error_msg("No matching side found within the specified tolerance");
        }

      side_id_map[boundary_elem->id()] = interior_parent_side_index;

      UniquePtr<const Elem> interior_parent_side = interior_parent->build_side_ptr(interior_parent_side_index);
      for(unsigned char local_node_index=0; local_node_index<boundary_elem->n_nodes(); local_node_index++)
        {
          dof_id_type boundary_node_id = boundary_elem->node_id(local_node_index);
          dof_id_type interior_node_id = interior_parent_side->node_id(local_node_index);

          node_id_map[interior_node_id] = boundary_node_id;
        }

    }
}


void BoundaryInfo::add_elements(const std::set<boundary_id_type> & requested_boundary_ids,
                                UnstructuredMesh & boundary_mesh)
{
  LOG_SCOPE("add_elements()", "BoundaryInfo");

  // We're not prepared to mix serial and distributed meshes in this
  // method, so make sure they match from the start.
  libmesh_assert_equal_to(_mesh.is_serial(),
                          boundary_mesh.is_serial());

  std::map<std::pair<dof_id_type, unsigned char>, dof_id_type> side_id_map;
  this->_find_id_maps(requested_boundary_ids,
                      0,
                      libmesh_nullptr,
                      boundary_mesh.max_elem_id(),
                      &side_id_map);

  // We have to add sides *outside* any element loop, because if
  // boundary_mesh and _mesh are the same then those additions can
  // invalidate our element iterators.  So we just use the element
  // loop to make a list of sides to add.
  typedef std::vector<std::pair<dof_id_type, unsigned char> >
    side_container;
  side_container sides_to_add;

  const MeshBase::const_element_iterator end_el = _mesh.elements_end();
  for (MeshBase::const_element_iterator el = _mesh.elements_begin();
       el != end_el; ++el)
    {
      const Elem * elem = *el;

      for (unsigned int s=0; s<elem->n_sides(); s++)
        if (elem->neighbor_ptr(s) == libmesh_nullptr) // on the boundary
          {
            // Get the top-level parent for this element
            const Elem * top_parent = elem->top_parent();

            // Find the right id number for that side
            std::pair<boundary_side_iter, boundary_side_iter> pos = _boundary_side_id.equal_range(top_parent);

            bool add_this_side = false;
            boundary_id_type this_bcid = invalid_id;

            for (; pos.first != pos.second; ++pos.first)
              {
                this_bcid = pos.first->second.second;

                // if this side is flagged with a boundary condition
                // and the user wants this id
                if ((pos.first->second.first == s) &&
                    (requested_boundary_ids.count(this_bcid)))
                  {
                    add_this_side = true;
                    break;
                  }
              }

            // if side s wasn't found or doesn't have a boundary
            // condition we may still want to add it
            if (pos.first == pos.second)
              {
                this_bcid = invalid_id;
                if (requested_boundary_ids.count(this_bcid))
                  add_this_side = true;
              }

            if (add_this_side)
              sides_to_add.push_back
                (std::make_pair(elem->id(), s));
          }
    }

#ifdef LIBMESH_ENABLE_UNIQUE_ID
  unique_id_type old_max_unique_id = boundary_mesh.parallel_max_unique_id();
#endif

  for (side_container::const_iterator it = sides_to_add.begin();
       it != sides_to_add.end(); ++it)
    {
      const dof_id_type elem_id = it->first;
      const unsigned char s = it->second;
      Elem * elem = _mesh.elem_ptr(elem_id);

      // Build the side - do not use a "proxy" element here:
      // This will be going into the boundary_mesh and needs to
      // stand on its own.
      UniquePtr<Elem> side (elem->build_side_ptr(s, false));

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
          for (unsigned int v=0; v != new_elem->n_vertices(); ++v)
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

      new_elem->set_interior_parent (const_cast<Elem *>(elem));

      // On non-local elements on DistributedMesh we might have
      // RemoteElem neighbor links to construct
      if (!_mesh.is_serial() &&
          (elem->processor_id() != this->processor_id()))
        {
          // Check every interior side for a RemoteElem
          for (unsigned int interior_side = 0;
               interior_side != elem->n_sides();
               ++interior_side)
            {
              // Might this interior side have a RemoteElem that
              // needs a corresponding Remote on a boundary side?
              if (elem->neighbor_ptr(interior_side) != remote_elem)
                continue;

              // Which boundary side?
              for (unsigned int boundary_side = 0;
                   boundary_side != new_elem->n_sides();
                   ++boundary_side)
                {
                  // Look for matching node points.  This is safe in
                  // *this* context.
                  bool found_all_nodes = true;
                  for (unsigned int boundary_node = 0;
                       boundary_node != new_elem->n_nodes();
                       ++boundary_node)
                    {
                      if (!new_elem->is_node_on_side(boundary_node,
                                                     boundary_side))
                        continue;

                      bool found_this_node = false;
                      for (unsigned int interior_node = 0;
                           interior_node != elem->n_nodes();
                           ++interior_node)
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
                         const_cast<RemoteElem *>(remote_elem));
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



void BoundaryInfo::add_node(const dof_id_type node,
                            const boundary_id_type id)
{
  this->add_node (_mesh.node_ptr(node), id);
}



void BoundaryInfo::add_node(const Node * node,
                            const boundary_id_type id)
{
  if (id == invalid_id)
    libmesh_error_msg("ERROR: You may not set a boundary ID of "   \
                      << invalid_id                                \
                      << "\n That is reserved for internal use.");

  // Don't add the same ID twice
  std::pair<boundary_node_iter, boundary_node_iter> pos = _boundary_node_id.equal_range(node);

  for (; pos.first != pos.second; ++pos.first)
    if (pos.first->second == id)
      return;

  _boundary_node_id.insert(std::make_pair(node, id));
  _boundary_ids.insert(id);
  _node_boundary_ids.insert(id); // Also add this ID to the set of node boundary IDs
}



void BoundaryInfo::add_node(const Node * node,
                            const std::vector<boundary_id_type> & ids)
{
  if (ids.empty())
    return;

  libmesh_assert(node);

  // Don't add the same ID twice
  std::pair<boundary_node_iter, boundary_node_iter> pos = _boundary_node_id.equal_range(node);

  // The entries in the ids vector may be non-unique.  If we expected
  // *lots* of ids, it might be fastest to construct a std::set from
  // the entries, but for a small number of entries, which is more
  // typical, it is probably faster to copy the vector and do sort+unique.
  // http://stackoverflow.com/questions/1041620/whats-the-most-efficient-way-to-erase-duplicates-and-sort-a-vector
  std::vector<boundary_id_type> unique_ids(ids.begin(), ids.end());
  std::sort(unique_ids.begin(), unique_ids.end());
  std::vector<boundary_id_type>::iterator new_end =
    std::unique(unique_ids.begin(), unique_ids.end());

  std::vector<boundary_id_type>::iterator it = unique_ids.begin();
  for (; it != new_end; ++it)
    {
      boundary_id_type id = *it;

      if (id == invalid_id)
        libmesh_error_msg("ERROR: You may not set a boundary ID of "    \
                          << invalid_id                                 \
                          << "\n That is reserved for internal use.");

      bool already_inserted = false;
      for (boundary_node_iter p = pos.first; p != pos.second; ++p)
        if (p->second == id)
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



void BoundaryInfo::clear_boundary_node_ids()
{
  _boundary_node_id.clear();
}

void BoundaryInfo::add_edge(const dof_id_type e,
                            const unsigned short int edge,
                            const boundary_id_type id)
{
  this->add_edge (_mesh.elem_ptr(e), edge, id);
}



void BoundaryInfo::add_edge(const Elem * elem,
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
  std::pair<boundary_edge_iter, boundary_edge_iter> pos = _boundary_edge_id.equal_range(elem);

  for (; pos.first != pos.second; ++pos.first)
    if (pos.first->second.first == edge &&
        pos.first->second.second == id)
      return;

  _boundary_edge_id.insert(std::make_pair(elem, std::make_pair(edge, id)));
  _boundary_ids.insert(id);
  _edge_boundary_ids.insert(id); // Also add this ID to the set of edge boundary IDs
}



void BoundaryInfo::add_edge(const Elem * elem,
                            const unsigned short int edge,
                            const std::vector<boundary_id_type> & ids)
{
  if (ids.empty())
    return;

  libmesh_assert(elem);

  // Only add BCs for level-0 elements.
  libmesh_assert_equal_to (elem->level(), 0);

  // Don't add the same ID twice
  std::pair<boundary_edge_iter, boundary_edge_iter> pos = _boundary_edge_id.equal_range(elem);

  // The entries in the ids vector may be non-unique.  If we expected
  // *lots* of ids, it might be fastest to construct a std::set from
  // the entries, but for a small number of entries, which is more
  // typical, it is probably faster to copy the vector and do sort+unique.
  // http://stackoverflow.com/questions/1041620/whats-the-most-efficient-way-to-erase-duplicates-and-sort-a-vector
  std::vector<boundary_id_type> unique_ids(ids.begin(), ids.end());
  std::sort(unique_ids.begin(), unique_ids.end());
  std::vector<boundary_id_type>::iterator new_end =
    std::unique(unique_ids.begin(), unique_ids.end());

  std::vector<boundary_id_type>::iterator it = unique_ids.begin();
  for (; it != new_end; ++it)
    {
      boundary_id_type id = *it;

      if (id == invalid_id)
        libmesh_error_msg("ERROR: You may not set a boundary ID of "   \
                          << invalid_id                                \
                          << "\n That is reserved for internal use.");

      bool already_inserted = false;
      for (boundary_edge_iter p = pos.first; p != pos.second; ++p)
        if (p->second.first == edge &&
            p->second.second == id)
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



void BoundaryInfo::add_shellface(const dof_id_type e,
                                 const unsigned short int shellface,
                                 const boundary_id_type id)
{
  this->add_shellface (_mesh.elem_ptr(e), shellface, id);
}



void BoundaryInfo::add_shellface(const Elem * elem,
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
  std::pair<boundary_shellface_iter, boundary_shellface_iter> pos = _boundary_shellface_id.equal_range(elem);

  for (; pos.first != pos.second; ++pos.first)
    if (pos.first->second.first == shellface &&
        pos.first->second.second == id)
      return;

  _boundary_shellface_id.insert(std::make_pair(elem, std::make_pair(shellface, id)));
  _boundary_ids.insert(id);
  _shellface_boundary_ids.insert(id); // Also add this ID to the set of shellface boundary IDs
}



void BoundaryInfo::add_shellface(const Elem * elem,
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
  std::pair<boundary_shellface_iter, boundary_shellface_iter> pos = _boundary_shellface_id.equal_range(elem);

  // The entries in the ids vector may be non-unique.  If we expected
  // *lots* of ids, it might be fastest to construct a std::set from
  // the entries, but for a small number of entries, which is more
  // typical, it is probably faster to copy the vector and do sort+unique.
  // http://stackoverflow.com/questions/1041620/whats-the-most-efficient-way-to-erase-duplicates-and-sort-a-vector
  std::vector<boundary_id_type> unique_ids(ids.begin(), ids.end());
  std::sort(unique_ids.begin(), unique_ids.end());
  std::vector<boundary_id_type>::iterator new_end =
    std::unique(unique_ids.begin(), unique_ids.end());

  std::vector<boundary_id_type>::iterator it = unique_ids.begin();
  for (; it != new_end; ++it)
    {
      boundary_id_type id = *it;

      if (id == invalid_id)
        libmesh_error_msg("ERROR: You may not set a boundary ID of "   \
                          << invalid_id                                \
                          << "\n That is reserved for internal use.");

      bool already_inserted = false;
      for (boundary_shellface_iter p = pos.first; p != pos.second; ++p)
        if (p->second.first == shellface &&
            p->second.second == id)
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


void BoundaryInfo::add_side(const dof_id_type e,
                            const unsigned short int side,
                            const boundary_id_type id)
{
  this->add_side (_mesh.elem_ptr(e), side, id);
}



void BoundaryInfo::add_side(const Elem * elem,
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
  std::pair<boundary_side_iter, boundary_side_iter> pos = _boundary_side_id.equal_range(elem);

  for (; pos.first != pos.second; ++pos.first)
    if (pos.first->second.first == side &&
        pos.first->second.second == id)
      return;

  _boundary_side_id.insert(std::make_pair(elem, std::make_pair(side, id)));
  _boundary_ids.insert(id);
  _side_boundary_ids.insert(id); // Also add this ID to the set of side boundary IDs
}



void BoundaryInfo::add_side(const Elem * elem,
                            const unsigned short int side,
                            const std::vector<boundary_id_type> & ids)
{
  if (ids.empty())
    return;

  libmesh_assert(elem);

  // Only add BCs for level-0 elements.
  libmesh_assert_equal_to (elem->level(), 0);

  // Don't add the same ID twice
  std::pair<boundary_side_iter, boundary_side_iter> pos = _boundary_side_id.equal_range(elem);

  // The entries in the ids vector may be non-unique.  If we expected
  // *lots* of ids, it might be fastest to construct a std::set from
  // the entries, but for a small number of entries, which is more
  // typical, it is probably faster to copy the vector and do sort+unique.
  // http://stackoverflow.com/questions/1041620/whats-the-most-efficient-way-to-erase-duplicates-and-sort-a-vector
  std::vector<boundary_id_type> unique_ids(ids.begin(), ids.end());
  std::sort(unique_ids.begin(), unique_ids.end());
  std::vector<boundary_id_type>::iterator new_end =
    std::unique(unique_ids.begin(), unique_ids.end());

  std::vector<boundary_id_type>::const_iterator it = unique_ids.begin();
  for (; it != new_end; ++it)
    {
      boundary_id_type id = *it;

      if (id == invalid_id)
        libmesh_error_msg("ERROR: You may not set a boundary ID of "    \
                          << invalid_id                                 \
                          << "\n That is reserved for internal use.");

      bool already_inserted = false;
      for (boundary_side_iter p = pos.first; p != pos.second; ++p)
        if (p->second.first == side && p->second.second == id)
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



bool BoundaryInfo::has_boundary_id(const Node * const node,
                                   const boundary_id_type id) const
{
  std::pair<boundary_node_iter, boundary_node_iter> pos = _boundary_node_id.equal_range(node);

  for (; pos.first != pos.second; ++pos.first)
    if (pos.first->second == id)
      return true;

  return false;
}



std::vector<boundary_id_type> BoundaryInfo::boundary_ids(const Node * node) const
{
  libmesh_deprecated();

  std::vector<boundary_id_type> ids;
  this->boundary_ids(node, ids);
  return ids;
}



void BoundaryInfo::boundary_ids (const Node * node,
                                 std::vector<boundary_id_type> & vec_to_fill) const
{
  // Clear out any previous contents
  vec_to_fill.clear();

  std::pair<boundary_node_iter, boundary_node_iter>
    pos = _boundary_node_id.equal_range(node);

  for (; pos.first != pos.second; ++pos.first)
    vec_to_fill.push_back(pos.first->second);
}



unsigned int BoundaryInfo::n_boundary_ids(const Node * node) const
{
  std::pair<boundary_node_iter, boundary_node_iter> pos = _boundary_node_id.equal_range(node);

  return cast_int<unsigned int>(std::distance(pos.first, pos.second));
}



std::vector<boundary_id_type> BoundaryInfo::edge_boundary_ids (const Elem * const elem,
                                                               const unsigned short int edge) const
{
  libmesh_deprecated();

  std::vector<boundary_id_type> ids;
  this->edge_boundary_ids(elem, edge, ids);
  return ids;
}



void BoundaryInfo::edge_boundary_ids (const Elem * const elem,
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
      for (unsigned int side=0; side<elem->n_sides(); side++)
        {
          if (elem->is_edge_on_side(edge,side))
            {
              if (elem->neighbor_ptr(side) == libmesh_nullptr)
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
          while (searched_elem->parent() != libmesh_nullptr)
            {
              const Elem * parent = searched_elem->parent();
              if (parent->is_child_on_edge(parent->which_child_am_i(searched_elem), edge) == false)
                return;
              searched_elem = parent;
            }
        }
    }
#endif

  std::pair<boundary_edge_iter, boundary_edge_iter>
    e = _boundary_edge_id.equal_range(searched_elem);

  // Check each element in the range to see if its edge matches the requested edge.
  for (; e.first != e.second; ++e.first)
    if (e.first->second.first == edge)
      vec_to_fill.push_back(e.first->second.second);
}



unsigned int BoundaryInfo::n_edge_boundary_ids (const Elem * const elem,
                                                const unsigned short int edge) const
{
  std::vector<boundary_id_type> ids;
  this->edge_boundary_ids(elem, edge, ids);
  return ids.size();
}



std::vector<boundary_id_type> BoundaryInfo::raw_edge_boundary_ids (const Elem * const elem,
                                                                   const unsigned short int edge) const
{
  libmesh_deprecated();

  std::vector<boundary_id_type> ids;
  this->raw_edge_boundary_ids(elem, edge, ids);
  return ids;
}



void BoundaryInfo::raw_edge_boundary_ids (const Elem * const elem,
                                          const unsigned short int edge,
                                          std::vector<boundary_id_type> & vec_to_fill) const
{
  libmesh_assert(elem);

  // Clear out any previous contents
  vec_to_fill.clear();

  // Only level-0 elements store BCs.
  if (elem->parent())
    return;

  std::pair<boundary_edge_iter, boundary_edge_iter>
    e = _boundary_edge_id.equal_range(elem);

  // Check each element in the range to see if its edge matches the requested edge.
  for (; e.first != e.second; ++e.first)
    if (e.first->second.first == edge)
      vec_to_fill.push_back(e.first->second.second);
}



void BoundaryInfo::shellface_boundary_ids (const Elem * const elem,
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
      while (searched_elem->parent() != libmesh_nullptr)
        {
          const Elem * parent = searched_elem->parent();
          searched_elem = parent;
        }
    }
#endif

  std::pair<boundary_shellface_iter, boundary_shellface_iter>
    e = _boundary_shellface_id.equal_range(searched_elem);

  // Check each element in the range to see if its shellface matches the requested shellface.
  for (; e.first != e.second; ++e.first)
    if (e.first->second.first == shellface)
      vec_to_fill.push_back(e.first->second.second);
}



unsigned int BoundaryInfo::n_shellface_boundary_ids (const Elem * const elem,
                                                     const unsigned short int shellface) const
{
  std::vector<boundary_id_type> ids;
  this->shellface_boundary_ids(elem, shellface, ids);
  return ids.size();
}



void BoundaryInfo::raw_shellface_boundary_ids (const Elem * const elem,
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

  std::pair<boundary_shellface_iter, boundary_shellface_iter>
    e = _boundary_shellface_id.equal_range(elem);

  // Check each element in the range to see if its shellface matches the requested shellface.
  for (; e.first != e.second; ++e.first)
    if (e.first->second.first == shellface)
      vec_to_fill.push_back(e.first->second.second);
}


boundary_id_type BoundaryInfo::boundary_id(const Elem * const elem,
                                           const unsigned short int side) const
{
  // Asking for just one boundary id means your code isn't safe to use
  // on meshes with overlapping boundary ids.  Try using
  // BoundaryInfo::boundary_ids or BoundaryInfo::has_boundary_id
  // instead.
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



bool BoundaryInfo::has_boundary_id(const Elem * const elem,
                                   const unsigned short int side,
                                   const boundary_id_type id) const
{
  std::vector<boundary_id_type> ids;
  this->boundary_ids(elem, side, ids);
  return (std::find(ids.begin(), ids.end(), id) != ids.end());
}



std::vector<boundary_id_type> BoundaryInfo::boundary_ids (const Elem * const elem,
                                                          const unsigned short int side) const
{
  libmesh_deprecated();

  std::vector<boundary_id_type> ids;
  this->boundary_ids(elem, side, ids);
  return ids;
}



void BoundaryInfo::boundary_ids (const Elem * const elem,
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
      if (elem->neighbor_ptr(side) == libmesh_nullptr)
        searched_elem = elem->top_parent ();
#ifdef LIBMESH_ENABLE_AMR
      else
        while (searched_elem->parent() != libmesh_nullptr)
          {
            const Elem * parent = searched_elem->parent();
            if (parent->is_child_on_side(parent->which_child_am_i(searched_elem), side) == false)
              return;
            searched_elem = parent;
          }
#endif
    }

  std::pair<boundary_side_iter, boundary_side_iter>
    e = _boundary_side_id.equal_range(searched_elem);

  // Check each element in the range to see if its side matches the requested side.
  for (; e.first != e.second; ++e.first)
    if (e.first->second.first == side)
      vec_to_fill.push_back(e.first->second.second);
}




unsigned int BoundaryInfo::n_boundary_ids (const Elem * const elem,
                                           const unsigned short int side) const
{
  std::vector<boundary_id_type> ids;
  this->boundary_ids(elem, side, ids);
  return ids.size();
}



std::vector<boundary_id_type> BoundaryInfo::raw_boundary_ids (const Elem * const elem,
                                                              const unsigned short int side) const
{
  libmesh_deprecated();

  std::vector<boundary_id_type> ids;
  this->raw_boundary_ids(elem, side, ids);
  return ids;
}



void BoundaryInfo::raw_boundary_ids (const Elem * const elem,
                                     const unsigned short int side,
                                     std::vector<boundary_id_type> & vec_to_fill) const
{
  libmesh_assert(elem);

  // Clear out any previous contents
  vec_to_fill.clear();

  // Only level-0 elements store BCs.
  if (elem->parent())
    return;

  std::pair<boundary_side_iter, boundary_side_iter>
    e = _boundary_side_id.equal_range(elem);

  // Check each element in the range to see if its side matches the requested side.
  for (; e.first != e.second; ++e.first)
    if (e.first->second.first == side)
      vec_to_fill.push_back(e.first->second.second);
}



void BoundaryInfo::copy_boundary_ids (const BoundaryInfo & old_boundary_info,
                                      const Elem * const old_elem,
                                      const Elem * const new_elem)
{
  libmesh_assert_equal_to (old_elem->n_sides(), new_elem->n_sides());
  libmesh_assert_equal_to (old_elem->n_edges(), new_elem->n_edges());

  std::vector<boundary_id_type> bndry_ids;

  for (unsigned short s=0; s<old_elem->n_sides(); s++)
    {
      old_boundary_info.raw_boundary_ids (old_elem, s, bndry_ids);
      this->add_side (new_elem, s, bndry_ids);
    }

  for (unsigned short e=0; e<old_elem->n_edges(); e++)
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



void BoundaryInfo::remove (const Node * node)
{
  libmesh_assert(node);

  // Erase everything associated with node
  _boundary_node_id.erase (node);
}



void BoundaryInfo::remove (const Elem * elem)
{
  libmesh_assert(elem);

  // Erase everything associated with elem
  _boundary_edge_id.erase (elem);
  _boundary_side_id.erase (elem);
  _boundary_shellface_id.erase (elem);
}



void BoundaryInfo::remove_edge (const Elem * elem,
                                const unsigned short int edge)
{
  libmesh_assert(elem);

  // The user shouldn't be trying to remove only one child's boundary
  // id
  libmesh_assert_equal_to (elem->level(), 0);

  // Some older compilers don't support erasing from a map with
  // const_iterators, so we explicitly use non-const iterators here.
  std::pair<erase_iter, erase_iter>
    e = _boundary_edge_id.equal_range(elem);

  // elem may be there, maybe multiple occurrences
  while (e.first != e.second)
    {
      // if this is true we found the requested edge
      // of the element and want to erase the id
      if (e.first->second.first == edge)
        {
          // (postfix++ - increment the iterator before it's invalid)
          _boundary_edge_id.erase(e.first++);
        }
      else
        ++e.first;
    }
}



void BoundaryInfo::remove_edge (const Elem * elem,
                                const unsigned short int edge,
                                const boundary_id_type id)
{
  libmesh_assert(elem);

  // The user shouldn't be trying to remove only one child's boundary
  // id
  libmesh_assert_equal_to (elem->level(), 0);

  // Some older compilers don't support erasing from a map with
  // const_iterators, so we explicitly use non-const iterators here.
  std::pair<erase_iter, erase_iter>
    e = _boundary_edge_id.equal_range(elem);

  // elem may be there, maybe multiple occurrences
  while (e.first != e.second)
    {
      // if this is true we found the requested edge
      // of the element and want to erase the requested id
      if (e.first->second.first == edge &&
          e.first->second.second == id)
        {
          // (postfix++ - increment the iterator before it's invalid)
          _boundary_edge_id.erase(e.first++);
        }
      else
        ++e.first;
    }
}


void BoundaryInfo::remove_shellface (const Elem * elem,
                                     const unsigned short int shellface)
{
  libmesh_assert(elem);

  // The user shouldn't be trying to remove only one child's boundary
  // id
  libmesh_assert_equal_to (elem->level(), 0);

  // Shells only have 2 faces
  libmesh_assert_less(shellface, 2);

  // Some older compilers don't support erasing from a map with
  // const_iterators, so we explicitly use non-const iterators here.
  std::pair<erase_iter, erase_iter>
    e = _boundary_shellface_id.equal_range(elem);

  // elem may be there, maybe multiple occurrences
  while (e.first != e.second)
    {
      // if this is true we found the requested shellface
      // of the element and want to erase the id
      if (e.first->second.first == shellface)
        {
          // (postfix++ - increment the iterator before it's invalid)
          _boundary_shellface_id.erase(e.first++);
        }
      else
        ++e.first;
    }
}



void BoundaryInfo::remove_shellface (const Elem * elem,
                                     const unsigned short int shellface,
                                     const boundary_id_type id)
{
  libmesh_assert(elem);

  // The user shouldn't be trying to remove only one child's boundary
  // id
  libmesh_assert_equal_to (elem->level(), 0);

  // Shells only have 2 faces
  libmesh_assert_less(shellface, 2);

  // Some older compilers don't support erasing from a map with
  // const_iterators, so we explicitly use non-const iterators here.
  std::pair<erase_iter, erase_iter>
    e = _boundary_shellface_id.equal_range(elem);

  // elem may be there, maybe multiple occurrences
  while (e.first != e.second)
    {
      // if this is true we found the requested shellface
      // of the element and want to erase the requested id
      if (e.first->second.first == shellface &&
          e.first->second.second == id)
        {
          // (postfix++ - increment the iterator before it's invalid)
          _boundary_shellface_id.erase(e.first++);
        }
      else
        ++e.first;
    }
}

void BoundaryInfo::remove_side (const Elem * elem,
                                const unsigned short int side)
{
  libmesh_assert(elem);

  // The user shouldn't be trying to remove only one child's boundary
  // id
  libmesh_assert_equal_to (elem->level(), 0);

  // Some older compilers don't support erasing from a map with
  // const_iterators, so we explicitly use non-const iterators here.
  std::pair<erase_iter, erase_iter>
    e = _boundary_side_id.equal_range(elem);

  // elem may be there, maybe multiple occurrences
  while (e.first != e.second)
    {
      // if this is true we found the requested side
      // of the element and want to erase the id
      if (e.first->second.first == side)
        {
          // (postfix++ - increment the iterator before it's invalid)
          _boundary_side_id.erase(e.first++);
        }
      else
        ++e.first;
    }
}



void BoundaryInfo::remove_side (const Elem * elem,
                                const unsigned short int side,
                                const boundary_id_type id)
{
  libmesh_assert(elem);

  // Some older compilers don't support erasing from a map with
  // const_iterators, so we explicitly use non-const iterators here.
  std::pair<erase_iter, erase_iter>
    e = _boundary_side_id.equal_range(elem);

  // elem may be there, maybe multiple occurrences
  while (e.first != e.second)
    {
      // if this is true we found the requested side
      // of the element and want to erase the requested id
      if (e.first->second.first == side &&
          e.first->second.second == id)
        {
          // (postfix++ - increment the iterator before it's invalid)
          _boundary_side_id.erase(e.first++);
        }
      else
        ++e.first;
    }
}



void BoundaryInfo::remove_id (boundary_id_type id)
{
  // Erase id from ids containers
  _boundary_ids.erase(id);
  _side_boundary_ids.erase(id);
  _edge_boundary_ids.erase(id);
  _shellface_boundary_ids.erase(id);
  _node_boundary_ids.erase(id);
  _ss_id_to_name.erase(id);
  _ns_id_to_name.erase(id);

  // Erase pointers to geometric entities with this id.
  for (boundary_node_erase_iter it = _boundary_node_id.begin(); it != _boundary_node_id.end(); /*below*/)
    {
      if (it->second == id)
        _boundary_node_id.erase(it++);
      else
        ++it;
    }

  for (erase_iter it = _boundary_edge_id.begin(); it != _boundary_edge_id.end(); /*below*/)
    {
      if (it->second.second == id)
        _boundary_edge_id.erase(it++);
      else
        ++it;
    }

  for (erase_iter it = _boundary_shellface_id.begin(); it != _boundary_shellface_id.end(); /*below*/)
    {
      if (it->second.second == id)
        _boundary_shellface_id.erase(it++);
      else
        ++it;
    }

  for (erase_iter it = _boundary_side_id.begin(); it != _boundary_side_id.end(); /*below*/)
    {
      if (it->second.second == id)
        _boundary_side_id.erase(it++);
      else
        ++it;
    }
}



unsigned int BoundaryInfo::side_with_boundary_id(const Elem * const elem,
                                                 const boundary_id_type boundary_id_in) const
{
  const Elem * searched_elem = elem;
  if (elem->level() != 0)
    searched_elem = elem->top_parent();

  std::pair<boundary_side_iter, boundary_side_iter>
    e = _boundary_side_id.equal_range(searched_elem);

  // elem may have zero or multiple occurrences
  for (; e.first != e.second; ++e.first)
    {
      // if this is true we found the requested boundary_id
      // of the element and want to return the side
      if (e.first->second.second == boundary_id_in)
        {
          unsigned int side = e.first->second.first;

          // If we're on this external boundary then we share this
          // external boundary id
          if (elem->neighbor_ptr(side) == libmesh_nullptr)
            return side;

          // If we're on an internal boundary then we need to be sure
          // it's the same internal boundary as our top_parent
          const Elem * p = elem;

#ifdef LIBMESH_ENABLE_AMR

          while (p != libmesh_nullptr)
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

void
BoundaryInfo::build_node_boundary_ids(std::vector<boundary_id_type> & b_ids) const
{
  b_ids.clear();

  boundary_node_iter pos = _boundary_node_id.begin();
  for (; pos != _boundary_node_id.end(); ++pos)
    {
      boundary_id_type id = pos->second;

      if (std::find(b_ids.begin(),b_ids.end(),id) == b_ids.end())
        b_ids.push_back(id);
    }
}

void
BoundaryInfo::build_side_boundary_ids(std::vector<boundary_id_type> & b_ids) const
{
  b_ids.clear();

  boundary_side_iter pos = _boundary_side_id.begin();
  for (; pos != _boundary_side_id.end(); ++pos)
    {
      boundary_id_type id = pos->second.second;

      if (std::find(b_ids.begin(),b_ids.end(),id) == b_ids.end())
        b_ids.push_back(id);
    }
}

void
BoundaryInfo::build_shellface_boundary_ids(std::vector<boundary_id_type> & b_ids) const
{
  b_ids.clear();

  boundary_side_iter pos = _boundary_shellface_id.begin();
  for (; pos != _boundary_shellface_id.end(); ++pos)
    {
      boundary_id_type id = pos->second.second;

      if (std::find(b_ids.begin(),b_ids.end(),id) == b_ids.end())
        b_ids.push_back(id);
    }
}

std::size_t BoundaryInfo::n_boundary_conds () const
{
  // in serial we know the number of bcs from the
  // size of the container
  if (_mesh.is_serial())
    return _boundary_side_id.size();

  // in parallel we need to sum the number of local bcs
  parallel_object_only();

  std::size_t nbcs=0;

  boundary_side_iter pos = _boundary_side_id.begin();
  for (; pos != _boundary_side_id.end(); ++pos)
    if (pos->first->processor_id() == this->processor_id())
      nbcs++;

  this->comm().sum (nbcs);

  return nbcs;
}

std::size_t BoundaryInfo::n_edge_conds () const
{
  // in serial we know the number of nodesets from the
  // size of the container
  if (_mesh.is_serial())
    return _boundary_edge_id.size();

  // in parallel we need to sum the number of local nodesets
  parallel_object_only();

  std::size_t n_edge_bcs=0;

  boundary_edge_iter pos = _boundary_edge_id.begin();
  for (; pos != _boundary_edge_id.end(); ++pos)
    if (pos->first->processor_id() == this->processor_id())
      n_edge_bcs++;

  this->comm().sum (n_edge_bcs);

  return n_edge_bcs;
}


std::size_t BoundaryInfo::n_shellface_conds () const
{
  // in serial we know the number of nodesets from the
  // size of the container
  if (_mesh.is_serial())
    return _boundary_shellface_id.size();

  // in parallel we need to sum the number of local nodesets
  parallel_object_only();

  std::size_t n_shellface_bcs=0;

  boundary_shellface_iter pos = _boundary_shellface_id.begin();
  for (; pos != _boundary_shellface_id.end(); ++pos)
    if (pos->first->processor_id() == this->processor_id())
      n_shellface_bcs++;

  this->comm().sum (n_shellface_bcs);

  return n_shellface_bcs;
}


std::size_t BoundaryInfo::n_nodeset_conds () const
{
  // in serial we know the number of nodesets from the
  // size of the container
  if (_mesh.is_serial())
    return _boundary_node_id.size();

  // in parallel we need to sum the number of local nodesets
  parallel_object_only();

  std::size_t n_nodesets=0;

  boundary_node_iter pos = _boundary_node_id.begin();
  for (; pos != _boundary_node_id.end(); ++pos)
    if (pos->first->processor_id() == this->processor_id())
      n_nodesets++;

  this->comm().sum (n_nodesets);

  return n_nodesets;
}



void BoundaryInfo::build_node_list (std::vector<dof_id_type> & nl,
                                    std::vector<boundary_id_type> & il) const
{
  // Clear the input vectors, just in case they were used for
  // something else recently...
  nl.clear();
  il.clear();

  // Reserve the size, then use push_back
  nl.reserve (_boundary_node_id.size());
  il.reserve (_boundary_node_id.size());

  boundary_node_iter pos = _boundary_node_id.begin();
  for (; pos != _boundary_node_id.end(); ++pos)
    {
      nl.push_back (pos->first->id());
      il.push_back (pos->second);
    }
}


void
BoundaryInfo::build_node_list_from_side_list()
{
  // Loop over the side list
  boundary_side_iter pos = _boundary_side_id.begin();
  for (; pos != _boundary_side_id.end(); ++pos)
    {
      // Don't add remote sides
      if (pos->first->is_remote())
        continue;

      // Need to loop over the sides of any possible children
      std::vector<const Elem *> family;
#ifdef LIBMESH_ENABLE_AMR
      pos->first->active_family_tree_by_side (family, pos->second.first);
#else
      family.push_back(pos->first);
#endif

      for (std::size_t elem_it=0; elem_it < family.size(); elem_it++)
        {
          const Elem * cur_elem = family[elem_it];

          UniquePtr<const Elem> side = cur_elem->build_side_ptr(pos->second.first);

          // Add each node node on the side with the side's boundary id
          for (unsigned int i=0; i<side->n_nodes(); i++)
            this->add_node(side->node_ptr(i), pos->second.second);
        }
    }
}




void BoundaryInfo::build_side_list_from_node_list()
{
  // Check for early return
  if (_boundary_node_id.empty())
    {
      libMesh::out << "No boundary node IDs have been added: cannot build side list!" << std::endl;
      return;
    }

  MeshBase::const_element_iterator el = _mesh.active_elements_begin();
  const MeshBase::const_element_iterator end_el = _mesh.active_elements_end();

  for (; el != end_el; ++el)
    {
      const Elem * elem = *el;

      for (unsigned short side=0; side<elem->n_sides(); ++side)
        {
          UniquePtr<const Elem> side_elem = elem->build_side_ptr(side);

          // map from nodeset_id to count for that ID
          std::map<boundary_id_type, unsigned> nodesets_node_count;
          for (unsigned node_num=0; node_num < side_elem->n_nodes(); ++node_num)
            {
              const Node * node = side_elem->node_ptr(node_num);
              std::pair<boundary_node_iter, boundary_node_iter>
                range = _boundary_node_id.equal_range(node);

              // For each nodeset that this node is a member of, increment the associated
              // nodeset ID count
              for (boundary_node_iter pos = range.first; pos != range.second; ++pos)
                nodesets_node_count[pos->second]++;
            }

          // Now check to see what nodeset_counts have the correct
          // number of nodes in them.  For any that do, add this side to
          // the sideset, making sure the sideset inherits the
          // nodeset's name, if there is one.
          std::map<boundary_id_type, unsigned>::const_iterator nodesets = nodesets_node_count.begin();
          for (; nodesets != nodesets_node_count.end(); ++nodesets)
            if (nodesets->second == side_elem->n_nodes())
              {
                add_side(elem, side, nodesets->first);

                // Let the sideset inherit any non-empty name from the nodeset
                std::string & nset_name = nodeset_name(nodesets->first);

                if (nset_name != "")
                  sideset_name(nodesets->first) = nset_name;
              }
        } // end for side
    } // end for el
}




void BoundaryInfo::build_side_list (std::vector<dof_id_type> & el,
                                    std::vector<unsigned short int> & sl,
                                    std::vector<boundary_id_type> & il) const
{
  // Clear the input vectors, just in case they were used for
  // something else recently...
  el.clear();
  sl.clear();
  il.clear();

  // Reserve the size, then use push_back
  el.reserve (_boundary_side_id.size());
  sl.reserve (_boundary_side_id.size());
  il.reserve (_boundary_side_id.size());

  boundary_side_iter pos = _boundary_side_id.begin();
  for (; pos != _boundary_side_id.end(); ++pos)
    {
      el.push_back (pos->first->id());
      sl.push_back (pos->second.first);
      il.push_back (pos->second.second);
    }
}

void BoundaryInfo::build_active_side_list (std::vector<dof_id_type> & el,
                                           std::vector<unsigned short int> & sl,
                                           std::vector<boundary_id_type> & il) const
{
  // Clear the input vectors, just in case they were used for
  // something else recently...
  el.clear();
  sl.clear();
  il.clear();

  boundary_side_iter pos = _boundary_side_id.begin();
  for (; pos != _boundary_side_id.end(); ++pos)
    {
      // Don't add remote sides
      if (pos->first->is_remote())
        continue;

      // Loop over the sides of possible children
      std::vector< const Elem * > family;
#ifdef LIBMESH_ENABLE_AMR
      pos->first->active_family_tree_by_side(family, pos->second.first);
#else
      family.push_back(pos->first);
#endif

      // Populate the list items
      for (std::vector<const Elem *>::iterator elem_it = family.begin(); elem_it != family.end(); elem_it++)
        {
          el.push_back ((*elem_it)->id());
          sl.push_back (pos->second.first);
          il.push_back (pos->second.second);
        }
    }
}


void BoundaryInfo::build_edge_list (std::vector<dof_id_type> & el,
                                    std::vector<unsigned short int> & sl,
                                    std::vector<boundary_id_type> & il) const
{
  // Clear the input vectors, just in case they were used for
  // something else recently...
  el.clear();
  sl.clear();
  il.clear();

  // Reserve the size, then use push_back
  el.reserve (_boundary_edge_id.size());
  sl.reserve (_boundary_edge_id.size());
  il.reserve (_boundary_edge_id.size());

  boundary_edge_iter pos = _boundary_edge_id.begin();
  for (; pos != _boundary_edge_id.end(); ++pos)
    {
      el.push_back (pos->first->id());
      sl.push_back (pos->second.first);
      il.push_back (pos->second.second);
    }
}


void BoundaryInfo::build_shellface_list (std::vector<dof_id_type> & el,
                                         std::vector<unsigned short int> & sl,
                                         std::vector<boundary_id_type> & il) const
{
  // Clear the input vectors, just in case they were used for
  // something else recently...
  el.clear();
  sl.clear();
  il.clear();

  // Reserve the size, then use push_back
  el.reserve (_boundary_shellface_id.size());
  sl.reserve (_boundary_shellface_id.size());
  il.reserve (_boundary_shellface_id.size());

  boundary_shellface_iter pos = _boundary_shellface_id.begin();
  for (; pos != _boundary_shellface_id.end(); ++pos)
    {
      el.push_back (pos->first->id());
      sl.push_back (pos->second.first);
      il.push_back (pos->second.second);
    }
}


void BoundaryInfo::print_info(std::ostream & out_stream) const
{
  // Print out the nodal BCs
  if (!_boundary_node_id.empty())
    {
      out_stream << "Nodal Boundary conditions:" << std::endl
                 << "--------------------------" << std::endl
                 << "  (Node No., ID)               " << std::endl;

      boundary_node_iter it        = _boundary_node_id.begin();
      const boundary_node_iter end = _boundary_node_id.end();

      for (; it != end; ++it)
        out_stream << "  (" << (*it).first->id()
                   << ", "  << (*it).second
                   << ")"  << std::endl;
    }

  // Print out the element edge BCs
  if (!_boundary_edge_id.empty())
    {
      out_stream << std::endl
                 << "Edge Boundary conditions:" << std::endl
                 << "-------------------------" << std::endl
                 << "  (Elem No., Edge No., ID)      " << std::endl;

      boundary_edge_iter it = _boundary_edge_id.begin();
      const boundary_edge_iter end = _boundary_edge_id.end();

      for (; it != end; ++it)
        out_stream << "  (" << (*it).first->id()
                   << ", "  << (*it).second.first
                   << ", "  << (*it).second.second
                   << ")"   << std::endl;
    }

  // Print out the element shellface BCs
  if (!_boundary_shellface_id.empty())
    {
      out_stream << std::endl
                 << "Shell-face Boundary conditions:" << std::endl
                 << "-------------------------" << std::endl
                 << "  (Elem No., Shell-face No., ID)      " << std::endl;

      boundary_shellface_iter it = _boundary_shellface_id.begin();
      const boundary_shellface_iter end = _boundary_shellface_id.end();

      for (; it != end; ++it)
        out_stream << "  (" << (*it).first->id()
                   << ", "  << (*it).second.first
                   << ", "  << (*it).second.second
                   << ")"   << std::endl;
    }

  // Print out the element side BCs
  if (!_boundary_side_id.empty())
    {
      out_stream << std::endl
                 << "Side Boundary conditions:" << std::endl
                 << "-------------------------" << std::endl
                 << "  (Elem No., Side No., ID)      " << std::endl;

      boundary_side_iter it = _boundary_side_id.begin();
      const boundary_side_iter end = _boundary_side_id.end();

      for (; it != end; ++it)
        out_stream << "  (" << (*it).first->id()
                   << ", "  << (*it).second.first
                   << ", "  << (*it).second.second
                   << ")"   << std::endl;
    }
}



void BoundaryInfo::print_summary(std::ostream & out_stream) const
{
  // Print out the nodal BCs
  if (!_boundary_node_id.empty())
    {
      out_stream << "Nodal Boundary conditions:" << std::endl
                 << "--------------------------" << std::endl
                 << "  (ID, number of nodes)   " << std::endl;

      std::map<boundary_id_type, std::size_t> ID_counts;

      boundary_node_iter it        = _boundary_node_id.begin();
      const boundary_node_iter end = _boundary_node_id.end();

      for (; it != end; ++it)
        ID_counts[(*it).second]++;

      std::map<boundary_id_type, std::size_t>::const_iterator ID_it        = ID_counts.begin();
      const std::map<boundary_id_type, std::size_t>::const_iterator ID_end = ID_counts.end();

      for (; ID_it != ID_end; ++ID_it)
        out_stream << "  (" << (*ID_it).first
                   << ", "  << (*ID_it).second
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

      boundary_edge_iter it = _boundary_edge_id.begin();
      const boundary_edge_iter end = _boundary_edge_id.end();

      for (; it != end; ++it)
        ID_counts[(*it).second.second]++;

      std::map<boundary_id_type, std::size_t>::const_iterator ID_it        = ID_counts.begin();
      const std::map<boundary_id_type, std::size_t>::const_iterator ID_end = ID_counts.end();

      for (; ID_it != ID_end; ++ID_it)
        out_stream << "  (" << (*ID_it).first
                   << ", "  << (*ID_it).second
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

      boundary_shellface_iter it = _boundary_shellface_id.begin();
      const boundary_shellface_iter end = _boundary_shellface_id.end();

      for (; it != end; ++it)
        ID_counts[(*it).second.second]++;

      std::map<boundary_id_type, std::size_t>::const_iterator ID_it        = ID_counts.begin();
      const std::map<boundary_id_type, std::size_t>::const_iterator ID_end = ID_counts.end();

      for (; ID_it != ID_end; ++ID_it)
        out_stream << "  (" << (*ID_it).first
                   << ", "  << (*ID_it).second
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

      boundary_side_iter it = _boundary_side_id.begin();
      const boundary_side_iter end = _boundary_side_id.end();

      for (; it != end; ++it)
        ID_counts[(*it).second.second]++;

      std::map<boundary_id_type, std::size_t>::const_iterator ID_it        = ID_counts.begin();
      const std::map<boundary_id_type, std::size_t>::const_iterator ID_end = ID_counts.end();

      for (; ID_it != ID_end; ++ID_it)
        out_stream << "  (" << (*ID_it).first
                   << ", "  << (*ID_it).second
                   << ")"  << std::endl;
    }
}


const std::string & BoundaryInfo::get_sideset_name(boundary_id_type id) const
{
  static const std::string empty_string;
  std::map<boundary_id_type, std::string>::const_iterator it =
    _ss_id_to_name.find(id);
  if (it == _ss_id_to_name.end())
    return empty_string;
  else
    return it->second;
}


std::string & BoundaryInfo::sideset_name(boundary_id_type id)
{
  return _ss_id_to_name[id];
}

const std::string & BoundaryInfo::get_nodeset_name(boundary_id_type id) const
{
  static const std::string empty_string;
  std::map<boundary_id_type, std::string>::const_iterator it =
    _ns_id_to_name.find(id);
  if (it == _ns_id_to_name.end())
    return empty_string;
  else
    return it->second;
}

std::string & BoundaryInfo::nodeset_name(boundary_id_type id)
{
  return _ns_id_to_name[id];
}

boundary_id_type BoundaryInfo::get_id_by_name(const std::string & name) const
{
  // Search sidesets
  std::map<boundary_id_type, std::string>::const_iterator
    iter = _ss_id_to_name.begin(),
    end_iter = _ss_id_to_name.end();

  for (; iter != end_iter; ++iter)
    if (iter->second == name)
      return iter->first;

  // Search nodesets
  iter = _ns_id_to_name.begin();
  end_iter = _ns_id_to_name.end();
  for (; iter != end_iter; ++iter)
    if (iter->second == name)
      return iter->first;

  // If we made it here without returning, we don't have a sideset or
  // nodeset by the requested name, so return invalid_id
  return invalid_id;
}



void BoundaryInfo::_find_id_maps(const std::set<boundary_id_type> & requested_boundary_ids,
                                 dof_id_type first_free_node_id,
                                 std::map<dof_id_type, dof_id_type> * node_id_map,
                                 dof_id_type first_free_elem_id,
                                 std::map<std::pair<dof_id_type, unsigned char>, dof_id_type> * side_id_map)
{
  // We'll do the same modulus trick that DistributedMesh uses to avoid
  // id conflicts between different processors
  dof_id_type
    next_node_id = first_free_node_id + this->processor_id(),
    next_elem_id = first_free_elem_id + this->processor_id();

  // We'll pass through the mesh once first to build
  // the maps and count boundary nodes and elements.
  // To find local boundary nodes, we have to examine all elements
  // here rather than just local elements, because it's possible to
  // have a local boundary node that's not on a local boundary
  // element, e.g. at the tip of a triangle.

  // We'll loop through two different ranges here: first all elements,
  // looking for local nodes, and second through unpartitioned
  // elements, looking for all remaining nodes.
  const MeshBase::const_element_iterator end_el = _mesh.elements_end();
  bool hit_end_el = false;
  const MeshBase::const_element_iterator end_unpartitioned_el =
    _mesh.pid_elements_end(DofObject::invalid_processor_id);

  for (MeshBase::const_element_iterator el = _mesh.elements_begin();
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

      for (unsigned char s=0; s<elem->n_sides(); s++)
        if (elem->neighbor_ptr(s) == libmesh_nullptr) // on the boundary
          {
            // Get the top-level parent for this element
            const Elem * top_parent = elem->top_parent();

            // Find the right id number for that side
            std::pair<boundary_side_iter, boundary_side_iter> pos = _boundary_side_id.equal_range(top_parent);

            bool add_this_side = false;
            boundary_id_type this_bcid = invalid_id;

            for (; pos.first != pos.second; ++pos.first)
              {
                this_bcid = pos.first->second.second;

                // if this side is flagged with a boundary condition
                // and the user wants this id
                if ((pos.first->second.first == s) &&
                    (requested_boundary_ids.count(this_bcid)))
                  {
                    add_this_side = true;
                    break;
                  }
              }

            // if side s wasn't found or doesn't have a boundary
            // condition we may still want to add it
            if (pos.first == pos.second)
              {
                this_bcid = invalid_id;
                if (requested_boundary_ids.count(this_bcid))
                  add_this_side = true;
              }

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

                // Use a proxy element for the side to query nodes
                UniquePtr<const Elem> side (elem->build_side_ptr(s));
                for (unsigned int n = 0; n != side->n_nodes(); ++n)
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

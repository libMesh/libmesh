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
#include <iterator>  // std::distance

// Local includes
#include "libmesh/libmesh_config.h"
#include "libmesh/boundary_info.h"
#include "libmesh/elem.h"
#include "libmesh/mesh_data.h"
#include "libmesh/mesh_serializer.h"
#include "libmesh/parallel.h"
#include "libmesh/partitioner.h"
#include "libmesh/unstructured_mesh.h"

namespace libMesh
{



//------------------------------------------------------
// BoundaryInfo static member initializations
const boundary_id_type BoundaryInfo::invalid_id = -123;



//------------------------------------------------------
// BoundaryInfo functions
BoundaryInfo::BoundaryInfo(const MeshBase& m) :
  ParallelObject(m.comm()),
  _mesh (m)
{
}

BoundaryInfo& BoundaryInfo::operator=(const BoundaryInfo& other_boundary_info)
{
  /**
   * A quick note: We're going to attempt to pull _new_ pointers out of the mesh assigned to this boundary info.
   * This will only work if the mesh assigned to this BoundaryInfo is the same mesh object as other_boundary_info
   * _or_ was constructed in exactly the same way (or constructed as a copy).
   */

  // Copy node boundary info
  {
    std::multimap<const Node*, boundary_id_type>::const_iterator it = other_boundary_info._boundary_node_id.begin();
    const std::multimap<const Node*, boundary_id_type>::const_iterator end = other_boundary_info._boundary_node_id.end();

    for(; it != end; ++it)
    {
      const Node * other_node = it->first;
      _boundary_node_id.insert
        (std::pair<const Node*, boundary_id_type>
          (_mesh.node_ptr(other_node->id()), it->second) );
    }
  }

  // Copy edge boundary info
  {
    std::multimap<const Elem*, std::pair<unsigned short int, boundary_id_type> >::
      const_iterator it = other_boundary_info._boundary_edge_id.begin();
    const std::multimap<const Elem*, std::pair<unsigned short int, boundary_id_type> >::
      const_iterator end = other_boundary_info._boundary_edge_id.end();

    for(; it != end; ++it)
    {
      const Elem * other_elem = it->first;
      _boundary_edge_id.insert
        (std::pair<const Elem*, std::pair<unsigned short int, boundary_id_type> >
          (_mesh.elem(other_elem->id()), it->second) );
    }
  }

  // Copy side boundary info
  {
    std::multimap<const Elem*, std::pair<unsigned short int, boundary_id_type> >::
      const_iterator it = other_boundary_info._boundary_side_id.begin();
    const std::multimap<const Elem*, std::pair<unsigned short int, boundary_id_type> >::
      const_iterator end = other_boundary_info._boundary_side_id.end();

    for(; it != end; ++it)
    {
      const Elem * other_elem = it->first;
      _boundary_side_id.insert
        (std::pair<const Elem*, std::pair<unsigned short int, boundary_id_type> >
          (_mesh.elem(other_elem->id()), it->second) );
    }
  }

  _boundary_ids = other_boundary_info._boundary_ids;
  _side_boundary_ids = other_boundary_info._side_boundary_ids;
  _node_boundary_ids = other_boundary_info._node_boundary_ids;

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
  _boundary_ids.clear();
  _side_boundary_ids.clear();
  _node_boundary_ids.clear();
}



void BoundaryInfo::sync (UnstructuredMesh& boundary_mesh,
			 MeshData*     boundary_mesh_data,
			 MeshData*     this_mesh_data)
{
  std::set<boundary_id_type> request_boundary_ids(_boundary_ids);
  request_boundary_ids.insert(invalid_id);
  if (!_mesh.is_serial())
    this->comm().set_union(request_boundary_ids);

  this->sync(request_boundary_ids, boundary_mesh,
             boundary_mesh_data, this_mesh_data);
}



void BoundaryInfo::sync (const std::set<boundary_id_type> &requested_boundary_ids,
                         UnstructuredMesh& boundary_mesh,
			 MeshData*     boundary_mesh_data,
			 MeshData*     this_mesh_data)
{
  START_LOG("sync()", "BoundaryInfo");

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
  MeshSerializer(const_cast<MeshBase&>(_mesh), boundary_mesh.is_serial());

  /**
   * The boundary mesh elements will be one lower dimension than the
   * interior mesh elements
   */
  boundary_mesh.set_mesh_dimension(_mesh.mesh_dimension() - 1);

  /**
   * Re-create the boundary mesh.
   */

  boundary_mesh.set_n_partitions() = _mesh.n_partitions();

  std::map<dof_id_type, dof_id_type> node_id_map;
  std::map<std::pair<dof_id_type, unsigned char>, dof_id_type> side_id_map;

  // We'll do the same modulus trick that ParallelMesh uses to avoid
  // id conflicts between different processors
  dof_id_type next_node_id = this->processor_id(),
              next_elem_id = this->processor_id();

  // We'll pass through the mesh once first to build
  // the maps and count boundary nodes and elements
  // We have to examine all elements here rather than just local
  // elements, because it's possible to have a local boundary node
  // that's not on a local boundary element, e.g. at the tip of a
  // triangle.
  const MeshBase::const_element_iterator end_el = _mesh.elements_end();
  for (MeshBase::const_element_iterator el = _mesh.elements_begin();
       el != end_el; ++el)
    {
      const Elem *elem = *el;

      for (unsigned char s=0; s<elem->n_sides(); s++)
        if (elem->neighbor(s) == NULL) // on the boundary
          {
            // Get the top-level parent for this element
            const Elem* top_parent = elem->top_parent();

            // A convenient typedef
            typedef
              std::multimap<const Elem*, std::pair<unsigned short int, boundary_id_type> >::
                const_iterator Iter;

            // Find the right id number for that side
            std::pair<Iter, Iter> pos = _boundary_side_id.equal_range(top_parent);

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
                std::pair<dof_id_type, unsigned char> side_pair(elem->id(), s);
                libmesh_assert (!side_id_map.count(side_pair));
                side_id_map[side_pair] = next_elem_id;
                next_elem_id += this->n_processors() + 1;

                // Use a proxy element for the side to query nodes
                AutoPtr<Elem> side (elem->build_side(s));
                for (unsigned int n = 0; n != side->n_nodes(); ++n)
                  {
                    Node *node = side->get_node(n);
                    libmesh_assert(node);

                    // In parallel we only know enough to number our own nodes.
                    if (node->processor_id() != this->processor_id())
                      continue;

                    dof_id_type node_id = node->id();
                    if (!node_id_map.count(node_id))
                      {
                        node_id_map[node_id] = next_node_id;
                        next_node_id += this->n_processors() + 1;
                      }
                  }
              }
          }
    }

  // Join up the results from other processors
  this->comm().set_union(side_id_map);
  this->comm().set_union(node_id_map);

  // Finally we'll pass through any unpartitioned elements to add them
  // to the maps and counts.
  next_node_id = this->n_processors();
  next_elem_id = this->n_processors();

  const MeshBase::const_element_iterator end_unpartitioned_el =
    _mesh.pid_elements_end(DofObject::invalid_processor_id);
  for (MeshBase::const_element_iterator el =
         _mesh.pid_elements_begin(DofObject::invalid_processor_id);
       el != end_unpartitioned_el; ++el)
    {
      const Elem *elem = *el;

      for (unsigned char s=0; s<elem->n_sides(); s++)
        if (elem->neighbor(s) == NULL) // on the boundary
          {
            // Get the top-level parent for this element
            const Elem* top_parent = elem->top_parent();

            // A convenient typedef
            typedef
              std::multimap<const Elem*, std::pair<unsigned short int, boundary_id_type> >::
                const_iterator Iter;

            // Find the right id number for that side
            std::pair<Iter, Iter> pos = _boundary_side_id.equal_range(top_parent);

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

            // if side s doesn't have a boundary condition we may
            // still want to add it
            if (pos.first == pos.second)
              {
                this_bcid = invalid_id;
                if (requested_boundary_ids.count(this_bcid))
                  add_this_side = true;
              }

            if (add_this_side)
              {
                std::pair<dof_id_type, unsigned char> side_pair(elem->id(), s);
                libmesh_assert (!side_id_map.count(side_pair));
                side_id_map[side_pair] = next_elem_id;
                next_elem_id += this->n_processors() + 1;

                // Use a proxy element for the side to query nodes
                AutoPtr<Elem> side (elem->build_side(s));
                for (unsigned int n = 0; n != side->n_nodes(); ++n)
                  {
                    Node *node = side->get_node(n);
                    libmesh_assert(node);
                    dof_id_type node_id = node->id();
                    if (!node_id_map.count(node_id))
                      {
                        node_id_map[node_id] = next_node_id;
                        next_node_id += this->n_processors() + 1;
                      }
                  }
              }
          }
    }

  // FIXME: ought to renumber side/node_id_map image to be contiguous
  // to save memory, also ought to reserve memory

  // Let's add all the nodes to the boundary mesh

  MeshBase::const_node_iterator n_end  = _mesh.nodes_end();

  for(MeshBase::const_node_iterator n_it = _mesh.nodes_begin();
      n_it != n_end; ++n_it)
    {
      const Node* node = *n_it;
      dof_id_type node_id = node->id();
      if (node_id_map.count(node_id))
        boundary_mesh.add_point(*node, node_id_map[node_id], node->processor_id());
    }


  // Finally let's add the elements


  for (MeshBase::const_element_iterator el = _mesh.elements_begin();
       el != end_el; ++el)
    {
      const Elem* elem = *el;

      for (unsigned int s=0; s<elem->n_sides(); s++)
        if (elem->neighbor(s) == NULL) // on the boundary
          {
            // Get the top-level parent for this element
            const Elem* top_parent = elem->top_parent();

            // A convenient typedef
            typedef
              std::multimap<const Elem*, std::pair<unsigned short int, boundary_id_type> >::
                const_iterator Iter;

            // Find the right id number for that side
            std::pair<Iter, Iter> pos = _boundary_side_id.equal_range(top_parent);

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
                // Build the side - do not use a "proxy" element here:
                // This will be going into the boundary_mesh and needs to
                // stand on its own.
                AutoPtr<Elem> side (elem->build_side(s, false));

                side->processor_id() = elem->processor_id();

                const std::pair<dof_id_type, unsigned char> side_pair(elem->id(), s);

                libmesh_assert(side_id_map.count(side_pair));

                side->set_id(side_id_map[side_pair]);

                // Add the side
                Elem* new_elem = boundary_mesh.add_elem(side.release());

                // This side's Node pointers still point to the nodes of the original mesh.
                // We need to re-point them to the boundary mesh's nodes!  Since we copied *ALL* of
                // the original mesh's nodes over, we should be guaranteed to have the same ordering.
                for (unsigned int nn=0; nn<new_elem->n_nodes(); ++nn)
                  {
                    // Get the correct node pointer, based on the id()
                    Node* new_node = boundary_mesh.node_ptr(node_id_map[new_elem->node(nn)]);

                    // sanity check: be sure that the new Node exists
                    // and its global id really matches
                    libmesh_assert (new_node);
                    libmesh_assert_equal_to (new_node->id(), node_id_map[new_elem->node(nn)]);

                    // Assign the new node pointer
                    new_elem->set_node(nn) = new_node;
                  }

#ifdef LIBMESH_ENABLE_AMR
                // Finally, set the parent and interior_parent links
                if (elem->parent())
                  {
                    const std::pair<dof_id_type, unsigned char> parent_side_pair(elem->parent()->id(), s);

                    libmesh_assert(side_id_map.count(parent_side_pair));

                    Elem* side_parent = boundary_mesh.elem(side_id_map[parent_side_pair]);

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
                      if (new_elem->get_node(v) == side_parent->get_node(v))
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

                new_elem->set_interior_parent (const_cast<Elem*>(elem));
              }
          }
    }

  // When desired, copy the MeshData
  // to the boundary_mesh
  if ((boundary_mesh_data != NULL) && (this_mesh_data != NULL))
    boundary_mesh_data->assign(*this_mesh_data);

  // Don't repartition this mesh; we want it to stay in sync with the
  // interior partitioning.
  boundary_mesh.partitioner().reset(NULL);

  // Make boundary_mesh nodes and elements contiguous
  boundary_mesh.prepare_for_use(/*skip_renumber =*/ false);

  // and finally distribute element partitioning to the nodes
  Partitioner::set_node_processor_ids(boundary_mesh);

  STOP_LOG("sync()", "BoundaryInfo");
}



void BoundaryInfo::add_node(const dof_id_type node,
			    const boundary_id_type id)
{
  this->add_node (_mesh.node_ptr(node), id);
}



void BoundaryInfo::add_node(const Node* node,
			    const boundary_id_type id)
{
  if (id == invalid_id)
    {
      libMesh::err << "ERROR: You may not set a boundary ID of "
		    << invalid_id << std::endl
		    << " That is reserved for internal use.\n"
		    << std::endl;

      libmesh_error();
    }

  // A convenient typedef
  typedef std::multimap<const Node*, boundary_id_type>::const_iterator Iter;

  // Don't add the same ID twice
  std::pair<Iter, Iter> pos = _boundary_node_id.equal_range(node);

  for (;pos.first != pos.second; ++pos.first)
    if (pos.first->second == id)
      return;

  std::pair<const Node*, boundary_id_type> kv (node, id);

  _boundary_node_id.insert(kv);
  _boundary_ids.insert(id);
  _node_boundary_ids.insert(id); // Also add this ID to the set of node boundary IDs
}

void BoundaryInfo::add_node(const Node* node,
			    const std::vector<boundary_id_type>& ids)
{
  if (ids.empty())
    return;

  libmesh_assert(node);

  // A convenient typedef
  typedef std::multimap<const Node*, boundary_id_type>::const_iterator Iter;

  // Don't add the same ID twice
  std::pair<Iter, Iter> pos = _boundary_node_id.equal_range(node);

  for (unsigned int i=0; i!= ids.size(); ++i)
    {
      boundary_id_type id=ids[i];

      if (id == invalid_id)
        {
          libMesh::err << "ERROR: You may not set a boundary ID of "
		        << invalid_id << std::endl
		        << " That is reserved for internal use.\n"
		        << std::endl;

          libmesh_error();
        }

      bool already_inserted = false;
      for (Iter p = pos.first;p != pos.second; ++p)
        if (p->second == id)
          {
            already_inserted = true;
            break;
          }
      if (already_inserted)
        continue;

      std::pair<const Node*, boundary_id_type> kv (node, id);

      _boundary_node_id.insert(kv);
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
  this->add_edge (_mesh.elem(e), edge, id);
}



void BoundaryInfo::add_edge(const Elem* elem,
			    const unsigned short int edge,
			    const boundary_id_type id)
{
  libmesh_assert(elem);

  // Only add BCs for level-0 elements.
  libmesh_assert_equal_to (elem->level(), 0);

  if (id == invalid_id)
    {
      libMesh::err << "ERROR: You may not set a boundary ID of "
		    << invalid_id << std::endl
		    << " That is reserved for internal use.\n"
		    << std::endl;

      libmesh_error();
    }

  // A convenient typedef
  typedef std::multimap<const Elem*, std::pair<unsigned short int, boundary_id_type> >::
    const_iterator Iter;

  // Don't add the same ID twice
  std::pair<Iter, Iter> pos = _boundary_edge_id.equal_range(elem);

  for (;pos.first != pos.second; ++pos.first)
    if (pos.first->second.first == edge &&
        pos.first->second.second == id)
      return;

  std::pair<unsigned short int, boundary_id_type> p(edge,id);
  std::pair<const Elem*, std::pair<unsigned short int, boundary_id_type> >
    kv (elem, p);

  _boundary_edge_id.insert(kv);
  _boundary_ids.insert(id);
  _edge_boundary_ids.insert(id); // Also add this ID to the set of edge boundary IDs
}



void BoundaryInfo::add_edge(const Elem* elem,
			    const unsigned short int edge,
			    const std::vector<boundary_id_type>& ids)
{
  if (ids.empty())
    return;

  libmesh_assert(elem);

  // Only add BCs for level-0 elements.
  libmesh_assert_equal_to (elem->level(), 0);

  // A convenient typedef
  typedef std::multimap<const Elem*, std::pair<unsigned short int, boundary_id_type> >::
    const_iterator Iter;

  // Don't add the same ID twice
  std::pair<Iter, Iter> pos = _boundary_edge_id.equal_range(elem);

  for (unsigned int i=0; i!= ids.size(); ++i)
    {
      boundary_id_type id=ids[i];

      if (id == invalid_id)
        {
          libMesh::err << "ERROR: You may not set a boundary ID of "
		        << invalid_id << std::endl
		        << " That is reserved for internal use.\n"
		        << std::endl;

          libmesh_error();
        }

      bool already_inserted = false;
      for (Iter p = pos.first;p != pos.second; ++p)
        if (p->second.first == edge &&
            p->second.second == id)
          {
            already_inserted = true;
            break;
          }
      if (already_inserted)
        continue;

      std::pair<unsigned short int, boundary_id_type> p(edge,id);
      std::pair<const Elem*, std::pair<unsigned short int, boundary_id_type> >
        kv (elem, p);

      _boundary_edge_id.insert(kv);
      _boundary_ids.insert(id);
      _edge_boundary_ids.insert(id); // Also add this ID to the set of edge boundary IDs
    }
}

void BoundaryInfo::add_side(const dof_id_type e,
			    const unsigned short int side,
			    const boundary_id_type id)
{
  this->add_side (_mesh.elem(e), side, id);
}



void BoundaryInfo::add_side(const Elem* elem,
			    const unsigned short int side,
			    const boundary_id_type id)
{
  libmesh_assert(elem);

  // Only add BCs for level-0 elements.
  libmesh_assert_equal_to (elem->level(), 0);

  if (id == invalid_id)
    {
      libMesh::err << "ERROR: You may not set a boundary ID of "
		    << invalid_id << std::endl
		    << " That is reserved for internal use.\n"
		    << std::endl;

      libmesh_error();
    }

  // A convenient typedef
  typedef std::multimap<const Elem*, std::pair<unsigned short int, boundary_id_type> >::
    const_iterator Iter;

  // Don't add the same ID twice
  std::pair<Iter, Iter> pos = _boundary_side_id.equal_range(elem);

  for (;pos.first != pos.second; ++pos.first)
    if (pos.first->second.first == side &&
        pos.first->second.second == id)
      return;

  std::pair<unsigned short int, boundary_id_type> p(side,id);
  std::pair<const Elem*, std::pair<unsigned short int, boundary_id_type> >
    kv (elem, p);

  _boundary_side_id.insert(kv);
  _boundary_ids.insert(id);
  _side_boundary_ids.insert(id); // Also add this ID to the set of side boundary IDs
}



void BoundaryInfo::add_side(const Elem* elem,
			    const unsigned short int side,
			    const std::vector<boundary_id_type>& ids)
{
  if (ids.empty())
    return;

  libmesh_assert(elem);

  // Only add BCs for level-0 elements.
  libmesh_assert_equal_to (elem->level(), 0);

  // A convenient typedef
  typedef std::multimap<const Elem*, std::pair<unsigned short int, boundary_id_type> >::
    const_iterator Iter;

  // Don't add the same ID twice
  std::pair<Iter, Iter> pos = _boundary_side_id.equal_range(elem);

  for (unsigned int i=0; i!= ids.size(); ++i)
    {
      boundary_id_type id=ids[i];

      if (id == invalid_id)
        {
          libMesh::err << "ERROR: You may not set a boundary ID of "
		        << invalid_id << std::endl
		        << " That is reserved for internal use.\n"
		        << std::endl;

          libmesh_error();
        }

      bool already_inserted = false;
      for (Iter p = pos.first;p != pos.second; ++p)
        if (p->second.first == side &&
            p->second.second == id)
          {
            already_inserted = true;
            break;
          }
      if (already_inserted)
        continue;

      std::pair<unsigned short int, boundary_id_type> p(side,id);
      std::pair<const Elem*, std::pair<unsigned short int, boundary_id_type> >
        kv (elem, p);

      _boundary_side_id.insert(kv);
      _boundary_ids.insert(id);
      _side_boundary_ids.insert(id); // Also add this ID to the set of side boundary IDs
    }
}



bool BoundaryInfo::has_boundary_id(const Node* const node,
                                   const boundary_id_type id) const
{
  // A convenient typedef
  typedef std::multimap<const Node*, boundary_id_type>::const_iterator Iter;

  std::pair<Iter, Iter> pos = _boundary_node_id.equal_range(node);

  for (;pos.first != pos.second; ++pos.first)
    if (pos.first->second == id)
      return true;

  return false;
}



std::vector<boundary_id_type> BoundaryInfo::boundary_ids(const Node* node) const
{
  std::vector<boundary_id_type> ids;

  // A convenient typedef
  typedef std::multimap<const Node*, boundary_id_type>::const_iterator Iter;

  std::pair<Iter, Iter> pos = _boundary_node_id.equal_range(node);

  for (;pos.first != pos.second; ++pos.first)
    ids.push_back(pos.first->second);

  return ids;
}



unsigned int BoundaryInfo::n_boundary_ids(const Node* node) const
{
  // A convenient typedef
  typedef std::multimap<const Node*, boundary_id_type>::const_iterator Iter;

  std::pair<Iter, Iter> pos = _boundary_node_id.equal_range(node);

  return libmesh_cast_int<unsigned int>
    (std::distance(pos.first, pos.second));
}



std::vector<boundary_id_type> BoundaryInfo::edge_boundary_ids (const Elem* const elem,
                                                               const unsigned short int edge) const
{
  libmesh_assert(elem);

  std::vector<boundary_id_type> ids;

  // Only level-0 elements store BCs.  If this is not a level-0
  // element get its level-0 parent and infer the BCs.
  const Elem* searched_elem = elem;
  if (elem->level() != 0)
  {
    // Find all the sides that contain edge. If one of those is a boundary
    // side, then this must be a boundary edge. In that case, we just use the
    // top-level parent.
    bool found_boundary_edge = false;
    for(unsigned int side=0; side<elem->n_sides(); side++)
    {
      if(elem->is_edge_on_side(edge,side))
      {
        if (elem->neighbor(side) == NULL)
        {
          searched_elem = elem->top_parent ();
          found_boundary_edge = true;
          break;
        }
      }
    }
  
#ifdef LIBMESH_ENABLE_AMR
    if(!found_boundary_edge)
    {
      // Child element is not on external edge, but it may have internal
      // "boundary" IDs.  We will walk up the tree, at each level checking that
      // the current child is actually on the same edge of the parent that is
      // currently being searched for (i.e. that was passed in as "edge").
      while (searched_elem->parent() != NULL)
      {
        const Elem * parent = searched_elem->parent();
        if (parent->is_child_on_edge(parent->which_child_am_i(searched_elem), edge) == false)
          return ids;
        searched_elem = parent;
      }
    }
#endif
  }

  std::pair<std::multimap<const Elem*,
                          std::pair<unsigned short int, boundary_id_type> >::const_iterator,
            std::multimap<const Elem*,
                          std::pair<unsigned short int, boundary_id_type> >::const_iterator >
    e = _boundary_edge_id.equal_range(searched_elem);

  // elem not in the data structure
  if (e.first == e.second)
    return ids;

  // elem is there, maybe multiple occurrences
  for (; e.first != e.second; ++e.first)
    // if this is true we found the requested edge of the element
    if (e.first->second.first == edge)
      ids.push_back(e.first->second.second);

  // Whether or not we found anything, return "ids".  If it's empty, it
  // means no valid bounary IDs were found for "edge"
  return ids;
}



unsigned int BoundaryInfo::n_edge_boundary_ids (const Elem* const elem,
                                                const unsigned short int edge) const
{
  libmesh_assert(elem);
  
  // Only level-0 elements store BCs.  If this is not a level-0
  // element get its level-0 parent and infer the BCs.
  const Elem* searched_elem = elem;
  if (elem->level() != 0)
  {
    // Find all the sides that contain edge. If one of those is a boundary
    // side, then this must be a boundary edge. In that case, we just use the
    // top-level parent.
    bool found_boundary_edge = false;
    for(unsigned int side=0; side<elem->n_sides(); side++)
    {
      if(elem->is_edge_on_side(edge,side))
      {
        if (elem->neighbor(side) == NULL)
        {
          searched_elem = elem->top_parent ();
          found_boundary_edge = true;
          break;
        }
      }
    }
  
#ifdef LIBMESH_ENABLE_AMR
    if(!found_boundary_edge)
    {
      // Child element is not on external edge, but it may have internal
      // "boundary" IDs.  We will walk up the tree, at each level checking that
      // the current child is actually on the same edge of the parent that is
      // currently being searched for (i.e. that was passed in as "edge").
      while (searched_elem->parent() != NULL)
      {
        const Elem * parent = searched_elem->parent();
        if (parent->is_child_on_edge(parent->which_child_am_i(searched_elem), edge) == false)
          return 0;
        searched_elem = parent;
      }
    }
#endif
  }

  std::pair<std::multimap<const Elem*,
                          std::pair<unsigned short int, boundary_id_type> >::const_iterator,
            std::multimap<const Elem*,
                          std::pair<unsigned short int, boundary_id_type> >::const_iterator >
    e = _boundary_edge_id.equal_range(searched_elem);

  unsigned int n_ids = 0;

  // elem is there, maybe multiple occurrences
  for (; e.first != e.second; ++e.first)
    // if this is true we found the requested edge of the element
    if (e.first->second.first == edge)
      n_ids++;

  return n_ids;
}



std::vector<boundary_id_type> BoundaryInfo::raw_edge_boundary_ids (const Elem* const elem,
                                                                   const unsigned short int edge) const
{
  libmesh_assert(elem);

  std::vector<boundary_id_type> ids;

  // Only level-0 elements store BCs.
  if (elem->parent())
    return ids;

  std::pair<std::multimap<const Elem*,
                          std::pair<unsigned short int, boundary_id_type> >::const_iterator,
            std::multimap<const Elem*,
                          std::pair<unsigned short int, boundary_id_type> >::const_iterator >
    e = _boundary_edge_id.equal_range(elem);

  // Check any occurrences
  for (; e.first != e.second; ++e.first)
    // if this is true we found the requested edge of the element
    if (e.first->second.first == edge)
      ids.push_back(e.first->second.second);

  // if nothing got pushed back, we didn't find elem in the data
  // structure with the requested edge, so return the default empty
  // vector
  return ids;
}


boundary_id_type BoundaryInfo::boundary_id(const Elem* const elem,
				           const unsigned short int side) const
{
  // Asking for just one boundary id means your code isn't safe to use
  // on meshes with overlapping boundary ids.  Try using
  // BoundaryInfo::boundary_ids or BoundaryInfo::has_boundary_id
  // instead.
  libmesh_deprecated();

  libmesh_assert(elem);

  // Only level-0 elements store BCs.  If this is not a level-0
  // element, one of its parent elements may have either internal
  // or external boundary IDs.  We find that parent now.
  const Elem*  searched_elem = elem;
  if (elem->level() != 0)
  {
    // Child element on external side: the top_parent will have the BCs
    if (elem->neighbor(side) == NULL)
      searched_elem = elem->top_parent ();

#ifdef LIBMESH_ENABLE_AMR
    // Child element is not on external side, but it may have internal
    // "boundary" IDs.  We will walk up the tree, at each level checking that
    // the current child is actually on the same side of the parent that is
    // currently being searched for (i.e. that was passed in as "side").
    else
      while (searched_elem->parent() != NULL)
	{
	  const Elem * parent = searched_elem->parent();
	  if (parent->is_child_on_side(parent->which_child_am_i(searched_elem), side) == false)
	    return invalid_id;
	  searched_elem = parent;
	}
#endif
  }

  std::pair<std::multimap<const Elem*,
                          std::pair<unsigned short int, boundary_id_type> >::const_iterator,
            std::multimap<const Elem*,
                          std::pair<unsigned short int, boundary_id_type> >::const_iterator >
    e = _boundary_side_id.equal_range(searched_elem);

  // elem not in the data structure
  if (e.first == e.second)
    return invalid_id;

  // elem is there, maybe multiple occurrences
  for (; e.first != e.second; ++e.first)
      // if this is true we found the requested side
      // of the element and want to return the id
      if (e.first->second.first == side)
	return e.first->second.second;

  // if we get here, we found elem in the data structure but not
  // the requested side, so return the default value
  return invalid_id;
}



bool BoundaryInfo::has_boundary_id(const Elem* const elem,
                                   const unsigned short int side,
                                   const boundary_id_type id) const
{
  libmesh_assert(elem);

  // Only level-0 elements store BCs.  If this is not a level-0
  // element get its level-0 parent and infer the BCs.
  const Elem*  searched_elem = elem;
  if (elem->level() != 0)
  {
    if (elem->neighbor(side) == NULL)
      searched_elem = elem->top_parent ();
#ifdef LIBMESH_ENABLE_AMR
    else
      while (searched_elem->parent() != NULL)
	{
	  const Elem * parent = searched_elem->parent();
	  if (parent->is_child_on_side(parent->which_child_am_i(searched_elem), side) == false)
	    return false;
	  searched_elem = parent;
	}
#endif
  }

  std::pair<std::multimap<const Elem*,
                          std::pair<unsigned short int, boundary_id_type> >::const_iterator,
            std::multimap<const Elem*,
                          std::pair<unsigned short int, boundary_id_type> >::const_iterator >
    e = _boundary_side_id.equal_range(searched_elem);

  // elem is there, maybe multiple occurrences
  for (; e.first != e.second; ++e.first)
    // if this is true we found the requested id on this side of the element
    if (e.first->second.first == side &&
        e.first->second.second == id)
      return true;

  return false;
}



std::vector<boundary_id_type> BoundaryInfo::boundary_ids (const Elem* const elem,
                                                          const unsigned short int side) const
{
  libmesh_assert(elem);

  std::vector<boundary_id_type> ids;

  // Only level-0 elements store BCs.  If this is not a level-0
  // element get its level-0 parent and infer the BCs.
  const Elem*  searched_elem = elem;
  if (elem->level() != 0)
  {
    if (elem->neighbor(side) == NULL)
      searched_elem = elem->top_parent ();
#ifdef LIBMESH_ENABLE_AMR
    else
      while (searched_elem->parent() != NULL)
	{
	  const Elem * parent = searched_elem->parent();
	  if (parent->is_child_on_side(parent->which_child_am_i(searched_elem), side) == false)
	    return ids;
	  searched_elem = parent;
	}
#endif
  }

  std::pair<std::multimap<const Elem*,
                          std::pair<unsigned short int, boundary_id_type> >::const_iterator,
            std::multimap<const Elem*,
                          std::pair<unsigned short int, boundary_id_type> >::const_iterator >
    e = _boundary_side_id.equal_range(searched_elem);

  // elem not in the data structure
  if (e.first == e.second)
    return ids;

  // elem is there, maybe multiple occurrences
  for (; e.first != e.second; ++e.first)
    // if this is true we found the requested side of the element
    if (e.first->second.first == side)
      ids.push_back(e.first->second.second);

  // Whether or not we found anything, return "ids".  If it's empty, it
  // means no valid bounary IDs were found for "side"
  return ids;
}



unsigned int BoundaryInfo::n_boundary_ids (const Elem* const elem,
                                           const unsigned short int side) const
{
  libmesh_assert(elem);

  // Only level-0 elements store BCs.  If this is not a level-0
  // element get its level-0 parent and infer the BCs.
  const Elem*  searched_elem = elem;
  if (elem->level() != 0)
  {
    if (elem->neighbor(side) == NULL)
      searched_elem = elem->top_parent ();
#ifdef LIBMESH_ENABLE_AMR
    else
      while (searched_elem->parent() != NULL)
	{
	  const Elem * parent = searched_elem->parent();
	  if (parent->is_child_on_side(parent->which_child_am_i(searched_elem), side) == false)
	    return 0;
	  searched_elem = parent;
	}
#endif
  }

  std::pair<std::multimap<const Elem*,
                          std::pair<unsigned short int, boundary_id_type> >::const_iterator,
            std::multimap<const Elem*,
                          std::pair<unsigned short int, boundary_id_type> >::const_iterator >
    e = _boundary_side_id.equal_range(searched_elem);

  unsigned int n_ids = 0;

  // elem is there, maybe multiple occurrences
  for (; e.first != e.second; ++e.first)
    // if this is true we found the requested side of the element
    if (e.first->second.first == side)
      n_ids++;

  return n_ids;
}



std::vector<boundary_id_type> BoundaryInfo::raw_boundary_ids (const Elem* const elem,
                                                              const unsigned short int side) const
{
  libmesh_assert(elem);

  std::vector<boundary_id_type> ids;

  // Only level-0 elements store BCs.
  if (elem->parent())
    return ids;

  std::pair<std::multimap<const Elem*,
                          std::pair<unsigned short int, boundary_id_type> >::const_iterator,
            std::multimap<const Elem*,
                          std::pair<unsigned short int, boundary_id_type> >::const_iterator >
    e = _boundary_side_id.equal_range(elem);

  // Check any occurrences
  for (; e.first != e.second; ++e.first)
    // if this is true we found the requested side of the element
    if (e.first->second.first == side)
      ids.push_back(e.first->second.second);

  // if nothing got pushed back, we didn't find elem in the data
  // structure with the requested side, so return the default empty
  // vector
  return ids;
}


void BoundaryInfo::remove_edge (const Elem* elem,
                                const unsigned short int edge)
{
  libmesh_assert(elem);

  // The user shouldn't be trying to remove only one child's boundary
  // id
  libmesh_assert_equal_to (elem->level(), 0);

  std::pair<std::multimap<const Elem*,
                          std::pair<unsigned short int, boundary_id_type> >::iterator,
            std::multimap<const Elem*,
                          std::pair<unsigned short int, boundary_id_type> >::iterator >
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



void BoundaryInfo::remove_edge (const Elem* elem,
                                const unsigned short int edge,
                                const boundary_id_type id)
{
  libmesh_assert(elem);

  // The user shouldn't be trying to remove only one child's boundary
  // id
  libmesh_assert_equal_to (elem->level(), 0);

  std::pair<std::multimap<const Elem*,
                          std::pair<unsigned short int, boundary_id_type> >::iterator,
            std::multimap<const Elem*,
                          std::pair<unsigned short int, boundary_id_type> >::iterator >
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

void BoundaryInfo::remove_side (const Elem* elem,
                                const unsigned short int side)
{
  libmesh_assert(elem);

  // The user shouldn't be trying to remove only one child's boundary
  // id
  libmesh_assert_equal_to (elem->level(), 0);

  std::pair<std::multimap<const Elem*,
                          std::pair<unsigned short int, boundary_id_type> >::iterator,
            std::multimap<const Elem*,
                          std::pair<unsigned short int, boundary_id_type> >::iterator >
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



void BoundaryInfo::remove_side (const Elem* elem,
                                const unsigned short int side,
                                const boundary_id_type id)
{
  libmesh_assert(elem);

  // The user shouldn't be trying to remove only one child's boundary
  // id
  libmesh_assert_equal_to (elem->level(), 0);

  std::pair<std::multimap<const Elem*,
                          std::pair<unsigned short int, boundary_id_type> >::iterator,
            std::multimap<const Elem*,
                          std::pair<unsigned short int, boundary_id_type> >::iterator >
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



unsigned int BoundaryInfo::side_with_boundary_id(const Elem* const elem,
                                                 const boundary_id_type boundary_id_in) const
{
  const Elem* searched_elem = elem;
  if (elem->level() != 0)
    searched_elem = elem->top_parent();

  std::pair<std::multimap<const Elem*,
                          std::pair<unsigned short int, boundary_id_type> >::const_iterator,
            std::multimap<const Elem*,
                          std::pair<unsigned short int, boundary_id_type> >::const_iterator >
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
         if (elem->neighbor(side) == NULL)
           return side;

         // If we're on an internal boundary then we need to be sure
         // it's the same internal boundary as our top_parent
         const Elem *p = elem;

#ifdef LIBMESH_ENABLE_AMR

         while (p != NULL)
           {
             const Elem *parent = p->parent();
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

void BoundaryInfo::build_node_boundary_ids(std::vector<boundary_id_type> &b_ids)
{
  b_ids.clear();

  std::multimap<const Node*, boundary_id_type>::const_iterator pos
    = _boundary_node_id.begin();

  for (; pos != _boundary_node_id.end(); ++pos)
    {
      boundary_id_type id = pos->second;

      if(std::find(b_ids.begin(),b_ids.end(),id) == b_ids.end())
        b_ids.push_back(id);
    }
}

void BoundaryInfo::build_side_boundary_ids(std::vector<boundary_id_type> &b_ids)
{
  b_ids.clear();

  std::multimap<const Elem*, std::pair<unsigned short int, boundary_id_type> >::const_iterator pos
    = _boundary_side_id.begin();

  for (; pos != _boundary_side_id.end(); ++pos)
    {
      boundary_id_type id = pos->second.second;

      if(std::find(b_ids.begin(),b_ids.end(),id) == b_ids.end())
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

  std::multimap<const Elem*,
                std::pair<unsigned short int,
                          boundary_id_type> >::const_iterator pos;

  for (pos=_boundary_side_id.begin(); pos != _boundary_side_id.end(); ++pos)
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

  std::multimap<const Elem*,
                std::pair<unsigned short int,
                          boundary_id_type> >::const_iterator pos;

  for (pos=_boundary_edge_id.begin(); pos != _boundary_edge_id.end(); ++pos)
    if (pos->first->processor_id() == this->processor_id())
      n_edge_bcs++;

  this->comm().sum (n_edge_bcs);

  return n_edge_bcs;
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

  std::multimap<const Node*, boundary_id_type>::const_iterator pos;

  for (pos=_boundary_node_id.begin(); pos != _boundary_node_id.end(); ++pos)
    if (pos->first->processor_id() == this->processor_id())
      n_nodesets++;

  this->comm().sum (n_nodesets);

  return n_nodesets;
}



void BoundaryInfo::build_node_list (std::vector<dof_id_type>& nl,
				    std::vector<boundary_id_type>&    il) const
{
  // Reserve the size, then use push_back
  nl.reserve (_boundary_node_id.size());
  il.reserve (_boundary_node_id.size());

  std::multimap<const Node*, boundary_id_type>::const_iterator pos
    = _boundary_node_id.begin();

  for (; pos != _boundary_node_id.end(); ++pos)
    {
      nl.push_back (pos->first->id());
      il.push_back (pos->second);
    }
}


void
BoundaryInfo::build_node_list_from_side_list()
{
  std::multimap<const Elem*,
                std::pair<unsigned short int,
                          boundary_id_type> >::const_iterator pos;

  //Loop over the side list
  for (pos=_boundary_side_id.begin(); pos != _boundary_side_id.end(); ++pos)
  {
    // Don't add remote sides
    if(pos->first->is_remote())
      continue;

    //Need to loop over the sides of any possible children
    std::vector< const Elem * > family;
#ifdef LIBMESH_ENABLE_AMR
    pos->first->active_family_tree_by_side (family, pos->second.first);
#else
    family.push_back(pos->first);
#endif

    for(std::size_t elem_it=0; elem_it < family.size(); elem_it++)
    {
      const Elem * cur_elem = family[elem_it];

      AutoPtr<Elem> side = cur_elem->build_side(pos->second.first);

      //Add each node node on the side with the side's boundary id
      for(unsigned int i=0; i<side->n_nodes(); i++)
      {
        Node * node = side->get_node(i);

        this->add_node(node, pos->second.second);
      }
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

  // typedef for less typing!
  typedef std::multimap<const Node*, boundary_id_type>::const_iterator iterator_t;

  // Return value and iterator for equal_range()
  iterator_t pos;
  std::pair<iterator_t, iterator_t> range;

  MeshBase::const_element_iterator el = _mesh.active_elements_begin();
  const MeshBase::const_element_iterator end_el = _mesh.active_elements_end();

  for (; el != end_el; ++el)
    {
      const Elem* elem = *el;

      for (unsigned side=0; side<elem->n_sides(); ++side)
	{
	  AutoPtr<Elem> side_elem = elem->build_side(side);

	  // map from nodeset_id to count for that ID
	  std::map<dof_id_type, unsigned> nodesets_node_count;
	  for (unsigned node_num=0; node_num < side_elem->n_nodes(); ++node_num)
	    {
	      Node* node = side_elem->get_node(node_num);
	      range = _boundary_node_id.equal_range(node);

	      // For each nodeset that this node is a member of, increment the associated
	      // nodeset ID count
	      for (pos = range.first; pos != range.second; ++pos)
		{
		  nodesets_node_count[pos->second]++;
		}
	    }

	  // Now check to see what nodeset_counts have the correct number of nodes in them
	  for (std::map<dof_id_type, unsigned>::const_iterator nodesets = nodesets_node_count.begin();
	       nodesets != nodesets_node_count.end(); ++nodesets)
	    {
	      if (nodesets->second == side_elem->n_nodes())
		{
		  // Add this side to the sideset
		  add_side(elem, side, nodesets->first);
		}
	    }
	} // end for side
    } // end for el
}




void BoundaryInfo::build_side_list (std::vector<dof_id_type>&        el,
				    std::vector<unsigned short int>& sl,
				    std::vector<boundary_id_type>&   il) const
{
  // Reserve the size, then use push_back
  el.reserve (_boundary_side_id.size());
  sl.reserve (_boundary_side_id.size());
  il.reserve (_boundary_side_id.size());

  std::multimap<const Elem*,
                std::pair<unsigned short int,
                          boundary_id_type> >::const_iterator pos;

  for (pos=_boundary_side_id.begin(); pos != _boundary_side_id.end();
       ++pos)
    {
      el.push_back (pos->first->id());
      sl.push_back (pos->second.first);
      il.push_back (pos->second.second);
    }
}



void BoundaryInfo::print_info(std::ostream& out_stream) const
{
  // Print out the nodal BCs
  if (!_boundary_node_id.empty())
    {
      out_stream << "Nodal Boundary conditions:" << std::endl
	  << "--------------------------" << std::endl
	  << "  (Node No., ID)               " << std::endl;

//       std::for_each(_boundary_node_id.begin(),
//		    _boundary_node_id.end(),
//		    PrintNodeInfo());

      std::multimap<const Node*, boundary_id_type>::const_iterator it        = _boundary_node_id.begin();
      const std::multimap<const Node*, boundary_id_type>::const_iterator end = _boundary_node_id.end();

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

//       std::for_each(_boundary_edge_id.begin(),
//		    _boundary_edge_id.end(),
//		    PrintSideInfo());

      std::multimap<const Elem*,
	std::pair<unsigned short int, boundary_id_type> >::const_iterator it = _boundary_edge_id.begin();
      const std::multimap<const Elem*,
	std::pair<unsigned short int, boundary_id_type> >::const_iterator end = _boundary_edge_id.end();

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

//       std::for_each(_boundary_side_id.begin(),
//		    _boundary_side_id.end(),
//		    PrintSideInfo());

      std::multimap<const Elem*,
	std::pair<unsigned short int, boundary_id_type> >::const_iterator it = _boundary_side_id.begin();
      const std::multimap<const Elem*,
	std::pair<unsigned short int, boundary_id_type> >::const_iterator end = _boundary_side_id.end();

      for (; it != end; ++it)
        out_stream << "  (" << (*it).first->id()
	    << ", "  << (*it).second.first
	    << ", "  << (*it).second.second
	    << ")"   << std::endl;
    }
}



void BoundaryInfo::print_summary(std::ostream& out_stream) const
{
  // Print out the nodal BCs
  if (!_boundary_node_id.empty())
    {
      out_stream << "Nodal Boundary conditions:" << std::endl
	  << "--------------------------" << std::endl
	  << "  (ID, number of nodes)   " << std::endl;

      std::map<boundary_id_type, std::size_t> ID_counts;

      std::multimap<const Node*, boundary_id_type>::const_iterator it        = _boundary_node_id.begin();
      const std::multimap<const Node*, boundary_id_type>::const_iterator end = _boundary_node_id.end();

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

      std::multimap<const Elem*,
	std::pair<unsigned short int, boundary_id_type> >::const_iterator it = _boundary_edge_id.begin();
      const std::multimap<const Elem*,
	std::pair<unsigned short int, boundary_id_type> >::const_iterator end = _boundary_edge_id.end();

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

      std::multimap<const Elem*,
	std::pair<unsigned short int, boundary_id_type> >::const_iterator it = _boundary_side_id.begin();
      const std::multimap<const Elem*,
	std::pair<unsigned short int, boundary_id_type> >::const_iterator end = _boundary_side_id.end();

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

std::string& BoundaryInfo::sideset_name(boundary_id_type id)
{
  return _ss_id_to_name[id];
}

std::string& BoundaryInfo::nodeset_name(boundary_id_type id)
{
  return _ns_id_to_name[id];
}

boundary_id_type BoundaryInfo::get_id_by_name(const std::string& name) const
{
  // This function is searching the keys of the map
  // We might want to make this more efficient
  std::map<boundary_id_type, std::string>::const_iterator iter = _ss_id_to_name.begin();
  std::map<boundary_id_type, std::string>::const_iterator end_iter = _ss_id_to_name.end();

  for ( ; iter != end_iter; ++iter)
  {
    if (iter->second == name)
      return iter->first;
  }

  // Loop over nodesets
  iter = _ns_id_to_name.begin();
  end_iter = _ns_id_to_name.end();
  for ( ; iter != end_iter; ++iter)
  {
    if (iter->second == name)
      return iter->first;
  }

  libMesh::err << "The sideset/nodeset named " << name << " does not exist in mesh!" << std::endl;
  libmesh_error();
}

} // namespace libMesh

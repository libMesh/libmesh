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



// Local includes
#include "libmesh/boundary_info.h"
#include "libmesh/elem.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/metis_partitioner.h"
#include "libmesh/serial_mesh.h"

#include LIBMESH_INCLUDE_UNORDERED_SET
LIBMESH_DEFINE_HASH_POINTERS

namespace libMesh
{

// ------------------------------------------------------------
// SerialMesh class member functions
SerialMesh::SerialMesh (unsigned int d,
			const Parallel::Communicator &comm) :
  UnstructuredMesh (d,comm)
{
  _partitioner = AutoPtr<Partitioner>(new MetisPartitioner());
}


SerialMesh::~SerialMesh ()
{
  this->clear();  // Free nodes and elements
}


// This might be specialized later, but right now it's just here to
// make sure the compiler doesn't give us a default (non-deep) copy
// constructor instead.
SerialMesh::SerialMesh (const SerialMesh &other_mesh) :
  UnstructuredMesh (other_mesh)
{
  this->copy_nodes_and_elements(other_mesh);
  *this->boundary_info = *other_mesh.boundary_info;
}


SerialMesh::SerialMesh (const UnstructuredMesh &other_mesh) :
  UnstructuredMesh (other_mesh)
{
  this->copy_nodes_and_elements(other_mesh);
  *this->boundary_info = *other_mesh.boundary_info;
}


const Point& SerialMesh::point (const dof_id_type i) const
{
  libmesh_assert_less (i, this->n_nodes());
  libmesh_assert(_nodes[i]);
  libmesh_assert_equal_to (_nodes[i]->id(), i); // This will change soon

  return (*_nodes[i]);
}





const Node& SerialMesh::node (const dof_id_type i) const
{
  libmesh_assert_less (i, this->n_nodes());
  libmesh_assert(_nodes[i]);
  libmesh_assert_equal_to (_nodes[i]->id(), i); // This will change soon

  return (*_nodes[i]);
}





Node& SerialMesh::node (const dof_id_type i)
{
  if (i >= this->n_nodes())
    {
      libMesh::out << " i=" << i
		    << ", n_nodes()=" << this->n_nodes()
		    << std::endl;
      libmesh_error();
    }

  libmesh_assert_less (i, this->n_nodes());
  libmesh_assert(_nodes[i]);
  libmesh_assert_equal_to (_nodes[i]->id(), i); // This will change soon

  return (*_nodes[i]);
}



const Node* SerialMesh::node_ptr (const dof_id_type i) const
{
  libmesh_assert_less (i, this->n_nodes());
  libmesh_assert(_nodes[i]);
  libmesh_assert_equal_to (_nodes[i]->id(), i); // This will change soon

  return _nodes[i];
}




Node* SerialMesh::node_ptr (const dof_id_type i)
{
  libmesh_assert_less (i, this->n_nodes());
  libmesh_assert(_nodes[i]);
  libmesh_assert_equal_to (_nodes[i]->id(), i); // This will change soon

  return _nodes[i];
}




const Node* SerialMesh::query_node_ptr (const dof_id_type i) const
{
  if (i >= this->n_nodes())
    return NULL;
  libmesh_assert (_nodes[i] == NULL ||
                  _nodes[i]->id() == i); // This will change soon

  return _nodes[i];
}




Node* SerialMesh::query_node_ptr (const dof_id_type i)
{
  if (i >= this->n_nodes())
    return NULL;
  libmesh_assert (_nodes[i] == NULL ||
                  _nodes[i]->id() == i); // This will change soon

  return _nodes[i];
}




const Elem* SerialMesh::elem (const dof_id_type i) const
{
  libmesh_assert_less (i, this->n_elem());
  libmesh_assert(_elements[i]);
  libmesh_assert_equal_to (_elements[i]->id(), i); // This will change soon

  return _elements[i];
}




Elem* SerialMesh::elem (const dof_id_type i)
{
  libmesh_assert_less (i, this->n_elem());
  libmesh_assert(_elements[i]);
  libmesh_assert_equal_to (_elements[i]->id(), i); // This will change soon

  return _elements[i];
}




const Elem* SerialMesh::query_elem (const dof_id_type i) const
{
  if (i >= this->n_elem())
    return NULL;
  libmesh_assert (_elements[i] == NULL ||
                  _elements[i]->id() == i); // This will change soon

  return _elements[i];
}




Elem* SerialMesh::query_elem (const dof_id_type i)
{
  if (i >= this->n_elem())
    return NULL;
  libmesh_assert (_elements[i] == NULL ||
                  _elements[i]->id() == i); // This will change soon

  return _elements[i];
}




Elem* SerialMesh::add_elem (Elem* e)
{
  libmesh_assert(e);

  // We no longer merely append elements with SerialMesh

  // If the user requests a valid id that doesn't correspond to an
  // existing element, let's give them that id, resizing the elements
  // container if necessary.
  if (!e->valid_id())
    e->set_id (_elements.size());

  const dof_id_type id = e->id();

  if (id < _elements.size())
    {
      // Overwriting existing elements is still probably a mistake.
      libmesh_assert(!_elements[id]);
    }
  else
    {
      _elements.resize(id+1, NULL);
    }

  _elements[id] = e;

  return e;
}



Elem* SerialMesh::insert_elem (Elem* e)
{
  dof_id_type eid = e->id();
  libmesh_assert_less (eid, _elements.size());
  Elem *oldelem = _elements[eid];

  if (oldelem)
    {
      libmesh_assert_equal_to (oldelem->id(), eid);
      this->delete_elem(oldelem);
    }

  _elements[e->id()] = e;

  return e;
}



void SerialMesh::delete_elem(Elem* e)
{
  libmesh_assert(e);

  // Initialize an iterator to eventually point to the element we want to delete
  std::vector<Elem*>::iterator pos = _elements.end();

  // In many cases, e->id() gives us a clue as to where e
  // is located in the _elements vector.  Try that first
  // before trying the O(n_elem) search.
  libmesh_assert_less (e->id(), _elements.size());

  if (_elements[e->id()] == e)
    {
      // We found it!
      pos = _elements.begin();
      std::advance(pos, e->id());
    }

  else
    {
      // This search is O(n_elem)
      pos = std::find (_elements.begin(),
		       _elements.end(),
		       e);
    }

  // Huh? Element not in the vector?
  libmesh_assert (pos != _elements.end());

  // Remove the element from the BoundaryInfo object
  this->boundary_info->remove(e);

  // delete the element
  delete e;

  // explicitly NULL the pointer
  *pos = NULL;
}



void SerialMesh::renumber_elem(const dof_id_type old_id,
                               const dof_id_type new_id)
{
  // This doesn't get used in serial yet
  Elem *el = _elements[old_id];
  libmesh_assert (el);

  el->set_id(new_id);
  libmesh_assert (!_elements[new_id]);
  _elements[new_id] = el;
  _elements[old_id] = NULL;
}



Node* SerialMesh::add_point (const Point& p,
			     const dof_id_type id,
			     const processor_id_type proc_id)
{
//   // We only append points with SerialMesh
//   libmesh_assert(id == DofObject::invalid_id || id == _nodes.size());
//   Node *n = Node::build(p, _nodes.size()).release();
//   n->processor_id() = proc_id;
//   _nodes.push_back (n);

  Node *n = NULL;

  // If the user requests a valid id, either
  // provide the existing node or resize the container
  // to fit the new node.
  if (id != DofObject::invalid_id)
    if (id < _nodes.size())
      n = _nodes[id];
    else
      _nodes.resize(id+1);
  else
    _nodes.push_back (static_cast<Node*>(NULL));

  // if the node already exists, then assign new (x,y,z) values
  if (n)
    *n = p;
  // otherwise build a new node, put it in the right spot, and return
  // a valid pointer.
  else
    {
      n = Node::build(p, (id == DofObject::invalid_id) ? _nodes.size()-1 : id).release();
      n->processor_id() = proc_id;

      if (id == DofObject::invalid_id)
	_nodes.back() = n;
      else
	_nodes[id] = n;
    }

  // better not pass back a NULL pointer.
  libmesh_assert (n);

  return n;
}



Node* SerialMesh::add_node (Node* n)
{
  libmesh_assert(n);
  // We only append points with SerialMesh
  libmesh_assert(!n->valid_id() || n->id() == _nodes.size());

  n->set_id (_nodes.size());

  _nodes.push_back(n);

  return n;
}



void SerialMesh::delete_node(Node* n)
{
  libmesh_assert(n);
  libmesh_assert_less (n->id(), _nodes.size());

  // Initialize an iterator to eventually point to the element we want
  // to delete
  std::vector<Node*>::iterator pos;

  // In many cases, e->id() gives us a clue as to where e
  // is located in the _elements vector.  Try that first
  // before trying the O(n_elem) search.
  if (_nodes[n->id()] == n)
    {
      pos = _nodes.begin();
      std::advance(pos, n->id());
    }
  else
    {
      pos = std::find (_nodes.begin(),
		       _nodes.end(),
		       n);
    }

  // Huh? Node not in the vector?
  libmesh_assert (pos != _nodes.end());

  // Delete the node from the BoundaryInfo object
  this->boundary_info->remove(n);

  // delete the node
  delete n;

  // explicitly NULL the pointer
  *pos = NULL;
}



void SerialMesh::renumber_node(const dof_id_type old_id,
                               const dof_id_type new_id)
{
  // This doesn't get used in serial yet
  Node *nd = _nodes[old_id];
  libmesh_assert (nd);

  nd->set_id(new_id);
  libmesh_assert (!_nodes[new_id]);
  _nodes[new_id] = nd;
  _nodes[old_id] = NULL;
}



void SerialMesh::clear ()
{
  // Call parent clear function
  MeshBase::clear();


  // Clear our elements and nodes
  {
    std::vector<Elem*>::iterator       it  = _elements.begin();
    const std::vector<Elem*>::iterator end = _elements.end();

    // There is no need to remove the elements from
    // the BoundaryInfo data structure since we
    // already cleared it.
    for (; it != end; ++it)
      delete *it;

    _elements.clear();
  }

  // clear the nodes data structure
  {
    std::vector<Node*>::iterator       it  = _nodes.begin();
    const std::vector<Node*>::iterator end = _nodes.end();

    // There is no need to remove the nodes from
    // the BoundaryInfo data structure since we
    // already cleared it.
    for (; it != end; ++it)
      delete *it;

    _nodes.clear();
  }
}



void SerialMesh::renumber_nodes_and_elements ()
{

  START_LOG("renumber_nodes_and_elem()", "Mesh");

  // node and element id counters
  dof_id_type next_free_elem = 0;
  dof_id_type next_free_node = 0;

  // Will hold the set of nodes that are currently connected to elements
  LIBMESH_BEST_UNORDERED_SET<Node*> connected_nodes;

  // Loop over the elements.  Note that there may
  // be NULLs in the _elements vector from the coarsening
  // process.  Pack the elements in to a contiguous array
  // and then trim any excess.
  {
    std::vector<Elem*>::iterator in        = _elements.begin();
    std::vector<Elem*>::iterator out_iter  = _elements.begin();
    const std::vector<Elem*>::iterator end = _elements.end();

    for (; in != end; ++in)
      if (*in != NULL)
	{
	  Elem* el = *in;

	  *out_iter = *in;
	  ++out_iter;

	  // Increment the element counter
	  el->set_id (next_free_elem++);

          if(_skip_renumber_nodes_and_elements)
          {
            // Add this elements nodes to the connected list
            for (unsigned int n=0; n<el->n_nodes(); n++)
              connected_nodes.insert(el->get_node(n));
          }
          else  // We DO want node renumbering
          {
            // Loop over this element's nodes.  Number them,
            // if they have not been numbered already.  Also,
            // position them in the _nodes vector so that they
            // are packed contiguously from the beginning.
            for (unsigned int n=0; n<el->n_nodes(); n++)
              if (el->node(n) == next_free_node)     // don't need to process
                next_free_node++;                      // [(src == dst) below]

              else if (el->node(n) > next_free_node) // need to process
	      {
		// The source and destination indices
		// for this node
		const dof_id_type src_idx = el->node(n);
		const dof_id_type dst_idx = next_free_node++;

		// ensure we want to swap a valid nodes
		libmesh_assert(_nodes[src_idx]);

		// Swap the source and destination nodes
                std::swap(_nodes[src_idx],
                          _nodes[dst_idx] );

		// Set proper indices where that makes sense
		if (_nodes[src_idx] != NULL)
		  _nodes[src_idx]->set_id (src_idx);
		_nodes[dst_idx]->set_id (dst_idx);
	      }
          }
	}

    // Erase any additional storage. These elements have been
    // copied into NULL voids by the procedure above, and are
    // thus repeated and unnecessary.
    _elements.erase (out_iter, end);
  }


  if(_skip_renumber_nodes_and_elements)
  {
    // Loop over the nodes.  Note that there may
    // be NULLs in the _nodes vector from the coarsening
    // process.  Pack the nodes in to a contiguous array
    // and then trim any excess.

    std::vector<Node*>::iterator in        = _nodes.begin();
    std::vector<Node*>::iterator out_iter  = _nodes.begin();
    const std::vector<Node*>::iterator end = _nodes.end();

    for (; in != end; ++in)
      if (*in != NULL)
      {
        // This is a reference so that if we change the pointer it will change in the vector
        Node* & nd = *in;

        // If this node is still connected to an elem, put it in the list
        if(connected_nodes.find(nd) != connected_nodes.end())
        {
          *out_iter = nd;
          ++out_iter;

          // Increment the node counter
          nd->set_id (next_free_node++);
        }
        else // This node is orphaned, delete it!
        {
          this->boundary_info->remove (nd);

          // delete the node
          delete nd;
          nd = NULL;
        }
      }

    // Erase any additional storage.  Whatever was
    _nodes.erase (out_iter, end);
  }
  else // We really DO want node renumbering
  {
    // Any nodes in the vector >= _nodes[next_free_node]
    // are not connected to any elements and may be deleted
    // if desired.

    // (This code block will erase the unused nodes)
    // Now, delete the unused nodes
    {
      std::vector<Node*>::iterator nd        = _nodes.begin();
      const std::vector<Node*>::iterator end = _nodes.end();

      std::advance (nd, next_free_node);

      for (std::vector<Node*>::iterator it=nd;
           it != end; ++it)
      {
        // Mesh modification code might have already deleted some
        // nodes
	if (*it == NULL)
          continue;

	// remove any boundary information associated with
	// this node
	this->boundary_info->remove (*it);

	// delete the node
	delete *it;
	*it = NULL;
      }

      _nodes.erase (nd, end);
    }
  }

  libmesh_assert_equal_to (next_free_elem, _elements.size());
  libmesh_assert_equal_to (next_free_node, _nodes.size());

  STOP_LOG("renumber_nodes_and_elem()", "Mesh");
}



void SerialMesh::fix_broken_node_and_element_numbering ()
{
   // Nodes first
  for (dof_id_type n=0; n<this->_nodes.size(); n++)
    if (this->_nodes[n] != NULL)
      this->_nodes[n]->set_id() = n;

  // Elements next
  for (dof_id_type e=0; e<this->_elements.size(); e++)
    if (this->_elements[e] != NULL)
      this->_elements[e]->set_id() = e;
}


void SerialMesh::stitch_meshes (SerialMesh& other_mesh,
                                boundary_id_type this_mesh_boundary_id,
                                boundary_id_type other_mesh_boundary_id,
                                Real tol,
                                bool clear_stitched_boundary_ids,
                                bool verbose)
{
  std::map<dof_id_type, dof_id_type> node_to_node_map;
  std::map<dof_id_type, std::vector<dof_id_type> > node_to_elems_map;

  if( (this_mesh_boundary_id  != BoundaryInfo::invalid_id) &&
      (other_mesh_boundary_id != BoundaryInfo::invalid_id) )
  {
    std::set<dof_id_type> this_boundary_node_ids;
    MeshBase::element_iterator elem_it  = this->elements_begin();
    MeshBase::element_iterator elem_end = this->elements_end();
    for ( ; elem_it != elem_end; ++elem_it)
    {
      Elem *el = *elem_it;

      // Now check whether elem has a face on the specified boundary
      for (unsigned int side_id=0; side_id<el->n_sides(); side_id++)
        if (el->neighbor(side_id) == NULL)
        {
          boundary_id_type bc_id = this->boundary_info->boundary_id (el, side_id);

          if(bc_id == this_mesh_boundary_id)
          {
            AutoPtr<Elem> side (el->build_side(side_id));
            for (unsigned int node_id=0; node_id<side->n_nodes(); node_id++)
            {
              this_boundary_node_ids.insert( side->node(node_id) );
            }
          }
        }
    }

    std::set<dof_id_type> other_boundary_node_ids;
    elem_it  = other_mesh.elements_begin();
    elem_end = other_mesh.elements_end();
    for ( ; elem_it != elem_end; ++elem_it)
    {
      Elem *el = *elem_it;

      // Now check whether elem has a face on the specified boundary
      for (unsigned int side_id=0; side_id<el->n_sides(); side_id++)
        if (el->neighbor(side_id) == NULL)
        {
          boundary_id_type bc_id = other_mesh.boundary_info->boundary_id (el, side_id);

          if(bc_id == other_mesh_boundary_id)
          {
            AutoPtr<Elem> side (el->build_side(side_id));
            for (unsigned int node_id=0; node_id<side->n_nodes(); node_id++)
            {
              other_boundary_node_ids.insert( side->node(node_id) );
            }
          }
        }
    }

    std::set<dof_id_type>::iterator set_it     = this_boundary_node_ids.begin();
    std::set<dof_id_type>::iterator set_it_end = this_boundary_node_ids.end();
    for( ; set_it != set_it_end; ++set_it)
    {
      dof_id_type this_node_id = *set_it;
      Node& this_node = this->node(this_node_id);

      bool found_matching_nodes = false;

      std::set<dof_id_type>::iterator other_set_it     = other_boundary_node_ids.begin();
      std::set<dof_id_type>::iterator other_set_it_end = other_boundary_node_ids.end();
      for( ; other_set_it != other_set_it_end; ++other_set_it)
      {
        dof_id_type other_node_id = *other_set_it;
        Node& other_node = other_mesh.node(other_node_id);

        Real node_distance = (this_node - other_node).size();

        if(node_distance < tol)
        {
          // Make sure we didn't already find a matching node!
          if(found_matching_nodes)
          {
            libMesh::out << "Error: Found multiple matching nodes in stitch_meshes" << std::endl;
            libmesh_error();
          }

          node_to_node_map[this_node_id] = other_node_id;

          // Build a vector of all the elements in other_mesh that contain other_node
          std::vector<dof_id_type> other_elem_ids;
          MeshBase::element_iterator other_elem_it  = other_mesh.elements_begin();
          MeshBase::element_iterator other_elem_end = other_mesh.elements_end();
          for (; other_elem_it != other_elem_end; ++other_elem_it)
          {
            Elem *el = *other_elem_it;

            if(el->contains_point(other_node))
              other_elem_ids.push_back(el->id());
          }

          node_to_elems_map[this_node_id] = other_elem_ids;

          found_matching_nodes = true;
        }
      }
    }

    if(verbose)
    {
      libMesh::out << "In SerialMesh::stitch_meshes:" << std::endl
                   << "This mesh has " << this_boundary_node_ids.size() << " nodes on specified boundary" << std::endl
                   << "Other mesh has " << other_boundary_node_ids.size() << " nodes on specified boundary" << std::endl
                   << "Found " << node_to_node_map.size() << " matching nodes." << std::endl << std::endl;
    }
  }
  else
  {
    if(verbose)
    {
      libMesh::out << "Skip node merging in SerialMesh::stitch_meshes:" << std::endl;
    }
  }



  dof_id_type node_delta = this->n_nodes();
  dof_id_type elem_delta = this->n_elem();

  // need to increment node and element IDs of other_mesh before copying to this mesh
  MeshBase::node_iterator node_it  = other_mesh.nodes_begin();
  MeshBase::node_iterator node_end = other_mesh.nodes_end();
  for (; node_it != node_end; ++node_it)
  {
    Node *nd = *node_it;
    dof_id_type new_id = nd->id() + node_delta;
    nd->set_id(new_id);
  }

  MeshBase::element_iterator elem_it  = other_mesh.elements_begin();
  MeshBase::element_iterator elem_end = other_mesh.elements_end();
  for (; elem_it != elem_end; ++elem_it)
  {
    Elem *el = *elem_it;
    dof_id_type new_id = el->id() + elem_delta;
    el->set_id(new_id);
  }

  // Also, increment the node_to_node_map and node_to_elems_map
  std::map<dof_id_type, dof_id_type>::iterator node_map_it     = node_to_node_map.begin();
  std::map<dof_id_type, dof_id_type>::iterator node_map_it_end = node_to_node_map.end();
  for( ; node_map_it != node_map_it_end; ++node_map_it)
  {
    node_map_it->second += node_delta;
  }
  std::map<dof_id_type, std::vector<dof_id_type> >::iterator elem_map_it     = node_to_elems_map.begin();
  std::map<dof_id_type, std::vector<dof_id_type> >::iterator elem_map_it_end = node_to_elems_map.end();
  for( ; elem_map_it != elem_map_it_end; ++elem_map_it)
  {
    dof_id_type n_elems = elem_map_it->second.size();
    for(dof_id_type i=0; i<n_elems; i++)
    {
      (elem_map_it->second)[i] += elem_delta;
    }
  }

  // Copy mesh data
  this->copy_nodes_and_elements(other_mesh);

  // then decrement node and element IDs of mesh_i to return to original state
  node_it  = other_mesh.nodes_begin();
  node_end = other_mesh.nodes_end();
  for (; node_it != node_end; ++node_it)
  {
    Node *nd = *node_it;
    dof_id_type new_id = nd->id() - node_delta;
    nd->set_id(new_id);
  }

  elem_it  = other_mesh.elements_begin();
  elem_end = other_mesh.elements_end();
  for (; elem_it != elem_end; ++elem_it)
  {
    Elem *el = *elem_it;

    // First copy boundary info to the stitched mesh
    for (unsigned int side_id=0; side_id<el->n_sides(); side_id++)
      if (el->neighbor(side_id) == NULL)
      {
        boundary_id_type bc_id = other_mesh.boundary_info->boundary_id (el, side_id);

        if(bc_id != BoundaryInfo::invalid_id)
        {
          this->boundary_info->add_side(el->id(), side_id, bc_id);
        }
      }

    // Then decrement
    dof_id_type new_id = el->id() - elem_delta;
    el->set_id(new_id);
  }

  // Finally, we need to "merge" the overlapping nodes
  // We do this by iterating over node_to_elems_map and updating
  // the elements so that they "point" to the nodes that came
  // from this mesh, rather than from other_mesh.
  // Then we iterate over node_to_node_map and delete the
  // duplicate nodes that came from other_mesh.
  elem_map_it     = node_to_elems_map.begin();
  elem_map_it_end = node_to_elems_map.end();
  for( ; elem_map_it != elem_map_it_end; ++elem_map_it)
  {
    dof_id_type target_node_id = elem_map_it->first;
    dof_id_type other_node_id = node_to_node_map[target_node_id];
    Node& target_node = this->node(target_node_id);

    dof_id_type n_elems = elem_map_it->second.size();
    for(unsigned int i=0; i<n_elems; i++)
    {
      dof_id_type elem_id = elem_map_it->second[i];
      Elem* el = this->elem(elem_id);

      // find the local node index that we want to update
      unsigned int local_node_index = el->local_node(other_node_id);

      el->set_node(local_node_index) = &target_node;
    }
  }

  node_map_it     = node_to_node_map.begin();
  node_map_it_end = node_to_node_map.end();
  for( ; node_map_it != node_map_it_end; ++node_map_it)
  {
    dof_id_type node_id = node_map_it->second;
    this->delete_node( this->node_ptr(node_id) );
  }

  this->prepare_for_use( /*skip_renumber_nodes_and_elements= */ false);

  // After the stitching, we may want to clear boundary IDs from element
  // faces that are now internal to the mesh
  if(clear_stitched_boundary_ids)
  {
    elem_it  = this->elements_begin();
    elem_end = this->elements_end();
    for (; elem_it != elem_end; ++elem_it)
    {
      Elem *el = *elem_it;

      for (unsigned int side_id=0; side_id<el->n_sides(); side_id++)
      {
        if (el->neighbor(side_id) != NULL)
        {
          boundary_id_type bc_id = this->boundary_info->boundary_id (el, side_id);

          if( (bc_id == this_mesh_boundary_id)  ||
              (bc_id == other_mesh_boundary_id) )
          {
            this->boundary_info->remove_side(el, side_id);
          }
        }
      }
    }
  }

}


dof_id_type SerialMesh::n_active_elem () const
{
  return static_cast<dof_id_type>(std::distance (this->active_elements_begin(),
						 this->active_elements_end()));
}

} // namespace libMesh

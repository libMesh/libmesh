// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2008 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "boundary_info.h"
#include "elem.h"
#include "libmesh_logging.h"
#include "metis_partitioner.h"
#include "serial_mesh.h"

namespace libMesh
{

// ------------------------------------------------------------
// SerialMesh class member functions
SerialMesh::SerialMesh (unsigned int d) :
  UnstructuredMesh (d)
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


const Point& SerialMesh::point (const unsigned int i) const
{
  libmesh_assert (i < this->n_nodes());
  libmesh_assert (_nodes[i] != NULL);
  libmesh_assert (_nodes[i]->id() == i); // This will change soon

  return (*_nodes[i]);
}





const Node& SerialMesh::node (const unsigned int i) const
{
  libmesh_assert (i < this->n_nodes());
  libmesh_assert (_nodes[i] != NULL);
  libmesh_assert (_nodes[i]->id() == i); // This will change soon

  return (*_nodes[i]);
}





Node& SerialMesh::node (const unsigned int i)
{
  if (i >= this->n_nodes())
    {
      libMesh::out << " i=" << i
		    << ", n_nodes()=" << this->n_nodes()
		    << std::endl;
      libmesh_error();
    }

  libmesh_assert (i < this->n_nodes());
  libmesh_assert (_nodes[i] != NULL);
  libmesh_assert (_nodes[i]->id() == i); // This will change soon

  return (*_nodes[i]);
}



const Node* SerialMesh::node_ptr (const unsigned int i) const
{
  libmesh_assert (i < this->n_nodes());
  libmesh_assert (_nodes[i] != NULL);
  libmesh_assert (_nodes[i]->id() == i); // This will change soon

  return _nodes[i];
}




Node* & SerialMesh::node_ptr (const unsigned int i)
{
  libmesh_assert (i < this->n_nodes());
  libmesh_assert (_nodes[i] != NULL);
  libmesh_assert (_nodes[i]->id() == i); // This will change soon

  return _nodes[i];
}




Elem* SerialMesh::elem (const unsigned int i) const
{
  libmesh_assert (i < this->n_elem());
  libmesh_assert (_elements[i] != NULL);
  libmesh_assert (_elements[i]->id() == i); // This will change soon

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

  const unsigned int id = e->id();

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
  unsigned int eid = e->id();
  libmesh_assert(eid < _elements.size());
  Elem *oldelem = _elements[eid];

  if (oldelem)
    {
      libmesh_assert(oldelem->id() == eid);
      this->delete_elem(oldelem);
    }

  _elements[e->id()] = e;

  return e;
}



void SerialMesh::delete_elem(Elem* e)
{
  libmesh_assert (e != NULL);

  // Initialize an iterator to eventually point to the element we want to delete
  std::vector<Elem*>::iterator pos = _elements.end();

  // In many cases, e->id() gives us a clue as to where e
  // is located in the _elements vector.  Try that first
  // before trying the O(n_elem) search.
  libmesh_assert (e->id() < _elements.size());

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



void SerialMesh::renumber_elem(const unsigned int old_id,
                               const unsigned int new_id)
{
  // This doesn't get used in serial yet
  Elem *elem = _elements[old_id];
  libmesh_assert (elem);

  elem->set_id(new_id);
  libmesh_assert (!_elements[new_id]);
  _elements[new_id] = elem;
  _elements[old_id] = NULL;
}



Node* SerialMesh::add_point (const Point& p,
			     const unsigned int id,
			     const unsigned int proc_id)
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
  libmesh_assert (n != NULL);
  libmesh_assert (n->id() < _nodes.size());

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



void SerialMesh::renumber_node(const unsigned int old_id,
                               const unsigned int new_id)
{
  // This doesn't get used in serial yet
  Node *node = _nodes[old_id];
  libmesh_assert (node);

  node->set_id(new_id);
  libmesh_assert (!_nodes[new_id]);
  _nodes[new_id] = node;
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
  unsigned int next_free_elem = 0;
  unsigned int next_free_node = 0;

  // Loop over the elements.  Note that there may
  // be NULLs in the _elements vector from the coarsening
  // process.  Pack the elements in to a contiguous array
  // and then trim any excess.
  {
    std::vector<Elem*>::iterator in        = _elements.begin();
    std::vector<Elem*>::iterator out       = _elements.begin();
    const std::vector<Elem*>::iterator end = _elements.end();

    for (; in != end; ++in)
      if (*in != NULL)
	{
	  Elem* elem = *in;

	  *out = *in;
	  ++out;

	  // Increment the element counter
	  elem->set_id (next_free_elem++);

	  // Loop over this element's nodes.  Number them,
	  // if they have not been numbered already.  Also,
	  // position them in the _nodes vector so that they
	  // are packed contiguously from the beginning.
	  for (unsigned int n=0; n<elem->n_nodes(); n++)
	    if (elem->node(n) == next_free_node)     // don't need to process
	      next_free_node++;                      // [(src == dst) below]

	    else if (elem->node(n) > next_free_node) // need to process
	      {
		// The source and destination indices
		// for this node
		const unsigned int src_idx = elem->node(n);
		const unsigned int dst_idx = next_free_node++;

		// ensure we want to swap a valid nodes
		libmesh_assert (_nodes[src_idx] != NULL);

		// Swap the source and destination nodes
                std::swap(_nodes[src_idx],
                          _nodes[dst_idx] );

		// Set proper indices where that makes sense
		if (_nodes[src_idx] != NULL)
		  _nodes[src_idx]->set_id (src_idx);
		_nodes[dst_idx]->set_id (dst_idx);
	      }
	}

    // Erase any additional storage. These elements have been
    // copied into NULL voids by the procedure above, and are
    // thus repeated and unnecessary.
    _elements.erase (out, end);
  }

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


  libmesh_assert (next_free_elem == _elements.size());
  libmesh_assert (next_free_node == _nodes.size());

  STOP_LOG("renumber_nodes_and_elem()", "Mesh");
}



void SerialMesh::fix_broken_node_and_element_numbering ()
{
   // Nodes first
  for (unsigned int n=0; n<this->_nodes.size(); n++)
    if (this->_nodes[n] != NULL)
      this->_nodes[n]->set_id() = n;

  // Elements next
  for (unsigned int e=0; e<this->_elements.size(); e++)
    if (this->_elements[e] != NULL)
      this->_elements[e]->set_id() = e;
}



unsigned int SerialMesh::n_active_elem () const
{
  return static_cast<unsigned int>(std::distance (this->active_elements_begin(),
						  this->active_elements_end()));
}

} // namespace libMesh

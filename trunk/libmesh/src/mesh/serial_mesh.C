// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2007  Benjamin S. Kirk, John W. Peterson
  
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
#include "serial_mesh.h"

// ------------------------------------------------------------
// SerialMesh class member functions
SerialMesh::SerialMesh (unsigned int d) :
  UnstructuredMesh (d)
{
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
}


SerialMesh::SerialMesh (const UnstructuredMesh &other_mesh) :
  UnstructuredMesh (other_mesh)
{
  this->copy_nodes_and_elements(other_mesh);
}


const Point& SerialMesh::point (const unsigned int i) const
{
  assert (i < this->n_nodes());
  assert (_nodes[i] != NULL);
  assert (_nodes[i]->id() == i); // This will change soon

  return (*_nodes[i]);
}





const Node& SerialMesh::node (const unsigned int i) const
{
  assert (i < this->n_nodes());
  assert (_nodes[i] != NULL);
  assert (_nodes[i]->id() == i); // This will change soon
  
  return (*_nodes[i]);
}





Node& SerialMesh::node (const unsigned int i)
{
  if (i >= this->n_nodes())
    {
      std::cout << " i=" << i
		<< ", n_nodes()=" << this->n_nodes()
		<< std::endl;
      error();
    }
  
  assert (i < this->n_nodes());
  assert (_nodes[i] != NULL);
  assert (_nodes[i]->id() == i); // This will change soon

  return (*_nodes[i]);
}



const Node* SerialMesh::node_ptr (const unsigned int i) const
{
  assert (i < this->n_nodes());
  assert (_nodes[i] != NULL);
  assert (_nodes[i]->id() == i); // This will change soon
  
  return _nodes[i];
}




Node* & SerialMesh::node_ptr (const unsigned int i)
{
  assert (i < this->n_nodes());
  assert (_nodes[i] != NULL);
  assert (_nodes[i]->id() == i); // This will change soon

  return _nodes[i];
}




Elem* SerialMesh::elem (const unsigned int i) const
{
  assert (i < this->n_elem());
  assert (_elements[i] != NULL);
  assert (_elements[i]->id() == i); // This will change soon
  
  return _elements[i];
}




Elem* SerialMesh::add_elem (Elem* e)
{
  if (e != NULL)
    e->set_id (_elements.size());
  
  _elements.push_back(e);

  return e;
}



Elem* SerialMesh::insert_elem (Elem* e)
{
  unsigned int eid = e->id();
  assert(eid < _elements.size());
  Elem *oldelem = _elements[eid];

  if (oldelem)
    {
      assert(oldelem->id() == eid);
      this->delete_elem(oldelem);
    }

  _elements[e->id()] = e;

  return e;
}



void SerialMesh::delete_elem(Elem* e)
{
  assert (e != NULL);

  // Initialize an iterator to eventually point to the element we want to delete
  std::vector<Elem*>::iterator pos = _elements.end();
  
  // In many cases, e->id() gives us a clue as to where e
  // is located in the _elements vector.  Try that first
  // before trying the O(n_elem) search.
  assert (e->id() < _elements.size());

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
  assert (pos != _elements.end());

  // Remove the element from the BoundaryInfo object
  this->boundary_info->remove(e);
  
  // delete the element
  delete e;
  
  // explicitly NULL the pointer
  *pos = NULL;
}



Node* SerialMesh::add_point (const Point& p)
{  
  _nodes.push_back (Node::build(p, this->n_nodes()).release());
  
  return _nodes.back();
}



void SerialMesh::delete_node(Node* n)
{
  assert (n != NULL);
  assert (n->id() < _nodes.size());

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
  assert (pos != _nodes.end());

  // Delete the node from the BoundaryInfo object
  this->boundary_info->remove(n);
  
  // delete the node
  delete n;
  
  // explicitly NULL the pointer
  *pos = NULL;
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

		// ensure we want to swap valid nodes
		assert (_nodes[src_idx] != NULL);
		assert (_nodes[dst_idx] != NULL);
		
		// Swap the source and destination nodes
                std::swap(_nodes[src_idx],
                          _nodes[dst_idx] );

		// Set proper indices
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
	assert (*it != NULL);

	// remove any boundary information associated with
	// this node
	this->boundary_info->remove (*it);
	
	// delete the node
	delete *it;
	*it = NULL;
      }

    _nodes.erase (nd, end);
  }
  

  assert (next_free_elem == _elements.size());
  assert (next_free_node == _nodes.size());
  
  STOP_LOG("renumber_nodes_and_elem()", "Mesh");
}

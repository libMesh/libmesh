// $Id: parallel_mesh.C,v 1.9 2007-10-22 23:06:31 roystgnr Exp $

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
#include "parallel_mesh.h"

// ------------------------------------------------------------
// ParallelMesh class member functions
ParallelMesh::ParallelMesh (unsigned int d) :
  UnstructuredMesh (d)
{
}


ParallelMesh::~ParallelMesh ()
{
  this->clear();  // Free nodes and elements
}


// This might be specialized later, but right now it's just here to
// make sure the compiler doesn't give us a default (non-deep) copy
// constructor instead.
ParallelMesh::ParallelMesh (const ParallelMesh &other_mesh) :
  UnstructuredMesh (other_mesh)
{
  this->copy_nodes_and_elements(other_mesh);
}



ParallelMesh::ParallelMesh (const UnstructuredMesh &other_mesh) :
  UnstructuredMesh (other_mesh)
{
  this->copy_nodes_and_elements(other_mesh);
}



const Point& ParallelMesh::point (const unsigned int i) const
{
  assert (i < this->max_node_id());
  assert (_nodes[i] != NULL);
  assert (_nodes[i]->id() != Node::invalid_id);  

  return (*_nodes[i]);
}





const Node& ParallelMesh::node (const unsigned int i) const
{
  assert (i < this->max_node_id());
  assert (_nodes[i] != NULL);
  assert (_nodes[i]->id() != Node::invalid_id);  
  
  return (*_nodes[i]);
}





Node& ParallelMesh::node (const unsigned int i)
{
  assert (i < this->max_node_id());
  assert (_nodes[i] != NULL);
  assert (_nodes[i]->id() != Node::invalid_id);  

  return (*_nodes[i]);
}



const Node* ParallelMesh::node_ptr (const unsigned int i) const
{
  assert (i < this->max_node_id());
  assert (_nodes[i] != NULL);
  assert (_nodes[i]->id() != Node::invalid_id);  
  
  return _nodes[i];
}




Node* & ParallelMesh::node_ptr (const unsigned int i)
{
  assert (i < this->max_node_id());

  return _nodes[i];
}




Elem* ParallelMesh::elem (const unsigned int i) const
{
  assert (i < this->max_elem_id());
  assert (_elements[i] != NULL);
  
  return _elements[i];
}




Elem* ParallelMesh::add_elem (Elem* e)
{
  if (e != NULL)
    e->set_id (this->max_elem_id());
  
  _elements[this->max_elem_id()] = e;

  return e;
}



Elem* ParallelMesh::insert_elem (Elem* e)
{
  if (_elements[e->id()])
    this->delete_elem(_elements[e->id()]);

  _elements[e->id()] = e;

  return e;
}



void ParallelMesh::delete_elem(Elem* e)
{
  assert (e != NULL);

  // Delete the element from the BoundaryInfo object
  this->boundary_info->remove(e);

  // But not yet from the container; we might invalidate
  // an iterator that way!

  //_elements.erase(e->id());

  // Instead, we set it to NULL for now

  _elements[e->id()] = NULL;
  
  // delete the element
  delete e;
}



Node* ParallelMesh::add_point (const Point& p)
{  
  Node* n = Node::build(p, this->max_node_id()).release();
  _nodes[this->max_node_id()] = n;
  
  return n;
}



void ParallelMesh::delete_node(Node* n)
{
  assert (n != NULL);
  assert (n->id() < this->max_node_id());

  // Delete the node from the BoundaryInfo object
  this->boundary_info->remove(n);

  // And from the container
  _nodes.erase(n->id());
  
  // delete the node
  delete n;
}



void ParallelMesh::clear ()
{
  // Call parent clear function
  MeshBase::clear();

  
  // Clear our elements and nodes
  {
    elem_iterator_imp        it = _elements.begin();
    const elem_iterator_imp end = _elements.end();

    // There is no need to remove the elements from
    // the BoundaryInfo data structure since we
    // already cleared it.
    for (; it != end; ++it)
      delete *it;

    _elements.clear();
  }

  // clear the nodes data structure
  {
    node_iterator_imp it  = _nodes.begin();
    node_iterator_imp end = _nodes.end();

    // There is no need to remove the nodes from
    // the BoundaryInfo data structure since we
    // already cleared it.
    for (; it != end; ++it)
      delete *it;
    
    _nodes.clear();
  }
}



void ParallelMesh::renumber_nodes_and_elements ()
{
  
  START_LOG("renumber_nodes_and_elem()", "Mesh");
  
  // In Parallel we do *not* want to renumber anything, just to delete
  // any unused nodes or NULL elements

  // Loop over the elements.  Remember any node ids we see.
  std::map<unsigned int, bool> used_nodes;

  {      
    elem_iterator_imp  in = _elements.begin();
    elem_iterator_imp end = _elements.end();

    for (; in != end;)
      {
        Elem* elem = *in;

        if (elem)
          {
            // Notice this element's nodes.
            for (unsigned int n=0; n<elem->n_nodes(); n++)
              used_nodes[elem->node(n)] = true;
            ++in;
          }
        else
          {
            // Remove this non-element
            _elements.erase(in++);
          }
      }
  }

  // Nodes not connected to any elements are deleted
  {
    node_iterator_imp  it = _nodes.begin();
    node_iterator_imp end = _nodes.end();

    for (; it != end;)
      {
	Node *node = *it;
        if (!used_nodes[node->id()])
          {
	    // remove any boundary information associated with
	    // this node
	    this->boundary_info->remove (node);

	    // delete the node
	    delete node;

            _nodes.erase(it++);
	  }
        else
          ++it;
      }
  }
  
  STOP_LOG("renumber_nodes_and_elem()", "Mesh");
}



void ParallelMesh::delete_nonlocal_elements()
{
  std::vector<bool> local_nodes(this->max_node_id(), false);
  std::vector<bool> semilocal_elems(this->max_elem_id(), false);

  const_element_iterator l_elem_it = this->local_elements_begin(),
                         l_end     = this->local_elements_end();
  for (; l_elem_it != l_end; ++l_elem_it)
    {
      const Elem *elem = *l_elem_it;
      for (unsigned int n=0; n != elem->n_nodes(); ++n)
        local_nodes[elem->node(n)] = true;
    }

  element_iterator nl_elem_it = this->not_local_elements_begin(),
                   nl_end     = this->not_local_elements_end();
  for (; nl_elem_it != nl_end; ++nl_elem_it)
    {
      Elem *elem = *nl_elem_it;
      for (unsigned int n=0; n != elem->n_nodes(); ++n)
        if (local_nodes[elem->node(n)])
          {
            semilocal_elems[elem->id()] = true;
            break;
          }
      if (!semilocal_elems[elem->id()])
        {
          // delete_elem doesn't currently invalidate element
          // iterators... that had better not change
          this->delete_elem(elem);
        }
    }
//  We eventually want to compact the _nodes and _elems vectors
//  But not do a repartition or anything
//  do not this->prepare_for_use();
}

void ParallelMesh::restore_nonlocal_elements()
{
  // Someday...
  error();
}

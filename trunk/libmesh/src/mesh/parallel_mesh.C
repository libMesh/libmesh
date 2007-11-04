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
#include "mesh_communication.h"
#include "parallel_mesh.h"
#include "parallel.h"

// ------------------------------------------------------------
// ParallelMesh class member functions
ParallelMesh::ParallelMesh (unsigned int d) :
  UnstructuredMesh (d), _is_serial(true)
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



Node* ParallelMesh::insert_node (Node* n)
{
  // If we already have this node we cannot
  // simply delete it, because we may have elements
  // which are attached to its address.
  //
  // Instead, call the Node copy constructor to
  // overwrite the current node (but keep its address),
  // delete the provided node, and return the address of
  // the one we already had.
  if (_nodes.count(n->id()))
    {
      Node *my_n = _nodes[n->id()];

      *my_n = *n;
      delete n;
      n = my_n;
    }
  else
    _nodes[n->id()] = n;

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



template <typename T>
void ParallelMesh::renumber_dof_objects (mapvector<T*> &objects)
{
  typedef typename mapvector<T*>::veclike_iterator object_iterator;

  // In parallel we may not know what objects other processors have.
  // Start by figuring out how many
  unsigned int unpartitioned_objects = 0;

  std::vector<unsigned int>
    ghost_objects_from_proc(libMesh::n_processors(), 0);

  object_iterator it  = objects.begin();
  object_iterator end = objects.end();

  for (; it != end;)
    {
      T *obj = *it;

      // Remove any NULL container entries while we're here,
      // being careful not to invalidate our iterator
      if (!*it)
        objects.erase(it++);
      else
        {
          unsigned int obj_procid = obj->processor_id();
          if (obj_procid == DofObject::invalid_processor_id)
            unpartitioned_objects++;
          else
            ghost_objects_from_proc[obj_procid]++;
          ++it;
        }
    }

  std::vector<unsigned int> objects_on_proc(libMesh::n_processors(), 0);
  Parallel::allgather(ghost_objects_from_proc[libMesh::processor_id()],
                      objects_on_proc);

#ifndef NDEBUG
  unsigned int global_unpartitioned_objects = unpartitioned_objects;
  Parallel::max(global_unpartitioned_objects);
  assert(global_unpartitioned_objects == unpartitioned_objects);
  for (unsigned int p=0; p != libMesh::n_processors(); ++p)
    assert(ghost_objects_from_proc[p] <= objects_on_proc[p]);
#endif

  // We'll renumber objects in blocks by processor id
  std::vector<unsigned int> first_object_on_proc(libMesh::n_processors());
  for (unsigned int i=1; i != libMesh::n_processors(); ++i)
    first_object_on_proc[i] = first_object_on_proc[i-1] +
                              objects_on_proc[i-1];
  unsigned int next_id = first_object_on_proc[libMesh::processor_id()];

  // First set new local object ids and build request sets 
  // for non-local object ids
  
  // Request sets to send to each processor
  std::vector<std::vector<unsigned int> >
    requested_ids(libMesh::n_processors());

  // We know how many objects live on each processor, so reseve() space for
  // each.
  for (unsigned int p=0; p != libMesh::n_processors(); ++p)
    if (p != libMesh::processor_id())
      requested_ids[p].reserve(ghost_objects_from_proc[p]);

  end = objects.end();
  for (it = objects.begin(); it != end; ++it)
    {
      T *obj = *it;
      if (obj->processor_id() == libMesh::processor_id())
        obj->set_id(next_id++);
      else if (obj->processor_id() != DofObject::invalid_processor_id)
        requested_ids[obj->processor_id()].push_back(obj->id());
    }

  // Next set ghost object ids from other processors
  if (libMesh::n_processors() > 1)
    {
      for (unsigned int p=1; p != libMesh::n_processors(); ++p)
        {
          // Trade my requests with processor procup and procdown
          unsigned int procup = (libMesh::processor_id() + p) %
                                 libMesh::n_processors();
          unsigned int procdown = (libMesh::n_processors() +
                                   libMesh::processor_id() - p) %
                                   libMesh::n_processors();
          std::vector<unsigned int> request_to_fill;
          Parallel::send_receive(procup, requested_ids[procup],
                                 procdown, request_to_fill);

          // Fill those requests
          std::vector<unsigned int> new_ids(request_to_fill.size());
          for (unsigned int i=0; i != request_to_fill.size(); ++i)
            {
              assert(objects[request_to_fill[i]]);
              assert(objects[request_to_fill[i]]->processor_id()
                     == libMesh::processor_id());
              new_ids[i] = objects[request_to_fill[i]]->id();
              assert(new_ids[i] >=
                     first_object_on_proc[libMesh::processor_id()]);
              assert(new_ids[i] <
                     first_object_on_proc[libMesh::processor_id()] +
                     objects_on_proc[libMesh::processor_id()]);
            }

          // Trade back the results
          std::vector<unsigned int> filled_request;
          Parallel::send_receive(procdown, new_ids,
                                 procup, filled_request);

          // And copy the id changes we've now been informed of
          for (unsigned int i=0; i != filled_request.size(); ++i)
            {
              assert (objects[requested_ids[procup][i]]->processor_id()
                      == procup);
              assert(filled_request[i] >=
                     first_object_on_proc[procup]);
              assert(filled_request[i] <
                     first_object_on_proc[procup] +
                     objects_on_proc[procup]);
              objects[requested_ids[procup][i]]->set_id(filled_request[i]);
            }
        }
    }

  // Next set unpartitioned object ids
  next_id = 0;
  for (unsigned int i=0; i != libMesh::n_processors(); ++i)
    next_id += objects_on_proc[i];
  for (it = objects.begin(); it != end; ++it)
    {
      T *obj = *it;
      if (obj->processor_id() == DofObject::invalid_processor_id)
        obj->set_id(next_id++);
    }

  // Finally shuffle around objects so that container indices
  // match ids
  end = objects.end();
  for (it = objects.begin(); it != end;)
    {
      T *obj = *it;
      if (obj) // don't try shuffling already-NULL entries
        {
          T *next = objects[obj->id()];
          // If we have to move this object
          if (next != obj)
            {
              // NULL out its original position for now
              // (our shuffling may put another object there shortly)
              *it = NULL;

              // There may already be another object with this id that
              // needs to be moved itself
              while (next)
                {
                  // We shouldn't be trying to give two objects the
                  // same id
                  assert (next->id() != obj->id());
                  objects[obj->id()] = obj;
                  obj = next;
                  next = objects[obj->id()];
                }
              objects[obj->id()] = obj;
            }
        }
      // Remove any container entries that were left as NULL,
      // being careful not to invalidate our iterator
      if (!*it)
        objects.erase(it++);
      else
        ++it;
    }
}


void ParallelMesh::renumber_nodes_and_elements ()
{
  START_LOG("renumber_nodes_and_elem()", "ParallelMesh");

  std::set<unsigned int> used_nodes;

  // flag the nodes we need
  {      
    element_iterator  it = elements_begin();
    element_iterator end = elements_end();

    for (; it != end; ++it)
      {
        Elem *elem = *it;

        for (unsigned int n=0; n != elem->n_nodes(); ++n)
          used_nodes.insert(elem->node(n));
      }
  }

  // Nodes not connected to any local elements are deleted
  {
    node_iterator_imp  it = _nodes.begin();
    node_iterator_imp end = _nodes.end();

    for (; it != end;)
      {
	Node *node = *it;
        if (!used_nodes.count(node->id()))
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

  // Finally renumber all the elements
  this->renumber_dof_objects (_elements);

  // and all the remaining nodes
  this->renumber_dof_objects (_nodes);

  STOP_LOG("renumber_nodes_and_elem()", "ParallelMesh");
}



void ParallelMesh::delete_remote_elements()
{
  _is_serial = false;
  MeshCommunication().delete_remote_elements(*this);
}



void ParallelMesh::allgather()
{
  _is_serial = true;
  MeshCommunication().allgather(*this);
}

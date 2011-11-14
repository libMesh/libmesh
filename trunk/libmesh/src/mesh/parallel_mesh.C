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
#include "mesh_communication.h"
#include "parallel_mesh.h"
#include "parallel.h"
#include "parmetis_partitioner.h"

namespace libMesh
{

// ------------------------------------------------------------
// ParallelMesh class member functions
ParallelMesh::ParallelMesh (unsigned int d) :
  UnstructuredMesh (d), _is_serial(true),
  _n_nodes(0), _n_elem(0), _max_node_id(0), _max_elem_id(0),
  _next_free_local_node_id(libMesh::processor_id()),
  _next_free_local_elem_id(libMesh::processor_id()),
  _next_free_unpartitioned_node_id(libMesh::n_processors()),
  _next_free_unpartitioned_elem_id(libMesh::n_processors())
{
  _partitioner = AutoPtr<Partitioner>(new ParmetisPartitioner());
}


ParallelMesh::~ParallelMesh ()
{
  this->clear();  // Free nodes and elements
}


// This might be specialized later, but right now it's just here to
// make sure the compiler doesn't give us a default (non-deep) copy
// constructor instead.
ParallelMesh::ParallelMesh (const ParallelMesh &other_mesh) :
  UnstructuredMesh (other_mesh), _is_serial(other_mesh._is_serial)
{
  this->copy_nodes_and_elements(other_mesh);
  _n_nodes = other_mesh.n_nodes();
  _n_elem  = other_mesh.n_elem();
  _max_node_id = other_mesh.max_node_id();
  _max_elem_id = other_mesh.max_elem_id();
  *this->boundary_info = *other_mesh.boundary_info;

  // Need to copy extra_ghost_elems
  for(std::set<Elem *>::iterator it = other_mesh._extra_ghost_elems.begin();
      it != other_mesh._extra_ghost_elems.end();
      ++it)
    _extra_ghost_elems.insert(elem((*it)->id()));
}



ParallelMesh::ParallelMesh (const UnstructuredMesh &other_mesh) :
  UnstructuredMesh (other_mesh)
{
  this->copy_nodes_and_elements(other_mesh);
  _n_nodes = other_mesh.n_nodes();
  _n_elem  = other_mesh.n_elem();
  _max_node_id = other_mesh.max_node_id();
  _max_elem_id = other_mesh.max_elem_id();
  *this->boundary_info = *other_mesh.boundary_info;
}


// We use cached values for these so they can be called
// from one processor without bothering the rest, but
// we may need to update those caches before doing a full
// renumbering
void ParallelMesh::update_parallel_id_counts()
{
  // This function must be run on all processors at once
  parallel_only();

  _n_elem = this->n_local_elem();
  Parallel::sum(_n_elem);
  _n_elem += this->n_unpartitioned_elem();

  _max_elem_id = _elements.empty() ?
    0 : _elements.rbegin()->first + 1;
  Parallel::max(_max_elem_id);

  _n_nodes = this->n_local_nodes();
  Parallel::sum(_n_nodes);
  _n_nodes += this->n_unpartitioned_nodes();

  _max_node_id = _nodes.empty() ?
    0 : _nodes.rbegin()->first + 1;
  Parallel::max(_max_node_id);
}


// Or in debug mode we may want to test the uncached values without
// changing the cache
unsigned int ParallelMesh::parallel_n_elem() const
{
  // This function must be run on all processors at once
  parallel_only();

  unsigned int n_local = this->n_local_elem();
  Parallel::sum(n_local);
  n_local += this->n_unpartitioned_elem();
  return n_local;
}



unsigned int ParallelMesh::parallel_max_elem_id() const
{
  // This function must be run on all processors at once
  parallel_only();

  unsigned int max_local = _elements.empty() ?
    0 : _elements.rbegin()->first + 1;
  Parallel::max(max_local);
  return max_local;
}



unsigned int ParallelMesh::parallel_n_nodes() const
{
  // This function must be run on all processors at once
  parallel_only();

  unsigned int n_local = this->n_local_nodes();
  Parallel::sum(n_local);
  n_local += this->n_unpartitioned_nodes();
  return n_local;
}



unsigned int ParallelMesh::parallel_max_node_id() const
{
  // This function must be run on all processors at once
  parallel_only();

  unsigned int max_local = _nodes.empty() ?
    0 : _nodes.rbegin()->first + 1;
  Parallel::max(max_local);
  return max_local;
}



const Point& ParallelMesh::point (const unsigned int i) const
{
  libmesh_assert (_nodes[i] != NULL);
  libmesh_assert (_nodes[i]->id() == i);  

  return (*_nodes[i]);
}





const Node& ParallelMesh::node (const unsigned int i) const
{
  libmesh_assert (_nodes[i] != NULL);
  libmesh_assert (_nodes[i]->id() == i);  
  
  return (*_nodes[i]);
}





Node& ParallelMesh::node (const unsigned int i)
{
  libmesh_assert (_nodes[i] != NULL);
  libmesh_assert (_nodes[i]->id() == i);  

  return (*_nodes[i]);
}



const Node* ParallelMesh::node_ptr (const unsigned int i) const
{
//  libmesh_assert (_nodes[i] != NULL);
  libmesh_assert (_nodes[i] == NULL || _nodes[i]->id() == i);  
  
  return _nodes[i];
}




Node* & ParallelMesh::node_ptr (const unsigned int i)
{
//  libmesh_assert (_nodes[i] != NULL);
  libmesh_assert (_nodes[i] == NULL || _nodes[i]->id() == i);

  return _nodes[i];
}




Elem* ParallelMesh::elem (const unsigned int i) const
{
//  libmesh_assert (_elements[i] != NULL);
  libmesh_assert (_elements[i] == NULL || _elements[i]->id() == i);
  
  return _elements[i];
}




Elem* ParallelMesh::add_elem (Elem *e)
{
  // Don't try to add NULLs!
  libmesh_assert(e);

  const unsigned int elem_procid = e->processor_id();

  if (!e->valid_id())
    {
      // Use the unpartitioned ids for unpartitioned elems,
      // in serial, and temporarily for ghost elems
      unsigned int *next_id = &_next_free_unpartitioned_elem_id;
      if (elem_procid == libMesh::processor_id() &&
          !this->is_serial())
        next_id = &_next_free_local_elem_id;
      e->set_id (*next_id);
      *next_id += libMesh::n_processors() + 1;
    }
  else
    {
      if (_next_free_unpartitioned_elem_id <= e->id())
        _next_free_unpartitioned_elem_id += libMesh::n_processors() + 1;
      if (_next_free_local_elem_id <= e->id())
        _next_free_local_elem_id += libMesh::n_processors() + 1;
    }

  // Don't try to overwrite existing elems
  libmesh_assert (!_elements[e->id()]);

  _elements[e->id()] = e;

  // Make the cached elem data more accurate
  _n_elem++;
  _max_elem_id = std::max(_max_elem_id, e->id()+1);
  
// Unpartitioned elems should be added on every processor
// And shouldn't be added in the same batch as ghost elems
// But we might be just adding on processor 0 to
// broadcast later
#if 0
#ifdef DEBUG
  if (elem_procid == DofObject::invalid_processor_id)
    {
      unsigned int elem_id = e->id();
      Parallel::max(elem_id);
      libmesh_assert(elem_id == e->id());
    }
#endif
#endif

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
  libmesh_assert (e);

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



void ParallelMesh::renumber_elem(const unsigned int old_id,
                                 const unsigned int new_id)
{
  Elem *elem = _elements[old_id];
  libmesh_assert (elem);
  libmesh_assert (elem->id() == old_id);

  elem->set_id(new_id);
  libmesh_assert (!_elements[new_id]);
  _elements[new_id] = elem;
  _elements.erase(old_id);
}



Node* ParallelMesh::add_point (const Point& p,
			       const unsigned int id,
			       const unsigned int proc_id)
{  
  if (_nodes.count(id))
    {
      Node *n = _nodes[id];
      libmesh_assert (n);
      libmesh_assert (n->id() == id);
      
      *n = p;
      n->processor_id() = proc_id;

      return n;
    }

  Node* n = Node::build(p, id).release();
  n->processor_id() = proc_id;

  return ParallelMesh::add_node(n);
}



Node* ParallelMesh::add_node (Node *n)
{
  // Don't try to add NULLs!
  libmesh_assert(n);

  const unsigned int node_procid = n->processor_id();

  if (!n->valid_id())
    {
      // Use the unpartitioned ids for unpartitioned nodes,
      // in serial, and temporarily for ghost nodes
      unsigned int *next_id = &_next_free_unpartitioned_node_id;
      if (node_procid == libMesh::processor_id() &&
          !this->is_serial())
        next_id = &_next_free_local_node_id;
      n->set_id (*next_id);
      *next_id += libMesh::n_processors() + 1;
    }
  else
    {
      if (_next_free_unpartitioned_node_id <= n->id())
        _next_free_unpartitioned_node_id += libMesh::n_processors() + 1;
      if (_next_free_local_node_id <= n->id())
        _next_free_local_node_id += libMesh::n_processors() + 1;
    }

  // Don't try to overwrite existing nodes
  libmesh_assert (!_nodes[n->id()]);

  _nodes[n->id()] = n;
  
  // Make the cached elem data more accurate
  _n_nodes++;
  _max_node_id = std::max(_max_node_id, n->id()+1);
  
// Unpartitioned nodes should be added on every processor
// And shouldn't be added in the same batch as ghost nodes
// But we might be just adding on processor 0 to
// broadcast later
#if 0
#ifdef DEBUG
  if (node_procid == DofObject::invalid_processor_id)
    {
      unsigned int node_id = n->id();
      Parallel::max(node_id);
      libmesh_assert(node_id == n->id());
    }
#endif
#endif

  return n;
}



Node* ParallelMesh::insert_node (Node* n)
{
  // If we already have this node we cannot
  // simply delete it, because we may have elements
  // which are attached to its address.
  //
  // Instead, call the Node = Point assignment operator
  // to overwrite the spatial coordinates (but keep its
  // address), delete the provided node, and return the
  // address of the one we already had.
  if (_nodes.count(n->id()))
    {
      Node *my_n = _nodes[n->id()];

      *my_n = static_cast<Point>(*n);
      delete n;
      n = my_n;
    }
  else
    _nodes[n->id()] = n;

  return n;
}



void ParallelMesh::delete_node(Node* n)
{
  libmesh_assert (n != NULL);
  libmesh_assert (_nodes[n->id()] != NULL);

  // Delete the node from the BoundaryInfo object
  this->boundary_info->remove(n);

  // And from the container
  _nodes.erase(n->id());
  
  // delete the node
  delete n;
}



void ParallelMesh::renumber_node(const unsigned int old_id,
                                 const unsigned int new_id)
{
  Node *node = _nodes[old_id];
  libmesh_assert (node);
  libmesh_assert (node->id() == old_id);

  node->set_id(new_id);
  libmesh_assert (!_nodes[new_id]);
  _nodes[new_id] = node;
  _nodes.erase(old_id);
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

  // We're no longer distributed if we were before
  _is_serial = true;

  // Correct our caches
  _n_nodes = 0;
  _n_elem = 0;
  _max_node_id = 0;
  _max_elem_id = 0;
  _next_free_local_node_id = libMesh::processor_id();
  _next_free_local_elem_id = libMesh::processor_id();
  _next_free_unpartitioned_node_id = libMesh::n_processors();
  _next_free_unpartitioned_elem_id = libMesh::n_processors();
}



void ParallelMesh::partition (const unsigned int n_parts)
{
  // FIXME: We still need to MeshCommunication::redistribute to handle
  // dof indices before repartitioning won't break AMR
  if (!this->is_serial())
    return;

  if(!skip_partitioning())
  { 
    // Call base class' partition() function.
    MeshBase::partition(n_parts);

    // Partitioning changes our numbers of unpartitioned objects
    this->update_parallel_id_counts();

    // If this is a truly parallel mesh, go through the redistribution/gather/delete remote steps
    if (!this->is_serial())
      {
	// Construct a MeshCommunication object to actually redistribute the nodes
	// and elements according to the partitioner, and then to re-gather the neighbors.
	MeshCommunication mc;
	mc.redistribute(*this);
	mc.gather_neighboring_elements(*this);

	// Is this necessary?  If we are called from prepare_for_use(), this will be called
	// anyway... but users can always call partition directly, in which case we do need
	// to call delete_remote_elements()...
	this->delete_remote_elements();
      }
  }
}



template <typename T>
void ParallelMesh::libmesh_assert_valid_parallel_object_ids 
  (const mapvector<T*> &objects) const
{
  // This function must be run on all processors at once
  parallel_only();

  unsigned int pmax_elem_id = this->parallel_max_elem_id();

  for (unsigned int i=0; i != pmax_elem_id; ++i)
    {
      T* obj = objects[i]; // Returns NULL if there's no map entry

      unsigned int dofid = obj && obj->valid_id() ?
        obj->id() : DofObject::invalid_id;
      // Local lookups by id should return the requested object
      libmesh_assert(!obj || obj->id() == i);

      unsigned int min_dofid = dofid;
      Parallel::min(min_dofid);
      // All processors with an object should agree on id
      libmesh_assert (!obj || dofid == min_dofid);

      unsigned int procid = obj && obj->valid_processor_id() ?
        obj->processor_id() : DofObject::invalid_processor_id;

      unsigned int min_procid = procid;
      Parallel::min(min_procid);

      // All processors with an object should agree on processor id
      libmesh_assert (!obj || procid == min_procid);

      // Either:
      // 1.) I own this elem (min_procid == libMesh::processor_id()) *and* I have a valid pointer to it (obj != NULL)
      // or
      // 2.) I don't own this elem (min_procid != libMesh::processor_id()).  (In this case I may or may not have a valid pointer to it.)

      // Original assert logic
      // libmesh_assert (min_procid != libMesh::processor_id() || obj);

      // More human-understandable logic...
      libmesh_assert (
		      ((min_procid == libMesh::processor_id()) && obj)
		      ||
		      (min_procid != libMesh::processor_id())
		       );
    }
}



void ParallelMesh::libmesh_assert_valid_parallel_ids () const
{
  this->libmesh_assert_valid_parallel_object_ids (this->_elements);
  this->libmesh_assert_valid_parallel_object_ids (this->_nodes);
}



void ParallelMesh::libmesh_assert_valid_parallel_flags () const
{
#ifdef LIBMESH_ENABLE_AMR
  // This function must be run on all processors at once
  parallel_only();

  unsigned int pmax_elem_id = this->parallel_max_elem_id();

  for (unsigned int i=0; i != pmax_elem_id; ++i)
    {
      Elem* elem = _elements[i]; // Returns NULL if there's no map entry

      unsigned int refinement_flag   = elem ?
        static_cast<unsigned int> (elem->refinement_flag()) : libMesh::invalid_uint;
#ifndef NDEBUG
      unsigned int p_refinement_flag = elem ?
        static_cast<unsigned int> (elem->p_refinement_flag()) : libMesh::invalid_uint;
#endif

      unsigned int min_rflag = refinement_flag;
      Parallel::min(min_rflag);
      // All processors with this element should agree on flag
      libmesh_assert (!elem || min_rflag == refinement_flag);

#ifndef NDEBUG
      unsigned int min_pflag = p_refinement_flag;
#endif
      // All processors with this element should agree on flag
      libmesh_assert (!elem || min_pflag == p_refinement_flag);
    }
#endif // LIBMESH_ENABLE_AMR
}



template <typename T>
unsigned int ParallelMesh::renumber_dof_objects (mapvector<T*> &objects)
{
  // This function must be run on all processors at once
  parallel_only();

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
  libmesh_assert(global_unpartitioned_objects == unpartitioned_objects);
  for (unsigned int p=0; p != libMesh::n_processors(); ++p)
    libmesh_assert(ghost_objects_from_proc[p] <= objects_on_proc[p]);
#endif

  // We'll renumber objects in blocks by processor id
  std::vector<unsigned int> first_object_on_proc(libMesh::n_processors());
  for (unsigned int i=1; i != libMesh::n_processors(); ++i)
    first_object_on_proc[i] = first_object_on_proc[i-1] +
                              objects_on_proc[i-1];
  unsigned int next_id = first_object_on_proc[libMesh::processor_id()];
  unsigned int first_free_id =
    first_object_on_proc[libMesh::n_processors()-1] +
    objects_on_proc[libMesh::n_processors()-1] + 
    unpartitioned_objects;

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
              T *obj = objects[request_to_fill[i]];
              libmesh_assert(obj);
              libmesh_assert(obj->processor_id() == libMesh::processor_id());
              new_ids[i] = obj->id();
              libmesh_assert(new_ids[i] >=
                     first_object_on_proc[libMesh::processor_id()]);
              libmesh_assert(new_ids[i] <
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
              T *obj = objects[requested_ids[procup][i]];
              libmesh_assert (obj);
              libmesh_assert (obj->processor_id() == procup);
              libmesh_assert(filled_request[i] >=
                     first_object_on_proc[procup]);
              libmesh_assert(filled_request[i] <
                     first_object_on_proc[procup] +
                     objects_on_proc[procup]);
              obj->set_id(filled_request[i]);
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
                  libmesh_assert (next->id() != obj->id());
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

  return first_free_id;
}


void ParallelMesh::renumber_nodes_and_elements ()
{
  START_LOG("renumber_nodes_and_elements()", "ParallelMesh");

#ifdef DEBUG
// Make sure our ids and flags are consistent
  this->libmesh_assert_valid_parallel_ids();
  this->libmesh_assert_valid_parallel_flags();
#endif

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
  _n_elem = this->renumber_dof_objects (this->_elements);
  _max_elem_id = _n_elem;
  _next_free_local_elem_id = _n_elem;

  // and all the remaining nodes
  _n_nodes = this->renumber_dof_objects (this->_nodes);
  _max_node_id = _n_nodes;
  _next_free_local_node_id = _n_nodes;

  // And figure out what IDs we should use when adding new nodes and
  // new elements
  unsigned int cycle = libMesh::n_processors()+1;
  unsigned int offset = _next_free_local_elem_id % cycle;
  if (offset)
    _next_free_local_elem_id += cycle - offset;
  _next_free_unpartitioned_elem_id = _next_free_local_elem_id +
                                     libMesh::n_processors();
  _next_free_local_elem_id += libMesh::processor_id();

  offset = _next_free_local_node_id % cycle;
  if (offset)
    _next_free_local_node_id += cycle - offset;
  _next_free_unpartitioned_node_id = _next_free_local_node_id +
                                     libMesh::n_processors();
  _next_free_local_node_id += libMesh::processor_id();

// Make sure our caches are up to date and our
// DofObjects are well packed
#ifdef DEBUG
  libmesh_assert(this->n_nodes() == this->parallel_n_nodes());
  libmesh_assert(this->n_elem() == this->parallel_n_elem());
  libmesh_assert(this->max_node_id() == this->parallel_max_node_id());
  libmesh_assert(this->max_elem_id() == this->parallel_max_elem_id());
  libmesh_assert(this->n_nodes() == this->max_node_id());
  libmesh_assert(this->n_elem() == this->max_elem_id());

  // Make sure our ids and flags are consistent
  this->libmesh_assert_valid_parallel_ids();
  this->libmesh_assert_valid_parallel_flags();

// And make sure we've made our numbering monotonic
  MeshTools::libmesh_assert_valid_elem_ids(*this);
#endif

  STOP_LOG("renumber_nodes_and_elements()", "ParallelMesh");
}



void ParallelMesh::fix_broken_node_and_element_numbering ()
{
  // We need access to iterators for the underlying containers,
  // not the mapvector<> reimplementations.
  mapvector<Node*>::maptype &nodes = this->_nodes;
  mapvector<Elem*>::maptype &elem  = this->_elements;
  
  // Nodes first
  {
    mapvector<Node*>::maptype::iterator
      it  = nodes.begin(),
      end = nodes.end();

    for (; it != end; ++it)
      if (it->second != NULL)
	it->second->set_id() = it->first;
  }

  // Elements next
  {
    mapvector<Elem*>::maptype::iterator
      it  = elem.begin(),
      end = elem.end();

    for (; it != end; ++it)
      if (it->second != NULL)
	it->second->set_id() = it->first;
  }
}



unsigned int ParallelMesh::n_active_elem () const
{
  parallel_only();

  // Get local active elements first
  unsigned int active_elements =
    static_cast<unsigned int>(std::distance (this->active_local_elements_begin(),
					     this->active_local_elements_end()));
  Parallel::sum(active_elements);

  // Then add unpartitioned active elements, which should exist on
  // every processor
  active_elements +=
    static_cast<unsigned int>(std::distance
      (this->active_pid_elements_begin(DofObject::invalid_processor_id),
       this->active_pid_elements_end(DofObject::invalid_processor_id)));
  return active_elements;
}



void ParallelMesh::delete_remote_elements()
{
#ifdef DEBUG
// Make sure our neighbor links are all fine
  MeshTools::libmesh_assert_valid_neighbors(*this);

// And our child/parent links, and our flags
  MeshTools::libmesh_assert_valid_refinement_tree(*this);
#endif

  _is_serial = false;
  MeshCommunication().delete_remote_elements(*this, _extra_ghost_elems);

#ifdef DEBUG
// Make sure our caches are up to date and our
// DofObjects are well packed
  libmesh_assert(this->n_nodes() == this->parallel_n_nodes());
  libmesh_assert(this->n_elem() == this->parallel_n_elem());
  libmesh_assert(this->max_node_id() == this->parallel_max_node_id());
  libmesh_assert(this->max_elem_id() == this->parallel_max_elem_id());
  libmesh_assert(this->n_nodes() == this->max_node_id());
  libmesh_assert(this->n_elem() == this->max_elem_id());

// Make sure our neighbor links are all fine
  MeshTools::libmesh_assert_valid_neighbors(*this);

// And our child/parent links, and our flags
  MeshTools::libmesh_assert_valid_refinement_tree(*this);

// Make sure our ids and flags are consistent
  this->libmesh_assert_valid_parallel_ids();
  this->libmesh_assert_valid_parallel_flags();

// And make sure our numbering is still monotonic
  MeshTools::libmesh_assert_valid_elem_ids(*this);
#endif
}


void ParallelMesh::insert_extra_ghost_elem(Elem* e)
{
  // First insert the elem like normal
  insert_elem(e);

  // Now add it to the set that won't be deleted when we call
  // delete_remote_elements()
  _extra_ghost_elems.insert(e);
}


void ParallelMesh::allgather()
{
  if (_is_serial)
    return;
  _is_serial = true;
  MeshCommunication().allgather(*this);

// Make sure our caches are up to date and our
// DofObjects are well packed
#ifdef DEBUG
  libmesh_assert(this->n_nodes() == this->parallel_n_nodes());
  libmesh_assert(this->n_elem() == this->parallel_n_elem());
  libmesh_assert(this->max_node_id() == this->parallel_max_node_id());
  libmesh_assert(this->max_elem_id() == this->parallel_max_elem_id());
  libmesh_assert(this->n_nodes() == this->max_node_id());
  libmesh_assert(this->n_elem() == this->max_elem_id());

// Make sure our neighbor links are all fine
  MeshTools::libmesh_assert_valid_neighbors(*this);

// Make sure our ids and flags are consistent
  this->libmesh_assert_valid_parallel_ids();
  this->libmesh_assert_valid_parallel_flags();
#endif
}

} // namespace libMesh

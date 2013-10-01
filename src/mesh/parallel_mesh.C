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
#include "libmesh/mesh_communication.h"
#include "libmesh/parallel_mesh.h"
#include "libmesh/parallel.h"
#include "libmesh/parmetis_partitioner.h"

namespace libMesh
{

// ------------------------------------------------------------
// ParallelMesh class member functions
ParallelMesh::ParallelMesh (const Parallel::Communicator &comm,
			    unsigned int d) :
  UnstructuredMesh (comm,d), _is_serial(true),
  _n_nodes(0), _n_elem(0), _max_node_id(0), _max_elem_id(0),
  _next_free_local_node_id(this->processor_id()),
  _next_free_local_elem_id(this->processor_id()),
  _next_free_unpartitioned_node_id(this->n_processors()),
  _next_free_unpartitioned_elem_id(this->n_processors())
{
#ifdef LIBMESH_ENABLE_UNIQUE_ID
  _next_unique_id = this->processor_id();
#endif

  // FIXME: give parmetis the communicator!
  _partitioner = AutoPtr<Partitioner>(new ParmetisPartitioner());
}


#ifndef LIBMESH_DISABLE_COMMWORLD
// ------------------------------------------------------------
// ParallelMesh class member functions
ParallelMesh::ParallelMesh (unsigned int d) :
  UnstructuredMesh (d), _is_serial(true),
  _n_nodes(0), _n_elem(0), _max_node_id(0), _max_elem_id(0),
  _next_free_local_node_id(this->processor_id()),
  _next_free_local_elem_id(this->processor_id()),
  _next_free_unpartitioned_node_id(this->n_processors()),
  _next_free_unpartitioned_elem_id(this->n_processors())
{
#ifdef LIBMESH_ENABLE_UNIQUE_ID
  _next_unique_id = this->processor_id();
#endif

  // FIXME: give parmetis the communicator!
  _partitioner = AutoPtr<Partitioner>(new ParmetisPartitioner());
}
#endif


ParallelMesh::~ParallelMesh ()
{
  this->clear();  // Free nodes and elements
}


// This might be specialized later, but right now it's just here to
// make sure the compiler doesn't give us a default (non-deep) copy
// constructor instead.
ParallelMesh::ParallelMesh (const ParallelMesh &other_mesh) :
  UnstructuredMesh (other_mesh), _is_serial(other_mesh._is_serial),
  _n_nodes(0), _n_elem(0), _max_node_id(0), _max_elem_id(0),
  _next_free_local_node_id(this->processor_id()),
  _next_free_local_elem_id(this->processor_id()),
  _next_free_unpartitioned_node_id(this->n_processors()),
  _next_free_unpartitioned_elem_id(this->n_processors())
{
  this->copy_nodes_and_elements(other_mesh);
  _n_nodes = other_mesh.n_nodes();
  _n_elem  = other_mesh.n_elem();
  _max_node_id = other_mesh.max_node_id();
  _max_elem_id = other_mesh.max_elem_id();
  _next_free_local_node_id =
    other_mesh._next_free_local_node_id;
  _next_free_local_elem_id =
    other_mesh._next_free_local_elem_id;
  _next_free_unpartitioned_node_id =
    other_mesh._next_free_unpartitioned_node_id;
  _next_free_unpartitioned_elem_id =
    other_mesh._next_free_unpartitioned_elem_id;
  *this->boundary_info = *other_mesh.boundary_info;

  // Need to copy extra_ghost_elems
  for(std::set<Elem *>::iterator it = other_mesh._extra_ghost_elems.begin();
      it != other_mesh._extra_ghost_elems.end();
      ++it)
    _extra_ghost_elems.insert(elem((*it)->id()));
}



ParallelMesh::ParallelMesh (const UnstructuredMesh &other_mesh) :
  UnstructuredMesh (other_mesh), _is_serial(other_mesh.is_serial()),
  _n_nodes(0), _n_elem(0), _max_node_id(0), _max_elem_id(0),
  _next_free_local_node_id(this->processor_id()),
  _next_free_local_elem_id(this->processor_id()),
  _next_free_unpartitioned_node_id(this->n_processors()),
  _next_free_unpartitioned_elem_id(this->n_processors())
{
  this->copy_nodes_and_elements(other_mesh);
  *this->boundary_info = *other_mesh.boundary_info;

  this->update_parallel_id_counts();
}


// We use cached values for these so they can be called
// from one processor without bothering the rest, but
// we may need to update those caches before doing a full
// renumbering
void ParallelMesh::update_parallel_id_counts()
{
  // This function must be run on all processors at once
  parallel_object_only();

  _n_elem  = this->parallel_n_elem();
  _n_nodes = this->parallel_n_nodes();
  _max_node_id = this->parallel_max_node_id();
  _max_elem_id = this->parallel_max_elem_id();

  if (_next_free_unpartitioned_elem_id < _max_elem_id)
    _next_free_unpartitioned_elem_id =
      ((_max_elem_id-1) / (this->n_processors() + 1) + 1) *
        (this->n_processors() + 1) + this->n_processors();
  if (_next_free_local_elem_id < _max_elem_id)
    _next_free_local_elem_id =
      ((_max_elem_id + this->n_processors() - 1) / (this->n_processors() + 1) + 1) *
        (this->n_processors() + 1) + this->processor_id();

  if (_next_free_unpartitioned_node_id < _max_node_id)
    _next_free_unpartitioned_node_id =
      ((_max_node_id-1) / (this->n_processors() + 1) + 1) *
        (this->n_processors() + 1) + this->n_processors();
  if (_next_free_local_node_id < _max_node_id)
    _next_free_local_node_id =
      ((_max_node_id + this->n_processors() - 1) / (this->n_processors() + 1) + 1) *
        (this->n_processors() + 1) + this->processor_id();
}


// Or in debug mode we may want to test the uncached values without
// changing the cache
dof_id_type ParallelMesh::parallel_n_elem() const
{
  // This function must be run on all processors at once
  parallel_object_only();

  dof_id_type n_local = this->n_local_elem();
  this->comm().sum(n_local);
  n_local += this->n_unpartitioned_elem();
  return n_local;
}



dof_id_type ParallelMesh::parallel_max_elem_id() const
{
  // This function must be run on all processors at once
  parallel_object_only();

  dof_id_type max_local = _elements.empty() ?
    0 : _elements.rbegin()->first + 1;
  this->comm().max(max_local);
  return max_local;
}



dof_id_type ParallelMesh::parallel_n_nodes() const
{
  // This function must be run on all processors at once
  parallel_object_only();

  dof_id_type n_local = this->n_local_nodes();
  this->comm().sum(n_local);
  n_local += this->n_unpartitioned_nodes();
  return n_local;
}



dof_id_type ParallelMesh::parallel_max_node_id() const
{
  // This function must be run on all processors at once
  parallel_object_only();

  dof_id_type max_local = _nodes.empty() ?
    0 : _nodes.rbegin()->first + 1;
  this->comm().max(max_local);
  return max_local;
}



const Point& ParallelMesh::point (const dof_id_type i) const
{
  libmesh_assert(_nodes[i]);
  libmesh_assert_equal_to (_nodes[i]->id(), i);

  return (*_nodes[i]);
}





const Node& ParallelMesh::node (const dof_id_type i) const
{
  libmesh_assert(_nodes[i]);
  libmesh_assert_equal_to (_nodes[i]->id(), i);

  return (*_nodes[i]);
}





Node& ParallelMesh::node (const dof_id_type i)
{
  libmesh_assert(_nodes[i]);
  libmesh_assert_equal_to (_nodes[i]->id(), i);

  return (*_nodes[i]);
}



const Node* ParallelMesh::node_ptr (const dof_id_type i) const
{
  libmesh_assert(_nodes[i]);
  libmesh_assert_equal_to (_nodes[i]->id(), i);

  return _nodes[i];
}




Node* ParallelMesh::node_ptr (const dof_id_type i)
{
  libmesh_assert(_nodes[i]);
  libmesh_assert_equal_to (_nodes[i]->id(), i);

  return _nodes[i];
}




const Node* ParallelMesh::query_node_ptr (const dof_id_type i) const
{
  std::map<dof_id_type, Node*>::const_iterator it = _nodes.find(i);
  if (it != _nodes.end().it)
    {
      const Node* n = it->second;
      libmesh_assert (!n || n->id() == i);
      return n;
    }

  return NULL;
}




Node* ParallelMesh::query_node_ptr (const dof_id_type i)
{
  std::map<dof_id_type, Node*>::const_iterator it = _nodes.find(i);
  if (it != _nodes.end().it)
    {
      Node* n = it->second;
      libmesh_assert (!n || n->id() == i);
      return n;
    }

  return NULL;
}




const Elem* ParallelMesh::elem (const dof_id_type i) const
{
  libmesh_assert(_elements[i]);
  libmesh_assert_equal_to (_elements[i]->id(), i);

  return _elements[i];
}




Elem* ParallelMesh::elem (const dof_id_type i)
{
  libmesh_assert(_elements[i]);
  libmesh_assert_equal_to (_elements[i]->id(), i);

  return _elements[i];
}




const Elem* ParallelMesh::query_elem (const dof_id_type i) const
{
  std::map<dof_id_type, Elem*>::const_iterator it = _elements.find(i);
  if (it != _elements.end().it)
    {
      const Elem* e = it->second;
      libmesh_assert (!e || e->id() == i);
      return e;
    }

  return NULL;
}




Elem* ParallelMesh::query_elem (const dof_id_type i)
{
  std::map<dof_id_type, Elem*>::const_iterator it = _elements.find(i);
  if (it != _elements.end().it)
    {
      Elem* e = _elements[i];
      libmesh_assert (!e || e->id() == i);
      return e;
    }

  return NULL;
}




Elem* ParallelMesh::add_elem (Elem *e)
{
  // Don't try to add NULLs!
  libmesh_assert(e);

  // Trying to add an existing element is a no-op
  if (e->valid_id() && _elements[e->id()] == e)
    return e;

  const processor_id_type elem_procid = e->processor_id();

  if (!e->valid_id())
    {
      // We should only be creating new ids past the end of the range
      // of existing ids
      libmesh_assert_greater_equal(_next_free_unpartitioned_elem_id,
                                   _max_elem_id);
      libmesh_assert_greater_equal(_next_free_local_elem_id, _max_elem_id);

      // Use the unpartitioned ids for unpartitioned elems,
      // in serial, and temporarily for ghost elems
      dof_id_type *next_id = &_next_free_unpartitioned_elem_id;
      if (elem_procid == this->processor_id() &&
          !this->is_serial())
        next_id = &_next_free_local_elem_id;
      e->set_id (*next_id);
    }

    {
      // Advance next_ids up high enough that each is pointing to an
      // unused id and any subsequent increments will still point us
      // to unused ids
      _max_elem_id = std::max(_max_elem_id,
			      static_cast<dof_id_type>(e->id()+1));

      if (_next_free_unpartitioned_elem_id < _max_elem_id)
        _next_free_unpartitioned_elem_id =
          ((_max_elem_id-1) / (this->n_processors() + 1) + 1) *
            (this->n_processors() + 1) + this->n_processors();
      if (_next_free_local_elem_id < _max_elem_id)
        _next_free_local_elem_id =
          ((_max_elem_id + this->n_processors() - 1) / (this->n_processors() + 1) + 1) *
            (this->n_processors() + 1) + this->processor_id();

#ifndef NDEBUG
      // We need a const mapvector so we don't inadvertently create
      // NULL entries when testing for non-NULL ones
      const mapvector<Elem*,dof_id_type>& const_elements = _elements;
#endif
      libmesh_assert(!const_elements[_next_free_unpartitioned_elem_id]);
      libmesh_assert(!const_elements[_next_free_local_elem_id]);
    }

  // Don't try to overwrite existing elems
  libmesh_assert (!_elements[e->id()]);

  _elements[e->id()] = e;

  // Try to make the cached elem data more accurate
  if (elem_procid == this->processor_id() ||
      elem_procid == DofObject::invalid_processor_id)
    _n_elem++;

// Unpartitioned elems should be added on every processor
// And shouldn't be added in the same batch as ghost elems
// But we might be just adding on processor 0 to
// broadcast later
#if 0
#ifdef DEBUG
  if (elem_procid == DofObject::invalid_processor_id)
    {
      dof_id_type elem_id = e->id();
      this->comm().max(elem_id);
      libmesh_assert_equal_to (elem_id, e->id());
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



void ParallelMesh::renumber_elem(const dof_id_type old_id,
                                 const dof_id_type new_id)
{
  Elem *el = _elements[old_id];
  libmesh_assert (el);
  libmesh_assert_equal_to (el->id(), old_id);

  el->set_id(new_id);
  libmesh_assert (!_elements[new_id]);
  _elements[new_id] = el;
  _elements.erase(old_id);
}



Node* ParallelMesh::add_point (const Point& p,
			       const dof_id_type id,
			       const processor_id_type proc_id)
{
  if (_nodes.count(id))
    {
      Node *n = _nodes[id];
      libmesh_assert (n);
      libmesh_assert_equal_to (n->id(), id);

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

  // Trying to add an existing node is a no-op
  if (n->valid_id() && _nodes[n->id()] == n)
    return n;

  const processor_id_type node_procid = n->processor_id();

  if (!n->valid_id())
    {
      // We should only be creating new ids past the end of the range
      // of existing ids
      libmesh_assert_greater_equal(_next_free_unpartitioned_node_id,
                                   _max_node_id);
      libmesh_assert_greater_equal(_next_free_local_node_id, _max_node_id);

      // Use the unpartitioned ids for unpartitioned nodes,
      // in serial, and temporarily for ghost nodes
      dof_id_type *next_id = &_next_free_unpartitioned_node_id;
      if (node_procid == this->processor_id() &&
          !this->is_serial())
        next_id = &_next_free_local_node_id;
      n->set_id (*next_id);
    }

    {
      // Advance next_ids up high enough that each is pointing to an
      // unused id and any subsequent increments will still point us
      // to unused ids
      _max_node_id = std::max(_max_node_id,
			      static_cast<dof_id_type>(n->id()+1));

      if (_next_free_unpartitioned_node_id < _max_node_id)
        _next_free_unpartitioned_node_id =
          ((_max_node_id-1) / (this->n_processors() + 1) + 1) *
            (this->n_processors() + 1) + this->n_processors();
      if (_next_free_local_node_id < _max_node_id)
        _next_free_local_node_id =
          ((_max_node_id + this->n_processors() - 1) / (this->n_processors() + 1) + 1) *
            (this->n_processors() + 1) + this->processor_id();

#ifndef NDEBUG
      // We need a const mapvector so we don't inadvertently create
      // NULL entries when testing for non-NULL ones
      const mapvector<Node*,dof_id_type>& const_nodes = _nodes;
#endif
      libmesh_assert(!const_nodes[_next_free_unpartitioned_node_id]);
      libmesh_assert(!const_nodes[_next_free_local_node_id]);
    }

  // Don't try to overwrite existing nodes
  libmesh_assert (!_nodes[n->id()]);

  _nodes[n->id()] = n;

  // Try to make the cached node data more accurate
  if (node_procid == this->processor_id() ||
      node_procid == DofObject::invalid_processor_id)
    _n_nodes++;

// Unpartitioned nodes should be added on every processor
// And shouldn't be added in the same batch as ghost nodes
// But we might be just adding on processor 0 to
// broadcast later
#if 0
#ifdef DEBUG
  if (node_procid == DofObject::invalid_processor_id)
    {
      dof_id_type node_id = n->id();
      this->comm().max(node_id);
      libmesh_assert_equal_to (node_id, n->id());
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
  libmesh_assert(n);
  libmesh_assert(_nodes[n->id()]);

  // Delete the node from the BoundaryInfo object
  this->boundary_info->remove(n);

  // But not yet from the container; we might invalidate
  // an iterator that way!

  //_nodes.erase(n->id());

  // Instead, we set it to NULL for now

  _nodes[n->id()] = NULL;

  // delete the node
  delete n;
}



void ParallelMesh::renumber_node(const dof_id_type old_id,
                                 const dof_id_type new_id)
{
  Node *nd = _nodes[old_id];
  libmesh_assert (nd);
  libmesh_assert_equal_to (nd->id(), old_id);

  nd->set_id(new_id);
  libmesh_assert (!_nodes[new_id]);
  _nodes[new_id] = nd;
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
  _next_free_local_node_id = this->processor_id();
  _next_free_local_elem_id = this->processor_id();
  _next_free_unpartitioned_node_id = this->n_processors();
  _next_free_unpartitioned_elem_id = this->n_processors();
}



void ParallelMesh::redistribute ()
{
  // If this is a truly parallel mesh, go through the redistribution/gather/delete remote steps
  if (!this->is_serial())
    {
      // Construct a MeshCommunication object to actually redistribute the nodes
      // and elements according to the partitioner, and then to re-gather the neighbors.
      MeshCommunication mc;
      mc.redistribute(*this);

      this->update_parallel_id_counts();

      // Is this necessary?  If we are called from prepare_for_use(), this will be called
      // anyway... but users can always call partition directly, in which case we do need
      // to call delete_remote_elements()...
      //
      // Regardless of whether it's necessary, it isn't safe.  We
      // haven't communicated new node processor_ids yet, and we can't
      // delete nodes until we do.
      // this->delete_remote_elements();
    }
}



void ParallelMesh::update_post_partitioning ()
{
  // this->recalculate_n_partitions();

  // Partitioning changes our numbers of unpartitioned objects
  this->update_parallel_id_counts();
}



template <typename T>
void ParallelMesh::libmesh_assert_valid_parallel_object_ids
  (const mapvector<T*,dof_id_type> &objects) const
{
  // This function must be run on all processors at once
  parallel_object_only();

  const dof_id_type pmax_node_id = this->parallel_max_node_id();
  const dof_id_type pmax_elem_id = this->parallel_max_elem_id();
  const dof_id_type pmax_id = std::max(pmax_node_id, pmax_elem_id);

  for (dof_id_type i=0; i != pmax_id; ++i)
    {
      T* obj = objects[i]; // Returns NULL if there's no map entry

      dof_id_type dofid = obj && obj->valid_id() ?
        obj->id() : DofObject::invalid_id;
      // Local lookups by id should return the requested object
      libmesh_assert(!obj || obj->id() == i);

      dof_id_type min_dofid = dofid;
      this->comm().min(min_dofid);
      // All processors with an object should agree on id
      libmesh_assert (!obj || dofid == min_dofid);

      dof_id_type procid = obj && obj->valid_processor_id() ?
        obj->processor_id() : DofObject::invalid_processor_id;

      dof_id_type min_procid = procid;
      this->comm().min(min_procid);

      // All processors with an object should agree on processor id
      libmesh_assert (!obj || procid == min_procid);

      // Either:
      // 1.) I own this elem (min_procid == this->processor_id()) *and* I have a valid pointer to it (obj != NULL)
      // or
      // 2.) I don't own this elem (min_procid != this->processor_id()).  (In this case I may or may not have a valid pointer to it.)

      // Original assert logic
      // libmesh_assert (min_procid != this->processor_id() || obj);

      // More human-understandable logic...
      libmesh_assert (
		      ((min_procid == this->processor_id()) && obj)
		      ||
		      (min_procid != this->processor_id())
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
  parallel_object_only();

  dof_id_type pmax_elem_id = this->parallel_max_elem_id();

  for (dof_id_type i=0; i != pmax_elem_id; ++i)
    {
      Elem* el = _elements[i]; // Returns NULL if there's no map entry

      unsigned int refinement_flag   = el ?
        static_cast<unsigned int> (el->refinement_flag()) : libMesh::invalid_uint;
#ifndef NDEBUG
      unsigned int p_refinement_flag = el ?
        static_cast<unsigned int> (el->p_refinement_flag()) : libMesh::invalid_uint;
#endif

      unsigned int min_rflag = refinement_flag;
      this->comm().min(min_rflag);
      // All processors with this element should agree on flag
      libmesh_assert (!el || min_rflag == refinement_flag);

#ifndef NDEBUG
      unsigned int min_pflag = p_refinement_flag;
#endif
      // All processors with this element should agree on flag
      libmesh_assert (!el || min_pflag == p_refinement_flag);
    }
#endif // LIBMESH_ENABLE_AMR
}



template <typename T>
dof_id_type ParallelMesh::renumber_dof_objects
  (mapvector<T*,dof_id_type> &objects)
{
  // This function must be run on all processors at once
  parallel_object_only();

  typedef typename mapvector<T*,dof_id_type>::veclike_iterator object_iterator;

  // In parallel we may not know what objects other processors have.
  // Start by figuring out how many
  dof_id_type unpartitioned_objects = 0;

  std::vector<dof_id_type>
    ghost_objects_from_proc(this->n_processors(), 0);

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
          processor_id_type obj_procid = obj->processor_id();
          if (obj_procid == DofObject::invalid_processor_id)
            unpartitioned_objects++;
          else
            ghost_objects_from_proc[obj_procid]++;
          ++it;
        }
    }

  std::vector<dof_id_type> objects_on_proc(this->n_processors(), 0);
  this->comm().allgather(ghost_objects_from_proc[this->processor_id()],
				 objects_on_proc);

#ifndef NDEBUG
  libmesh_assert(this->comm().verify(unpartitioned_objects));
  for (processor_id_type p=0; p != this->n_processors(); ++p)
    libmesh_assert_less_equal (ghost_objects_from_proc[p], objects_on_proc[p]);
#endif

  // We'll renumber objects in blocks by processor id
  std::vector<dof_id_type> first_object_on_proc(this->n_processors());
  for (processor_id_type i=1; i != this->n_processors(); ++i)
    first_object_on_proc[i] = first_object_on_proc[i-1] +
                              objects_on_proc[i-1];
  dof_id_type next_id = first_object_on_proc[this->processor_id()];
  dof_id_type first_free_id =
    first_object_on_proc[this->n_processors()-1] +
    objects_on_proc[this->n_processors()-1] +
    unpartitioned_objects;

  // First set new local object ids and build request sets
  // for non-local object ids

  // Request sets to send to each processor
  std::vector<std::vector<dof_id_type> >
    requested_ids(this->n_processors());

  // We know how many objects live on each processor, so reseve() space for
  // each.
  for (processor_id_type p=0; p != this->n_processors(); ++p)
    if (p != this->processor_id())
      requested_ids[p].reserve(ghost_objects_from_proc[p]);

  end = objects.end();
  for (it = objects.begin(); it != end; ++it)
    {
      T *obj = *it;
      if (obj->processor_id() == this->processor_id())
        obj->set_id(next_id++);
      else if (obj->processor_id() != DofObject::invalid_processor_id)
        requested_ids[obj->processor_id()].push_back(obj->id());
    }

  // Next set ghost object ids from other processors
  if (this->n_processors() > 1)
    {
      for (processor_id_type p=1; p != this->n_processors(); ++p)
        {
          // Trade my requests with processor procup and procdown
          processor_id_type procup = (this->processor_id() + p) %
                                      this->n_processors();
          processor_id_type procdown = (this->n_processors() +
                                        this->processor_id() - p) %
                                        this->n_processors();
          std::vector<dof_id_type> request_to_fill;
          this->comm().send_receive(procup, requested_ids[procup],
					    procdown, request_to_fill);

          // Fill those requests
          std::vector<dof_id_type> new_ids(request_to_fill.size());
          for (std::size_t i=0; i != request_to_fill.size(); ++i)
            {
              T *obj = objects[request_to_fill[i]];
              libmesh_assert(obj);
              libmesh_assert_equal_to (obj->processor_id(), this->processor_id());
              new_ids[i] = obj->id();
              libmesh_assert_greater_equal (new_ids[i],
                     first_object_on_proc[this->processor_id()]);
              libmesh_assert_less (new_ids[i],
                     first_object_on_proc[this->processor_id()] +
                     objects_on_proc[this->processor_id()]);
            }

          // Trade back the results
          std::vector<dof_id_type> filled_request;
          this->comm().send_receive(procdown, new_ids,
					    procup, filled_request);

          // And copy the id changes we've now been informed of
          for (std::size_t i=0; i != filled_request.size(); ++i)
            {
              T *obj = objects[requested_ids[procup][i]];
              libmesh_assert (obj);
              libmesh_assert_equal_to (obj->processor_id(), procup);
              libmesh_assert_greater_equal (filled_request[i],
                     first_object_on_proc[procup]);
              libmesh_assert_less (filled_request[i],
                     first_object_on_proc[procup] +
                     objects_on_proc[procup]);
              obj->set_id(filled_request[i]);
            }
        }
    }

  // Next set unpartitioned object ids
  next_id = 0;
  for (processor_id_type i=0; i != this->n_processors(); ++i)
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
                  libmesh_assert_not_equal_to (next->id(), obj->id());
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
  parallel_object_only();

  if (_skip_renumber_nodes_and_elements)
    {
      this->update_parallel_id_counts();
      return;
    }

  START_LOG("renumber_nodes_and_elements()", "ParallelMesh");

#ifdef DEBUG
// Make sure our ids and flags are consistent
  this->libmesh_assert_valid_parallel_ids();
  this->libmesh_assert_valid_parallel_flags();
#endif

  std::set<dof_id_type> used_nodes;

  // flag the nodes we need
  {
    element_iterator  it = elements_begin();
    element_iterator end = elements_end();

    for (; it != end; ++it)
      {
        Elem *el = *it;

        for (unsigned int n=0; n != el->n_nodes(); ++n)
          used_nodes.insert(el->node(n));
      }
  }

  // Nodes not connected to any local elements, and NULL node entries
  // in our container, are deleted
  {
    node_iterator_imp  it = _nodes.begin();
    node_iterator_imp end = _nodes.end();

    for (; it != end;)
      {
	Node *nd = *it;
        if (!nd)
          _nodes.erase(it++);
        else if (!used_nodes.count(nd->id()))
          {
	    // remove any boundary information associated with
	    // this node
	    this->boundary_info->remove (nd);

	    // delete the node
	    delete nd;

            _nodes.erase(it++);
	  }
        else
          ++it;
      }
  }

  // Finally renumber all the elements
  _n_elem = this->renumber_dof_objects (this->_elements);

  // and all the remaining nodes
  _n_nodes = this->renumber_dof_objects (this->_nodes);

  // And figure out what IDs we should use when adding new nodes and
  // new elements
  this->update_parallel_id_counts();

// Make sure our caches are up to date and our
// DofObjects are well packed
#ifdef DEBUG
  libmesh_assert_equal_to (this->n_nodes(), this->parallel_n_nodes());
  libmesh_assert_equal_to (this->n_elem(), this->parallel_n_elem());
  const dof_id_type pmax_node_id = this->parallel_max_node_id();
  const dof_id_type pmax_elem_id = this->parallel_max_elem_id();
  libmesh_assert_equal_to (this->max_node_id(), pmax_node_id);
  libmesh_assert_equal_to (this->max_elem_id(), pmax_elem_id);
  libmesh_assert_equal_to (this->n_nodes(), this->max_node_id());
  libmesh_assert_equal_to (this->n_elem(), this->max_elem_id());

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
  mapvector<Node*,dof_id_type>::maptype &nodes = this->_nodes;
  mapvector<Elem*,dof_id_type>::maptype &elems = this->_elements;

  // Nodes first
  {
    mapvector<Node*,dof_id_type>::maptype::iterator
      it  = nodes.begin(),
      end = nodes.end();

    for (; it != end; ++it)
      if (it->second != NULL)
	it->second->set_id() = it->first;
  }

  // Elements next
  {
    mapvector<Elem*,dof_id_type>::maptype::iterator
      it  = elems.begin(),
      end = elems.end();

    for (; it != end; ++it)
      if (it->second != NULL)
	it->second->set_id() = it->first;
  }
}



dof_id_type ParallelMesh::n_active_elem () const
{
  parallel_object_only();

  // Get local active elements first
  dof_id_type active_elements =
    static_cast<dof_id_type>(std::distance (this->active_local_elements_begin(),
					    this->active_local_elements_end()));
  this->comm().sum(active_elements);

  // Then add unpartitioned active elements, which should exist on
  // every processor
  active_elements +=
    static_cast<dof_id_type>(std::distance
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

// Make sure our ids and flags are consistent
  this->libmesh_assert_valid_parallel_ids();
  this->libmesh_assert_valid_parallel_flags();

  libmesh_assert_equal_to (this->n_nodes(), this->parallel_n_nodes());
  libmesh_assert_equal_to (this->n_elem(), this->parallel_n_elem());
  const dof_id_type pmax_node_id = this->parallel_max_node_id();
  const dof_id_type pmax_elem_id = this->parallel_max_elem_id();
  libmesh_assert_equal_to (this->max_node_id(), pmax_node_id);
  libmesh_assert_equal_to (this->max_elem_id(), pmax_elem_id);
#endif

  _is_serial = false;
  MeshCommunication().delete_remote_elements(*this, _extra_ghost_elems);

  libmesh_assert_equal_to (this->max_elem_id(), this->parallel_max_elem_id());

  // Now make sure the containers actually shrink - strip
  // any newly-created NULL voids out of the element array
  mapvector<Elem*,dof_id_type>::veclike_iterator e_it        = _elements.begin();
  const mapvector<Elem*,dof_id_type>::veclike_iterator e_end = _elements.end();
  for (; e_it != e_end;)
    if (!*e_it)
      _elements.erase(e_it++);
    else
      ++e_it;

  mapvector<Node*,dof_id_type>::veclike_iterator n_it        = _nodes.begin();
  const mapvector<Node*,dof_id_type>::veclike_iterator n_end = _nodes.end();
  for (; n_it != n_end;)
    if (!*n_it)
      _nodes.erase(n_it++);
    else
      ++n_it;

  // We may have deleted no-longer-connected nodes or coarsened-away
  // elements; let's update our caches.
  this->update_parallel_id_counts();

#ifdef DEBUG
  // We might not have well-packed objects if the user didn't allow us
  // to renumber
  // libmesh_assert_equal_to (this->n_nodes(), this->max_node_id());
  // libmesh_assert_equal_to (this->n_elem(), this->max_elem_id());

// Make sure our neighbor links are all fine
  MeshTools::libmesh_assert_valid_neighbors(*this);

// And our child/parent links, and our flags
  MeshTools::libmesh_assert_valid_refinement_tree(*this);

// Make sure our ids and flags are consistent
  this->libmesh_assert_valid_parallel_ids();
  this->libmesh_assert_valid_parallel_flags();
#endif
}


void ParallelMesh::add_extra_ghost_elem(Elem* e)
{
  // First add the elem like normal
  add_elem(e);

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
  libmesh_assert_equal_to (this->n_nodes(), this->parallel_n_nodes());
  libmesh_assert_equal_to (this->n_elem(), this->parallel_n_elem());
  const dof_id_type pmax_node_id = this->parallel_max_node_id();
  const dof_id_type pmax_elem_id = this->parallel_max_elem_id();
  libmesh_assert_equal_to (this->max_node_id(), pmax_node_id);
  libmesh_assert_equal_to (this->max_elem_id(), pmax_elem_id);

  // If we've disabled renumbering we can't be sure we're contiguous
  // libmesh_assert_equal_to (this->n_nodes(), this->max_node_id());
  // libmesh_assert_equal_to (this->n_elem(), this->max_elem_id());

// Make sure our neighbor links are all fine
  MeshTools::libmesh_assert_valid_neighbors(*this);

// Make sure our ids and flags are consistent
  this->libmesh_assert_valid_parallel_ids();
  this->libmesh_assert_valid_parallel_flags();
#endif
}

#ifdef LIBMESH_ENABLE_UNIQUE_ID
void ParallelMesh::assign_unique_ids()
{
  {
    elem_iterator_imp        it = _elements.begin();
    const elem_iterator_imp end = _elements.end();

    for (; it != end; ++it)
      if ((*it) && ! (*it)->valid_unique_id() && libMesh::processor_id() == (*it)->processor_id())
      {
        (*it)->set_unique_id() = _next_unique_id;
        _next_unique_id += this->n_processors();
      }
  }

  {
    node_iterator_imp it  = _nodes.begin();
    node_iterator_imp end = _nodes.end();

    for (; it != end; ++it)
      if ((*it) && ! (*it)->valid_unique_id() && libMesh::processor_id() == (*it)->processor_id())
      {
        (*it)->set_unique_id() = _next_unique_id;
        _next_unique_id += this->n_processors();
      }
  }
}
#endif


} // namespace libMesh

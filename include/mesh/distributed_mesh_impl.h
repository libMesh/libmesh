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



// Local includes
#include "libmesh/distributed_mesh.h"

// libMesh includes
#include "libmesh/boundary_info.h"
#include "libmesh/elem.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/mesh_communication.h"
#include "libmesh/parmetis_partitioner.h"

// TIMPI includes
#include "timpi/parallel_implementation.h"
#include "timpi/parallel_sync.h"


namespace libMesh
{

// ------------------------------------------------------------
// DistributedMesh class member functions
template <typename RealType>
DistributedMeshTempl<RealType>::DistributedMeshTempl (const Parallel::Communicator & comm_in,
                                  unsigned char d) :
  UnstructuredMesh (comm_in,d), _is_serial(true),
  _is_serial_on_proc_0(true),
  _n_nodes(0), _n_elem(0), _max_node_id(0), _max_elem_id(0),
  _next_free_local_node_id(this->processor_id()),
  _next_free_local_elem_id(this->processor_id()),
  _next_free_unpartitioned_node_id(this->n_processors()),
  _next_free_unpartitioned_elem_id(this->n_processors())
#ifdef LIBMESH_ENABLE_UNIQUE_ID
  , _next_unpartitioned_unique_id(this->n_processors())
#endif
{
#ifdef LIBMESH_ENABLE_UNIQUE_ID
  this->_next_unique_id = this->processor_id();
#endif

  // FIXME: give parmetis the communicator!
  this->_partitioner = libmesh_make_unique<ParmetisPartitioner>();
}



template <typename RealType>
DistributedMeshTempl<RealType>::~DistributedMeshTempl ()
{
  this->clear();  // Free nodes and elements
}


// This might be specialized later, but right now it's just here to
// make sure the compiler doesn't give us a default (non-deep) copy
// constructor instead.
template <typename RealType>
DistributedMeshTempl<RealType>::DistributedMeshTempl (const DistributedMesh & other_mesh) :
  UnstructuredMesh (other_mesh), _is_serial(other_mesh._is_serial),
  _is_serial_on_proc_0(other_mesh._is_serial_on_proc_0),
  _n_nodes(0), _n_elem(0), _max_node_id(0), _max_elem_id(0),
  _next_free_local_node_id(this->processor_id()),
  _next_free_local_elem_id(this->processor_id()),
  _next_free_unpartitioned_node_id(this->n_processors()),
  _next_free_unpartitioned_elem_id(this->n_processors())
{
  this->copy_nodes_and_elements(other_mesh, true);
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
#ifdef LIBMESH_ENABLE_UNIQUE_ID
  this->_next_unique_id =
    other_mesh._next_unique_id;
  _next_unpartitioned_unique_id =
    other_mesh._next_unpartitioned_unique_id;
#endif

  auto & this_boundary_info = this->get_boundary_info();
  const auto & other_boundary_info = other_mesh.get_boundary_info();

  this_boundary_info = other_boundary_info;

  this->set_subdomain_name_map() = other_mesh.get_subdomain_name_map();

  // Use the first BoundaryInfo object to build the list of side boundary ids
  std::vector<boundary_id_type> side_boundaries;
  other_boundary_info.build_side_boundary_ids(side_boundaries);

  // Assign those boundary ids in our BoundaryInfo object
  for (const auto & side_bnd_id : side_boundaries)
    this_boundary_info.sideset_name(side_bnd_id) =
      other_boundary_info.get_sideset_name(side_bnd_id);

  // Do the same thing for node boundary ids
  std::vector<boundary_id_type> node_boundaries;
  other_boundary_info.build_node_boundary_ids(node_boundaries);

  for (const auto & node_bnd_id : node_boundaries)
    this_boundary_info.nodeset_name(node_bnd_id) =
      other_boundary_info.get_nodeset_name(node_bnd_id);

  // Need to copy extra_ghost_elems
  for (auto & elem : other_mesh._extra_ghost_elems)
    _extra_ghost_elems.insert(this->elem_ptr(elem->id()));
}



template <typename RealType>
DistributedMeshTempl<RealType>::DistributedMeshTempl (const UnstructuredMesh & other_mesh) :
  UnstructuredMesh (other_mesh), _is_serial(other_mesh.is_serial()),
  _is_serial_on_proc_0(other_mesh.is_serial()),
  _n_nodes(0), _n_elem(0), _max_node_id(0), _max_elem_id(0),
  _next_free_local_node_id(this->processor_id()),
  _next_free_local_elem_id(this->processor_id()),
  _next_free_unpartitioned_node_id(this->n_processors()),
  _next_free_unpartitioned_elem_id(this->n_processors())
{
  this->copy_nodes_and_elements(other_mesh, true);

  auto & this_boundary_info = this->get_boundary_info();
  const auto & other_boundary_info = other_mesh.get_boundary_info();

  this_boundary_info = other_boundary_info;

  this->set_subdomain_name_map() = other_mesh.get_subdomain_name_map();

  // Use the first BoundaryInfo object to build the list of side boundary ids
  std::vector<boundary_id_type> side_boundaries;
  other_boundary_info.build_side_boundary_ids(side_boundaries);

  // Assign those boundary ids in our BoundaryInfo object
  for (const auto & side_bnd_id : side_boundaries)
    this_boundary_info.sideset_name(side_bnd_id) =
      other_boundary_info.get_sideset_name(side_bnd_id);

  // Do the same thing for node boundary ids
  std::vector<boundary_id_type> node_boundaries;
  other_boundary_info.build_node_boundary_ids(node_boundaries);

  for (const auto & node_bnd_id : node_boundaries)
    this_boundary_info.nodeset_name(node_bnd_id) =
      other_boundary_info.get_nodeset_name(node_bnd_id);

#ifdef LIBMESH_ENABLE_UNIQUE_ID
  this->_next_unique_id = other_mesh.parallel_max_unique_id();
#endif
  this->update_parallel_id_counts();
}


// We use cached values for these so they can be called
// from one processor without bothering the rest, but
// we may need to update those caches before doing a full
// renumbering
template <typename RealType>
void DistributedMeshTempl<RealType>::update_parallel_id_counts()
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

#ifdef LIBMESH_ENABLE_UNIQUE_ID
  this->_next_unique_id = this->parallel_max_unique_id();
  _next_unpartitioned_unique_id =
    ((this->_next_unique_id-1) / (this->n_processors() + 1) + 1) *
    (this->n_processors() + 1) + this->n_processors();
  this->_next_unique_id =
    ((this->_next_unique_id + this->n_processors() - 1) / (this->n_processors() + 1) + 1) *
    (this->n_processors() + 1) + this->processor_id();
#endif
}


// Or in debug mode we may want to test the uncached values without
// changing the cache
template <typename RealType>
dof_id_type DistributedMeshTempl<RealType>::parallel_n_elem() const
{
  // This function must be run on all processors at once
  parallel_object_only();

  dof_id_type n_local = this->n_local_elem();
  this->comm().sum(n_local);
  n_local += this->n_unpartitioned_elem();
  return n_local;
}



template <typename RealType>
dof_id_type DistributedMeshTempl<RealType>::parallel_max_elem_id() const
{
  // This function must be run on all processors at once
  parallel_object_only();

  dof_id_type max_local = 0;

  typename mapvector<Elem *,dof_id_type>::maptype::const_reverse_iterator
    rit = _elements.rbegin();

  const typename mapvector<Elem *,dof_id_type>::maptype::const_reverse_iterator
    rend = _elements.rend();

  // Look for the maximum element id.  Search backwards through
  // elements so we can break out early.  Beware of nullptr entries that
  // haven't yet been cleared from _elements.
  for (; rit != rend; ++rit)
    if (rit->second)
      {
        libmesh_assert_equal_to(rit->second->id(), rit->first);
        max_local = rit->first + 1;
        break;
      }

  this->comm().max(max_local);
  return max_local;
}



#ifdef LIBMESH_ENABLE_UNIQUE_ID
template <typename RealType>
unique_id_type DistributedMeshTempl<RealType>::parallel_max_unique_id() const
{
  // This function must be run on all processors at once
  parallel_object_only();

  unique_id_type max_local = std::max(this->_next_unique_id,
                                      _next_unpartitioned_unique_id);
  this->comm().max(max_local);
  return max_local;
}
#endif



template <typename RealType>
dof_id_type DistributedMeshTempl<RealType>::parallel_n_nodes() const
{
  // This function must be run on all processors at once
  parallel_object_only();

  dof_id_type n_local = this->n_local_nodes();
  this->comm().sum(n_local);
  n_local += this->n_unpartitioned_nodes();
  return n_local;
}



template <typename RealType>
dof_id_type DistributedMeshTempl<RealType>::parallel_max_node_id() const
{
  // This function must be run on all processors at once
  parallel_object_only();

  dof_id_type max_local = 0;

  typename mapvector<Node *,dof_id_type>::maptype::const_reverse_iterator
    rit = _nodes.rbegin();

  const typename mapvector<Node *,dof_id_type>::maptype::const_reverse_iterator
    rend = _nodes.rend();

  // Look for the maximum element id.  Search backwards through
  // elements so we can break out early.  Beware of nullptr entries that
  // haven't yet been cleared from _elements.
  for (; rit != rend; ++rit)
    if (rit->second)
      {
        libmesh_assert_equal_to(rit->second->id(), rit->first);
        max_local = rit->first + 1;
        break;
      }

  this->comm().max(max_local);
  return max_local;
}



template <typename RealType>
const PointTempl<RealType> & DistributedMeshTempl<RealType>::point (const dof_id_type i) const
{
  return this->node_ref(i);
}



template <typename RealType>
const NodeTempl<RealType> * DistributedMeshTempl<RealType>::node_ptr (const dof_id_type i) const
{
  libmesh_assert(_nodes[i]);
  libmesh_assert_equal_to (_nodes[i]->id(), i);

  return _nodes[i];
}




template <typename RealType>
NodeTempl<RealType> * DistributedMeshTempl<RealType>::node_ptr (const dof_id_type i)
{
  libmesh_assert(_nodes[i]);
  libmesh_assert_equal_to (_nodes[i]->id(), i);

  return _nodes[i];
}




template <typename RealType>
const NodeTempl<RealType> * DistributedMeshTempl<RealType>::query_node_ptr (const dof_id_type i) const
{
  typename std::map<dof_id_type, Node *>::const_iterator it = _nodes.find(i);
  if (it != _nodes.end().it)
    {
      const Node * n = it->second;
      libmesh_assert (!n || n->id() == i);
      return n;
    }

  return nullptr;
}




template <typename RealType>
NodeTempl<RealType> * DistributedMeshTempl<RealType>::query_node_ptr (const dof_id_type i)
{
  typename std::map<dof_id_type, Node *>::const_iterator it = _nodes.find(i);
  if (it != _nodes.end().it)
    {
      Node * n = it->second;
      libmesh_assert (!n || n->id() == i);
      return n;
    }

  return nullptr;
}




template <typename RealType>
const ElemTempl<RealType> * DistributedMeshTempl<RealType>::elem_ptr (const dof_id_type i) const
{
  libmesh_assert(_elements[i]);
  libmesh_assert_equal_to (_elements[i]->id(), i);

  return _elements[i];
}




template <typename RealType>
ElemTempl<RealType> * DistributedMeshTempl<RealType>::elem_ptr (const dof_id_type i)
{
  libmesh_assert(_elements[i]);
  libmesh_assert_equal_to (_elements[i]->id(), i);

  return _elements[i];
}




template <typename RealType>
const ElemTempl<RealType> * DistributedMeshTempl<RealType>::query_elem_ptr (const dof_id_type i) const
{
  typename std::map<dof_id_type, Elem *>::const_iterator it = _elements.find(i);
  if (it != _elements.end().it)
    {
      const Elem * e = it->second;
      libmesh_assert (!e || e->id() == i);
      return e;
    }

  return nullptr;
}




template <typename RealType>
ElemTempl<RealType> * DistributedMeshTempl<RealType>::query_elem_ptr (const dof_id_type i)
{
  typename std::map<dof_id_type, Elem *>::const_iterator it = _elements.find(i);
  if (it != _elements.end().it)
    {
      Elem * e = _elements[i];
      libmesh_assert (!e || e->id() == i);
      return e;
    }

  return nullptr;
}




template <typename RealType>
ElemTempl<RealType> * DistributedMeshTempl<RealType>::add_elem (Elem * e)
{
  // Don't try to add nullptrs!
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

      // Use the unpartitioned ids for unpartitioned elems, and
      // temporarily for ghost elems
      dof_id_type * next_id = &_next_free_unpartitioned_elem_id;
      if (elem_procid == this->processor_id())
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
    // nullptr entries when testing for non-nullptr ones
    const mapvector<Elem *, dof_id_type> & const_elements = _elements;
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

#ifdef LIBMESH_ENABLE_UNIQUE_ID
  if (!e->valid_unique_id())
    {
      if (this->processor_id() == e->processor_id())
        {
          e->set_unique_id() = this->_next_unique_id;
          this->_next_unique_id += this->n_processors() + 1;
        }
      else
        {
          e->set_unique_id() = _next_unpartitioned_unique_id;
          _next_unpartitioned_unique_id += this->n_processors() + 1;
        }
    }
#endif

  // Unpartitioned elems should be added on every processor
  // And shouldn't be added in the same batch as ghost elems
  // But we might be just adding on processor 0 to
  // broadcast later
  // #ifdef DEBUG
  //   if (elem_procid == DofObject::invalid_processor_id)
  //     {
  //       dof_id_type elem_id = e->id();
  //       this->comm().max(elem_id);
  //       libmesh_assert_equal_to (elem_id, e->id());
  //     }
  // #endif

  // Make sure any new element is given space for any extra integers
  // we've requested
  e->add_extra_integers(this->_elem_integer_names.size());

  // And set mapping type and data on any new element
  e->set_mapping_type(this->default_mapping_type());
  e->set_mapping_data(this->default_mapping_data());

  return e;
}



template <typename RealType>
ElemTempl<RealType> * DistributedMeshTempl<RealType>::insert_elem (Elem * e)
{
  if (_elements[e->id()])
    this->delete_elem(_elements[e->id()]);

#ifdef LIBMESH_ENABLE_UNIQUE_ID
  if (!e->valid_unique_id())
    {
      if (this->processor_id() == e->processor_id())
        {
          e->set_unique_id() = this->_next_unique_id;
          this->_next_unique_id += this->n_processors() + 1;
        }
      else
        {
          e->set_unique_id() = _next_unpartitioned_unique_id;
          _next_unpartitioned_unique_id += this->n_processors() + 1;
        }
    }
#endif

  // Try to make the cached elem data more accurate
  processor_id_type elem_procid = e->processor_id();
  if (elem_procid == this->processor_id() ||
      elem_procid == DofObject::invalid_processor_id)
    _n_elem++;

  _elements[e->id()] = e;

  // Make sure any new element is given space for any extra integers
  // we've requested
  e->add_extra_integers(this->_elem_integer_names.size());

  // And set mapping type and data on any new element
  e->set_mapping_type(this->default_mapping_type());
  e->set_mapping_data(this->default_mapping_data());

  return e;
}



template <typename RealType>
void DistributedMeshTempl<RealType>::delete_elem(Elem * e)
{
  libmesh_assert (e);

  // Try to make the cached elem data more accurate
  processor_id_type elem_procid = e->processor_id();
  if (elem_procid == this->processor_id() ||
      elem_procid == DofObject::invalid_processor_id)
    _n_elem--;

  // Delete the element from the BoundaryInfo object
  this->get_boundary_info().remove(e);

  // But not yet from the container; we might invalidate
  // an iterator that way!

  //_elements.erase(e->id());

  // Instead, we set it to nullptr for now

  _elements[e->id()] = nullptr;

  // delete the element
  delete e;
}



template <typename RealType>
void DistributedMeshTempl<RealType>::renumber_elem(const dof_id_type old_id,
                                    const dof_id_type new_id)
{
  Elem * el = _elements[old_id];
  libmesh_assert (el);
  libmesh_assert_equal_to (el->id(), old_id);

  el->set_id(new_id);
  libmesh_assert (!_elements[new_id]);
  _elements[new_id] = el;
  _elements.erase(old_id);
}



template <typename RealType>
NodeTempl<RealType> * DistributedMeshTempl<RealType>::add_point (const Point & p,
                                   const dof_id_type id,
                                   const processor_id_type proc_id)
{
  auto n_it = _nodes.find(id);
  if (n_it != _nodes.end().it)
    {
      Node * n = n_it->second;
      libmesh_assert (n);
      libmesh_assert_equal_to (n->id(), id);

      *n = p;
      n->processor_id() = proc_id;

      return n;
    }

  Node * n = Node::build(p, id).release();
  n->processor_id() = proc_id;

  return DistributedMesh::add_node(n);
}


template <typename RealType>
void DistributedMeshTempl<RealType>::own_node (Node & n)
{
  // This had better be a node in our mesh
  libmesh_assert(_nodes[n.id()] == &n);

  _nodes[n.id()] = nullptr;

  n.set_id(DofObject::invalid_id);
  n.processor_id() = this->processor_id();

  this->add_node(&n);
}


template <typename RealType>
NodeTempl<RealType> * DistributedMeshTempl<RealType>::add_node (Node * n)
{
  // Don't try to add nullptrs!
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
      // and temporarily for ghost nodes
      dof_id_type * next_id = &_next_free_unpartitioned_node_id;
      if (node_procid == this->processor_id())
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
    // nullptr entries when testing for non-nullptr ones
    const mapvector<Node *,dof_id_type> & const_nodes = _nodes;
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

#ifdef LIBMESH_ENABLE_UNIQUE_ID
  if (!n->valid_unique_id())
    {
      if (this->processor_id() == n->processor_id())
        {
          n->set_unique_id() = this->_next_unique_id;
          this->_next_unique_id += this->n_processors();
        }
      else
        {
          n->set_unique_id() = _next_unpartitioned_unique_id;
          _next_unpartitioned_unique_id += this->n_processors() + 1;
        }
    }
#endif

  n->add_extra_integers(this->_node_integer_names.size());

  // Unpartitioned nodes should be added on every processor
  // And shouldn't be added in the same batch as ghost nodes
  // But we might be just adding on processor 0 to
  // broadcast later
  // #ifdef DEBUG
  //   if (node_procid == DofObject::invalid_processor_id)
  //     {
  //       dof_id_type node_id = n->id();
  //       this->comm().max(node_id);
  //       libmesh_assert_equal_to (node_id, n->id());
  //     }
  // #endif

  return n;
}



template <typename RealType>
NodeTempl<RealType> * DistributedMeshTempl<RealType>::insert_node(Node * n)
{
  return DistributedMesh::add_node(n);
}



template <typename RealType>
void DistributedMeshTempl<RealType>::delete_node(Node * n)
{
  libmesh_assert(n);
  libmesh_assert(_nodes[n->id()]);

  // Try to make the cached elem data more accurate
  processor_id_type node_procid = n->processor_id();
  if (node_procid == this->processor_id() ||
      node_procid == DofObject::invalid_processor_id)
    _n_nodes--;

  // Delete the node from the BoundaryInfo object
  this->get_boundary_info().remove(n);

  // But not yet from the container; we might invalidate
  // an iterator that way!

  //_nodes.erase(n->id());

  // Instead, we set it to nullptr for now

  _nodes[n->id()] = nullptr;

  // delete the node
  delete n;
}



template <typename RealType>
void DistributedMeshTempl<RealType>::renumber_node(const dof_id_type old_id,
                                    const dof_id_type new_id)
{
  Node * nd = _nodes[old_id];
  libmesh_assert (nd);
  libmesh_assert_equal_to (nd->id(), old_id);

  nd->set_id(new_id);

  // If we have nodes shipped to this processor for NodeConstraints
  // use, then those nodes will exist in _nodes, but may not be
  // locatable via a TopologyMap due to the insufficiency of elements
  // connecting to them.  If local refinement then wants to create a
  // *new* node in the same location, it will initially get a temporary
  // id, and then make_node_ids_parallel_consistent() will try to move
  // it to the canonical id.  We need to account for this case to
  // avoid false positives and memory leaks.
#ifdef LIBMESH_ENABLE_NODE_CONSTRAINTS
  if (_nodes[new_id])
    {
      libmesh_assert_equal_to (*(Point *)_nodes[new_id],
                               *(Point *)_nodes[old_id]);
      _nodes.erase(new_id);
    }
#else
  // If we aren't shipping nodes for NodeConstraints, there should be
  // no reason for renumbering one node onto another.
  libmesh_assert (!_nodes[new_id]);
#endif
  _nodes[new_id] = nd;
  _nodes.erase(old_id);
}



template <typename RealType>
void DistributedMeshTempl<RealType>::clear ()
{
  // Call parent clear function
  MeshBase::clear();

  // Clear our elements and nodes
  // There is no need to remove the elements from
  // the BoundaryInfo data structure since we
  // already cleared it.
  for (auto & elem : _elements)
    delete elem;

  // clear the nodes data structure
  // There is no need to remove the nodes from
  // the BoundaryInfo data structure since we
  // already cleared it.
  for (auto & node : _nodes)
    delete node;

  _elements.clear();
  _nodes.clear();

  // We're no longer distributed if we were before
  _is_serial = true;
  _is_serial_on_proc_0 = true;

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



template <typename RealType>
void DistributedMeshTempl<RealType>::redistribute ()
{
  // If this is a truly parallel mesh, go through the redistribution/gather/delete remote steps
  if (!this->is_serial())
    {
      // Construct a MeshCommunication object to actually redistribute the nodes
      // and elements according to the partitioner, and then to re-gather the neighbors.
      MeshCommunication mc;
      mc.redistribute(*this);

      this->update_parallel_id_counts();

      // We probably had valid neighbors previously, so that a quality
      // new partitioning could be found, but we might not have valid
      // neighbor links on the newly-redistributed elements
      this->find_neighbors();

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



template <typename RealType>
void DistributedMeshTempl<RealType>::update_post_partitioning ()
{
  // this->recalculate_n_partitions();

  // Partitioning changes our numbers of unpartitioned objects
  this->update_parallel_id_counts();
}



template <typename RealType>
template <typename T>
void DistributedMeshTempl<RealType>::libmesh_assert_valid_parallel_object_ids(const mapvector<T *, dof_id_type> & objects) const
{
  // This function must be run on all processors at once
  parallel_object_only();

  const dof_id_type pmax_node_id = this->parallel_max_node_id();
  const dof_id_type pmax_elem_id = this->parallel_max_elem_id();
  const dof_id_type pmax_id = std::max(pmax_node_id, pmax_elem_id);

  for (dof_id_type i=0; i != pmax_id; ++i)
    {
      T * obj = objects[i]; // Returns nullptr if there's no map entry

      // Local lookups by id should return the requested object
      libmesh_assert(!obj || obj->id() == i);

      // All processors with an object should agree on id
#ifndef NDEBUG
      const dof_id_type dofid = obj && obj->valid_id() ?
        obj->id() : DofObject::invalid_id;
      libmesh_assert(this->comm().semiverify(obj ? &dofid : nullptr));
#endif

      // All processors with an object should agree on processor id
      const dof_id_type procid = obj && obj->valid_processor_id() ?
        obj->processor_id() : DofObject::invalid_processor_id;
      libmesh_assert(this->comm().semiverify(obj ? &procid : nullptr));

      dof_id_type min_procid = procid;
      this->comm().min(min_procid);

      // Either:
      // 1.) I own this elem (min_procid == this->processor_id()) *and* I have a valid pointer to it (obj != nullptr)
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

#if defined(LIBMESH_ENABLE_UNIQUE_ID) && !defined(NDEBUG)
      // All processors with an object should agree on unique id
      const unique_id_type uniqueid = obj ? obj->unique_id() : 0;
      libmesh_assert(this->comm().semiverify(obj ? &uniqueid : nullptr));
#endif
    }
}



template <typename RealType>
void DistributedMeshTempl<RealType>::libmesh_assert_valid_parallel_ids () const
{
  this->libmesh_assert_valid_parallel_object_ids (this->_elements);
  this->libmesh_assert_valid_parallel_object_ids (this->_nodes);
}



template <typename RealType>
void DistributedMeshTempl<RealType>::libmesh_assert_valid_parallel_p_levels () const
{
#ifndef NDEBUG
  // This function must be run on all processors at once
  parallel_object_only();

  dof_id_type pmax_elem_id = this->parallel_max_elem_id();

  for (dof_id_type i=0; i != pmax_elem_id; ++i)
    {
      Elem * el = _elements[i]; // Returns nullptr if there's no map entry

      unsigned int p_level = el ?  (el->p_level()) : libMesh::invalid_uint;

      // All processors with an active element should agree on p level
      libmesh_assert(this->comm().semiverify((el && el->active()) ? &p_level : nullptr));
    }
#endif
}




template <typename RealType>
void DistributedMeshTempl<RealType>::libmesh_assert_valid_parallel_flags () const
{
#if defined(LIBMESH_ENABLE_AMR) && !defined(NDEBUG)
  // This function must be run on all processors at once
  parallel_object_only();

  dof_id_type pmax_elem_id = this->parallel_max_elem_id();

  for (dof_id_type i=0; i != pmax_elem_id; ++i)
    {
      Elem * el = _elements[i]; // Returns nullptr if there's no map entry

      unsigned int refinement_flag   = el ?
        static_cast<unsigned int> (el->refinement_flag()) : libMesh::invalid_uint;
      unsigned int p_refinement_flag = el ?
        static_cast<unsigned int> (el->p_refinement_flag()) : libMesh::invalid_uint;

      libmesh_assert(this->comm().semiverify(el ? &refinement_flag : nullptr));

      // p refinement flags aren't always kept correct on inactive
      // ghost elements
      libmesh_assert(this->comm().semiverify((el && el->active()) ? &p_refinement_flag : nullptr));
    }
#endif // LIBMESH_ENABLE_AMR
}



template <typename RealType>
template <typename T>
dof_id_type
DistributedMeshTempl<RealType>::renumber_dof_objects(mapvector<T *, dof_id_type> & objects)
{
  // This function must be run on all processors at once
  parallel_object_only();

  typedef typename mapvector<T *,dof_id_type>::veclike_iterator object_iterator;

  // In parallel we may not know what objects other processors have.
  // Start by figuring out how many
  dof_id_type unpartitioned_objects = 0;

  std::unordered_map<processor_id_type, dof_id_type>
    ghost_objects_from_proc;

  object_iterator it  = objects.begin();
  object_iterator end = objects.end();

  while (it != end)
    {
      T * obj = *it;

      // Remove any nullptr container entries while we're here.
      if (!obj)
        it = objects.erase(it);
      else
        {
          processor_id_type obj_procid = obj->processor_id();
          if (obj_procid == DofObject::invalid_processor_id)
            unpartitioned_objects++;
          else
            ghost_objects_from_proc[obj_procid]++;

          // Finally, increment the iterator
          ++it;
        }
    }

  std::vector<dof_id_type> objects_on_proc(this->n_processors(), 0);
  auto this_it = ghost_objects_from_proc.find(this->processor_id());
  this->comm().allgather
    ((this_it == ghost_objects_from_proc.end()) ?
     dof_id_type(0) : this_it->second, objects_on_proc);

#ifndef NDEBUG
  libmesh_assert(this->comm().verify(unpartitioned_objects));
  for (processor_id_type p=0, np=this->n_processors(); p != np; ++p)
    if (ghost_objects_from_proc.count(p))
      libmesh_assert_less_equal (ghost_objects_from_proc[p], objects_on_proc[p]);
    else
      libmesh_assert_less_equal (0, objects_on_proc[p]);
#endif

  // We'll renumber objects in blocks by processor id
  std::vector<dof_id_type> first_object_on_proc(this->n_processors());
  for (processor_id_type i=1, np=this->n_processors(); i != np; ++i)
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
  std::map<processor_id_type, std::vector<dof_id_type>>
    requested_ids;

  // We know how many objects live on each processor, so reserve() space for
  // each.
  auto ghost_end = ghost_objects_from_proc.end();
  for (auto p : IntRange<processor_id_type>(0, this->n_processors()))
    if (p != this->processor_id())
      {
        const auto p_it = ghost_objects_from_proc.find(p);
        if (p_it != ghost_end)
          requested_ids[p].reserve(p_it->second);
      }

  end = objects.end();
  for (it = objects.begin(); it != end; ++it)
    {
      T * obj = *it;
      if (obj->processor_id() == this->processor_id())
        obj->set_id(next_id++);
      else if (obj->processor_id() != DofObject::invalid_processor_id)
        requested_ids[obj->processor_id()].push_back(obj->id());
    }

  // Next set ghost object ids from other processors

  auto gather_functor =
    [
#ifndef NDEBUG
     this,
     &first_object_on_proc,
     &objects_on_proc,
#endif
     &objects]
    (processor_id_type, const std::vector<dof_id_type> & ids,
     std::vector<dof_id_type> & new_ids)
    {
      std::size_t ids_size = ids.size();
      new_ids.resize(ids_size);

      for (std::size_t i=0; i != ids_size; ++i)
        {
          T * obj = objects[ids[i]];
          libmesh_assert(obj);
          libmesh_assert_equal_to (obj->processor_id(), this->processor_id());
          new_ids[i] = obj->id();

          libmesh_assert_greater_equal (new_ids[i],
                                        first_object_on_proc[this->processor_id()]);
          libmesh_assert_less (new_ids[i],
                               first_object_on_proc[this->processor_id()] +
                               objects_on_proc[this->processor_id()]);
        }
    };

  auto action_functor =
    [
#ifndef NDEBUG
     &first_object_on_proc,
     &objects_on_proc,
#endif
     &objects]
    (processor_id_type libmesh_dbg_var(pid),
     const std::vector<dof_id_type> & ids,
     const std::vector<dof_id_type> & data)
    {
      // Copy the id changes we've now been informed of
      for (auto i : index_range(ids))
        {
          T * obj = objects[ids[i]];
          libmesh_assert (obj);
          libmesh_assert_equal_to (obj->processor_id(), pid);
          libmesh_assert_greater_equal (data[i],
                                        first_object_on_proc[pid]);
          libmesh_assert_less (data[i],
                               first_object_on_proc[pid] +
                               objects_on_proc[pid]);
          obj->set_id(data[i]);
        }
    };

  const dof_id_type * ex = nullptr;
  Parallel::pull_parallel_vector_data
    (this->comm(), requested_ids, gather_functor, action_functor, ex);

#ifdef LIBMESH_ENABLE_UNIQUE_ID
  auto unique_gather_functor =
    [
#ifndef NDEBUG
     this,
#endif
     &objects]
    (processor_id_type, const std::vector<dof_id_type> & ids,
     std::vector<unique_id_type> & data)
    {
      std::size_t ids_size = ids.size();
      data.resize(ids_size);

      for (std::size_t i=0; i != ids_size; ++i)
        {
          T * obj = objects[ids[i]];
          libmesh_assert(obj);
          libmesh_assert_equal_to (obj->processor_id(), this->processor_id());
          data[i] = obj->valid_unique_id() ? obj->unique_id() : DofObject::invalid_unique_id;
        }
    };

  auto unique_action_functor =
    [&objects]
    (processor_id_type libmesh_dbg_var(pid),
     const std::vector<dof_id_type> & ids,
     const std::vector<unique_id_type> & data)
    {
      for (auto i : index_range(ids))
        {
          T * obj = objects[ids[i]];
          libmesh_assert (obj);
          libmesh_assert_equal_to (obj->processor_id(), pid);
          if (!obj->valid_unique_id() && data[i] != DofObject::invalid_unique_id)
            obj->set_unique_id() = (data[i]);
        }
    };

  const unique_id_type * unique_ex = nullptr;
  Parallel::pull_parallel_vector_data
    (this->comm(), requested_ids, unique_gather_functor,
     unique_action_functor, unique_ex);
#endif

  // Next set unpartitioned object ids
  next_id = 0;
  for (auto i : IntRange<processor_id_type>(0, this->n_processors()))
    next_id += objects_on_proc[i];
  for (it = objects.begin(); it != end; ++it)
    {
      T * obj = *it;
      if (obj->processor_id() == DofObject::invalid_processor_id)
        obj->set_id(next_id++);
    }

  // Finally shuffle around objects so that container indices
  // match ids
  it = objects.begin();
  end = objects.end();
  while (it != end)
    {
      T * obj = *it;
      if (obj) // don't try shuffling already-nullptr entries
        {
          T * next = objects[obj->id()];
          // If we have to move this object
          if (next != obj)
            {
              // nullptr out its original position for now
              // (our shuffling may put another object there shortly)
              *it = nullptr;

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

      // Remove any container entries that were left as nullptr.
      if (!obj)
        it = objects.erase(it);
      else
        ++it;
    }

  return first_free_id;
}


template <typename RealType>
void DistributedMeshTempl<RealType>::renumber_nodes_and_elements ()
{
  parallel_object_only();

#ifdef DEBUG
  // Make sure our ids and flags are consistent
  this->libmesh_assert_valid_parallel_ids();
  this->libmesh_assert_valid_parallel_flags();
  this->libmesh_assert_valid_parallel_p_levels();
#endif

  LOG_SCOPE("renumber_nodes_and_elements()", "DistributedMesh");

  std::set<dof_id_type> used_nodes;

  // flag the nodes we need
  for (auto & elem : this->element_ptr_range())
    for (const Node & node : elem->node_ref_range())
      used_nodes.insert(node.id());

  // Nodes not connected to any local elements, and nullptr node entries
  // in our container, are deleted
  {
    node_iterator_imp  it = _nodes.begin();
    node_iterator_imp end = _nodes.end();

    while (it != end)
      {
        Node * nd = *it;
        if (!nd)
          it = _nodes.erase(it);
        else if (!used_nodes.count(nd->id()))
          {
            // remove any boundary information associated with
            // this node
            this->get_boundary_info().remove (nd);

            // delete the node
            delete nd;

            it = _nodes.erase(it);
          }
        else
          ++it;
      }
  }

  if (this->_skip_renumber_nodes_and_elements)
    {
      this->update_parallel_id_counts();
      return;
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
}



template <typename RealType>
void DistributedMeshTempl<RealType>::fix_broken_node_and_element_numbering ()
{
  // We need access to iterators for the underlying containers,
  // not the mapvector<> reimplementations.
  typename mapvector<Node *,dof_id_type>::maptype & nodes = this->_nodes;
  typename mapvector<Elem *,dof_id_type>::maptype & elems = this->_elements;

  // Nodes first
  for (auto & pr : nodes)
    if (pr.second != nullptr)
      pr.second->set_id() = pr.first;

  // Elements next
  for (const auto & pr : elems)
    if (pr.second != nullptr)
      pr.second->set_id() = pr.first;
}



template <typename RealType>
dof_id_type DistributedMeshTempl<RealType>::n_active_elem () const
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



template <typename RealType>
void DistributedMeshTempl<RealType>::delete_remote_elements()
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
  _is_serial_on_proc_0 = false;

  MeshCommunication().delete_remote_elements(*this, _extra_ghost_elems);

  libmesh_assert_equal_to (this->max_elem_id(), this->parallel_max_elem_id());

  // Now make sure the containers actually shrink - strip
  // any newly-created nullptr voids out of the element array
  typename mapvector<Elem *,dof_id_type>::veclike_iterator e_it        = _elements.begin();
  const typename mapvector<Elem *,dof_id_type>::veclike_iterator e_end = _elements.end();
  while (e_it != e_end)
    if (!*e_it)
      e_it = _elements.erase(e_it);
    else
      ++e_it;

  typename mapvector<Node *,dof_id_type>::veclike_iterator n_it        = _nodes.begin();
  const typename mapvector<Node *,dof_id_type>::veclike_iterator n_end = _nodes.end();
  while (n_it != n_end)
    if (!*n_it)
      n_it = _nodes.erase(n_it);
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


template <typename RealType>
void DistributedMeshTempl<RealType>::add_extra_ghost_elem(Elem * e)
{
  // First add the elem like normal
  add_elem(e);

  // Now add it to the set that won't be deleted when we call
  // delete_remote_elements()
  _extra_ghost_elems.insert(e);
}

template <typename RealType>
void
DistributedMeshTempl<RealType>::clear_extra_ghost_elems(const std::set<Elem *> & extra_ghost_elems)
{
  std::set<Elem *> tmp;
  std::set_difference(_extra_ghost_elems.begin(), _extra_ghost_elems.end(),
                      extra_ghost_elems.begin(), extra_ghost_elems.end(),
                      std::inserter(tmp, tmp.begin()));
  _extra_ghost_elems = tmp;
}

template <typename RealType>
void DistributedMeshTempl<RealType>::allgather()
{
  if (_is_serial)
    return;
  MeshCommunication().allgather(*this);
  _is_serial = true;
  _is_serial_on_proc_0 = true;

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

template <typename RealType>
void DistributedMeshTempl<RealType>::gather_to_zero()
{
  if (_is_serial_on_proc_0)
    return;

  _is_serial_on_proc_0 = true;
  MeshCommunication().gather(0, *this);
}


} // namespace libMesh

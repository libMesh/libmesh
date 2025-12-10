// The libMesh Finite Element Library.
// Copyright (C) 2002-2025 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/mesh_tools.h"
#include "libmesh/partitioner.h"
#include "libmesh/string_to_enum.h"

// TIMPI includes
#include "timpi/parallel_implementation.h"
#include "timpi/parallel_sync.h"


namespace libMesh
{

// ------------------------------------------------------------
// DistributedMesh class member functions
DistributedMesh::DistributedMesh (const Parallel::Communicator & comm_in,
                                  unsigned char d) :
  UnstructuredMesh (comm_in,d), _is_serial(true),
  _is_serial_on_proc_0(true),
  _deleted_coarse_elements(false),
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
  _next_unique_id = this->processor_id();
#endif

  const std::string default_partitioner = "parmetis";
  const std::string my_partitioner =
    libMesh::command_line_value("--default-partitioner",
                                default_partitioner);
  _partitioner = Partitioner::build
    (Utility::string_to_enum<PartitionerType>(my_partitioner));
}

DistributedMesh & DistributedMesh::operator= (DistributedMesh && other_mesh)
{
  LOG_SCOPE("operator=(&&)", "DistributedMesh");

  // Move assign as an UnstructuredMesh.
  this->UnstructuredMesh::operator=(std::move(other_mesh));

  // Nodes and elements belong to DistributedMesh and have to be
  // moved before we can move arbitrary GhostingFunctor, Partitioner,
  // etc. subclasses.
  this->move_nodes_and_elements(std::move(other_mesh));

  // But move_nodes_and_elems misses (or guesses about) some of our
  // subclass values, and we want more precision than a guess.
  _deleted_coarse_elements = other_mesh._deleted_coarse_elements;
  _extra_ghost_elems = std::move(other_mesh._extra_ghost_elems);

  // Handle remaining MeshBase moves.
  this->post_dofobject_moves(std::move(other_mesh));

  return *this;
}

MeshBase & DistributedMesh::assign(MeshBase && other_mesh)
{
  *this = std::move(cast_ref<DistributedMesh&>(other_mesh));

  return *this;
}

std::string_view DistributedMesh::subclass_first_difference_from (const MeshBase & other_mesh_base) const
{
  const DistributedMesh * dist_mesh_ptr =
    dynamic_cast<const DistributedMesh *>(&other_mesh_base);
  if (!dist_mesh_ptr)
    return "DistributedMesh class";
  const DistributedMesh & other_mesh = *dist_mesh_ptr;

#define CHECK_MEMBER(member_name) \
  if (member_name != other_mesh.member_name) \
    return #member_name;

  CHECK_MEMBER(_is_serial);
  CHECK_MEMBER(_is_serial_on_proc_0);
  CHECK_MEMBER(_deleted_coarse_elements);
  CHECK_MEMBER(_n_nodes);
  CHECK_MEMBER(_n_elem);
  CHECK_MEMBER(_max_node_id);
  CHECK_MEMBER(_max_elem_id);
  // We expect these things to change in a prepare_for_use();
  // they're conceptually "mutable"...
/*
  CHECK_MEMBER(_next_free_local_node_id);
  CHECK_MEMBER(_next_free_local_elem_id);
  CHECK_MEMBER(_next_free_unpartitioned_node_id);
  CHECK_MEMBER(_next_free_unpartitioned_elem_id);
#ifdef LIBMESH_ENABLE_UNIQUE_ID
  CHECK_MEMBER(_next_unpartitioned_unique_id);
#endif
*/
  if (!this->nodes_and_elements_equal(other_mesh))
    return "nodes and/or elements";

  if (_extra_ghost_elems.size() !=
      other_mesh._extra_ghost_elems.size())
    return "_extra_ghost_elems size";
  for (auto & elem : _extra_ghost_elems)
    {
      libmesh_assert(this->query_elem_ptr(elem->id()) == elem);
      const Elem * other_elem = other_mesh.query_elem_ptr(elem->id());
      if (!other_elem ||
          !other_mesh._extra_ghost_elems.count(const_cast<Elem *>(other_elem)))
        return "_extra_ghost_elems entry";
    }

  return "";
}

DistributedMesh::~DistributedMesh ()
{
  this->DistributedMesh::clear();  // Free nodes and elements
}


// This might be specialized later, but right now it's just here to
// make sure the compiler doesn't give us a default (non-deep) copy
// constructor instead.
DistributedMesh::DistributedMesh (const DistributedMesh & other_mesh) :
  DistributedMesh(static_cast<const MeshBase &>(other_mesh))
{
  _is_serial = other_mesh._is_serial;
  _is_serial_on_proc_0 = other_mesh._is_serial_on_proc_0;
  _deleted_coarse_elements = other_mesh._deleted_coarse_elements;

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
  _next_unique_id =
    other_mesh._next_unique_id;
  _next_unpartitioned_unique_id =
    other_mesh._next_unpartitioned_unique_id;
#endif

  // Need to copy extra_ghost_elems
  for (auto & elem : other_mesh._extra_ghost_elems)
    _extra_ghost_elems.insert(this->elem_ptr(elem->id()));
}



DistributedMesh::DistributedMesh (const MeshBase & other_mesh) :
  UnstructuredMesh (other_mesh), _is_serial(other_mesh.is_serial()),
  _is_serial_on_proc_0(other_mesh.is_serial()),
  _deleted_coarse_elements(true), // better safe than sorry...
  _n_nodes(0), _n_elem(0), _max_node_id(0), _max_elem_id(0),
  _next_free_local_node_id(this->processor_id()),
  _next_free_local_elem_id(this->processor_id()),
  _next_free_unpartitioned_node_id(this->n_processors()),
  _next_free_unpartitioned_elem_id(this->n_processors())
{
  // Just copy, skipping preparation
  this->copy_nodes_and_elements(other_mesh, true, 0, 0, 0, nullptr, true);

  this->allow_find_neighbors(other_mesh.allow_find_neighbors());
  this->allow_renumbering(other_mesh.allow_renumbering());
  this->allow_remote_element_removal(other_mesh.allow_remote_element_removal());
  this->skip_partitioning(other_mesh.skip_partitioning());

  this->copy_constraint_rows(other_mesh);

  this->_preparation = other_mesh.preparation();

  auto & this_boundary_info = this->get_boundary_info();
  const auto & other_boundary_info = other_mesh.get_boundary_info();

  this_boundary_info = other_boundary_info;

  this->set_subdomain_name_map() = other_mesh.get_subdomain_name_map();

#ifdef LIBMESH_ENABLE_UNIQUE_ID
  _next_unique_id = other_mesh.parallel_max_unique_id() +
                    this->processor_id();
  _next_unpartitioned_unique_id = _next_unique_id +
    (this->n_processors() - this->processor_id());
#endif
  this->update_parallel_id_counts();
}

void DistributedMesh::move_nodes_and_elements(MeshBase && other_meshbase)
{
  DistributedMesh & other_mesh = cast_ref<DistributedMesh&>(other_meshbase);

  this->_nodes = std::move(other_mesh._nodes);
  this->_n_nodes = other_mesh.n_nodes();

  this->_elements = std::move(other_mesh._elements);
  this->_n_elem = other_mesh.n_elem();

  _is_serial = other_mesh.is_serial();
  _is_serial_on_proc_0 = other_mesh.is_serial_on_zero();
  _deleted_coarse_elements = true; // Better safe than sorry

  _max_node_id = other_mesh.max_node_id();
  _max_elem_id = other_mesh.max_elem_id();

  _next_free_local_node_id = other_mesh._next_free_local_node_id;
  _next_free_local_elem_id = other_mesh._next_free_local_elem_id;
  _next_free_unpartitioned_node_id = other_mesh._next_free_unpartitioned_node_id;
  _next_free_unpartitioned_elem_id = other_mesh._next_free_unpartitioned_elem_id;

  #ifdef LIBMESH_ENABLE_UNIQUE_ID
  _next_unpartitioned_unique_id = other_mesh._next_unpartitioned_unique_id;
  #endif
}

// We use cached values for these so they can be called
// from one processor without bothering the rest, but
// we may need to update those caches before doing a full
// renumbering
void DistributedMesh::update_parallel_id_counts()
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
  _next_unique_id = this->parallel_max_unique_id();
  _next_unpartitioned_unique_id =
    ((_next_unique_id-1) / (this->n_processors() + 1) + 1) *
    (this->n_processors() + 1) + this->n_processors();
  _next_unique_id =
    ((_next_unique_id + this->n_processors() - 1) / (this->n_processors() + 1) + 1) *
    (this->n_processors() + 1) + this->processor_id();
#endif

  this->_preparation.has_synched_id_counts = true;
}


// Or in debug mode we may want to test the uncached values without
// changing the cache
dof_id_type DistributedMesh::parallel_n_elem() const
{
  // This function must be run on all processors at once
  parallel_object_only();

  dof_id_type n_local = this->n_local_elem();
  this->comm().sum(n_local);
  n_local += this->n_unpartitioned_elem();
  return n_local;
}



dof_id_type DistributedMesh::parallel_max_elem_id() const
{
  // This function must be run on all processors at once
  parallel_object_only();

  dof_id_type max_local = 0;

  dofobject_container<Elem>::const_reverse_veclike_iterator
    rit = _elements.rbegin();

  const dofobject_container<Elem>::const_reverse_veclike_iterator
    rend = _elements.rend();

  // Look for the maximum element id.  Search backwards through
  // elements so we can break out early.  Beware of nullptr entries that
  // haven't yet been cleared from _elements.
  for (; rit != rend; ++rit)
    {
      const DofObject *d = *rit;
      if (d)
        {
          libmesh_assert(_elements[d->id()] == d);
          max_local = d->id() + 1;
          break;
        }
    }

  this->comm().max(max_local);
  return max_local;
}



#ifdef LIBMESH_ENABLE_UNIQUE_ID
unique_id_type DistributedMesh::parallel_max_unique_id() const
{
  // This function must be run on all processors at once
  parallel_object_only();

  unique_id_type max_local = std::max(_next_unique_id,
                                      _next_unpartitioned_unique_id);
  this->comm().max(max_local);
  return max_local;
}



void DistributedMesh::set_next_unique_id(unique_id_type id)
{
  _next_unique_id = id;
  _next_unpartitioned_unique_id =
    ((_next_unique_id-1) / (this->n_processors() + 1) + 1) *
    (this->n_processors() + 1) + this->n_processors();
  _next_unique_id =
    ((_next_unique_id + this->n_processors() - 1) / (this->n_processors() + 1) + 1) *
    (this->n_processors() + 1) + this->processor_id();
}
#endif



dof_id_type DistributedMesh::parallel_n_nodes() const
{
  // This function must be run on all processors at once
  parallel_object_only();

  dof_id_type n_local = this->n_local_nodes();
  this->comm().sum(n_local);
  n_local += this->n_unpartitioned_nodes();
  return n_local;
}



dof_id_type DistributedMesh::parallel_max_node_id() const
{
  // This function must be run on all processors at once
  parallel_object_only();

  dof_id_type max_local = 0;

  dofobject_container<Node>::const_reverse_veclike_iterator
    rit = _nodes.rbegin();

  const dofobject_container<Node>::const_reverse_veclike_iterator
    rend = _nodes.rend();

  // Look for the maximum node id.  Search backwards through
  // nodes so we can break out early.  Beware of nullptr entries that
  // haven't yet been cleared from _nodes
  for (; rit != rend; ++rit)
    {
      const DofObject *d = *rit;
      if (d)
        {
          libmesh_assert(_nodes[d->id()] == d);
          max_local = d->id() + 1;
          break;
        }
    }

  this->comm().max(max_local);
  return max_local;
}



const Point & DistributedMesh::point (const dof_id_type i) const
{
  return this->node_ref(i);
}



const Node * DistributedMesh::node_ptr (const dof_id_type i) const
{
  libmesh_assert(_nodes[i]);
  libmesh_assert_equal_to (_nodes[i]->id(), i);

  return _nodes[i];
}




Node * DistributedMesh::node_ptr (const dof_id_type i)
{
  libmesh_assert(_nodes[i]);
  libmesh_assert_equal_to (_nodes[i]->id(), i);

  return _nodes[i];
}




const Node * DistributedMesh::query_node_ptr (const dof_id_type i) const
{
  if (const auto it = _nodes.find(i);
      it != _nodes.end())
    {
      const Node * n = *it;
      libmesh_assert (!n || n->id() == i);
      return n;
    }

  return nullptr;
}




Node * DistributedMesh::query_node_ptr (const dof_id_type i)
{
  if (auto it = _nodes.find(i);
      it != _nodes.end())
    {
      Node * n = *it;
      libmesh_assert (!n || n->id() == i);
      return n;
    }

  return nullptr;
}




const Elem * DistributedMesh::elem_ptr (const dof_id_type i) const
{
  libmesh_assert(_elements[i]);
  libmesh_assert_equal_to (_elements[i]->id(), i);

  return _elements[i];
}




Elem * DistributedMesh::elem_ptr (const dof_id_type i)
{
  libmesh_assert(_elements[i]);
  libmesh_assert_equal_to (_elements[i]->id(), i);

  return _elements[i];
}




const Elem * DistributedMesh::query_elem_ptr (const dof_id_type i) const
{
  if (const auto it = _elements.find(i);
      it != _elements.end())
    {
      const Elem * e = *it;
      libmesh_assert (!e || e->id() == i);
      return e;
    }

  return nullptr;
}




Elem * DistributedMesh::query_elem_ptr (const dof_id_type i)
{
  if (auto it = _elements.find(i);
      it != _elements.end())
    {
      Elem * e = *it;
      libmesh_assert (!e || e->id() == i);
      return e;
    }

  return nullptr;
}




Elem * DistributedMesh::add_elem (Elem * e)
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
    // We need a const dofobject_container so we don't inadvertently create
    // nullptr entries when testing for non-nullptr ones
    const dofobject_container<Elem> & const_elements = _elements;
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
      if (processor_id() == e->processor_id())
        {
          e->set_unique_id(_next_unique_id);
          _next_unique_id += this->n_processors() + 1;
        }
      else
        {
          e->set_unique_id(_next_unpartitioned_unique_id);
          _next_unpartitioned_unique_id += this->n_processors() + 1;
        }
    }
  else
    {
      _next_unique_id = std::max(_next_unique_id, e->unique_id()+1);
      _next_unique_id =
        ((_next_unique_id + this->n_processors() - 1) / (this->n_processors() + 1) + 1) *
        (this->n_processors() + 1) + this->processor_id();
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
  e->add_extra_integers(_elem_integer_names.size(),
                        _elem_integer_default_values);

  // And set mapping type and data on any new element
  e->set_mapping_type(this->default_mapping_type());
  e->set_mapping_data(this->default_mapping_data());

  return e;
}



Elem * DistributedMesh::add_elem (std::unique_ptr<Elem> e)
{
  // The mesh now takes ownership of the Elem. Eventually the guts of
  // add_elem() will get moved to a private helper function, and
  // calling add_elem() directly will be deprecated.
  return add_elem(e.release());
}



Elem * DistributedMesh::insert_elem (Elem * e)
{
  if (_elements[e->id()])
    this->delete_elem(_elements[e->id()]);

#ifdef LIBMESH_ENABLE_UNIQUE_ID
  if (!e->valid_unique_id())
    {
      if (processor_id() == e->processor_id())
        {
          e->set_unique_id(_next_unique_id);
          _next_unique_id += this->n_processors() + 1;
        }
      else
        {
          e->set_unique_id(_next_unpartitioned_unique_id);
          _next_unpartitioned_unique_id += this->n_processors() + 1;
        }
    }
  else
    {
      _next_unique_id = std::max(_next_unique_id, e->unique_id()+1);
      _next_unique_id =
        ((_next_unique_id + this->n_processors() - 1) / (this->n_processors() + 1) + 1) *
        (this->n_processors() + 1) + this->processor_id();
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
  e->add_extra_integers(_elem_integer_names.size(),
                        _elem_integer_default_values);

  // And set mapping type and data on any new element
  e->set_mapping_type(this->default_mapping_type());
  e->set_mapping_data(this->default_mapping_data());

  return e;
}

Elem * DistributedMesh::insert_elem (std::unique_ptr<Elem> e)
{
  // The mesh now takes ownership of the Elem. Eventually the guts of
  // insert_elem(Elem*) will get moved to a private helper function, and
  // calling insert_elem(Elem*) directly will be deprecated.
  return insert_elem(e.release());
}


void DistributedMesh::delete_elem(Elem * e)
{
  libmesh_assert (e);

  // Try to make the cached elem data more accurate
  _n_elem--;

  // Was this a coarse element, not just a coarsening where we still
  // have some ancestor structure?  Was it a *local* element, that we
  // might have been depending on as an owner of local nodes?  We'll
  // have to be more careful with our nodes in contract() later; no
  // telling if we just locally orphaned a node that should be
  // globally retained.
  if (e->processor_id() == this->processor_id() &&
      !e->parent())
    _deleted_coarse_elements = true;

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



void DistributedMesh::renumber_elem(const dof_id_type old_id,
                                    const dof_id_type new_id)
{
  // This could be a no-op
  if (old_id == new_id)
    return;

  Elem * el = _elements[old_id];
  libmesh_assert (el);
  libmesh_assert_equal_to (el->id(), old_id);

  el->set_id(new_id);
  libmesh_assert (!_elements[new_id]);
  _elements[new_id] = el;
  _elements.erase(old_id);
}



Node * DistributedMesh::add_point (const Point & p,
                                   const dof_id_type id,
                                   const processor_id_type proc_id)
{
  Node * old_n = this->query_node_ptr(id);

  if (old_n)
    {
      *old_n = p;
      old_n->processor_id() = proc_id;

      return old_n;
    }

  Node * n = Node::build(p, id).release();
  n->processor_id() = proc_id;

  return DistributedMesh::add_node(n);
}


void DistributedMesh::own_node (Node & n)
{
  // This had better be a node in our mesh
  libmesh_assert(_nodes[n.id()] == &n);

  _nodes[n.id()] = nullptr;
  _n_nodes--;

  n.set_id(DofObject::invalid_id);
  n.processor_id() = this->processor_id();

  this->add_node(&n);
}


Node * DistributedMesh::add_node (Node * n)
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
    // We need a const dofobject_container so we don't inadvertently create
    // nullptr entries when testing for non-nullptr ones
    const dofobject_container<Node> & const_nodes = _nodes;
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
      if (processor_id() == n->processor_id())
        {
          n->set_unique_id(_next_unique_id);
          _next_unique_id += this->n_processors() + 1;
        }
      else
        {
          n->set_unique_id(_next_unpartitioned_unique_id);
          _next_unpartitioned_unique_id += this->n_processors() + 1;
        }
    }
  else
    {
      _next_unique_id = std::max(_next_unique_id, n->unique_id()+1);
      _next_unique_id =
        ((_next_unique_id + this->n_processors() - 1) / (this->n_processors() + 1) + 1) *
        (this->n_processors() + 1) + this->processor_id();
    }
#endif

  n->add_extra_integers(_node_integer_names.size(),
                        _node_integer_default_values);

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

Node * DistributedMesh::add_node (std::unique_ptr<Node> n)
{
  // The mesh now takes ownership of the Node. Eventually the guts of
  // add_node() will get moved to a private helper function, and
  // calling add_node() directly will be deprecated.
  return add_node(n.release());
}

void DistributedMesh::delete_node(Node * n)
{
  libmesh_assert(n);
  libmesh_assert(_nodes[n->id()]);

  // Try to make the cached elem data more accurate
  _n_nodes--;

  // Delete the node from the BoundaryInfo object
  this->get_boundary_info().remove(n);
  _constraint_rows.erase(n);

  // But not yet from the container; we might invalidate
  // an iterator that way!

  //_nodes.erase(n->id());

  // Instead, we set it to nullptr for now

  _nodes[n->id()] = nullptr;

  // delete the node
  delete n;
}



void DistributedMesh::renumber_node(const dof_id_type old_id,
                                    const dof_id_type new_id)
{
  // This could be a no-op
  if (old_id == new_id)
    return;

  Node * nd = _nodes[old_id];
  libmesh_assert (nd);
  libmesh_assert_equal_to (nd->id(), old_id);

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
  nd->set_id(new_id);

  _nodes.erase(old_id);
}



void DistributedMesh::clear ()
{
  // Call parent clear function
  MeshBase::clear();

  // Clear our elements and nodes
  // There is no need to remove them from
  // the BoundaryInfo data structure since we
  // already cleared it.
  this->DistributedMesh::clear_elems();

  for (auto & node : _nodes)
    delete node;

  _nodes.clear();

  // We're no longer distributed if we were before
  _is_serial = true;
  _is_serial_on_proc_0 = true;

  // We deleted a ton of coarse elements, but their nodes got deleted too so
  // all is copacetic.
  _deleted_coarse_elements = false;

  // Correct our caches
  _n_nodes = 0;
  _max_node_id = 0;
  _next_free_local_node_id = this->processor_id();
  _next_free_unpartitioned_node_id = this->n_processors();
}



void DistributedMesh::clear_elems ()
{
  for (auto & elem : _elements)
    delete elem;

  _elements.clear();

  // Correct our caches
  _n_elem = 0;
  _max_elem_id = 0;
  _next_free_local_elem_id = this->processor_id();
  _next_free_unpartitioned_elem_id = this->n_processors();
}



void DistributedMesh::redistribute ()
{
  // If this is a truly parallel mesh, go through the redistribution/gather/delete remote steps
  if (!this->is_serial())
    {
      // Construct a MeshCommunication object to actually redistribute the nodes
      // and elements according to the partitioner, and then to re-gather the neighbors.
      MeshCommunication mc;
      mc.redistribute(*this);

      this->update_parallel_id_counts();

      // We communicate valid neighbor links for newly-redistributed
      // elements, but we may still have remote_elem links that should
      // be corrected but whose correction wasn't communicated.
      this->find_neighbors(/*reset_remote_elements*/ false,
                           /*reset_current_list*/ false /*non-default*/,
                           /*assert_valid*/ false,
                           /*check_non_remote*/ false /*non-default!*/);

      // Is this necessary?  If we are called from prepare_for_use(), this will be called
      // anyway... but users can always call partition directly, in which case we do need
      // to call delete_remote_elements()...
      //
      // Regardless of whether it's necessary, it isn't safe.  We
      // haven't communicated new node processor_ids yet, and we can't
      // delete nodes until we do.
      // this->delete_remote_elements();
    }
  else
    // The base class can handle non-distributed things, like
    // notifying any GhostingFunctors of changes
    MeshBase::redistribute();
}



void DistributedMesh::update_post_partitioning ()
{
  // this->recalculate_n_partitions();

  // Partitioning changes our numbers of unpartitioned objects
  this->update_parallel_id_counts();
}



template <typename T>
void DistributedMesh::libmesh_assert_valid_parallel_object_ids(const dofobject_container<T> & objects) const
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



void DistributedMesh::libmesh_assert_valid_parallel_ids () const
{
  this->libmesh_assert_valid_parallel_object_ids (this->_elements);
  this->libmesh_assert_valid_parallel_object_ids (this->_nodes);
}



void DistributedMesh::libmesh_assert_valid_parallel_p_levels () const
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




void DistributedMesh::libmesh_assert_valid_parallel_flags () const
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



template <typename T>
dof_id_type
DistributedMesh::renumber_dof_objects(dofobject_container<T> & objects)
{
  // This function must be run on all processors at once
  parallel_object_only();

  typedef typename dofobject_container<T>::veclike_iterator object_iterator;

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
  for (auto p : make_range(this->n_processors()))
    if (p != this->processor_id())
      {
        if (const auto p_it = ghost_objects_from_proc.find(p);
            p_it != ghost_end)
          requested_ids[p].reserve(p_it->second);
      }

  end = objects.end();
  for (it = objects.begin(); it != end; ++it)
    {
      T * obj = *it;
      if (!obj)
        continue;
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
            obj->set_unique_id(data[i]);
        }
    };

  const unique_id_type * unique_ex = nullptr;
  Parallel::pull_parallel_vector_data
    (this->comm(), requested_ids, unique_gather_functor,
     unique_action_functor, unique_ex);
#endif

  // Next set unpartitioned object ids
  next_id = 0;
  for (auto i : make_range(this->n_processors()))
    next_id += objects_on_proc[i];
  for (it = objects.begin(); it != end; ++it)
    {
      T * obj = *it;
      if (!obj)
        continue;
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


void DistributedMesh::renumber_nodes_and_elements ()
{
  parallel_object_only();

#ifdef DEBUG
  // Make sure our ids and flags are consistent
  this->libmesh_assert_valid_parallel_ids();
  this->libmesh_assert_valid_parallel_flags();
  this->libmesh_assert_valid_parallel_p_levels();
#endif

  LOG_SCOPE("renumber_nodes_and_elements()", "DistributedMesh");

  // Nodes not connected to any elements, and nullptr node entries
  // in our container, should be deleted.  But wait!  If we've deleted coarse
  // local elements on some processor, other processors might have ghosted
  // nodes from it that are now no longer connected to any elements on it, but
  // that are connected to their own semilocal elements.  We'll have to
  // communicate to ascertain if that's the case.
  this->comm().max(_deleted_coarse_elements);

  // What used nodes do we see on our proc?
  std::set<dof_id_type> used_nodes;

  // What used node info should we send from our proc?  Could we take ownership
  // of each node if we needed to?
  std::map<processor_id_type, std::map<dof_id_type, bool>>
    used_nodes_on_proc;

  // flag the nodes we need
  for (auto & elem : this->element_ptr_range())
    for (const Node & node : elem->node_ref_range())
      {
        const dof_id_type n = node.id();
        used_nodes.insert(n);
        if (_deleted_coarse_elements)
          {
            const processor_id_type p = node.processor_id();
            if (p != this->processor_id())
              {
                auto & used_nodes_on_p = used_nodes_on_proc[p];
                if (elem->processor_id() == this->processor_id())
                  used_nodes_on_p[n] = true;
                else
                  if (!used_nodes_on_p.count(n))
                    used_nodes_on_p[n] = false;
              }
          }
      }

  if (_deleted_coarse_elements)
    {
      // "unsigned char" == "bool, but MPI::BOOL is iffy to use"
      typedef unsigned char boolish;
      std::map<processor_id_type, std::vector<std::pair<dof_id_type, boolish>>>
        used_nodes_on_proc_vecs;
      for (auto & [pid, nodemap] : used_nodes_on_proc)
        used_nodes_on_proc_vecs[pid].assign(nodemap.begin(), nodemap.end());

      std::map<dof_id_type,processor_id_type> repartitioned_node_pids;
      std::map<processor_id_type, std::set<dof_id_type>>
        repartitioned_node_sets_to_push;

      auto ids_action_functor =
        [&used_nodes, &repartitioned_node_pids,
         &repartitioned_node_sets_to_push]
        (processor_id_type pid,
         const std::vector<std::pair<dof_id_type, boolish>> & ids_and_bools)
        {
          for (auto [n, sender_could_become_owner] : ids_and_bools)
            {
              // If we don't see a use for our own node, but someone
              // else does, better figure out who should own it next.
              if (!used_nodes.count(n))
                {
                  if (auto it = repartitioned_node_pids.find(n);
                      sender_could_become_owner)
                    {
                      if (it != repartitioned_node_pids.end() &&
                          pid < it->second)
                        it->second = pid;
                      else
                        repartitioned_node_pids[n] = pid;
                    }
                  else
                    if (it == repartitioned_node_pids.end())
                      repartitioned_node_pids[n] =
                        DofObject::invalid_processor_id;

                  repartitioned_node_sets_to_push[pid].insert(n);
                }
            }
        };

      // We need two pushes instead of a pull here because we need to
      // know *all* the queries for a particular node before we can
      // respond to *any* of them.
      Parallel::push_parallel_vector_data
      (this->comm(), used_nodes_on_proc_vecs, ids_action_functor);

      // Repartition (what used to be) our own nodes first
      for (auto & [n, p] : repartitioned_node_pids)
        {
          Node & node = this->node_ref(n);
          libmesh_assert_equal_to(node.processor_id(), this->processor_id());
          libmesh_assert_not_equal_to_msg(p, DofObject::invalid_processor_id, "Node " << n << " is lost?");
          node.processor_id() = p;
        }

      // Then push to repartition others' ghosted copies.

      std::map<processor_id_type, std::vector<std::pair<dof_id_type,processor_id_type>>>
        repartitioned_node_vecs;

      for (auto & [p, nodeset] : repartitioned_node_sets_to_push)
        {
          auto & rn_vec = repartitioned_node_vecs[p];
          for (auto n : nodeset)
            rn_vec.emplace_back(n, repartitioned_node_pids[n]);
        }

      auto repartition_node_functor =
        [this]
        (processor_id_type libmesh_dbg_var(pid),
         const std::vector<std::pair<dof_id_type, processor_id_type>> & ids_and_pids)
        {
          for (auto [n, p] : ids_and_pids)
            {
              libmesh_assert_not_equal_to(p, DofObject::invalid_processor_id);
              Node & node = this->node_ref(n);
              libmesh_assert_equal_to(node.processor_id(), pid);
              node.processor_id() = p;
            }
        };

      Parallel::push_parallel_vector_data
      (this->comm(), repartitioned_node_vecs, repartition_node_functor);
    }

  _deleted_coarse_elements = false;

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
            _constraint_rows.erase(nd);

            // delete the node
            delete nd;

            it = _nodes.erase(it);
          }
        else
          ++it;
      }
  }

  this->_preparation.has_removed_orphaned_nodes = true;

  if (_skip_renumber_nodes_and_elements)
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



void DistributedMesh::fix_broken_node_and_element_numbering ()
{
  // We can't use range-for here because we need access to the special
  // iterators' methods, not just to their dereferenced values.

  // Nodes first
  for (auto pr = this->_nodes.begin(),
           end = this->_nodes.end(); pr != end; ++pr)
    {
      Node * n = *pr;
      if (n != nullptr)
        {
          const dof_id_type id = pr.index();
          n->set_id() = id;
          libmesh_assert_equal_to(this->node_ptr(id), n);
        }
    }

  // Elements next
  for (auto pr = this->_elements.begin(),
           end = this->_elements.end(); pr != end; ++pr)
    {
      Elem * e = *pr;
      if (e != nullptr)
        {
          const dof_id_type id = pr.index();
          e->set_id() = id;
          libmesh_assert_equal_to(this->elem_ptr(id), e);
        }
    }
}



dof_id_type DistributedMesh::n_active_elem () const
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



void DistributedMesh::delete_remote_elements()
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
  dofobject_container<Elem>::veclike_iterator e_it        = _elements.begin();
  const dofobject_container<Elem>::veclike_iterator e_end = _elements.end();
  while (e_it != e_end)
    if (!*e_it)
      e_it = _elements.erase(e_it);
    else
      ++e_it;

  dofobject_container<Node>::veclike_iterator n_it        = _nodes.begin();
  const dofobject_container<Node>::veclike_iterator n_end = _nodes.end();
  while (n_it != n_end)
    if (!*n_it)
      n_it = _nodes.erase(n_it);
    else
      ++n_it;

  // We may have deleted no-longer-connected nodes or coarsened-away
  // elements; let's update our caches.
  this->update_parallel_id_counts();

  // We may have deleted nodes or elements that were the only local
  // representatives of some particular boundary id(s); let's update
  // those caches.
  this->get_boundary_info().regenerate_id_sets();

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

  this->_preparation.has_removed_remote_elements = true;
}


void DistributedMesh::add_extra_ghost_elem(Elem * e)
{
  // First add the elem like normal
  add_elem(e);

  // Now add it to the set that won't be deleted when we call
  // delete_remote_elements()
  _extra_ghost_elems.insert(e);
}

void
DistributedMesh::clear_extra_ghost_elems(const std::set<Elem *> & extra_ghost_elems)
{
  std::set<Elem *> tmp;
  std::set_difference(_extra_ghost_elems.begin(), _extra_ghost_elems.end(),
                      extra_ghost_elems.begin(), extra_ghost_elems.end(),
                      std::inserter(tmp, tmp.begin()));
  _extra_ghost_elems = tmp;
}

void DistributedMesh::allgather()
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

void DistributedMesh::gather_to_zero()
{
  if (_is_serial_on_proc_0)
    return;

  _is_serial_on_proc_0 = true;
  MeshCommunication().gather(0, *this);
}


} // namespace libMesh

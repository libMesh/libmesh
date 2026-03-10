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
#include "libmesh/replicated_mesh.h"

#include "libmesh/boundary_info.h"
#include "libmesh/elem.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/mesh_communication.h"
#include "libmesh/parallel_implementation.h"
#include "libmesh/partitioner.h"
#include "libmesh/point.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/utility.h"

// C++ includes
#include <unordered_map>
#include <unordered_set>

namespace libMesh
{

// ------------------------------------------------------------
// ReplicatedMesh class member functions
ReplicatedMesh::ReplicatedMesh (const Parallel::Communicator & comm_in,
                                unsigned char d) :
  UnstructuredMesh (comm_in,d),
  _n_nodes(0), _n_elem(0)
{
#ifdef LIBMESH_ENABLE_UNIQUE_ID
  // In serial we just need to reset the next unique id to zero
  // here in the constructor.
  _next_unique_id = 0;
#endif

  const std::string default_partitioner = "metis";
  const std::string my_partitioner =
    libMesh::command_line_value("--default-partitioner",
                                default_partitioner);
  _partitioner = Partitioner::build
    (Utility::string_to_enum<PartitionerType>(my_partitioner));
}


std::string_view ReplicatedMesh::subclass_first_difference_from (const MeshBase & other_mesh_base) const
{
  const ReplicatedMesh * rep_mesh_ptr =
    dynamic_cast<const ReplicatedMesh *>(&other_mesh_base);
  if (!rep_mesh_ptr)
    return "ReplicatedMesh class";
  const ReplicatedMesh & other_mesh = *rep_mesh_ptr;

#define CHECK_MEMBER(member_name) \
  if (member_name != other_mesh.member_name) \
    return #member_name;

  CHECK_MEMBER(_n_nodes);
  CHECK_MEMBER(_n_elem);
#ifdef LIBMESH_ENABLE_UNIQUE_ID
  CHECK_MEMBER(_next_unique_id);
#endif
  if (!this->nodes_and_elements_equal(other_mesh))
    return "nodes and/or elements";

  return "";
}


ReplicatedMesh::~ReplicatedMesh ()
{
  this->ReplicatedMesh::clear();  // Free nodes and elements
}


// This might be specialized later, but right now it's just here to
// make sure the compiler doesn't give us a default (non-deep) copy
// constructor instead.
ReplicatedMesh::ReplicatedMesh (const ReplicatedMesh & other_mesh) :
  ReplicatedMesh(static_cast<const MeshBase&>(other_mesh))
{
#ifdef LIBMESH_ENABLE_UNIQUE_ID
  this->_next_unique_id = other_mesh._next_unique_id;
#endif
}


ReplicatedMesh::ReplicatedMesh (const MeshBase & other_mesh) :
  UnstructuredMesh (other_mesh),
  _n_nodes(0), _n_elem(0) // copy_* will increment this
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

  // If other_mesh is distributed, then we've got parts of it on each
  // processor but we're not replicated yet; fix that.
  if (!other_mesh.is_serial())
    MeshCommunication().allgather(*this);
}

ReplicatedMesh & ReplicatedMesh::operator= (ReplicatedMesh && other_mesh)
{
  LOG_SCOPE("operator=(&&)", "ReplicatedMesh");

  // Move assign as an UnstructuredMesh
  this->UnstructuredMesh::operator=(std::move(other_mesh));

  // Nodes and elements belong to ReplicatedMesh and have to be
  // moved before we can move arbitrary GhostingFunctor, Partitioner,
  // etc. subclasses.
  this->move_nodes_and_elements(std::move(other_mesh));

  // Handle those remaining moves.
  this->post_dofobject_moves(std::move(other_mesh));

  return *this;
}

MeshBase & ReplicatedMesh::assign(MeshBase && other_mesh)
{
  *this = std::move(cast_ref<ReplicatedMesh&>(other_mesh));

  return *this;
}

void ReplicatedMesh::move_nodes_and_elements(MeshBase && other_meshbase)
{
  ReplicatedMesh & other_mesh = cast_ref<ReplicatedMesh&>(other_meshbase);

  this->_nodes = std::move(other_mesh._nodes);
  this->_n_nodes = other_mesh.n_nodes();

  this->_elements = std::move(other_mesh._elements);
  this->_n_elem = other_mesh.n_elem();
}


const Point & ReplicatedMesh::point (const dof_id_type i) const
{
  return this->node_ref(i);
}




const Node * ReplicatedMesh::node_ptr (const dof_id_type i) const
{
  libmesh_assert_less (i, this->max_node_id());
  libmesh_assert(_nodes[i]);
  libmesh_assert_equal_to (_nodes[i]->id(), i); // This will change soon

  return _nodes[i];
}




Node * ReplicatedMesh::node_ptr (const dof_id_type i)
{
  libmesh_assert_less (i, this->max_node_id());
  libmesh_assert(_nodes[i]);
  libmesh_assert_equal_to (_nodes[i]->id(), i); // This will change soon

  return _nodes[i];
}




const Node * ReplicatedMesh::query_node_ptr (const dof_id_type i) const
{
  if (i >= this->max_node_id())
    return nullptr;
  libmesh_assert (_nodes[i] == nullptr ||
                  _nodes[i]->id() == i); // This will change soon

  return _nodes[i];
}




Node * ReplicatedMesh::query_node_ptr (const dof_id_type i)
{
  if (i >= this->max_node_id())
    return nullptr;
  libmesh_assert (_nodes[i] == nullptr ||
                  _nodes[i]->id() == i); // This will change soon

  return _nodes[i];
}




const Elem * ReplicatedMesh::elem_ptr (const dof_id_type i) const
{
  libmesh_assert_less (i, this->max_elem_id());
  libmesh_assert(_elements[i]);
  libmesh_assert_equal_to (_elements[i]->id(), i); // This will change soon

  return _elements[i];
}




Elem * ReplicatedMesh::elem_ptr (const dof_id_type i)
{
  libmesh_assert_less (i, this->max_elem_id());
  libmesh_assert(_elements[i]);
  libmesh_assert_equal_to (_elements[i]->id(), i); // This will change soon

  return _elements[i];
}




const Elem * ReplicatedMesh::query_elem_ptr (const dof_id_type i) const
{
  if (i >= this->max_elem_id())
    return nullptr;
  libmesh_assert (_elements[i] == nullptr ||
                  _elements[i]->id() == i); // This will change soon

  return _elements[i];
}




Elem * ReplicatedMesh::query_elem_ptr (const dof_id_type i)
{
  if (i >= this->max_elem_id())
    return nullptr;
  libmesh_assert (_elements[i] == nullptr ||
                  _elements[i]->id() == i); // This will change soon

  return _elements[i];
}




Elem * ReplicatedMesh::add_elem (Elem * e)
{
  libmesh_assert(e);

  // We no longer merely append elements with ReplicatedMesh

  // If the user requests a valid id that doesn't correspond to an
  // existing element, let's give them that id, resizing the elements
  // container if necessary.
  if (!e->valid_id())
    e->set_id (cast_int<dof_id_type>(_elements.size()));

#ifdef LIBMESH_ENABLE_UNIQUE_ID
  if (!e->valid_unique_id())
    e->set_unique_id(_next_unique_id++);
  else
   _next_unique_id = std::max(_next_unique_id, e->unique_id()+1);
#endif

  const dof_id_type id = e->id();

  if (id < _elements.size())
    {
      // This should *almost* never happen, but we rely on it when
      // using allgather to replicate a not-yet-actually-replicated
      // ReplicatedMesh under construction in parallel.
      if (e == _elements[id])
        return e;

      // Overwriting existing elements is still probably a mistake.
      libmesh_assert(!_elements[id]);
    }
  else
    {
      _elements.resize(id+1, nullptr);
    }

  ++_n_elem;
  _elements[id] = e;

  // Make sure any new element is given space for any extra integers
  // we've requested
  e->add_extra_integers(_elem_integer_names.size(),
                        _elem_integer_default_values);

  // And set mapping type and data on any new element
  e->set_mapping_type(this->default_mapping_type());
  e->set_mapping_data(this->default_mapping_data());

  return e;
}

Elem * ReplicatedMesh::add_elem (std::unique_ptr<Elem> e)
{
  // The mesh now takes ownership of the Elem. Eventually the guts of
  // add_elem() will get moved to a private helper function, and
  // calling add_elem() directly will be deprecated.
  return add_elem(e.release());
}



Elem * ReplicatedMesh::insert_elem (Elem * e)
{
#ifdef LIBMESH_ENABLE_UNIQUE_ID
  if (!e->valid_unique_id())
    e->set_unique_id(_next_unique_id++);
  else
   _next_unique_id = std::max(_next_unique_id, e->unique_id()+1);
#endif

  dof_id_type eid = e->id();
  libmesh_assert_less (eid, _elements.size());
  Elem * oldelem = _elements[eid];

  if (oldelem)
    {
      libmesh_assert_equal_to (oldelem->id(), eid);
      this->delete_elem(oldelem);
    }

  ++_n_elem;
  _elements[eid] = e;

  // Make sure any new element is given space for any extra integers
  // we've requested
  e->add_extra_integers(_elem_integer_names.size(),
                        _elem_integer_default_values);

  // And set mapping type and data on any new element
  e->set_mapping_type(this->default_mapping_type());
  e->set_mapping_data(this->default_mapping_data());

  return e;
}

Elem * ReplicatedMesh::insert_elem (std::unique_ptr<Elem> e)
{
  // The mesh now takes ownership of the Elem. Eventually the guts of
  // insert_elem(Elem*) will get moved to a private helper function, and
  // calling insert_elem(Elem*) directly will be deprecated.
  return insert_elem(e.release());
}



void ReplicatedMesh::delete_elem(Elem * e)
{
  libmesh_assert(e);

  // Initialize an iterator to eventually point to the element we want to delete
  std::vector<Elem *>::iterator pos = _elements.end();

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
  this->get_boundary_info().remove(e);

  // delete the element
  --_n_elem;
  delete e;

  // explicitly zero the pointer
  *pos = nullptr;
}



void ReplicatedMesh::renumber_elem(const dof_id_type old_id,
                                   const dof_id_type new_id)
{
  // This could be a no-op
  if (old_id == new_id)
    return;

  // This doesn't get used in serial yet
  Elem * el = _elements[old_id];
  libmesh_assert (el);

  if (new_id >= _elements.size())
    _elements.resize(new_id+1, nullptr);

  el->set_id(new_id);
  libmesh_assert (!_elements[new_id]);
  _elements[new_id] = el;
  _elements[old_id] = nullptr;
}



Node * ReplicatedMesh::add_point (const Point & p,
                                  const dof_id_type id,
                                  const processor_id_type proc_id)
{
  Node * n = nullptr;

  // If the user requests a valid id, either
  // provide the existing node or resize the container
  // to fit the new node.
  if (id != DofObject::invalid_id)
    if (id < _nodes.size())
      n = _nodes[id];
    else
      _nodes.resize(id+1);
  else
    _nodes.push_back (static_cast<Node *>(nullptr));

  // if the node already exists, then assign new (x,y,z) values
  if (n)
    *n = p;
  // otherwise build a new node, put it in the right spot, and return
  // a valid pointer.
  else
    {
      n = Node::build(p, (id == DofObject::invalid_id) ?
                      cast_int<dof_id_type>(_nodes.size()-1) : id).release();
      n->processor_id() = proc_id;

      n->add_extra_integers(_node_integer_names.size(),
                            _node_integer_default_values);

#ifdef LIBMESH_ENABLE_UNIQUE_ID
      if (!n->valid_unique_id())
        n->set_unique_id(_next_unique_id++);
      else
       _next_unique_id = std::max(_next_unique_id, n->unique_id()+1);
#endif

      ++_n_nodes;
      if (id == DofObject::invalid_id)
        _nodes.back() = n;
      else
        _nodes[id] = n;
    }

  // better not pass back a nullptr.
  libmesh_assert (n);

  return n;
}



Node * ReplicatedMesh::add_node (Node * n)
{
  libmesh_assert(n);

  // If the user requests a valid id, either set the existing
  // container entry or resize the container to fit the new node.
  if (n->valid_id())
    {
      const dof_id_type id = n->id();
      if (id < _nodes.size())
        libmesh_assert(!_nodes[id]);
      else
        _nodes.resize(id+1); // default nullptr

      _nodes[id] = n;
    }
  else
    {
      n->set_id (cast_int<dof_id_type>(_nodes.size()));
      _nodes.push_back(n);
    }

  ++_n_nodes;

#ifdef LIBMESH_ENABLE_UNIQUE_ID
  if (!n->valid_unique_id())
    n->set_unique_id(_next_unique_id++);
  else
   _next_unique_id = std::max(_next_unique_id, n->unique_id()+1);
#endif

  n->add_extra_integers(_node_integer_names.size(),
                        _node_integer_default_values);

  return n;
}

Node * ReplicatedMesh::add_node (std::unique_ptr<Node> n)
{
  // The mesh now takes ownership of the Node. Eventually the guts of
  // add_node() will get moved to a private helper function, and
  // calling add_node() directly will be deprecated.
  return add_node(n.release());
}

void ReplicatedMesh::delete_node(Node * n)
{
  libmesh_assert(n);
  libmesh_assert_less (n->id(), _nodes.size());

  // Initialize an iterator to eventually point to the element we want
  // to delete
  std::vector<Node *>::iterator pos;

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
  this->get_boundary_info().remove(n);
  _constraint_rows.erase(n);

  // delete the node
  --_n_nodes;
  delete n;

  // explicitly zero the pointer
  *pos = nullptr;
}



void ReplicatedMesh::renumber_node(const dof_id_type old_id,
                                   const dof_id_type new_id)
{
  // This could be a no-op
  if (old_id == new_id)
    return;

  // This doesn't get used in serial yet
  Node * nd = _nodes[old_id];
  libmesh_assert (nd);

  if (new_id >= _nodes.size())
    _nodes.resize(new_id+1, nullptr);

  nd->set_id(new_id);
  libmesh_assert (!_nodes[new_id]);
  _nodes[new_id] = nd;
  _nodes[old_id] = nullptr;
}



void ReplicatedMesh::clear ()
{
  // Call parent clear function
  MeshBase::clear();

  // Clear our elements and nodes
  // There is no need to remove them from
  // the BoundaryInfo data structure since we
  // already cleared it.
  this->ReplicatedMesh::clear_elems();

  for (auto & node : _nodes)
    delete node;

  _n_nodes = 0;
  _nodes.clear();
}



void ReplicatedMesh::clear_elems ()
{
  for (auto & elem : _elements)
    delete elem;

  _n_elem = 0;
  _elements.clear();
}



void ReplicatedMesh::update_parallel_id_counts()
{
#ifdef LIBMESH_ENABLE_UNIQUE_ID
  _next_unique_id = this->parallel_max_unique_id();
#endif

  this->_preparation.has_synched_id_counts = true;
}



#ifdef LIBMESH_ENABLE_UNIQUE_ID
unique_id_type ReplicatedMesh::parallel_max_unique_id() const
{
  // This function must be run on all processors at once
  parallel_object_only();

  unique_id_type max_local = _next_unique_id;
  this->comm().max(max_local);
  return max_local;
}



void ReplicatedMesh::set_next_unique_id(unique_id_type id)
{
  _next_unique_id = id;
}
#endif



void ReplicatedMesh::renumber_nodes_and_elements ()
{
  LOG_SCOPE("renumber_nodes_and_elem()", "Mesh");

  // node and element id counters
  dof_id_type next_free_elem = 0;
  dof_id_type next_free_node = 0;

  // Will hold the set of nodes that are currently connected to elements
  std::unordered_set<Node *> connected_nodes;

  // Loop over the elements.  Note that there may
  // be nullptrs in the _elements vector from the coarsening
  // process.  Pack the elements in to a contiguous array
  // and then trim any excess.
  {
    std::vector<Elem *>::iterator in        = _elements.begin();
    std::vector<Elem *>::iterator out_iter  = _elements.begin();
    const std::vector<Elem *>::iterator end = _elements.end();

    for (; in != end; ++in)
      if (*in != nullptr)
        {
          Elem * el = *in;

          *out_iter = *in;
          ++out_iter;

          // Increment the element counter
          el->set_id (next_free_elem++);

          if (_skip_renumber_nodes_and_elements)
            {
              // Add this elements nodes to the connected list
              for (auto & n : el->node_ref_range())
                connected_nodes.insert(&n);
            }
          else  // We DO want node renumbering
            {
              // Loop over this element's nodes.  Number them,
              // if they have not been numbered already.  Also,
              // position them in the _nodes vector so that they
              // are packed contiguously from the beginning.
              for (auto & n : el->node_ref_range())
                if (n.id() == next_free_node)    // don't need to process
                  next_free_node++;                      // [(src == dst) below]

                else if (n.id() > next_free_node) // need to process
                  {
                    // The source and destination indices
                    // for this node
                    const dof_id_type src_idx = n.id();
                    const dof_id_type dst_idx = next_free_node++;

                    // ensure we want to swap a valid nodes
                    libmesh_assert(_nodes[src_idx]);

                    // Swap the source and destination nodes
                    std::swap(_nodes[src_idx],
                              _nodes[dst_idx] );

                    // Set proper indices where that makes sense
                    if (_nodes[src_idx] != nullptr)
                      _nodes[src_idx]->set_id (src_idx);
                    _nodes[dst_idx]->set_id (dst_idx);
                  }
            }
        }

    // Erase any additional storage. These elements have been
    // copied into nullptr voids by the procedure above, and are
    // thus repeated and unnecessary.
    _elements.erase (out_iter, end);
  }


  if (_skip_renumber_nodes_and_elements)
    {
      // Loop over the nodes.  Note that there may
      // be nullptrs in the _nodes vector from the coarsening
      // process.  Pack the nodes in to a contiguous array
      // and then trim any excess.

      std::vector<Node *>::iterator in        = _nodes.begin();
      std::vector<Node *>::iterator out_iter  = _nodes.begin();
      const std::vector<Node *>::iterator end = _nodes.end();

      for (; in != end; ++in)
        if (*in != nullptr)
          {
            // This is a reference so that if we change the pointer it will change in the vector
            Node * & nd = *in;

            // If this node is still connected to an elem, put it in the list
            if (connected_nodes.count(nd))
              {
                *out_iter = nd;
                ++out_iter;

                // Increment the node counter
                nd->set_id (next_free_node++);
              }
            else // This node is orphaned, delete it!
              {
                this->get_boundary_info().remove (nd);
                _constraint_rows.erase(nd);

                // delete the node
                --_n_nodes;
                delete nd;
                nd = nullptr;
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

      // Now, delete the unused nodes
      {
        std::vector<Node *>::iterator nd        = _nodes.begin();
        const std::vector<Node *>::iterator end = _nodes.end();

        std::advance (nd, next_free_node);

        for (auto & node : as_range(nd, end))
          {
            // Mesh modification code might have already deleted some
            // nodes
            if (node == nullptr)
              continue;

            // remove any boundary information associated with
            // this node
            this->get_boundary_info().remove (node);
            _constraint_rows.erase(node);

            // delete the node
            --_n_nodes;
            delete node;
            node = nullptr;
          }

        _nodes.erase (nd, end);
      }
    }

  this->_preparation.has_removed_orphaned_nodes = true;

  libmesh_assert_equal_to (next_free_elem, _elements.size());
  libmesh_assert_equal_to (next_free_node, _nodes.size());

  this->update_parallel_id_counts();
}



void ReplicatedMesh::fix_broken_node_and_element_numbering ()
{
  // Nodes first
  for (auto n : index_range(_nodes))
    if (this->_nodes[n] != nullptr)
      this->_nodes[n]->set_id() = cast_int<dof_id_type>(n);

  // Elements next
  for (auto e : index_range(_elements))
    if (this->_elements[e] != nullptr)
      this->_elements[e]->set_id() = cast_int<dof_id_type>(e);
}


dof_id_type ReplicatedMesh::n_active_elem () const
{
  return static_cast<dof_id_type>(std::distance (this->active_elements_begin(),
                                                 this->active_elements_end()));
}

std::vector<dof_id_type>
ReplicatedMesh::get_disconnected_subdomains(std::vector<subdomain_id_type> * subdomain_ids) const
{
  // find number of disconnected subdomains
  std::vector<dof_id_type> representative_elem_ids;

  // use subdomain_ids as markers for all elements to indicate if the elements
  // have been visited. Note: here subdomain ID is unrelated with element
  // subdomain_id().
  std::vector<subdomain_id_type> subdomains;
  if (!subdomain_ids)
    subdomain_ids = &subdomains;
  subdomain_ids->clear();
  subdomain_ids->resize(max_elem_id() + 1, Elem::invalid_subdomain_id);

  // counter of disconnected subdomains
  subdomain_id_type subdomain_counter = 0;

  // a stack for visiting elements, make its capacity sufficiently large to avoid
  // memory allocation and deallocation when the vector size changes
  std::vector<const Elem *> list;
  list.reserve(n_elem());

  // counter of visited elements
  dof_id_type visited = 0;
  dof_id_type n_active = n_active_elem();
  do
  {
    for (const auto & elem : active_element_ptr_range())
      if ((*subdomain_ids)[elem->id()] == Elem::invalid_subdomain_id)
      {
        list.push_back(elem);
        (*subdomain_ids)[elem->id()] = subdomain_counter;
        break;
      }
    // we should be able to find a seed here
    libmesh_assert(list.size() > 0);

    dof_id_type min_id = std::numeric_limits<dof_id_type>::max();
    while (list.size() > 0)
    {
      // pop up an element
      const Elem * elem = list.back(); list.pop_back(); ++visited;

      min_id = std::min(elem->id(), min_id);

      for (auto s : elem->side_index_range())
      {
        const Elem * neighbor = elem->neighbor_ptr(s);
        if (neighbor != nullptr && (*subdomain_ids)[neighbor->id()] == Elem::invalid_subdomain_id)
        {
          // neighbor must be active
          libmesh_assert(neighbor->active());
          list.push_back(neighbor);
          (*subdomain_ids)[neighbor->id()] = subdomain_counter;
        }
      }
    }

    representative_elem_ids.push_back(min_id);
    subdomain_counter++;
  }
  while (visited != n_active);

  return representative_elem_ids;
}

std::unordered_map<dof_id_type, std::vector<std::vector<Point>>>
ReplicatedMesh::get_boundary_points() const
{
  libmesh_error_msg_if(mesh_dimension() != 2,
                       "Error: get_boundary_points only works for 2D now");

  // find number of disconnected subdomains
  // subdomains will hold the IDs of disconnected subdomains for all elements.
  std::vector<subdomain_id_type> subdomains;
  std::vector<dof_id_type> elem_ids = get_disconnected_subdomains(&subdomains);

  std::unordered_map<dof_id_type, std::vector<std::vector<Point>>> boundary_points;

  // get all boundary sides that are to be erased later during visiting
  // use a comparison functor to avoid run-time randomness due to pointers
  struct boundary_side_compare
  {
    bool operator()(const std::pair<const Elem *, unsigned int> & lhs,
                    const std::pair<const Elem *, unsigned int> & rhs) const
      {
        if (lhs.first->id() < rhs.first->id())
          return true;
        else if (lhs.first->id() == rhs.first->id())
        {
          if (lhs.second < rhs.second)
            return true;
        }
        return false;
      }
  };
  std::set<std::pair<const Elem *, unsigned int>, boundary_side_compare> boundary_elements;
  for (const auto & elem : active_element_ptr_range())
    for (auto s : elem->side_index_range())
      if (elem->neighbor_ptr(s) == nullptr)
        boundary_elements.insert(std::pair<const Elem *, unsigned int>(elem, s));

  while (!boundary_elements.empty())
  {
    // get the first entry as the seed
    const Elem * eseed = boundary_elements.begin()->first;
    unsigned int sseed = boundary_elements.begin()->second;

    // get the subdomain ID that these boundary sides attached to
    subdomain_id_type subdomain_id = subdomains[eseed->id()];

    // start visiting the mesh to find all boundary nodes with the seed
    std::vector<Point> bpoints;
    const Elem * elem = eseed;
    unsigned int s = sseed;
    std::vector<unsigned int> local_side_nodes = elem->nodes_on_side(s);
    while (true)
    {
      std::pair<const Elem *, unsigned int> side(elem, s);
      libmesh_assert(boundary_elements.count(side));
      boundary_elements.erase(side);

      // push all nodes on the side except the node on the other end of the side (index 1)
      for (auto i : index_range(local_side_nodes))
        if (i != 1)
          bpoints.push_back(*static_cast<const Point *>(elem->node_ptr(local_side_nodes[i])));

      // use the last node to find next element and side
      const Node * node = elem->node_ptr(local_side_nodes[1]);
      std::set<const Elem *> neighbors;
      elem->find_point_neighbors(*node, neighbors);

      // if only one neighbor is found (itself), this node is a cornor node on boundary
      if (neighbors.size() != 1)
        neighbors.erase(elem);

      // find the connecting side
      bool found = false;
      for (const auto & neighbor : neighbors)
      {
        for (auto ss : neighbor->side_index_range())
          if (neighbor->neighbor_ptr(ss) == nullptr && !(elem == neighbor && s == ss))
          {
            local_side_nodes = neighbor->nodes_on_side(ss);
            // we expect the starting point of the side to be the same as the end of the previous side
            if (neighbor->node_ptr(local_side_nodes[0]) == node)
            {
              elem = neighbor;
              s = ss;
              found = true;
              break;
            }
            else if (neighbor->node_ptr(local_side_nodes[1]) == node)
            {
              elem = neighbor;
              s = ss;
              found = true;
              // flip nodes in local_side_nodes because the side is in an opposite direction
              auto temp(local_side_nodes);
              local_side_nodes[0] = temp[1];
              local_side_nodes[1] = temp[0];
              for (unsigned int i = 2; i < temp.size(); ++i)
                local_side_nodes[temp.size() + 1 - i] = temp[i];
              break;
            }
          }
        if (found)
          break;
      }

      libmesh_error_msg_if(!found, "ERROR: mesh topology error on visiting boundary sides");

      // exit if we reach the starting point
      if (elem == eseed && s == sseed)
        break;
    }
    boundary_points[elem_ids[subdomain_id]].push_back(bpoints);
  }

  return boundary_points;
}

} // namespace libMesh

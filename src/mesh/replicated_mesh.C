// The libMesh Finite Element Library.
// Copyright (C) 2002-2021 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/replicated_mesh.h"
#include "libmesh/utility.h"
#include "libmesh/parallel.h"
#include "libmesh/point.h"
#ifdef LIBMESH_HAVE_NANOFLANN
#include "libmesh/nanoflann.hpp"
#endif

// C++ includes
#include <unordered_map>
#include <unordered_set>

namespace libMesh
{

// This class adapts a vector of Nodes (represented by a pair of a Point and a dof_id_type)
// for use in a nanoflann KD-Tree

class VectorOfNodesAdaptor
{
private:
  const std::vector<std::pair<Point, dof_id_type>> _nodes;

public:
  VectorOfNodesAdaptor(const std::vector<std::pair<Point, dof_id_type>> & nodes) :
    _nodes(nodes)
  {}

  /**
   * Must return the number of data points
   */
  inline size_t kdtree_get_point_count() const { return _nodes.size(); }

  /**
   * \returns The dim'th component of the idx'th point in the class:
   * Since this is inlined and the "dim" argument is typically an immediate value, the
   *  "if's" are actually solved at compile time.
   */
  inline Real kdtree_get_pt(const size_t idx, int dim) const
    {
      libmesh_assert_less (idx, _nodes.size());
      libmesh_assert_less (dim, 3);

      const Point & p(_nodes[idx].first);

      if (dim==0) return p(0);
      if (dim==1) return p(1);
      return p(2);
    }

  /*
   * Optional bounding-box computation
   */
  template <class BBOX>
  bool kdtree_get_bbox(BBOX & /* bb */) const { return false; }
};


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
  _partitioner = libmesh_make_unique<MetisPartitioner>();
}



ReplicatedMesh::~ReplicatedMesh ()
{
  this->clear();  // Free nodes and elements
}


// This might be specialized later, but right now it's just here to
// make sure the compiler doesn't give us a default (non-deep) copy
// constructor instead.
ReplicatedMesh::ReplicatedMesh (const ReplicatedMesh & other_mesh) :
  UnstructuredMesh (other_mesh),
  _n_nodes(0), _n_elem(0) // copy_* will increment this
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
  this->_next_unique_id = other_mesh._next_unique_id;
#endif
}


ReplicatedMesh::ReplicatedMesh (const UnstructuredMesh & other_mesh) :
  UnstructuredMesh (other_mesh),
  _n_nodes(0), _n_elem(0) // copy_* will increment this
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
}

ReplicatedMesh & ReplicatedMesh::operator= (ReplicatedMesh && other_mesh)
{
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
  *this = std::move(*(cast_ptr<ReplicatedMesh*>(&other_mesh)));

  return *this;
}

void ReplicatedMesh::move_nodes_and_elements(MeshBase && other_mesh)
{
  this->_nodes = std::move((cast_ptr<ReplicatedMesh*>(&other_mesh))->_nodes);
  this->_n_nodes = other_mesh.n_nodes();

  this->_elements = std::move((cast_ptr<ReplicatedMesh*>(&other_mesh))->_elements);
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
  // This doesn't get used in serial yet
  Elem * el = _elements[old_id];
  libmesh_assert (el);

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



Node * ReplicatedMesh::insert_node(Node * n)
{
  libmesh_error_msg_if(!n, "Error, attempting to insert nullptr node.");
  libmesh_error_msg_if(n->id() == DofObject::invalid_id, "Error, cannot insert node with invalid id.");

  if (n->id() < _nodes.size())
    {
      // Trying to add an existing node is a no-op.  This lets us use
      // a straightforward MeshCommunication::allgather() to fix a
      // not-actually-replicated-yet ReplicatedMesh, as occurs when
      // reading Nemesis files.
      if (n->valid_id() && _nodes[n->id()] == n)
        return n;

      // Don't allow anything else to replace an existing Node.
      libmesh_error_msg_if(_nodes[ n->id() ] != nullptr,
                           "Error, cannot insert new node on top of existing node.");
    }
  else
    {
      // Allocate just enough space to store the new node.  This will
      // cause highly non-ideal memory allocation behavior if called
      // repeatedly...
      _nodes.resize(n->id() + 1);
    }

#ifdef LIBMESH_ENABLE_UNIQUE_ID
  if (!n->valid_unique_id())
    n->set_unique_id(_next_unique_id++);
  else
   _next_unique_id = std::max(_next_unique_id, n->unique_id()+1);
#endif

  n->add_extra_integers(_node_integer_names.size(),
                        _node_integer_default_values);

  // We have enough space and this spot isn't already occupied by
  // another node, so go ahead and add it.
  ++_n_nodes;
  _nodes[ n->id() ] = n;

  // If we made it this far, we just inserted the node the user handed
  // us, so we can give it right back.
  return n;
}

Node * ReplicatedMesh::insert_node(std::unique_ptr<Node> n)
{
  // The mesh now takes ownership of the Node. Eventually the guts of
  // insert_node(Node*) will get moved to a private helper function, and
  // calling insert_node(Node*) directly will be deprecated.
  return insert_node(n.release());
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
  // This doesn't get used in serial yet
  Node * nd = _nodes[old_id];
  libmesh_assert (nd);

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
  // There is no need to remove the elements from
  // the BoundaryInfo data structure since we
  // already cleared it.
  for (auto & elem : _elements)
    delete elem;

  _n_elem = 0;
  _elements.clear();

  // clear the nodes data structure
  // There is no need to remove the nodes from
  // the BoundaryInfo data structure since we
  // already cleared it.
  for (auto & node : _nodes)
    delete node;

  _n_nodes = 0;
  _nodes.clear();
}



void ReplicatedMesh::update_parallel_id_counts()
{
#ifdef LIBMESH_ENABLE_UNIQUE_ID
  _next_unique_id = this->parallel_max_unique_id();
#endif
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
            if (connected_nodes.find(nd) != connected_nodes.end())
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


void ReplicatedMesh::stitch_meshes (const ReplicatedMesh & other_mesh,
                                    boundary_id_type this_mesh_boundary_id,
                                    boundary_id_type other_mesh_boundary_id,
                                    Real tol,
                                    bool clear_stitched_boundary_ids,
                                    bool verbose,
                                    bool use_binary_search,
                                    bool enforce_all_nodes_match_on_boundaries)
{
  LOG_SCOPE("stitch_meshes()", "ReplicatedMesh");
  stitching_helper(&other_mesh,
                   this_mesh_boundary_id,
                   other_mesh_boundary_id,
                   tol,
                   clear_stitched_boundary_ids,
                   verbose,
                   use_binary_search,
                   enforce_all_nodes_match_on_boundaries,
                   true);
}

void ReplicatedMesh::stitch_surfaces (boundary_id_type boundary_id_1,
                                      boundary_id_type boundary_id_2,
                                      Real tol,
                                      bool clear_stitched_boundary_ids,
                                      bool verbose,
                                      bool use_binary_search,
                                      bool enforce_all_nodes_match_on_boundaries)
{
  stitching_helper(nullptr,
                   boundary_id_1,
                   boundary_id_2,
                   tol,
                   clear_stitched_boundary_ids,
                   verbose,
                   use_binary_search,
                   enforce_all_nodes_match_on_boundaries,
                   true);
}

void ReplicatedMesh::stitching_helper (const ReplicatedMesh * other_mesh,
                                       boundary_id_type this_mesh_boundary_id,
                                       boundary_id_type other_mesh_boundary_id,
                                       Real tol,
                                       bool clear_stitched_boundary_ids,
                                       bool verbose,
                                       bool use_binary_search,
                                       bool enforce_all_nodes_match_on_boundaries,
                                       bool skip_find_neighbors)
{
  std::map<dof_id_type, dof_id_type> node_to_node_map, other_to_this_node_map; // The second is the inverse map of the first
  std::map<dof_id_type, std::vector<dof_id_type>> node_to_elems_map;

  typedef dof_id_type                     key_type;
  typedef std::pair<Elem *, unsigned char> val_type;
  typedef std::pair<key_type, val_type>   key_val_pair;
  typedef std::unordered_multimap<key_type, val_type> map_type;
  // Mapping between all side keys in this mesh and elements+side numbers relevant to the boundary in this mesh as well.
  map_type side_to_elem_map;

  // If there is only one mesh (i.e. other_mesh == nullptr), then loop over this mesh twice
  if (!other_mesh)
    {
      other_mesh = this;
    }

  if ((this_mesh_boundary_id  != BoundaryInfo::invalid_id) &&
      (other_mesh_boundary_id != BoundaryInfo::invalid_id))
    {
      LOG_SCOPE("stitch_meshes node merging", "ReplicatedMesh");

      // While finding nodes on the boundary, also find the minimum edge length
      // of all faces on both boundaries.  This will later be used in relative
      // distance checks when stitching nodes.
      Real h_min = std::numeric_limits<Real>::max();
      bool h_min_updated = false;

      // Loop below fills in these sets for the two meshes.
      std::set<dof_id_type> this_boundary_node_ids, other_boundary_node_ids;

      // Pull objects out of the loop to reduce heap operations
      std::unique_ptr<Elem> side;

      {
        // Make temporary fixed-size arrays for loop
        boundary_id_type id_array[2]         = {this_mesh_boundary_id, other_mesh_boundary_id};
        std::set<dof_id_type> * set_array[2] = {&this_boundary_node_ids, &other_boundary_node_ids};
        const ReplicatedMesh * mesh_array[2] = {this, other_mesh};

        for (unsigned i=0; i<2; ++i)
          {
            // First we deal with node boundary IDs.
            // We only enter this loop if we have at least one
            // nodeset.
            if (mesh_array[i]->get_boundary_info().n_nodeset_conds() > 0)
              {
                // build_node_list() returns a vector of (node-id, bc-id) tuples
                for (const auto & t : mesh_array[i]->get_boundary_info().build_node_list())
                  {
                    boundary_id_type node_bc_id = std::get<1>(t);
                    if (node_bc_id == id_array[i])
                      {
                        dof_id_type this_node_id = std::get<0>(t);
                        set_array[i]->insert( this_node_id );

                        // We need to set h_min to some value. It's too expensive to
                        // search for the element that actually contains this node,
                        // since that would require a PointLocator. As a result, we
                        // just use the first (non-NodeElem!) element in the mesh to
                        // give us hmin if it's never been set before.
                        if (!h_min_updated)
                          {
                            for (const auto & elem : mesh_array[i]->active_element_ptr_range())
                              {
                                Real current_h_min = elem->hmin();
                                if (current_h_min > 0.)
                                  {
                                    h_min = current_h_min;
                                    h_min_updated = true;
                                    break;
                                  }
                              }

                            // If, after searching all the active elements, we did not update
                            // h_min, give up and set h_min to 1 so that we don't repeat this
                            // fruitless search
                            if (!h_min_updated)
                              {
                                h_min_updated = true;
                                h_min = 1.0;
                              }
                          }
                      }
                  }
              }

            // Container to catch boundary IDs passed back from BoundaryInfo.
            std::vector<boundary_id_type> bc_ids;

            for (auto & el : mesh_array[i]->element_ptr_range())
              {
                // Now check whether elem has a face on the specified boundary
                for (auto side_id : el->side_index_range())
                  if (el->neighbor_ptr(side_id) == nullptr)
                    {
                      // Get *all* boundary IDs on this side, not just the first one!
                      mesh_array[i]->get_boundary_info().boundary_ids (el, side_id, bc_ids);

                      if (std::find(bc_ids.begin(), bc_ids.end(), id_array[i]) != bc_ids.end())
                        {
                          el->build_side_ptr(side, side_id);
                          for (auto & n : side->node_ref_range())
                            set_array[i]->insert(n.id());

                          h_min = std::min(h_min, side->hmin());
                          h_min_updated = true;

                          // This side is on the boundary, add its information to side_to_elem
                          if (skip_find_neighbors && (i==0))
                            {
                              key_type key = el->key(side_id);
                              val_type val;
                              val.first = el;
                              val.second = cast_int<unsigned char>(side_id);

                              key_val_pair kvp;
                              kvp.first = key;
                              kvp.second = val;
                              side_to_elem_map.insert (kvp);
                            }
                        }

                      // Also, check the edges on this side. We don't have to worry about
                      // updating neighbor info in this case since elements don't store
                      // neighbor info on edges.
                      for (auto edge_id : el->edge_index_range())
                        {
                          if (el->is_edge_on_side(edge_id, side_id))
                            {
                              // Get *all* boundary IDs on this edge, not just the first one!
                              mesh_array[i]->get_boundary_info().edge_boundary_ids (el, edge_id, bc_ids);

                              if (std::find(bc_ids.begin(), bc_ids.end(), id_array[i]) != bc_ids.end())
                                {
                                  std::unique_ptr<Elem> edge (el->build_edge_ptr(edge_id));
                                  for (auto & n : edge->node_ref_range())
                                    set_array[i]->insert( n.id() );

                                  h_min = std::min(h_min, edge->hmin());
                                  h_min_updated = true;
                                }
                            }
                        }
                    }
              }
          }
      }

      if (verbose)
        {
          libMesh::out << "In ReplicatedMesh::stitch_meshes:\n"
                       << "This mesh has "  << this_boundary_node_ids.size()
                       << " nodes on boundary " << this_mesh_boundary_id  << ".\n"
                       << "Other mesh has " << other_boundary_node_ids.size()
                       << " nodes on boundary " << other_mesh_boundary_id << ".\n";

          if (h_min_updated)
            {
              libMesh::out << "Minimum edge length on both surfaces is " << h_min << ".\n";
            }
          else
            {
              libMesh::out << "No elements on specified surfaces." << std::endl;
            }
        }

      // At this point, if h_min==0 it means that there were at least two coincident
      // nodes on the surfaces being stitched, and we don't currently support that case.
      // (It might be possible to support, but getting it exactly right would be tricky
      // and probably not worth the extra complications to the "normal" case.)
      libmesh_error_msg_if(h_min < std::numeric_limits<Real>::epsilon(),
                           "Coincident nodes detected on source and/or target "
                           "surface, stitching meshes is not possible.");

      // We require nanoflann for the "binary search" (really kd-tree)
      // option to work. If it's not available, turn that option off,
      // warn the user, and fall back on the N^2 search algorithm.
      if (use_binary_search)
        {
#ifndef LIBMESH_HAVE_NANOFLANN
          use_binary_search = false;
          libmesh_warning("The use_binary_search option in the "
                          "ReplicatedMesh stitching algorithms requires nanoflann "
                          "support. Falling back on N^2 search algorithm.");
#endif
        }

      if (this_boundary_node_ids.size())
      {
        if (use_binary_search)
        {
#ifdef LIBMESH_HAVE_NANOFLANN
          typedef nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<Real, VectorOfNodesAdaptor>, VectorOfNodesAdaptor, 3> kd_tree_t;

          // Create the dataset needed to build the kd tree with nanoflann
          std::vector<std::pair<Point, dof_id_type>> this_mesh_nodes(this_boundary_node_ids.size());
          std::set<dof_id_type>::iterator current_node = this_boundary_node_ids.begin(),
                                          node_ids_end = this_boundary_node_ids.end();
          for (unsigned int ctr = 0; current_node != node_ids_end; ++current_node, ++ctr)
          {
            this_mesh_nodes[ctr].first = this->point(*current_node);
            this_mesh_nodes[ctr].second = *current_node;
          }

          VectorOfNodesAdaptor vec_nodes_adaptor(this_mesh_nodes);

          kd_tree_t this_kd_tree(3, vec_nodes_adaptor, 10);
          this_kd_tree.buildIndex();

          // Storage for nearest neighbor in the loop below
          std::vector<size_t> ret_index(1);
          std::vector<Real> ret_dist_sqr(1);

          // Loop over other mesh. For each node, find its nearest neighbor in this mesh, and fill in the maps.
          for (auto node : other_boundary_node_ids)
          {
            const Real query_pt[] = {other_mesh->point(node)(0), other_mesh->point(node)(1), other_mesh->point(node)(2)};
            this_kd_tree.knnSearch(&query_pt[0], 1, &ret_index[0], &ret_dist_sqr[0]);
            if (ret_dist_sqr[0] < TOLERANCE*TOLERANCE)
            {
              node_to_node_map[this_mesh_nodes[ret_index[0]].second] = node;
              other_to_this_node_map[node] = this_mesh_nodes[ret_index[0]].second;
            }
          }

          // If the 2 maps don't have the same size, it means we have overwritten a value in node_to_node_map
          // It means one node in this mesh is the nearest neighbor of several nodes in other mesh.
          // Not possible !
          libmesh_error_msg_if(node_to_node_map.size() != other_to_this_node_map.size(),
                               "Error: Found multiple matching nodes in stitch_meshes");
#endif
        }
        else
        {
          // In the unlikely event that two meshes composed entirely of
          // NodeElems are being stitched together, we will not have
          // selected a valid h_min value yet, and the distance
          // comparison below will be true for essentially any two
          // nodes. In this case we simply fall back on an absolute
          // distance check.
          if (!h_min_updated)
            {
              libmesh_warning("No valid h_min value was found, falling back on "
                              "absolute distance check in the N^2 search algorithm.");
              h_min = 1.;
            }

          // Otherwise, use a simple N^2 search to find the closest matching points. This can be helpful
          // in the case that we have tolerance issues which cause mismatch between the two surfaces
          // that are being stitched.
          for (const auto & this_node_id : this_boundary_node_ids)
          {
            Node & this_node = this->node_ref(this_node_id);

            bool found_matching_nodes = false;

            for (const auto & other_node_id : other_boundary_node_ids)
            {
              const Node & other_node = other_mesh->node_ref(other_node_id);

              Real node_distance = (this_node - other_node).norm();

              if (node_distance < tol*h_min)
              {
                // Make sure we didn't already find a matching node!
                libmesh_error_msg_if(found_matching_nodes,
                                     "Error: Found multiple matching nodes in stitch_meshes");

                node_to_node_map[this_node_id] = other_node_id;
                other_to_this_node_map[other_node_id] = this_node_id;

                found_matching_nodes = true;
              }
            }
          }
        }
      }

      // Build up the node_to_elems_map, using only one loop over other_mesh
      for (auto & el : other_mesh->element_ptr_range())
        {
          // For each node on the element, find the corresponding node
          // on "this" Mesh, 'this_node_id', if it exists, and push
          // the current element ID back onto node_to_elems_map[this_node_id].
          // For that we will use the reverse mapping we created at
          // the same time as the forward mapping.
          for (auto & n : el->node_ref_range())
            {
              dof_id_type other_node_id = n.id();
              std::map<dof_id_type, dof_id_type>::iterator it =
                other_to_this_node_map.find(other_node_id);

              if (it != other_to_this_node_map.end())
                {
                  dof_id_type this_node_id = it->second;
                  node_to_elems_map[this_node_id].push_back( el->id() );
                }
            }
        }

      if (verbose)
        {
          libMesh::out << "In ReplicatedMesh::stitch_meshes:\n"
                       << "Found " << node_to_node_map.size()
                       << " matching nodes.\n"
                       << std::endl;
        }

      if (enforce_all_nodes_match_on_boundaries)
        {
          std::size_t n_matching_nodes = node_to_node_map.size();
          std::size_t this_mesh_n_nodes = this_boundary_node_ids.size();
          std::size_t other_mesh_n_nodes = other_boundary_node_ids.size();
          libmesh_error_msg_if((n_matching_nodes != this_mesh_n_nodes) || (n_matching_nodes != other_mesh_n_nodes),
                               "Error: We expected the number of nodes to match.");
        }
    }
  else
    {
      if (verbose)
        {
          libMesh::out << "Skip node merging in ReplicatedMesh::stitch_meshes:" << std::endl;
        }
    }

  dof_id_type node_delta = this->max_node_id();
  dof_id_type elem_delta = this->max_elem_id();

  unique_id_type unique_delta =
#ifdef LIBMESH_ENABLE_UNIQUE_ID
    this->parallel_max_unique_id();
#else
    0;
#endif

  // If other_mesh != nullptr, then we have to do a bunch of work
  // in order to copy it to this mesh
  if (this!=other_mesh)
    {
      LOG_SCOPE("stitch_meshes copying", "ReplicatedMesh");

      // Increment the node_to_node_map and node_to_elems_map
      // to account for id offsets
      for (auto & pr : node_to_node_map)
        pr.second += node_delta;

      for (auto & pr : node_to_elems_map)
        for (auto & entry : pr.second)
          entry += elem_delta;

      // Copy mesh data. If we skip the call to find_neighbors(), the lists
      // of neighbors will be copied verbatim from the other mesh
      this->copy_nodes_and_elements(*other_mesh, skip_find_neighbors,
                                    elem_delta, node_delta,
                                    unique_delta);

      // Copy BoundaryInfo from other_mesh too.  We do this via the
      // list APIs rather than element-by-element for speed.
      BoundaryInfo & boundary = this->get_boundary_info();
      const BoundaryInfo & other_boundary = other_mesh->get_boundary_info();

      for (const auto & t : other_boundary.build_node_list())
        boundary.add_node(std::get<0>(t) + node_delta,
                          std::get<1>(t));

      for (const auto & t : other_boundary.build_side_list())
        boundary.add_side(std::get<0>(t) + elem_delta,
                          std::get<1>(t),
                          std::get<2>(t));

      for (const auto & t : other_boundary.build_edge_list())
        boundary.add_edge(std::get<0>(t) + elem_delta,
                          std::get<1>(t),
                          std::get<2>(t));

      for (const auto & t : other_boundary.build_shellface_list())
        boundary.add_shellface(std::get<0>(t) + elem_delta,
                               std::get<1>(t),
                               std::get<2>(t));

      const auto & other_ns_id_to_name = other_boundary.get_nodeset_name_map();
      auto & ns_id_to_name = boundary.set_nodeset_name_map();
      ns_id_to_name.insert(other_ns_id_to_name.begin(), other_ns_id_to_name.end());

      const auto & other_ss_id_to_name = other_boundary.get_sideset_name_map();
      auto & ss_id_to_name = boundary.set_sideset_name_map();
      ss_id_to_name.insert(other_ss_id_to_name.begin(), other_ss_id_to_name.end());

      const auto & other_es_id_to_name = other_boundary.get_edgeset_name_map();
      auto & es_id_to_name = boundary.set_edgeset_name_map();
      es_id_to_name.insert(other_es_id_to_name.begin(), other_es_id_to_name.end());
    } // end if (other_mesh)

  // Finally, we need to "merge" the overlapping nodes
  // We do this by iterating over node_to_elems_map and updating
  // the elements so that they "point" to the nodes that came
  // from this mesh, rather than from other_mesh.
  // Then we iterate over node_to_node_map and delete the
  // duplicate nodes that came from other_mesh.

  {
    LOG_SCOPE("stitch_meshes node updates", "ReplicatedMesh");

    // Container to catch boundary IDs passed back from BoundaryInfo.
    std::vector<boundary_id_type> bc_ids;

    for (const auto & pr : node_to_elems_map)
      {
        dof_id_type target_node_id = pr.first;
        dof_id_type other_node_id = node_to_node_map[target_node_id];
        Node & target_node = this->node_ref(target_node_id);

        std::size_t n_elems = pr.second.size();
        for (std::size_t i=0; i<n_elems; i++)
          {
            dof_id_type elem_id = pr.second[i];
            Elem * el = this->elem_ptr(elem_id);

            // find the local node index that we want to update
            unsigned int local_node_index = el->local_node(other_node_id);
            libmesh_assert_not_equal_to(local_node_index, libMesh::invalid_uint);

            // We also need to copy over the nodeset info here,
            // because the node will get deleted below
            this->get_boundary_info().boundary_ids(el->node_ptr(local_node_index), bc_ids);
            el->set_node(local_node_index) = &target_node;
            this->get_boundary_info().add_node(&target_node, bc_ids);
          }
      }
  }

  {
    LOG_SCOPE("stitch_meshes node deletion", "ReplicatedMesh");
    for (const auto & pr : node_to_node_map)
      {
        // In the case that this==other_mesh, the two nodes might be the same (e.g. if
        // we're stitching a "sliver"), hence we need to skip node deletion in that case.
        if ((this == other_mesh) && (pr.second == pr.first))
          continue;

        dof_id_type this_node_id = pr.second;
        this->delete_node( this->node_ptr(this_node_id) );
      }
  }

  // If find_neighbors() wasn't called in prepare_for_use(), we need to
  // manually loop once more over all elements adjacent to the stitched boundary
  // and fix their lists of neighbors.
  // This is done according to the following steps:
  //   1. Loop over all copied elements adjacent to the boundary using node_to_elems_map (trying to avoid duplicates)
  //   2. Look at all their sides with a nullptr neighbor and update them using side_to_elem_map if necessary
  //   3. Update the corresponding side in side_to_elem_map as well
  if (skip_find_neighbors)
    {
      LOG_SCOPE("stitch_meshes neighbor fixes", "ReplicatedMesh");

      // Pull objects out of the loop to reduce heap operations
      std::unique_ptr<Elem> my_side, their_side;

      std::set<dof_id_type> fixed_elems;
      for (const auto & pr : node_to_elems_map)
        {
          std::size_t n_elems = pr.second.size();
          for (std::size_t i=0; i<n_elems; i++)
            {
              dof_id_type elem_id = pr.second[i];
              if (fixed_elems.find(elem_id) == fixed_elems.end())
                {
                  Elem * el = this->elem_ptr(elem_id);
                  fixed_elems.insert(elem_id);
                  for (auto s : el->side_index_range())
                    {
                      if (el->neighbor_ptr(s) == nullptr)
                        {
                          key_type key = el->key(s);
                          auto bounds = side_to_elem_map.equal_range(key);

                          if (bounds.first != bounds.second)
                            {
                              // Get the side for this element
                              el->side_ptr(my_side, s);

                              // Look at all the entries with an equivalent key
                              while (bounds.first != bounds.second)
                                {
                                  // Get the potential element
                                  Elem * neighbor = bounds.first->second.first;

                                  // Get the side for the neighboring element
                                  const unsigned int ns = bounds.first->second.second;
                                  neighbor->side_ptr(their_side, ns);
                                  //libmesh_assert(my_side.get());
                                  //libmesh_assert(their_side.get());

                                  // If found a match with my side
                                  //
                                  // We need special tests here for 1D:
                                  // since parents and children have an equal
                                  // side (i.e. a node), we need to check
                                  // ns != ms, and we also check level() to
                                  // avoid setting our neighbor pointer to
                                  // any of our neighbor's descendants
                                  if ((*my_side == *their_side) &&
                                      (el->level() == neighbor->level()) &&
                                      ((el->dim() != 1) || (ns != s)))
                                    {
                                      // So share a side.  Is this a mixed pair
                                      // of subactive and active/ancestor
                                      // elements?
                                      // If not, then we're neighbors.
                                      // If so, then the subactive's neighbor is

                                      if (el->subactive() ==
                                          neighbor->subactive())
                                        {
                                          // an element is only subactive if it has
                                          // been coarsened but not deleted
                                          el->set_neighbor (s,neighbor);
                                          neighbor->set_neighbor(ns,el);
                                        }
                                      else if (el->subactive())
                                        {
                                          el->set_neighbor(s,neighbor);
                                        }
                                      else if (neighbor->subactive())
                                        {
                                          neighbor->set_neighbor(ns,el);
                                        }
                                      // It's OK to invalidate the
                                      // bounds.first iterator here,
                                      // as we are immediately going
                                      // to break out of this while
                                      // loop. bounds.first will
                                      // therefore not be used for
                                      // anything else.
                                      side_to_elem_map.erase (bounds.first);
                                      break;
                                    }

                                  ++bounds.first;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

  const bool old_allow_find_neighbors = this->allow_find_neighbors();
  this->allow_find_neighbors(!skip_find_neighbors);
  this->prepare_for_use();
  this->allow_find_neighbors(old_allow_find_neighbors);

  // After the stitching, we may want to clear boundary IDs from element
  // faces that are now internal to the mesh
  if (clear_stitched_boundary_ids)
    {
      LOG_SCOPE("stitch_meshes clear bcids", "ReplicatedMesh");

      this->get_boundary_info().clear_stitched_boundary_side_ids(
          this_mesh_boundary_id, other_mesh_boundary_id, /*clear_nodeset_data=*/true);
    }
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
  std::vector<Elem *> list;
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
      Elem * elem = list.back(); list.pop_back(); ++visited;

      min_id = std::min(elem->id(), min_id);

      for (auto s : elem->side_index_range())
      {
        Elem * neighbor = elem->neighbor_ptr(s);
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
      libmesh_assert(boundary_elements.find(side) != boundary_elements.end());
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

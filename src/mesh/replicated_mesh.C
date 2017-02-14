// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#include LIBMESH_INCLUDE_UNORDERED_MAP
#include LIBMESH_INCLUDE_UNORDERED_SET
#include LIBMESH_INCLUDE_HASH
LIBMESH_DEFINE_HASH_POINTERS


namespace
{
using namespace libMesh;

// A custom comparison function, based on Point::operator<,
// that tries to ignore floating point differences in components
// of the point
class FuzzyPointCompare
{
private:
  Real _tol;

public:
  // Constructor takes the tolerance to be used in fuzzy comparisons
  FuzzyPointCompare(Real tol) : _tol(tol) {}

  // This is inspired directly by Point::operator<
  bool operator()(const Point & lhs, const Point & rhs)
  {
    for (unsigned i=0; i<LIBMESH_DIM; ++i)
      {
        // If the current components are within some tolerance
        // of one another, then don't attempt the less-than comparison.
        // Note that this may cause something strange to happen, as Roy
        // believes he can prove it is not a total ordering...
        Real rel_size = std::max(std::abs(lhs(i)), std::abs(rhs(i)));

        // Don't use relative tolerance if both numbers are already small.
        // How small?  Some possible options are:
        // * std::numeric_limits<Real>::epsilon()
        // * TOLERANCE
        // * 1.0
        // If we use std::numeric_limits<Real>::epsilon(), we'll
        // do more relative comparisons for small numbers, but
        // increase the chance for false positives?  If we pick 1.0,
        // we'll never "increase" the difference between small numbers
        // in the test below.
        if (rel_size < 1.)
          rel_size = 1.;

        // Don't attempt the comparison if lhs(i) and rhs(i) are too close
        // together.
        if ( std::abs(lhs(i) - rhs(i)) / rel_size < _tol)
          continue;

        if (lhs(i) < rhs(i))
          return true;
        if (lhs(i) > rhs(i))
          return false;
      }

    // We compared all the components without returning yet, so
    // each component was neither greater than nor less than they other.
    // They might be equal, so return false.
    return false;
  }

  // Needed by std::sort on vector< pair<Point,id> >
  bool operator()(const std::pair<Point, dof_id_type> & lhs,
                  const std::pair<Point, dof_id_type> & rhs)
  {
    return (*this)(lhs.first, rhs.first);
  }

  // Comparsion function where lhs is a Point and rhs is a pair<Point,dof_id_type>.
  // This is used in routines like lower_bound, where a specific value is being
  // searched for.
  bool operator()(const Point & lhs, std::pair<Point, dof_id_type> & rhs)
  {
    return (*this)(lhs, rhs.first);
  }

  // And the other way around...
  bool operator()(std::pair<Point, dof_id_type> & lhs, const Point & rhs)
  {
    return (*this)(lhs.first, rhs);
  }
};
}



namespace libMesh
{

// ------------------------------------------------------------
// ReplicatedMesh class member functions
ReplicatedMesh::ReplicatedMesh (const Parallel::Communicator & comm_in,
                                unsigned char d) :
  UnstructuredMesh (comm_in,d)
{
#ifdef LIBMESH_ENABLE_UNIQUE_ID
  // In serial we just need to reset the next unique id to zero
  // here in the constructor.
  _next_unique_id = 0;
#endif
  _partitioner = UniquePtr<Partitioner>(new MetisPartitioner());
}



#ifndef LIBMESH_DISABLE_COMMWORLD
ReplicatedMesh::ReplicatedMesh (unsigned char d) :
  UnstructuredMesh (d)
{
#ifdef LIBMESH_ENABLE_UNIQUE_ID
  // In serial we just need to reset the next unique id to zero
  // here in the constructor.
  _next_unique_id = 0;
#endif
  _partitioner = UniquePtr<Partitioner>(new MetisPartitioner());
}
#endif


ReplicatedMesh::~ReplicatedMesh ()
{
  this->clear();  // Free nodes and elements
}


// This might be specialized later, but right now it's just here to
// make sure the compiler doesn't give us a default (non-deep) copy
// constructor instead.
ReplicatedMesh::ReplicatedMesh (const ReplicatedMesh & other_mesh) :
  UnstructuredMesh (other_mesh)
{
  this->copy_nodes_and_elements(other_mesh);
  this->get_boundary_info() = other_mesh.get_boundary_info();
#ifdef LIBMESH_ENABLE_UNIQUE_ID
  this->_next_unique_id = other_mesh._next_unique_id;
#endif
}


ReplicatedMesh::ReplicatedMesh (const UnstructuredMesh & other_mesh) :
  UnstructuredMesh (other_mesh)
{
  this->copy_nodes_and_elements(other_mesh);
  this->get_boundary_info() = other_mesh.get_boundary_info();
}


const Point & ReplicatedMesh::point (const dof_id_type i) const
{
  return this->node_ref(i);
}




const Node * ReplicatedMesh::node_ptr (const dof_id_type i) const
{
  libmesh_assert_less (i, this->n_nodes());
  libmesh_assert(_nodes[i]);
  libmesh_assert_equal_to (_nodes[i]->id(), i); // This will change soon

  return _nodes[i];
}




Node * ReplicatedMesh::node_ptr (const dof_id_type i)
{
  libmesh_assert_less (i, this->n_nodes());
  libmesh_assert(_nodes[i]);
  libmesh_assert_equal_to (_nodes[i]->id(), i); // This will change soon

  return _nodes[i];
}




const Node * ReplicatedMesh::query_node_ptr (const dof_id_type i) const
{
  if (i >= this->n_nodes())
    return libmesh_nullptr;
  libmesh_assert (_nodes[i] == libmesh_nullptr ||
                  _nodes[i]->id() == i); // This will change soon

  return _nodes[i];
}




Node * ReplicatedMesh::query_node_ptr (const dof_id_type i)
{
  if (i >= this->n_nodes())
    return libmesh_nullptr;
  libmesh_assert (_nodes[i] == libmesh_nullptr ||
                  _nodes[i]->id() == i); // This will change soon

  return _nodes[i];
}




const Elem * ReplicatedMesh::elem_ptr (const dof_id_type i) const
{
  libmesh_assert_less (i, this->n_elem());
  libmesh_assert(_elements[i]);
  libmesh_assert_equal_to (_elements[i]->id(), i); // This will change soon

  return _elements[i];
}




Elem * ReplicatedMesh::elem_ptr (const dof_id_type i)
{
  libmesh_assert_less (i, this->n_elem());
  libmesh_assert(_elements[i]);
  libmesh_assert_equal_to (_elements[i]->id(), i); // This will change soon

  return _elements[i];
}




const Elem * ReplicatedMesh::query_elem_ptr (const dof_id_type i) const
{
  if (i >= this->n_elem())
    return libmesh_nullptr;
  libmesh_assert (_elements[i] == libmesh_nullptr ||
                  _elements[i]->id() == i); // This will change soon

  return _elements[i];
}




Elem * ReplicatedMesh::query_elem_ptr (const dof_id_type i)
{
  if (i >= this->n_elem())
    return libmesh_nullptr;
  libmesh_assert (_elements[i] == libmesh_nullptr ||
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
    e->set_unique_id() = _next_unique_id++;
#endif

  const dof_id_type id = e->id();

  if (id < _elements.size())
    {
      // Overwriting existing elements is still probably a mistake.
      libmesh_assert(!_elements[id]);
    }
  else
    {
      _elements.resize(id+1, libmesh_nullptr);
    }

  _elements[id] = e;

  return e;
}



Elem * ReplicatedMesh::insert_elem (Elem * e)
{
#ifdef LIBMESH_ENABLE_UNIQUE_ID
  if (!e->valid_unique_id())
    e->set_unique_id() = _next_unique_id++;
#endif

  dof_id_type eid = e->id();
  libmesh_assert_less (eid, _elements.size());
  Elem * oldelem = _elements[eid];

  if (oldelem)
    {
      libmesh_assert_equal_to (oldelem->id(), eid);
      this->delete_elem(oldelem);
    }

  _elements[e->id()] = e;

  return e;
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
  delete e;

  // explicitly NULL the pointer
  *pos = libmesh_nullptr;
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
  _elements[old_id] = libmesh_nullptr;
}



Node * ReplicatedMesh::add_point (const Point & p,
                                  const dof_id_type id,
                                  const processor_id_type proc_id)
{
  //   // We only append points with ReplicatedMesh
  //   libmesh_assert(id == DofObject::invalid_id || id == _nodes.size());
  //   Node *n = Node::build(p, _nodes.size()).release();
  //   n->processor_id() = proc_id;
  //   _nodes.push_back (n);

  Node * n = libmesh_nullptr;

  // If the user requests a valid id, either
  // provide the existing node or resize the container
  // to fit the new node.
  if (id != DofObject::invalid_id)
    if (id < _nodes.size())
      n = _nodes[id];
    else
      _nodes.resize(id+1);
  else
    _nodes.push_back (static_cast<Node *>(libmesh_nullptr));

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

#ifdef LIBMESH_ENABLE_UNIQUE_ID
      if (!n->valid_unique_id())
        n->set_unique_id() = _next_unique_id++;
#endif

      if (id == DofObject::invalid_id)
        _nodes.back() = n;
      else
        _nodes[id] = n;
    }

  // better not pass back a NULL pointer.
  libmesh_assert (n);

  return n;
}



Node * ReplicatedMesh::add_node (Node * n)
{
  libmesh_assert(n);
  // We only append points with ReplicatedMesh
  libmesh_assert(!n->valid_id() || n->id() == _nodes.size());

  n->set_id (cast_int<dof_id_type>(_nodes.size()));

#ifdef LIBMESH_ENABLE_UNIQUE_ID
  if (!n->valid_unique_id())
    n->set_unique_id() = _next_unique_id++;
#endif

  _nodes.push_back(n);

  return n;
}



Node * ReplicatedMesh::insert_node(Node * n)
{
  if (!n)
    libmesh_error_msg("Error, attempting to insert NULL node.");

  if (n->id() == DofObject::invalid_id)
    libmesh_error_msg("Error, cannot insert node with invalid id.");

  if (n->id() < _nodes.size())
    {
      // Don't allow inserting on top of an existing Node.

      // Doing so doesn't have to be *error*, in the case where a
      // redundant insert is done, but when that happens we ought to
      // always be able to make the code more efficient by avoiding
      // the redundant insert, so let's keep screaming "Error" here.
      if (_nodes[ n->id() ] != libmesh_nullptr)
        libmesh_error_msg("Error, cannot insert node on top of existing node.");
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
    n->set_unique_id() = _next_unique_id++;
#endif

  // We have enough space and this spot isn't already occupied by
  // another node, so go ahead and add it.
  _nodes[ n->id() ] = n;

  // If we made it this far, we just inserted the node the user handed
  // us, so we can give it right back.
  return n;
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

  // delete the node
  delete n;

  // explicitly NULL the pointer
  *pos = libmesh_nullptr;
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
  _nodes[old_id] = libmesh_nullptr;
}



void ReplicatedMesh::clear ()
{
  // Call parent clear function
  MeshBase::clear();


  // Clear our elements and nodes
  {
    std::vector<Elem *>::iterator       it  = _elements.begin();
    const std::vector<Elem *>::iterator end = _elements.end();

    // There is no need to remove the elements from
    // the BoundaryInfo data structure since we
    // already cleared it.
    for (; it != end; ++it)
      delete *it;

    _elements.clear();
  }

  // clear the nodes data structure
  {
    std::vector<Node *>::iterator       it  = _nodes.begin();
    const std::vector<Node *>::iterator end = _nodes.end();

    // There is no need to remove the nodes from
    // the BoundaryInfo data structure since we
    // already cleared it.
    for (; it != end; ++it)
      delete *it;

    _nodes.clear();
  }
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
#endif



void ReplicatedMesh::renumber_nodes_and_elements ()
{
  LOG_SCOPE("renumber_nodes_and_elem()", "Mesh");

  // Even if we're not allowed to renumber, we can still delete nodes
  // which have been "orphaned" due to coarsening, and remove NULL
  // entries off the end of the _elements and _nodes vectors, since
  // that does not require renumbering.
  if (_skip_renumber_nodes_and_elements)
    {
      // Build set of nodes that are currently connected to elements.
      LIBMESH_BEST_UNORDERED_SET<Node *> connected_nodes;

      std::vector<Elem *>::iterator elem_it = _elements.begin();
      const std::vector<Elem *>::iterator elem_end = _elements.end();
      for (; elem_it != elem_end; ++elem_it)
        {
          Elem * elem = *elem_it;

          if (elem)
            {
              // Add this element's nodes to the connected list.
              for (unsigned int n=0; n<elem->n_nodes(); n++)
                connected_nodes.insert(elem->node_ptr(n));
            }
        }

      // Delete (and leave NULL entries for) the unconnected nodes.
      std::vector<Node *>::iterator node_it = _nodes.begin();
      const std::vector<Node *>::iterator node_end = _nodes.end();
      for (; node_it != node_end; ++node_it)
        {
          // Get a reference to the pointer in the actual vector, so
          // we can potentially delete it and reassign its value to
          // NULL, rather than setting a _copy_ of the pointer's value
          // to NULL.
          Node *& node = *node_it;

          if (node)
            {
              // If this node is not connected, delete it.  Note that
              // it is *not* being erased from the _nodes vector.
              if (connected_nodes.find(node) == connected_nodes.end())
                {
                  this->get_boundary_info().remove (node);
                  delete node;
                  node = libmesh_nullptr;
                }
            }
        }

      // Now, actually erase any NULL entries at the end of the
      // _elements and _nodes vectors.

      // Find the first non-NULL Elem, starting from the end.
      std::vector<Elem *>::reverse_iterator last_elem = _elements.rbegin();
      while (last_elem != _elements.rend() && *last_elem == libmesh_nullptr)
        ++last_elem;

      // Remove trailing NULL entries off the end of the _elements vector.
      _elements.erase(last_elem.base(), _elements.end());

      // Find the first non-NULL Node, starting from the end.
      std::vector<Node *>::reverse_iterator last_node = _nodes.rbegin();
      while (last_node != _nodes.rend() && *last_node == libmesh_nullptr)
        ++last_node;

      // Remove trailing NULL entries off the end of the _nodes vector.
      _nodes.erase(last_node.base(), _nodes.end());

      return;
    }

  // node and element id counters
  dof_id_type next_free_elem = 0;
  dof_id_type next_free_node = 0;

  // Loop over the elements.  Note that there may
  // be NULLs in the _elements vector from the coarsening
  // process.  Pack the elements in to a contiguous array
  // and then trim any excess.
  {
    std::vector<Elem *>::iterator in        = _elements.begin();
    std::vector<Elem *>::iterator out_iter  = _elements.begin();
    const std::vector<Elem *>::iterator end = _elements.end();

    for (; in != end; ++in)
      if (*in != libmesh_nullptr)
        {
          Elem * el = *in;

          *out_iter = *in;
          ++out_iter;

          // Increment the element counter
          el->set_id (next_free_elem++);

          // Loop over this element's nodes.  Number them,
          // if they have not been numbered already.  Also,
          // position them in the _nodes vector so that they
          // are packed contiguously from the beginning.
          for (unsigned int n=0; n<el->n_nodes(); n++)
            if (el->node_id(n) == next_free_node)    // don't need to process
              next_free_node++;                      // [(src == dst) below]

            else if (el->node_id(n) > next_free_node) // need to process
              {
                // The source and destination indices
                // for this node
                const dof_id_type src_idx = el->node_id(n);
                const dof_id_type dst_idx = next_free_node++;

                // ensure we want to swap a valid nodes
                libmesh_assert(_nodes[src_idx]);

                // Swap the source and destination nodes
                std::swap(_nodes[src_idx],
                          _nodes[dst_idx] );

                // Set proper indices where that makes sense
                if (_nodes[src_idx] != libmesh_nullptr)
                  _nodes[src_idx]->set_id (src_idx);
                _nodes[dst_idx]->set_id (dst_idx);
              }
        }

    // Erase any additional storage. These elements have been
    // copied into NULL voids by the procedure above, and are
    // thus repeated and unnecessary.
    _elements.erase (out_iter, end);
  }

  // Any nodes in the vector >= _nodes[next_free_node]
  // are not connected to any elements will now be deleted.
  {
    std::vector<Node *>::iterator nd        = _nodes.begin();
    const std::vector<Node *>::iterator end = _nodes.end();

    std::advance (nd, next_free_node);

    for (std::vector<Node *>::iterator it=nd;
         it != end; ++it)
      {
        // Mesh modification code might have already deleted some
        // nodes
        if (*it == libmesh_nullptr)
          continue;

        // remove any boundary information associated with
        // this node
        this->get_boundary_info().remove (*it);

        // delete the node
        delete *it;
        *it = libmesh_nullptr;
      }

    _nodes.erase (nd, end);
  }

  libmesh_assert_equal_to (next_free_elem, _elements.size());
  libmesh_assert_equal_to (next_free_node, _nodes.size());

  this->update_parallel_id_counts();
}



void ReplicatedMesh::fix_broken_node_and_element_numbering ()
{
  // Nodes first
  for (std::size_t n=0; n<this->_nodes.size(); n++)
    if (this->_nodes[n] != libmesh_nullptr)
      this->_nodes[n]->set_id() = cast_int<dof_id_type>(n);

  // Elements next
  for (std::size_t e=0; e<this->_elements.size(); e++)
    if (this->_elements[e] != libmesh_nullptr)
      this->_elements[e]->set_id() = cast_int<dof_id_type>(e);
}


void ReplicatedMesh::stitch_meshes (ReplicatedMesh & other_mesh,
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
  stitching_helper(libmesh_nullptr,
                   boundary_id_1,
                   boundary_id_2,
                   tol,
                   clear_stitched_boundary_ids,
                   verbose,
                   use_binary_search,
                   enforce_all_nodes_match_on_boundaries,
                   true);
}

void ReplicatedMesh::stitching_helper (ReplicatedMesh * other_mesh,
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
  std::map<dof_id_type, std::vector<dof_id_type> > node_to_elems_map;

  typedef dof_id_type                     key_type;
  typedef std::pair<Elem *, unsigned char> val_type;
  typedef std::pair<key_type, val_type>   key_val_pair;
  typedef LIBMESH_BEST_UNORDERED_MULTIMAP<key_type, val_type> map_type;
  // Mapping between all side keys in this mesh and elements+side numbers relevant to the boundary in this mesh as well.
  map_type side_to_elem_map;

  // If there is only one mesh (i.e. other_mesh==NULL), then loop over this mesh twice
  if(!other_mesh)
    {
      other_mesh = this;
    }

  if( (this_mesh_boundary_id  != BoundaryInfo::invalid_id) &&
      (other_mesh_boundary_id != BoundaryInfo::invalid_id) )
    {
      // While finding nodes on the boundary, also find the minimum edge length
      // of all faces on both boundaries.  This will later be used in relative
      // distance checks when stitching nodes.
      Real h_min = std::numeric_limits<Real>::max();
      bool h_min_updated = false;

      // Loop below fills in these sets for the two meshes.
      std::set<dof_id_type> this_boundary_node_ids, other_boundary_node_ids;
      {
        // Make temporary fixed-size arrays for loop
        boundary_id_type id_array[2]        = {this_mesh_boundary_id, other_mesh_boundary_id};
        std::set<dof_id_type> * set_array[2] = {&this_boundary_node_ids, &other_boundary_node_ids};
        ReplicatedMesh * mesh_array[2]           = {this, other_mesh};

        for (unsigned i=0; i<2; ++i)
          {
            // First we deal with node boundary IDs.
            // We only enter this loop if we have at least one
            // nodeset.
            if(mesh_array[i]->get_boundary_info().n_nodeset_conds() > 0)
              {
                std::vector<numeric_index_type> node_id_list;
                std::vector<boundary_id_type> bc_id_list;

                // Get the list of nodes with associated boundary IDs
                mesh_array[i]->get_boundary_info().build_node_list(node_id_list, bc_id_list);

                for (std::size_t node_index=0; node_index<bc_id_list.size(); node_index++)
                  {
                    boundary_id_type node_bc_id = bc_id_list[node_index];
                    if (node_bc_id == id_array[i])
                      {
                        dof_id_type this_node_id = node_id_list[node_index];
                        set_array[i]->insert( this_node_id );

                        // We need to set h_min to some value. It's too expensive to
                        // search for the element that actually contains this node,
                        // since that would require a PointLocator. As a result, we
                        // just use the first element in the mesh to give us hmin.
                        const Elem * first_active_elem = *mesh_array[i]->active_elements_begin();
                        h_min = first_active_elem->hmin();
                        h_min_updated = true;
                      }
                  }
              }

            // Container to catch boundary IDs passed back from BoundaryInfo.
            std::vector<boundary_id_type> bc_ids;

            MeshBase::element_iterator elem_it  = mesh_array[i]->elements_begin();
            MeshBase::element_iterator elem_end = mesh_array[i]->elements_end();
            for ( ; elem_it != elem_end; ++elem_it)
              {
                Elem * el = *elem_it;

                // Now check whether elem has a face on the specified boundary
                for (unsigned char side_id=0; side_id<el->n_sides(); ++side_id)
                  if (el->neighbor_ptr(side_id) == libmesh_nullptr)
                    {
                      // Get *all* boundary IDs on this side, not just the first one!
                      mesh_array[i]->get_boundary_info().boundary_ids (el, side_id, bc_ids);

                      if (std::find(bc_ids.begin(), bc_ids.end(), id_array[i]) != bc_ids.end())
                        {
                          UniquePtr<Elem> side (el->build_side_ptr(side_id));
                          for (unsigned int node_i=0; node_i<side->n_nodes(); ++node_i)
                            set_array[i]->insert( side->node_id(node_i) );

                          h_min = std::min(h_min, side->hmin());
                          h_min_updated = true;

                          // This side is on the boundary, add its information to side_to_elem
                          if(skip_find_neighbors && (i==0))
                            {
                              key_type key = el->key(side_id);
                              val_type val;
                              val.first = el;
                              val.second = side_id;

                              key_val_pair kvp;
                              kvp.first = key;
                              kvp.second = val;
                              // side_to_elem_map[key] = val;
#if defined(LIBMESH_HAVE_UNORDERED_MAP) || defined(LIBMESH_HAVE_TR1_UNORDERED_MAP) || defined(LIBMESH_HAVE_HASH_MAP) || defined(LIBMESH_HAVE_EXT_HASH_MAP)
                              side_to_elem_map.insert (kvp);
#else
                              side_to_elem_map.insert (side_to_elem_map.begin(),kvp);
#endif
                            }
                        }

                      // Also, check the edges on this side. We don't have to worry about
                      // updating neighbor info in this case since elements don't store
                      // neighbor info on edges.
                      for (unsigned short edge_id=0; edge_id<el->n_edges(); ++edge_id)
                        {
                          if(el->is_edge_on_side(edge_id, side_id))
                            {
                              // Get *all* boundary IDs on this edge, not just the first one!
                              mesh_array[i]->get_boundary_info().edge_boundary_ids (el, edge_id, bc_ids);

                              if (std::find(bc_ids.begin(), bc_ids.end(), id_array[i]) != bc_ids.end())
                                {
                                  UniquePtr<Elem> edge (el->build_edge_ptr(edge_id));
                                  for (unsigned int node_i=0; node_i<edge->n_nodes(); ++node_i)
                                    set_array[i]->insert( edge->node_id(node_i) );

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

          if(h_min_updated)
            {
              libMesh::out << "Minimum edge length on both surfaces is " << h_min << ".\n";
            }
          else
            {
              libMesh::out << "No elements on specified surfaces." << std::endl;
            }
        }


      if(use_binary_search)
        {
          // Store points from both stitched faces in sorted vectors for faster
          // searching later.
          typedef std::vector< std::pair<Point, dof_id_type> > PointVector;
          PointVector
            this_sorted_bndry_nodes(this_boundary_node_ids.size()),
            other_sorted_bndry_nodes(other_boundary_node_ids.size());

          // Comparison object that will be used later. So far, I've had reasonable success
          // with TOLERANCE...
          FuzzyPointCompare mein_comp(TOLERANCE);

          // Create and sort the vectors we will use to do the geometric searching
          {
            std::set<dof_id_type> * set_array[2] = {&this_boundary_node_ids, &other_boundary_node_ids};
            ReplicatedMesh * mesh_array[2]           = {this, other_mesh};
            PointVector * vec_array[2]           = {&this_sorted_bndry_nodes, &other_sorted_bndry_nodes};

            for (unsigned i=0; i<2; ++i)
              {
                std::set<dof_id_type>::iterator
                  set_it     = set_array[i]->begin(),
                  set_it_end = set_array[i]->end();

                // Fill up the vector with the contents of the set...
                for (unsigned ctr=0; set_it != set_it_end; ++set_it, ++ctr)
                  {
                    (*vec_array[i])[ctr] = std::make_pair(mesh_array[i]->point(*set_it), // The geometric point
                                                          *set_it);                      // Its ID
                  }

                // Sort the vectors based on the FuzzyPointCompare struct op()
                std::sort(vec_array[i]->begin(), vec_array[i]->end(), mein_comp);
              }
          }

          // Build up the node_to_node_map and node_to_elems_map using the sorted vectors of Points.
          for (std::size_t i=0; i<this_sorted_bndry_nodes.size(); ++i)
            {
              // Current point we're working on
              Point this_point = this_sorted_bndry_nodes[i].first;

              // FuzzyPointCompare does a fuzzy equality comparison internally to handle
              // slight differences between the list of nodes on each mesh.
              PointVector::iterator other_iter = Utility::binary_find(other_sorted_bndry_nodes.begin(),
                                                                      other_sorted_bndry_nodes.end(),
                                                                      this_point,
                                                                      mein_comp);

              // Not every node on this_sorted_bndry_nodes will necessarily be stitched, so
              // if its pair is not found on other_mesh, just continue.
              if (other_iter != other_sorted_bndry_nodes.end())
                {
                  // Check that the points do indeed match - should not be necessary unless something
                  // is wrong with binary_find.  To be on the safe side, we'll check.
                  {
                    // Grab the other point from the iterator
                    Point other_point = other_iter->first;

                    if (!this_point.absolute_fuzzy_equals(other_point, tol*h_min))
                      libmesh_error_msg("Error: mismatched points: " << this_point << " and " << other_point);
                  }


                  // Associate these two nodes in both the node_to_node_map and the other_to_this_node_map
                  dof_id_type
                    this_node_id = this_sorted_bndry_nodes[i].second,
                    other_node_id = other_iter->second;
                  node_to_node_map[this_node_id] = other_node_id;
                  other_to_this_node_map[other_node_id] = this_node_id;
                }

            }
        }
      else
        {
          // Otherwise, use a simple N^2 search to find the closest matching points. This can be helpful
          // in the case that we have tolerance issues which cause mismatch between the two surfaces
          // that are being stitched.

          std::set<dof_id_type>::iterator set_it     = this_boundary_node_ids.begin();
          std::set<dof_id_type>::iterator set_it_end = this_boundary_node_ids.end();
          for( ; set_it != set_it_end; ++set_it)
            {
              dof_id_type this_node_id = *set_it;
              Node & this_node = this->node_ref(this_node_id);

              bool found_matching_nodes = false;

              std::set<dof_id_type>::iterator other_set_it     = other_boundary_node_ids.begin();
              std::set<dof_id_type>::iterator other_set_it_end = other_boundary_node_ids.end();
              for( ; other_set_it != other_set_it_end; ++other_set_it)
                {
                  dof_id_type other_node_id = *other_set_it;
                  Node & other_node = other_mesh->node_ref(other_node_id);

                  Real node_distance = (this_node - other_node).norm();

                  if(node_distance < tol*h_min)
                    {
                      // Make sure we didn't already find a matching node!
                      if(found_matching_nodes)
                        libmesh_error_msg("Error: Found multiple matching nodes in stitch_meshes");

                      node_to_node_map[this_node_id] = other_node_id;
                      other_to_this_node_map[other_node_id] = this_node_id;

                      found_matching_nodes = true;
                    }
                }
            }
        }

      // Build up the node_to_elems_map, using only one loop over other_mesh
      {
        MeshBase::element_iterator other_elem_it  = other_mesh->elements_begin();
        MeshBase::element_iterator other_elem_end = other_mesh->elements_end();
        for (; other_elem_it != other_elem_end; ++other_elem_it)
          {
            Elem * el = *other_elem_it;

            // For each node on the element, find the corresponding node
            // on "this" Mesh, 'this_node_id', if it exists, and push
            // the current element ID back onto node_to_elems_map[this_node_id].
            // For that we will use the reverse mapping we created at
            // the same time as the forward mapping.
            for (unsigned n=0; n<el->n_nodes(); ++n)
              {
                dof_id_type other_node_id = el->node_id(n);
                std::map<dof_id_type, dof_id_type>::iterator it =
                  other_to_this_node_map.find(other_node_id);

                if (it != other_to_this_node_map.end())
                  {
                    dof_id_type this_node_id = it->second;
                    node_to_elems_map[this_node_id].push_back( el->id() );
                  }
              }
          }
      }

      if(verbose)
        {
          libMesh::out << "In ReplicatedMesh::stitch_meshes:\n"
                       << "Found " << node_to_node_map.size()
                       << " matching nodes.\n"
                       << std::endl;
        }

      if(enforce_all_nodes_match_on_boundaries)
        {
          std::size_t n_matching_nodes = node_to_node_map.size();
          std::size_t this_mesh_n_nodes = this_boundary_node_ids.size();
          std::size_t other_mesh_n_nodes = other_boundary_node_ids.size();
          if( (n_matching_nodes != this_mesh_n_nodes) ||
              (n_matching_nodes != other_mesh_n_nodes) )
            libmesh_error_msg("Error: We expected the number of nodes to match.");
        }
    }
  else
    {
      if(verbose)
        {
          libMesh::out << "Skip node merging in ReplicatedMesh::stitch_meshes:" << std::endl;
        }
    }



  dof_id_type node_delta = this->n_nodes();
  dof_id_type elem_delta = this->n_elem();

  // If other_mesh!=NULL, then we have to do a bunch of work
  // in order to copy it to this mesh
  if(this!=other_mesh)
    {
      // need to increment node and element IDs of other_mesh before copying to this mesh
      MeshBase::node_iterator node_it  = other_mesh->nodes_begin();
      MeshBase::node_iterator node_end = other_mesh->nodes_end();
      for (; node_it != node_end; ++node_it)
        {
          Node * nd = *node_it;
          dof_id_type new_id = nd->id() + node_delta;
          nd->set_id(new_id);
        }

      MeshBase::element_iterator elem_it  = other_mesh->elements_begin();
      MeshBase::element_iterator elem_end = other_mesh->elements_end();
      for (; elem_it != elem_end; ++elem_it)
        {
          Elem * el = *elem_it;
          dof_id_type new_id = el->id() + elem_delta;
          el->set_id(new_id);
        }

      // Also, increment the node_to_node_map and node_to_elems_map
      std::map<dof_id_type, dof_id_type>::iterator node_map_it     = node_to_node_map.begin();
      std::map<dof_id_type, dof_id_type>::iterator node_map_it_end = node_to_node_map.end();
      for( ; node_map_it != node_map_it_end; ++node_map_it)
        {
          node_map_it->second += node_delta;
        }
      std::map<dof_id_type, std::vector<dof_id_type> >::iterator elem_map_it     = node_to_elems_map.begin();
      std::map<dof_id_type, std::vector<dof_id_type> >::iterator elem_map_it_end = node_to_elems_map.end();
      for( ; elem_map_it != elem_map_it_end; ++elem_map_it)
        {
          std::size_t n_elems = elem_map_it->second.size();
          for(std::size_t i=0; i<n_elems; i++)
            {
              (elem_map_it->second)[i] += elem_delta;
            }
        }

      // Copy mesh data. If we skip the call to find_neighbors(), the lists
      // of neighbors will be copied verbatim from the other mesh
      this->copy_nodes_and_elements(*other_mesh, skip_find_neighbors);

      // Decrement node IDs of mesh to return to original state
      node_it  = other_mesh->nodes_begin();
      node_end = other_mesh->nodes_end();
      for (; node_it != node_end; ++node_it)
        {
          Node * nd = *node_it;
          dof_id_type new_id = nd->id() - node_delta;
          nd->set_id(new_id);
        }

      // Container to catch boundary IDs passed back from BoundaryInfo.
      std::vector<boundary_id_type> bc_ids;

      elem_it  = other_mesh->elements_begin();
      elem_end = other_mesh->elements_end();
      for (; elem_it != elem_end; ++elem_it)
        {
          Elem * other_elem = *elem_it;

          // Find the corresponding element on this mesh
          Elem * this_elem = this->elem_ptr(other_elem->id());

          // Decrement elem IDs of other_mesh to return it to original state
          dof_id_type new_id = other_elem->id() - elem_delta;
          other_elem->set_id(new_id);

          unsigned int other_n_nodes = other_elem->n_nodes();
          for (unsigned int n=0; n != other_n_nodes; ++n)
            {
              other_mesh->get_boundary_info().boundary_ids(other_elem->node_ptr(n), bc_ids);
              this->get_boundary_info().add_node(this_elem->node_ptr(n), bc_ids);
            }

          // Copy edge boundary info
          unsigned int n_edges = other_elem->n_edges();
          for (unsigned short edge=0; edge != n_edges; ++edge)
            {
              other_mesh->get_boundary_info().edge_boundary_ids(other_elem, edge, bc_ids);
              this->get_boundary_info().add_edge(this_elem, edge, bc_ids);
            }

          unsigned int n_sides = other_elem->n_sides();
          for (unsigned short s=0; s != n_sides; ++s)
            {
              other_mesh->get_boundary_info().boundary_ids(other_elem, s, bc_ids);
              this->get_boundary_info().add_side(this_elem, s, bc_ids);
            }

          // Copy shellface boundary ids
          unsigned int n_shellfaces = 2;
          for (unsigned short shellface=0; shellface != n_shellfaces; ++shellface)
            {
              other_mesh->get_boundary_info().shellface_boundary_ids(other_elem, shellface, bc_ids);
              this->get_boundary_info().add_shellface(this_elem, shellface, bc_ids);
            }
        }

    } // end if(other_mesh)

  // Finally, we need to "merge" the overlapping nodes
  // We do this by iterating over node_to_elems_map and updating
  // the elements so that they "point" to the nodes that came
  // from this mesh, rather than from other_mesh.
  // Then we iterate over node_to_node_map and delete the
  // duplicate nodes that came from other_mesh.

  // Container to catch boundary IDs passed back from BoundaryInfo.
  std::vector<boundary_id_type> bc_ids;

  std::map<dof_id_type, std::vector<dof_id_type> >::iterator elem_map_it     = node_to_elems_map.begin();
  std::map<dof_id_type, std::vector<dof_id_type> >::iterator elem_map_it_end = node_to_elems_map.end();
  for( ; elem_map_it != elem_map_it_end; ++elem_map_it)
    {
      dof_id_type target_node_id = elem_map_it->first;
      dof_id_type other_node_id = node_to_node_map[target_node_id];
      Node & target_node = this->node_ref(target_node_id);

      std::size_t n_elems = elem_map_it->second.size();
      for(std::size_t i=0; i<n_elems; i++)
        {
          dof_id_type elem_id = elem_map_it->second[i];
          Elem * el = this->elem_ptr(elem_id);

          // find the local node index that we want to update
          unsigned int local_node_index = el->local_node(other_node_id);

          // We also need to copy over the nodeset info here,
          // because the node will get deleted below
          this->get_boundary_info().boundary_ids(el->node_ptr(local_node_index), bc_ids);
          el->set_node(local_node_index) = &target_node;
          this->get_boundary_info().add_node(&target_node, bc_ids);
        }
    }

  std::map<dof_id_type, dof_id_type>::iterator node_map_it     = node_to_node_map.begin();
  std::map<dof_id_type, dof_id_type>::iterator node_map_it_end = node_to_node_map.end();
  for( ; node_map_it != node_map_it_end; ++node_map_it)
    {
      // In the case that this==other_mesh, the two nodes might be the same (e.g. if
      // we're stitching a "sliver"), hence we need to skip node deletion in that case.
      if ((this == other_mesh) && (node_map_it->second == node_map_it->first))
        continue;

      dof_id_type this_node_id = node_map_it->second;
      this->delete_node( this->node_ptr(this_node_id) );
    }

  // If find_neighbors() wasn't called in prepare_for_use(), we need to
  // manually loop once more over all elements adjacent to the stitched boundary
  // and fix their lists of neighbors.
  // This is done according to the following steps:
  //   1. Loop over all copied elements adjacent to the boundary using node_to_elems_map (trying to avoid duplicates)
  //   2. Look at all their sides with a NULL neighbor and update them using side_to_elem_map if necessary
  //   3. Update the corresponding side in side_to_elem_map as well
  if(skip_find_neighbors)
    {
      elem_map_it     = node_to_elems_map.begin();
      elem_map_it_end = node_to_elems_map.end();
      std::set<dof_id_type> fixed_elems;
      for( ; elem_map_it != elem_map_it_end; ++elem_map_it)
        {
          std::size_t n_elems = elem_map_it->second.size();
          for(std::size_t i=0; i<n_elems; i++)
            {
              dof_id_type elem_id = elem_map_it->second[i];
              if(fixed_elems.find(elem_id) == fixed_elems.end())
                {
                  Elem * el = this->elem_ptr(elem_id);
                  fixed_elems.insert(elem_id);
                  for(unsigned int s = 0; s < el->n_neighbors(); ++s)
                    {
                      if (el->neighbor_ptr(s) == libmesh_nullptr)
                        {
                          key_type key = el->key(s);
                          typedef map_type::iterator key_val_it_type;
                          std::pair<key_val_it_type, key_val_it_type>
                            bounds = side_to_elem_map.equal_range(key);

                          if(bounds.first != bounds.second)
                            {
                              // Get the side for this element
                              const UniquePtr<Elem> my_side(el->side_ptr(s));

                              // Look at all the entries with an equivalent key
                              while (bounds.first != bounds.second)
                                {
                                  // Get the potential element
                                  Elem * neighbor = bounds.first->second.first;

                                  // Get the side for the neighboring element
                                  const unsigned int ns = bounds.first->second.second;
                                  const UniquePtr<Elem> their_side(neighbor->side_ptr(ns));
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
                                  if( (*my_side == *their_side) &&
                                      (el->level() == neighbor->level()) &&
                                      ((el->dim() != 1) || (ns != s)) )
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

  this->prepare_for_use( /*skip_renumber_nodes_and_elements= */ false, skip_find_neighbors);

  // After the stitching, we may want to clear boundary IDs from element
  // faces that are now internal to the mesh
  if(clear_stitched_boundary_ids)
    {
      // Container to catch boundary IDs passed back from BoundaryInfo.
      std::vector<boundary_id_type> bc_ids;

      MeshBase::element_iterator elem_it  = this->elements_begin();
      MeshBase::element_iterator elem_end = this->elements_end();
      for (; elem_it != elem_end; ++elem_it)
        {
          Elem * el = *elem_it;

          for (unsigned short side_id=0; side_id<el->n_sides(); side_id++)
            {
              if (el->neighbor_ptr(side_id) != libmesh_nullptr)
                {
                  // Completely remove the side from the boundary_info object if it has either
                  // this_mesh_boundary_id or other_mesh_boundary_id.
                  this->get_boundary_info().boundary_ids (el, side_id, bc_ids);

                  if (std::find(bc_ids.begin(), bc_ids.end(), this_mesh_boundary_id) != bc_ids.end() ||
                      std::find(bc_ids.begin(), bc_ids.end(), other_mesh_boundary_id) != bc_ids.end())
                    this->get_boundary_info().remove_side(el, side_id);
                }
            }
        }
    }
}


dof_id_type ReplicatedMesh::n_active_elem () const
{
  return static_cast<dof_id_type>(std::distance (this->active_elements_begin(),
                                                 this->active_elements_end()));
}


} // namespace libMesh

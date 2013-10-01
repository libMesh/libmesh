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
#include "libmesh/metis_partitioner.h"
#include "libmesh/serial_mesh.h"
#include "libmesh/utility.h"

#include LIBMESH_INCLUDE_UNORDERED_SET
LIBMESH_DEFINE_HASH_POINTERS


namespace
{
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
    bool operator()(const Point& lhs, const Point& rhs)
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
    bool operator()(const std::pair<Point, dof_id_type>& lhs,
                    const std::pair<Point, dof_id_type>& rhs)
    {
      return (*this)(lhs.first, rhs.first);
    }

    // Comparsion function where lhs is a Point and rhs is a pair<Point,dof_id_type>.
    // This is used in routines like lower_bound, where a specific value is being
    // searched for.
    bool operator()(const Point& lhs, std::pair<Point, dof_id_type>& rhs)
    {
      return (*this)(lhs, rhs.first);
    }

    // And the other way around...
    bool operator()(std::pair<Point, dof_id_type>& lhs, const Point& rhs)
    {
      return (*this)(lhs.first, rhs);
    }
  };
}



namespace libMesh
{

// ------------------------------------------------------------
// SerialMesh class member functions
SerialMesh::SerialMesh (const Parallel::Communicator &comm,
			unsigned int d) :
  UnstructuredMesh (comm,d)
{
#ifdef LIBMESH_ENABLE_UNIQUE_ID
  // In serial we just need to reset the next unique id to zero
  // here in the constructor.
  _next_unique_id = 0;
#endif
  _partitioner = AutoPtr<Partitioner>(new MetisPartitioner());
}



#ifndef LIBMESH_DISABLE_COMMWORLD
SerialMesh::SerialMesh (unsigned int d) :
  UnstructuredMesh (d)
{
#ifdef LIBMESH_ENABLE_UNIQUE_ID
  // In serial we just need to reset the next unique id to zero
  // here in the constructor.
  _next_unique_id = 0;
#endif
  _partitioner = AutoPtr<Partitioner>(new MetisPartitioner());
}
#endif


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
  *this->boundary_info = *other_mesh.boundary_info;
}


SerialMesh::SerialMesh (const UnstructuredMesh &other_mesh) :
  UnstructuredMesh (other_mesh)
{
  this->copy_nodes_and_elements(other_mesh);
  *this->boundary_info = *other_mesh.boundary_info;
}


const Point& SerialMesh::point (const dof_id_type i) const
{
  libmesh_assert_less (i, this->n_nodes());
  libmesh_assert(_nodes[i]);
  libmesh_assert_equal_to (_nodes[i]->id(), i); // This will change soon

  return (*_nodes[i]);
}





const Node& SerialMesh::node (const dof_id_type i) const
{
  libmesh_assert_less (i, this->n_nodes());
  libmesh_assert(_nodes[i]);
  libmesh_assert_equal_to (_nodes[i]->id(), i); // This will change soon

  return (*_nodes[i]);
}





Node& SerialMesh::node (const dof_id_type i)
{
  if (i >= this->n_nodes())
    {
      libMesh::out << " i=" << i
		    << ", n_nodes()=" << this->n_nodes()
		    << std::endl;
      libmesh_error();
    }

  libmesh_assert_less (i, this->n_nodes());
  libmesh_assert(_nodes[i]);
  libmesh_assert_equal_to (_nodes[i]->id(), i); // This will change soon

  return (*_nodes[i]);
}



const Node* SerialMesh::node_ptr (const dof_id_type i) const
{
  libmesh_assert_less (i, this->n_nodes());
  libmesh_assert(_nodes[i]);
  libmesh_assert_equal_to (_nodes[i]->id(), i); // This will change soon

  return _nodes[i];
}




Node* SerialMesh::node_ptr (const dof_id_type i)
{
  libmesh_assert_less (i, this->n_nodes());
  libmesh_assert(_nodes[i]);
  libmesh_assert_equal_to (_nodes[i]->id(), i); // This will change soon

  return _nodes[i];
}




const Node* SerialMesh::query_node_ptr (const dof_id_type i) const
{
  if (i >= this->n_nodes())
    return NULL;
  libmesh_assert (_nodes[i] == NULL ||
                  _nodes[i]->id() == i); // This will change soon

  return _nodes[i];
}




Node* SerialMesh::query_node_ptr (const dof_id_type i)
{
  if (i >= this->n_nodes())
    return NULL;
  libmesh_assert (_nodes[i] == NULL ||
                  _nodes[i]->id() == i); // This will change soon

  return _nodes[i];
}




const Elem* SerialMesh::elem (const dof_id_type i) const
{
  libmesh_assert_less (i, this->n_elem());
  libmesh_assert(_elements[i]);
  libmesh_assert_equal_to (_elements[i]->id(), i); // This will change soon

  return _elements[i];
}




Elem* SerialMesh::elem (const dof_id_type i)
{
  libmesh_assert_less (i, this->n_elem());
  libmesh_assert(_elements[i]);
  libmesh_assert_equal_to (_elements[i]->id(), i); // This will change soon

  return _elements[i];
}




const Elem* SerialMesh::query_elem (const dof_id_type i) const
{
  if (i >= this->n_elem())
    return NULL;
  libmesh_assert (_elements[i] == NULL ||
                  _elements[i]->id() == i); // This will change soon

  return _elements[i];
}




Elem* SerialMesh::query_elem (const dof_id_type i)
{
  if (i >= this->n_elem())
    return NULL;
  libmesh_assert (_elements[i] == NULL ||
                  _elements[i]->id() == i); // This will change soon

  return _elements[i];
}




Elem* SerialMesh::add_elem (Elem* e)
{
  libmesh_assert(e);

  // We no longer merely append elements with SerialMesh

  // If the user requests a valid id that doesn't correspond to an
  // existing element, let's give them that id, resizing the elements
  // container if necessary.
  if (!e->valid_id())
    e->set_id (_elements.size());

  const dof_id_type id = e->id();

  if (id < _elements.size())
    {
      // Overwriting existing elements is still probably a mistake.
      libmesh_assert(!_elements[id]);
    }
  else
    {
      _elements.resize(id+1, NULL);
    }

  _elements[id] = e;

  return e;
}



Elem* SerialMesh::insert_elem (Elem* e)
{
  dof_id_type eid = e->id();
  libmesh_assert_less (eid, _elements.size());
  Elem *oldelem = _elements[eid];

  if (oldelem)
    {
      libmesh_assert_equal_to (oldelem->id(), eid);
      this->delete_elem(oldelem);
    }

  _elements[e->id()] = e;

  return e;
}



void SerialMesh::delete_elem(Elem* e)
{
  libmesh_assert(e);

  // Initialize an iterator to eventually point to the element we want to delete
  std::vector<Elem*>::iterator pos = _elements.end();

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
  this->boundary_info->remove(e);

  // delete the element
  delete e;

  // explicitly NULL the pointer
  *pos = NULL;
}



void SerialMesh::renumber_elem(const dof_id_type old_id,
                               const dof_id_type new_id)
{
  // This doesn't get used in serial yet
  Elem *el = _elements[old_id];
  libmesh_assert (el);

  el->set_id(new_id);
  libmesh_assert (!_elements[new_id]);
  _elements[new_id] = el;
  _elements[old_id] = NULL;
}



Node* SerialMesh::add_point (const Point& p,
			     const dof_id_type id,
			     const processor_id_type proc_id)
{
//   // We only append points with SerialMesh
//   libmesh_assert(id == DofObject::invalid_id || id == _nodes.size());
//   Node *n = Node::build(p, _nodes.size()).release();
//   n->processor_id() = proc_id;
//   _nodes.push_back (n);

  Node *n = NULL;

  // If the user requests a valid id, either
  // provide the existing node or resize the container
  // to fit the new node.
  if (id != DofObject::invalid_id)
    if (id < _nodes.size())
      n = _nodes[id];
    else
      _nodes.resize(id+1);
  else
    _nodes.push_back (static_cast<Node*>(NULL));

  // if the node already exists, then assign new (x,y,z) values
  if (n)
    *n = p;
  // otherwise build a new node, put it in the right spot, and return
  // a valid pointer.
  else
    {
      n = Node::build(p, (id == DofObject::invalid_id) ? _nodes.size()-1 : id).release();
      n->processor_id() = proc_id;

      if (id == DofObject::invalid_id)
	_nodes.back() = n;
      else
	_nodes[id] = n;
    }

  // better not pass back a NULL pointer.
  libmesh_assert (n);

  return n;
}



Node* SerialMesh::add_node (Node* n)
{
  libmesh_assert(n);
  // We only append points with SerialMesh
  libmesh_assert(!n->valid_id() || n->id() == _nodes.size());

  n->set_id (_nodes.size());

  _nodes.push_back(n);

  return n;
}



void SerialMesh::delete_node(Node* n)
{
  libmesh_assert(n);
  libmesh_assert_less (n->id(), _nodes.size());

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
  libmesh_assert (pos != _nodes.end());

  // Delete the node from the BoundaryInfo object
  this->boundary_info->remove(n);

  // delete the node
  delete n;

  // explicitly NULL the pointer
  *pos = NULL;
}



void SerialMesh::renumber_node(const dof_id_type old_id,
                               const dof_id_type new_id)
{
  // This doesn't get used in serial yet
  Node *nd = _nodes[old_id];
  libmesh_assert (nd);

  nd->set_id(new_id);
  libmesh_assert (!_nodes[new_id]);
  _nodes[new_id] = nd;
  _nodes[old_id] = NULL;
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
  dof_id_type next_free_elem = 0;
  dof_id_type next_free_node = 0;

  // Will hold the set of nodes that are currently connected to elements
  LIBMESH_BEST_UNORDERED_SET<Node*> connected_nodes;

  // Loop over the elements.  Note that there may
  // be NULLs in the _elements vector from the coarsening
  // process.  Pack the elements in to a contiguous array
  // and then trim any excess.
  {
    std::vector<Elem*>::iterator in        = _elements.begin();
    std::vector<Elem*>::iterator out_iter  = _elements.begin();
    const std::vector<Elem*>::iterator end = _elements.end();

    for (; in != end; ++in)
      if (*in != NULL)
	{
	  Elem* el = *in;

	  *out_iter = *in;
	  ++out_iter;

	  // Increment the element counter
	  el->set_id (next_free_elem++);

          if(_skip_renumber_nodes_and_elements)
          {
            // Add this elements nodes to the connected list
            for (unsigned int n=0; n<el->n_nodes(); n++)
              connected_nodes.insert(el->get_node(n));
          }
          else  // We DO want node renumbering
          {
            // Loop over this element's nodes.  Number them,
            // if they have not been numbered already.  Also,
            // position them in the _nodes vector so that they
            // are packed contiguously from the beginning.
            for (unsigned int n=0; n<el->n_nodes(); n++)
              if (el->node(n) == next_free_node)     // don't need to process
                next_free_node++;                      // [(src == dst) below]

              else if (el->node(n) > next_free_node) // need to process
	      {
		// The source and destination indices
		// for this node
		const dof_id_type src_idx = el->node(n);
		const dof_id_type dst_idx = next_free_node++;

		// ensure we want to swap a valid nodes
		libmesh_assert(_nodes[src_idx]);

		// Swap the source and destination nodes
                std::swap(_nodes[src_idx],
                          _nodes[dst_idx] );

		// Set proper indices where that makes sense
		if (_nodes[src_idx] != NULL)
		  _nodes[src_idx]->set_id (src_idx);
		_nodes[dst_idx]->set_id (dst_idx);
	      }
          }
	}

    // Erase any additional storage. These elements have been
    // copied into NULL voids by the procedure above, and are
    // thus repeated and unnecessary.
    _elements.erase (out_iter, end);
  }


  if(_skip_renumber_nodes_and_elements)
  {
    // Loop over the nodes.  Note that there may
    // be NULLs in the _nodes vector from the coarsening
    // process.  Pack the nodes in to a contiguous array
    // and then trim any excess.

    std::vector<Node*>::iterator in        = _nodes.begin();
    std::vector<Node*>::iterator out_iter  = _nodes.begin();
    const std::vector<Node*>::iterator end = _nodes.end();

    for (; in != end; ++in)
      if (*in != NULL)
      {
        // This is a reference so that if we change the pointer it will change in the vector
        Node* & nd = *in;

        // If this node is still connected to an elem, put it in the list
        if(connected_nodes.find(nd) != connected_nodes.end())
        {
          *out_iter = nd;
          ++out_iter;

          // Increment the node counter
          nd->set_id (next_free_node++);
        }
        else // This node is orphaned, delete it!
        {
          this->boundary_info->remove (nd);

          // delete the node
          delete nd;
          nd = NULL;
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

    // (This code block will erase the unused nodes)
    // Now, delete the unused nodes
    {
      std::vector<Node*>::iterator nd        = _nodes.begin();
      const std::vector<Node*>::iterator end = _nodes.end();

      std::advance (nd, next_free_node);

      for (std::vector<Node*>::iterator it=nd;
           it != end; ++it)
      {
        // Mesh modification code might have already deleted some
        // nodes
	if (*it == NULL)
          continue;

	// remove any boundary information associated with
	// this node
	this->boundary_info->remove (*it);

	// delete the node
	delete *it;
	*it = NULL;
      }

      _nodes.erase (nd, end);
    }
  }

  libmesh_assert_equal_to (next_free_elem, _elements.size());
  libmesh_assert_equal_to (next_free_node, _nodes.size());

  STOP_LOG("renumber_nodes_and_elem()", "Mesh");
}



void SerialMesh::fix_broken_node_and_element_numbering ()
{
   // Nodes first
  for (dof_id_type n=0; n<this->_nodes.size(); n++)
    if (this->_nodes[n] != NULL)
      this->_nodes[n]->set_id() = n;

  // Elements next
  for (dof_id_type e=0; e<this->_elements.size(); e++)
    if (this->_elements[e] != NULL)
      this->_elements[e]->set_id() = e;
}


void SerialMesh::stitch_meshes (SerialMesh& other_mesh,
                                boundary_id_type this_mesh_boundary_id,
                                boundary_id_type other_mesh_boundary_id,
                                Real tol,
                                bool clear_stitched_boundary_ids,
                                bool verbose,
                                bool use_binary_search,
                                bool enforce_all_nodes_match_on_boundaries)
{
  stitching_helper(&other_mesh,
                   this_mesh_boundary_id,
                   other_mesh_boundary_id,
                   tol,
                   clear_stitched_boundary_ids,
                   verbose,
                   use_binary_search,
                   enforce_all_nodes_match_on_boundaries);
}

void SerialMesh::stitch_surfaces (boundary_id_type boundary_id_1,
                                  boundary_id_type boundary_id_2,
                                  Real tol,
                                  bool clear_stitched_boundary_ids,
                                  bool verbose,
                                  bool use_binary_search,
                                  bool enforce_all_nodes_match_on_boundaries)
{
  stitching_helper(NULL,
                   boundary_id_1,
                   boundary_id_2,
                   tol,
                   clear_stitched_boundary_ids,
                   verbose,
                   use_binary_search,
                   enforce_all_nodes_match_on_boundaries);
}

void SerialMesh::stitching_helper (SerialMesh* other_mesh,
                                   boundary_id_type this_mesh_boundary_id,
                                   boundary_id_type other_mesh_boundary_id,
                                   Real tol,
                                   bool clear_stitched_boundary_ids,
                                   bool verbose,
                                   bool use_binary_search,
                                   bool enforce_all_nodes_match_on_boundaries)
{
  std::map<dof_id_type, dof_id_type> node_to_node_map, other_to_this_node_map; // The second is the inverse map of the first
  std::map<dof_id_type, std::vector<dof_id_type> > node_to_elems_map;
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

    // Loop below fills in these sets for the two meshes.
    std::set<dof_id_type> this_boundary_node_ids, other_boundary_node_ids;
    {
      // Make temporary fixed-size arrays for loop
      boundary_id_type id_array[2]        = {this_mesh_boundary_id, other_mesh_boundary_id};
      std::set<dof_id_type>* set_array[2] = {&this_boundary_node_ids, &other_boundary_node_ids};
      SerialMesh* mesh_array[2]           = {this, other_mesh};

      for (unsigned i=0; i<2; ++i)
        {
          MeshBase::element_iterator elem_it  = mesh_array[i]->elements_begin();
          MeshBase::element_iterator elem_end = mesh_array[i]->elements_end();
          for ( ; elem_it != elem_end; ++elem_it)
            {
              Elem *el = *elem_it;

              // Now check whether elem has a face on the specified boundary
              for (unsigned int side_id=0; side_id<el->n_sides(); ++side_id)
                if (el->neighbor(side_id) == NULL)
                  {
                    // Get *all* boundary IDs, not just the first one!
                    std::vector<boundary_id_type> bc_ids = mesh_array[i]->boundary_info->boundary_ids (el, side_id);

                    if (std::count(bc_ids.begin(), bc_ids.end(), id_array[i]))
                      {
                        AutoPtr<Elem> side (el->build_side(side_id));
                        for (unsigned int node_id=0; node_id<side->n_nodes(); ++node_id)
                          set_array[i]->insert( side->node(node_id) );

                        h_min = std::min(h_min, side->hmin());
                      }
                  }
            }
        }
    }

    if (verbose)
      {
        libMesh::out << "In SerialMesh::stitch_meshes:\n"
                     << "This mesh has "  << this_boundary_node_ids.size()
                     << " nodes on boundary " << this_mesh_boundary_id  << ".\n"
                     << "Other mesh has " << other_boundary_node_ids.size()
                     << " nodes on boundary " << other_mesh_boundary_id << ".\n"
                     << "Minimum edge length on both surfaces is " << h_min << ".\n"
                     << std::endl;
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
        std::set<dof_id_type>* set_array[2] = {&this_boundary_node_ids, &other_boundary_node_ids};
        SerialMesh* mesh_array[2]           = {this, other_mesh};
        PointVector* vec_array[2]           = {&this_sorted_bndry_nodes, &other_sorted_bndry_nodes};

        for (unsigned i=0; i<2; ++i)
        {
          std::set<dof_id_type>::iterator
            set_it     = set_array[i]->begin(),
            set_it_end = set_array[i]->end();

          // Fill up the vector with the contents of the set...
          for (unsigned ctr=0; set_it != set_it_end; ++set_it, ++ctr)
            {
              (*vec_array[i])[ctr] = std::make_pair( mesh_array[i]->point(*set_it), // The geometric point
                                                     *set_it );                     // Its ID
            }

          // Sort the vectors based on the FuzzyPointCompare struct op()
          std::sort(vec_array[i]->begin(), vec_array[i]->end(), mein_comp);
        }
      }

      // Build up the node_to_node_map and node_to_elems_map using the sorted vectors of Points.
      for (unsigned i=0; i<this_sorted_bndry_nodes.size(); ++i)
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
		{
		  libMesh::out << "Error: mismatched points: " << this_point << " and " << other_point << std::endl;
		  libmesh_error();
		}
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
        Node& this_node = this->node(this_node_id);

        bool found_matching_nodes = false;

        std::set<dof_id_type>::iterator other_set_it     = other_boundary_node_ids.begin();
        std::set<dof_id_type>::iterator other_set_it_end = other_boundary_node_ids.end();
        for( ; other_set_it != other_set_it_end; ++other_set_it)
        {
          dof_id_type other_node_id = *other_set_it;
          Node& other_node = other_mesh->node(other_node_id);

          Real node_distance = (this_node - other_node).size();

          if(node_distance < tol*h_min)
          {
            // Make sure we didn't already find a matching node!
            if(found_matching_nodes)
            {
              libMesh::out << "Error: Found multiple matching nodes in stitch_meshes" << std::endl;
              libmesh_error();
            }

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
        Elem *el = *other_elem_it;

        // For each node on the element, find the corresponding node
        // on "this" Mesh, 'this_node_id', if it exists, and push
        // the current element ID back onto node_to_elems_map[this_node_id].
        // For that we will use the reverse mapping we created at
        // the same time as the forward mapping.
        for (unsigned n=0; n<el->n_nodes(); ++n)
        {
          dof_id_type other_node_id = el->node(n);
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
      libMesh::out << "In SerialMesh::stitch_meshes:\n"
                   << "Found " << node_to_node_map.size()
                   << " matching nodes.\n"
                   << std::endl;
    }

    if(enforce_all_nodes_match_on_boundaries)
    {
      unsigned int n_matching_nodes = node_to_node_map.size();
      unsigned int this_mesh_n_nodes = this_boundary_node_ids.size();
      unsigned int other_mesh_n_nodes = other_boundary_node_ids.size();
      if( (n_matching_nodes != this_mesh_n_nodes) ||
          (n_matching_nodes != other_mesh_n_nodes) )
      {
        libMesh::out << "Error: We expected the number of nodes to match."
                     << std::endl;
        libmesh_error();
      }
    }
  }
  else
  {
    if(verbose)
    {
      libMesh::out << "Skip node merging in SerialMesh::stitch_meshes:" << std::endl;
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
      Node *nd = *node_it;
      dof_id_type new_id = nd->id() + node_delta;
      nd->set_id(new_id);
    }

    MeshBase::element_iterator elem_it  = other_mesh->elements_begin();
    MeshBase::element_iterator elem_end = other_mesh->elements_end();
    for (; elem_it != elem_end; ++elem_it)
    {
      Elem *el = *elem_it;
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
      dof_id_type n_elems = elem_map_it->second.size();
      for(dof_id_type i=0; i<n_elems; i++)
      {
        (elem_map_it->second)[i] += elem_delta;
      }
    }

    // Copy mesh data
    this->copy_nodes_and_elements(*other_mesh);

    // then decrement node and element IDs of mesh_i to return to original state
    node_it  = other_mesh->nodes_begin();
    node_end = other_mesh->nodes_end();
    for (; node_it != node_end; ++node_it)
    {
      Node *nd = *node_it;
      dof_id_type new_id = nd->id() - node_delta;
      nd->set_id(new_id);
    }

    elem_it  = other_mesh->elements_begin();
    elem_end = other_mesh->elements_end();
    for (; elem_it != elem_end; ++elem_it)
    {
      Elem *el = *elem_it;

      // First copy boundary info to the stitched mesh.  Note that other_mesh may
      // contain interior sidesets as well, so don't *just* copy boundary info for
      // elements on the boundary, copy it for all elements!
      for (unsigned int side_id=0; side_id<el->n_sides(); side_id++)
      {
        // There could be multiple boundary IDs on this side, so add them all.
        std::vector<boundary_id_type> bc_ids = other_mesh->boundary_info->boundary_ids (el, side_id);
        for (unsigned i=0; i<bc_ids.size(); ++i)
          {
            if (bc_ids[i] != BoundaryInfo::invalid_id)
              this->boundary_info->add_side(el->id(), side_id, bc_ids[i]);
          }
      }

      // Then decrement
      dof_id_type new_id = el->id() - elem_delta;
      el->set_id(new_id);
    }
  } // end if(other_mesh)

  // Finally, we need to "merge" the overlapping nodes
  // We do this by iterating over node_to_elems_map and updating
  // the elements so that they "point" to the nodes that came
  // from this mesh, rather than from other_mesh.
  // Then we iterate over node_to_node_map and delete the
  // duplicate nodes that came from other_mesh.
  std::map<dof_id_type, std::vector<dof_id_type> >::iterator elem_map_it     = node_to_elems_map.begin();
  std::map<dof_id_type, std::vector<dof_id_type> >::iterator elem_map_it_end = node_to_elems_map.end();
  for( ; elem_map_it != elem_map_it_end; ++elem_map_it)
  {
    dof_id_type target_node_id = elem_map_it->first;
    dof_id_type other_node_id = node_to_node_map[target_node_id];
    Node& target_node = this->node(target_node_id);

    dof_id_type n_elems = elem_map_it->second.size();
    for(unsigned int i=0; i<n_elems; i++)
    {
      dof_id_type elem_id = elem_map_it->second[i];
      Elem* el = this->elem(elem_id);

      // find the local node index that we want to update
      unsigned int local_node_index = el->local_node(other_node_id);

      el->set_node(local_node_index) = &target_node;
    }
  }

  std::map<dof_id_type, dof_id_type>::iterator node_map_it     = node_to_node_map.begin();
  std::map<dof_id_type, dof_id_type>::iterator node_map_it_end = node_to_node_map.end();
  for( ; node_map_it != node_map_it_end; ++node_map_it)
  {
    dof_id_type node_id = node_map_it->second;
    this->delete_node( this->node_ptr(node_id) );
  }

  this->prepare_for_use( /*skip_renumber_nodes_and_elements= */ false);

  // After the stitching, we may want to clear boundary IDs from element
  // faces that are now internal to the mesh
  if(clear_stitched_boundary_ids)
  {
    MeshBase::element_iterator elem_it  = this->elements_begin();
    MeshBase::element_iterator elem_end = this->elements_end();
    for (; elem_it != elem_end; ++elem_it)
    {
      Elem *el = *elem_it;

      for (unsigned int side_id=0; side_id<el->n_sides(); side_id++)
      {
        if (el->neighbor(side_id) != NULL)
        {
          // Completely remove the side from the boundary_info object if it has either
          // this_mesh_boundary_id or other_mesh_boundary_id.
          std::vector<boundary_id_type> bc_ids = this->boundary_info->boundary_ids (el, side_id);

          if (std::count(bc_ids.begin(), bc_ids.end(), this_mesh_boundary_id) ||
              std::count(bc_ids.begin(), bc_ids.end(), other_mesh_boundary_id))
            this->boundary_info->remove_side(el, side_id);
        }
      }
    }
  }

}


dof_id_type SerialMesh::n_active_elem () const
{
  return static_cast<dof_id_type>(std::distance (this->active_elements_begin(),
						 this->active_elements_end()));
}


#ifdef LIBMESH_ENABLE_UNIQUE_ID
void SerialMesh::assign_unique_ids()
{
  for (dof_id_type i=0; i<_elements.size(); ++i)
    if (_elements[i] && ! _elements[i]->valid_unique_id())
      _elements[i]->set_unique_id() = _next_unique_id++;

  for (dof_id_type i=0; i<_nodes.size(); ++i)
    if (_nodes[i] && ! _nodes[i]->valid_unique_id())
      _nodes[i]->set_unique_id() = _next_unique_id++;
}
#endif


} // namespace libMesh

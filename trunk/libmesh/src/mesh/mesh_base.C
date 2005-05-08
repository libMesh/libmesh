// $Id: mesh_base.C,v 1.92 2005-05-08 14:11:17 benkirk Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2005  Benjamin S. Kirk, John W. Peterson
  
// This library is free software; you can redistribute and/or
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



// library configuration
#include "libmesh_config.h"

// C++ includes
#include <algorithm> // for std::min
#include <sstream>   // for std::ostringstream
#include <map>       // for std::multimap

#if   defined(HAVE_HASH_MAP)
# include <hash_map>
#elif defined(HAVE_EXT_HASH_MAP)
# include <ext/hash_map>
#endif

// Local includes
#include "mesh_base.h"
#include "libmesh.h"
#include "face_tri3.h"
#include "face_tri6.h"
#include "libmesh_logging.h"
#include "gmv_io.h"
#include "metis_partitioner.h" // for default partitioning




// ------------------------------------------------------------
// MeshBase class member functions
MeshBase::MeshBase (unsigned int d) :
  boundary_info (*this),
  _n_sbd        (1),
  _n_parts      (1),
  _dim          (d),
  _is_prepared  (false)
{
  assert (DIM <= 3);
  assert (DIM >= _dim);
  assert (libMesh::initialized());
}



MeshBase::MeshBase (const MeshBase& other_mesh) :
  boundary_info      (*this),
  _nodes             (other_mesh._nodes),
  _elements          (other_mesh._elements),
  _n_sbd             (other_mesh._n_sbd),
  _n_parts           (other_mesh._n_parts),
  _dim               (other_mesh._dim),
  _is_prepared       (other_mesh._is_prepared)

{
}



MeshBase::~MeshBase()
{
  this->clear();

  assert (!libMesh::closed());
}



void MeshBase::prepare_for_use ()
{
  // Renumber the nodes and elements so that they in
  // contiguous blocks.
  this->renumber_nodes_and_elements();
  
  // Let all the elements find their neighbors
  this->find_neighbors();

  // Partition the mesh.
  this->partition();
  
  // The mesh is now prepared for use.
  _is_prepared = true;
}



Node* MeshBase::add_point (const Point& p)
{  
  _nodes.push_back (Node::build(p, this->n_nodes()).release());
  
  return _nodes.back();
}



Elem* MeshBase::add_elem (Elem* e)
{
  if (e != NULL)
    e->set_id (_elements.size());
  
  _elements.push_back(e);

  return e;
}



unsigned int MeshBase::n_active_elem () const
{
  unsigned int num=0;

  const_element_iterator       el  = this->active_elements_begin();
  const const_element_iterator end = this->active_elements_end(); 
  
  for (; el!=end; ++el)
    num++;

  return num;
}

 

void MeshBase::clear ()
{
  // Clear other data structures
  this->boundary_info.clear();

  
  // Reset the number of subdomains
  _n_sbd  = 1;

  // Reset the number of partitions
  _n_parts = 1;

  // Clear the elements data structure
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

  // Reset the _is_prepared flag
  _is_prepared = false;
}



unsigned int MeshBase::n_elem_on_proc (const unsigned int proc_id) const
{
  assert (proc_id < libMesh::n_processors());

  unsigned int ne=0;

  const_element_iterator       el  = this->pid_elements_begin(proc_id);
  const const_element_iterator end = this->pid_elements_end(proc_id);

  for (; el!=end; ++el)
    ne++;

  return ne;
}



unsigned int MeshBase::n_active_elem_on_proc (const unsigned int proc_id) const
{
  assert (proc_id < libMesh::n_processors());

  unsigned int ne=0;

  const_element_iterator       el  = this->active_pid_elements_begin(proc_id);
  const const_element_iterator end = this->active_pid_elements_end(proc_id);

  
  for (; el!=end; ++el)
    ne++;

  return ne;
}
  


unsigned int MeshBase::n_sub_elem () const
{
  unsigned int ne=0;

  const_element_iterator       el  = this->elements_begin();
  const const_element_iterator end = this->elements_end(); 

  for (; el!=end; ++el)
    ne += (*el)->n_sub_elem(); 

  return ne;
}



unsigned int MeshBase::n_active_sub_elem () const
{
  unsigned int ne=0;

  const_element_iterator       el  = this->active_elements_begin();
  const const_element_iterator end = this->active_elements_end(); 

  for (; el!=end; ++el)
    ne += (*el)->n_sub_elem(); 

  return ne;
}



std::vector<ElemType> MeshBase::elem_types() const
{
  std::vector<ElemType> et;

  assert (n_elem());

  const_element_iterator       el  = this->elements_begin();
  const const_element_iterator end = this->elements_end(); 

  // Automatically get the first type
  et.push_back((*el)->type());  ++el;

  // Loop over the rest of the elements.
  // If the current element type isn't in the
  // vector, insert it.
  for (; el != end; ++el)
    if (!std::count(et.begin(), et.end(), (*el)->type()))
      et.push_back((*el)->type());
  
  return et;
}



unsigned int MeshBase::n_elem_of_type(const ElemType type) const
{
  unsigned int cnt=0;

  const_element_iterator       el  = this->type_elements_begin(type);
  const const_element_iterator end = this->type_elements_end(type);

  for (; el!=end; ++el)
    cnt++;
  
  return cnt;
}



unsigned int MeshBase::n_active_elem_of_type(const ElemType type) const
{
  unsigned int cnt=0;

  const_element_iterator       el  = this->active_type_elements_begin(type);
  const const_element_iterator end = this->active_type_elements_end(type);
  
  for (; el!=end; ++el)
    cnt++;
    
  return cnt;
}



std::string MeshBase::get_info() const
{
  std::ostringstream out;

  out << " Mesh Information:"                                  << '\n'
      << "  mesh_dimension()="    << this->mesh_dimension()    << '\n'
      << "  spatial_dimension()=" << this->spatial_dimension() << '\n'
      << "  n_nodes()="           << this->n_nodes()           << '\n'
      << "  n_elem()="            << this->n_elem()            << '\n'
      << "   n_local_elem()="     << this->n_local_elem()      << '\n'
#ifdef ENABLE_AMR
      << "   n_active_elem()="    << this->n_active_elem()     << '\n'
#endif
      << "  n_subdomains()="      << this->n_subdomains()      << '\n'
      << "  n_processors()="      << this->n_processors()      << '\n'
      << "  processor_id()="      << this->processor_id()      << '\n';

  return out.str();
}


void MeshBase::print_info(std::ostream& os) const
{
  os << this->get_info()
     << std::endl;
}


std::ostream& operator << (std::ostream& os, const MeshBase& m)
{
  m.print_info(os);
  return os;
}



void MeshBase::find_neighbors()
{
  assert(this->n_nodes() != 0);
  assert(this->n_elem()  != 0);

  
  if (_dim == 1)
    error();


  START_LOG("find_neighbors()", "MeshBase");
  
  //TODO:[BSK] This should be removed later?!
  for (unsigned int e=0; e<this->n_elem(); e++)
    for (unsigned int s=0; s<this->elem(e)->n_neighbors(); s++)
      this->elem(e)->set_neighbor(s,NULL);

  // Find neighboring elements by first finding elements
  // with identical side keys and then check to see if they
  // are neighbors
  {
    // data structures -- Use the hash_multimap if available
    typedef unsigned int                    key_type;
    typedef std::pair<Elem*, unsigned char> val_type;
    typedef std::pair<key_type, val_type>   key_val_pair;
    
#if   defined(HAVE_HASH_MAP)    
    typedef std::hash_multimap<key_type, val_type> map_type;    
#elif defined(HAVE_EXT_HASH_MAP)
# if    (__GNUC__ == 3) && (__GNUC_MINOR__ == 0) // gcc 3.0   
    typedef std::hash_multimap<key_type, val_type> map_type;
# elif (__GNUC__ >= 3)                          // gcc 3.1 & newer
    typedef __gnu_cxx::hash_multimap<key_type, val_type> map_type;
# else
# error     DIE A HORRIBLE DEATH
# endif
#else
    typedef std::multimap<key_type, val_type>  map_type;
#endif
    
    // A map from side keys to corresponding elements & side numbers  
    map_type side_to_elem_map;
  
//     elem_iterator       el (this->elements_begin());
//     const elem_iterator end(this->elements_end());

    element_iterator       el  = this->elements_begin();
    const element_iterator end = this->elements_end();
    
    for (; el != end; ++el)
      {
	Elem* element = *el;
	
	for (unsigned int ms=0; ms<element->n_neighbors(); ms++)
	  {
	  next_side:
	    
	    if (element->neighbor(ms) == NULL)
	      {
		// Get the key for the side of this element
		const unsigned int key = element->key(ms);
		
		// Look for elements that have an identical side key
		std::pair <map_type::iterator, map_type::iterator>
		  bounds = side_to_elem_map.equal_range(key);
		
		// May be multiple keys, check all the possible
		// elements which _might_ be neighbors.
		if (bounds.first != bounds.second)
		  {
		    // Get the side for this element
		    const AutoPtr<Elem> my_side(element->side(ms));

		    // Look at all the entries with an equivalent key
		    while (bounds.first != bounds.second)
		      {
			// Get the potential element
			Elem* neighbor = bounds.first->second.first;
			
			// Get the side for the neighboring element
			const unsigned int ns = bounds.first->second.second;
			const AutoPtr<Elem> their_side(neighbor->side(ns));
			
			// If found a match wbounds.firsth my side
			if (*my_side == *their_side) 
			  {
			    // So share a side.  Is either of us
			    // the descendant of an active element?
                            // If not, then we're neighbors.
                            if (!element->subactive() &&
                                !neighbor->subactive())
                              {
			        element->set_neighbor (ms,neighbor);
			        neighbor->set_neighbor(ns,element);
                              }
			    
			    side_to_elem_map.erase (bounds.first);
			    
			    // get out of this nested crap
			    goto next_side; 
			  }

			++bounds.first;
		      }
		  }
		    
		// didn't find a match...
		// Build the map entry for this element
		key_val_pair kvp;
		
		kvp.first         = key;
		kvp.second.first  = element;
		kvp.second.second = ms;
		
		// use the lower bound as a hint for
		// where to put it.
#if defined(HAVE_HASH_MAP) || defined(HAVE_EXT_HASH_MAP)
		side_to_elem_map.insert (kvp);
#else
		side_to_elem_map.insert (bounds.first,kvp);
#endif
	      }
	  }
      }
  }

  
  
#ifdef ENABLE_AMR

  /**
   * Here we look at all of the child elements.
   * If a non-subactive child element has a NULL 
   * neighbor it is either because it is on the
   * boundary or because its neighbor is at a
   * different level.  In the latter case we must
   * get the neighbor from the parent.
   *
   * Furthermore, that neighbor better be active,
   * otherwise we missed a child somewhere.
   */
  element_iterator el  = this->not_level_elements_begin(0);
  element_iterator end = this->not_level_elements_end(0);

  for (; el != end; ++el)
    {
      Elem* elem = *el;
      
      assert (elem->parent() != NULL);

      if (elem->subactive()) 
        continue;
      
      for (unsigned int s=0; s < elem->n_neighbors(); s++)
	if (elem->neighbor(s) == NULL)
	  {	    
	    elem->set_neighbor(s, elem->parent()->neighbor(s));
	    
#ifdef DEBUG	    
            Elem *neigh = elem->neighbor(s);
	    if (neigh != NULL)
	      if (!neigh->active())
		{
		  std::cerr << "ERROR: " 
                    << (elem->active()?"Active":"Ancestor")
                    << "Element at level "
                    << elem->level() << " found "
		    << (neigh->subactive()?"subactive":"ancestor")
                    << " neighbor at level " << neigh->level()
                    << std::endl;
		  GMVIO(*dynamic_cast<Mesh*>(this)).write ("bad_mesh.gmv");
		  error();
		}
#endif // DEBUG
	  }
    }
  
#endif // AMR

  STOP_LOG("find_neighbors()", "MeshBase");
}



void MeshBase::partition (const unsigned int n_parts)
{
  MetisPartitioner partitioner;
  partitioner.partition (*this, n_parts); 
}



void MeshBase::renumber_nodes_and_elements ()
{
  // Only renumber the nodes & trim the element vector
  // if we have performed some local coarsening.  This
  // will be the case if the element vector is not the
  // same size as the number of active elements.
  //
  // In this case simply number the elements (in case this
  // was not already done) and return.
  if (_elements.size() == this->n_active_elem())
    { 
      START_LOG("renumber_nodes_and_elem()", "MeshBase");
      
      // Number the elements
      {
	element_iterator       it  = this->elements_begin();
	const element_iterator end = this->elements_end();
	
	unsigned int id=0;

	for (; it != end; ++it)
	  (*it)->set_id() = id++;
      }

      // Number the nodes
      {
	node_iterator       it  = this->nodes_begin();
	const node_iterator end = this->nodes_end();
	
	unsigned int id=0;

	for (; it != end; ++it)
	  (*it)->set_id() = id++;
      }      
      
      STOP_LOG("renumber_nodes_and_elem()", "MeshBase");
      return;
    }
  
  START_LOG("renumber_nodes_and_elem()", "MeshBase");
  
  // node and element id counters
  unsigned int next_free_elem = 0;
  unsigned int next_free_node = 0;

  // Loop over the elements.  Note that there may
  // be NULLs in the _elements vector from the coarsening
  // process.  Pack the elements in to a contiguous array
  // and then trim any excess.
  {      
    std::vector<Elem*>::iterator in        = _elements.begin();
    std::vector<Elem*>::iterator out       = _elements.begin();
    const std::vector<Elem*>::iterator end = _elements.end();

    for (; in != end; ++in)
      if (*in != NULL)
	{
	  Elem* elem = *in;
	  
	  *out = *in;
	  ++out;
	  
	  // Increment the element counter
	  elem->set_id (next_free_elem++);
	  
	  // Loop over this element's nodes.  Number them,
	  // if they have not been numbered already.  Also,
	  // position them in the _nodes vector so that they
	  // are packed contiguously from the beginning.
	  for (unsigned int n=0; n<elem->n_nodes(); n++)
	    if (elem->node(n) == next_free_node)     // don't need to process
	      next_free_node++;                      // [(src == dst) below]

	    else if (elem->node(n) > next_free_node) // need to process
	      {
		// The source and destination indices
		// for this node
		const unsigned int src_idx = elem->node(n);
		const unsigned int dst_idx = next_free_node++;

		// ensure we want to swap valid nodes
		assert (_nodes[src_idx] != NULL);
		assert (_nodes[dst_idx] != NULL);
		
		// Swap the source and destination nodes
		std::swap (_nodes[src_idx],
			   _nodes[dst_idx] );

		// Set proper indices
		_nodes[src_idx]->set_id (src_idx);
		_nodes[dst_idx]->set_id (dst_idx);
	      }
	}

    // Erase any additional storage. These elements have been
    // copied into NULL voids by the procedure above, and are
    // thus repeated and unnecessary.
    _elements.erase (out, end);
  }

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
	assert (*it != NULL);

	// remove any boundary information associated with
	// this node
	this->boundary_info.remove (*it);
	
	// delete the node
	delete *it;
	*it = NULL;
      }
    
    _nodes.erase (nd, end);
  }
  

//   // (This code block will keep the unused nodes)
//   // Now, number the unused nodes so that they are
//   // contiguous
//   {
//     std::vector<Node*>::iterator nd        = _nodes.begin();
//     const std::vector<Node*>::iterator end = _nodes.end();

//     std::advance (nd, next_free_node);
    
//     for (; nd != end; ++nd)
//       (*nd)->set_id(next_free_node++);
//   }
  

  assert (next_free_elem == _elements.size());
  assert (next_free_node == _nodes.size());
  
  STOP_LOG("renumber_nodes_and_elem()", "MeshBase");
}



void MeshBase::delete_elem(Elem* e)
{
  assert (e != NULL);

  std::vector<Elem*>::iterator pos = std::find (_elements.begin(),
						_elements.end(),
						e);

  // Huh? Element not in the vector?
  assert (pos != _elements.end());

  // delete the element
  delete e;
  
  // explicitly NULL the pointer
  e    = NULL;
  *pos = NULL;
  
  //_elements.erase(pos);
}



bool MeshBase::contract ()
{
  std::cerr << "ERROR:  This method must be re-implemented in derived classes!"
	    << std::endl;
  error();

  return false;
}






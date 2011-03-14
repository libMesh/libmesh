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



// C++ includes


// Local includes
#include "libmesh_config.h"
#include "boundary_info.h"
#include "boundary_mesh.h"
#include "elem.h"
#include "mesh_data.h"
#include "parallel.h"
#include "partitioner.h"

namespace libMesh
{



//------------------------------------------------------
// BoundaryInfo static member initializations
const short int BoundaryInfo::invalid_id = -1234;



//------------------------------------------------------
// BoundaryInfo functions
BoundaryInfo::BoundaryInfo(const MeshBase& m) :
  _mesh (m)
{
}

BoundaryInfo& BoundaryInfo::operator=(const BoundaryInfo& other_boundary_info)
{
  /**
   * A quick note: We're going to attempt to pull _new_ pointers out of the mesh assigned to this boundary info.
   * This will only work if the mesh assigned to this BoundaryInfo is the same mesh object as other_boundary_info
   * _or_ was constructed in exactly the same way (or constructed as a copy).
   */

  {
    std::multimap<const Node*, short int>::const_iterator it = other_boundary_info._boundary_node_id.begin();
    const std::multimap<const Node*, short int>::const_iterator end = other_boundary_info._boundary_node_id.end();

    for(; it != end; ++it)
    {
      const Node * other_node = it->first;
      _boundary_node_id.insert( std::pair<const Node*, short int>(_mesh.node_ptr(other_node->id()), it->second) );
    }
  }
  

  {
    std::multimap<const Elem*, std::pair<unsigned short int, short int> >::const_iterator it = other_boundary_info._boundary_side_id.begin();
    const std::multimap<const Elem*, std::pair<unsigned short int, short int> >::const_iterator end = other_boundary_info._boundary_side_id.end();

    for(; it != end; ++it)
    {
      const Elem * other_elem = it->first;
      _boundary_side_id.insert( std::pair<const Elem*, std::pair<unsigned short int, short int> >(_mesh.elem(other_elem->id()), it->second) );
    }
  }

  _boundary_ids = other_boundary_info._boundary_ids;

  return *this;
}


BoundaryInfo::~BoundaryInfo()
{
  this->clear();
}



void BoundaryInfo::clear()
{
  _boundary_node_id.clear();
  _boundary_side_id.clear();
  _boundary_ids.clear();
}



void BoundaryInfo::sync (BoundaryMesh& boundary_mesh,
			 MeshData*     boundary_mesh_data,
			 MeshData*     this_mesh_data)
{
  boundary_mesh.clear();

  /**
   * Re-create the boundary mesh.
   */

  // Map boundary ids to side subdomain/partition ids
  std::map<short int, unsigned int> id_map;

  // Original Code
  //     unsigned int cnt = 0;
  //     for (std::set<short int>::iterator pos = boundary_ids.begin();
  // 	 pos != boundary_ids.end(); ++pos)
  //       id_map[*pos] = cnt++;
  
  //     id_map[invalid_id] = cnt;

    
  // New code
  // Here we need to use iota() once it is in the
  // Utility namespace.
  std::for_each(_boundary_ids.begin(),
		_boundary_ids.end(),
		Fill(id_map));
    
  boundary_mesh.set_n_partitions() = id_map.size();


  // Make individual copies of all the nodes in the current mesh
  // and add them to the boundary mesh.  Yes, this is overkill because
  // all of the current mesh nodes will not end up in the the boundary
  // mesh.  These nodes can be trimmed later via a call to prepare_for_use().
  {
    libmesh_assert (boundary_mesh.n_nodes() == 0);
    boundary_mesh.reserve_nodes(_mesh.n_nodes());
    
    MeshBase::const_node_iterator it  = _mesh.nodes_begin();
    MeshBase::const_node_iterator end = _mesh.nodes_end();
    
    for(; it != end; ++it)
      {
	const Node* node = *it;
	boundary_mesh.add_point(*node); // calls Node::build(Point, id)
      }
  }

  // Add additional sides that aren't flagged with boundary conditions
  MeshBase::const_element_iterator       el     = _mesh.active_elements_begin();
  const MeshBase::const_element_iterator end_el = _mesh.active_elements_end(); 

  for ( ; el != end_el; ++el)
    {
      const Elem* elem = *el;
      
      for (unsigned int s=0; s<elem->n_sides(); s++)
	if (elem->neighbor(s) == NULL) // on the boundary
	  {

	    // Build the side - do not use a "proxy" element here:
	    // This will be going into the BoundaryMesh and needs to
	    // stand on its own.
	    AutoPtr<Elem> side (elem->build_side(s, false));
	    
	    // Get the top-level parent for this element
	    const Elem* top_parent = elem->top_parent();

	    // A convenient typedef
	    typedef
	      std::multimap<const Elem*, std::pair<unsigned short int, short int> >::const_iterator
	      Iter;
	      
	    // Find the right id number for that side
	    std::pair<Iter, Iter> pos = _boundary_side_id.equal_range(top_parent);

	    while (pos.first != pos.second)
	      {
		if (pos.first->second.first == s) // already flagged with a boundary condition
		  {
		    side->subdomain_id() =
		      id_map[pos.first->second.second];
		    break;
		  }
		
		++pos.first;
	      }

	    // either the element wasn't found or side s
	    // doesn't have a boundary condition
	    if (pos.first == pos.second)
	      {
		side->subdomain_id() = id_map[invalid_id];
	      }

	    side->processor_id() = side->subdomain_id(); //elem->processor_id();
	    
	    // Add the side
	    Elem* new_elem = boundary_mesh.add_elem(side.release());

	    // This side's Node pointers still point to the nodes of the  original mesh.
	    // We need to re-point them to the boundary mesh's nodes!  Since we copied *ALL* of
	    // the original mesh's nodes over, we should be guaranteed to have the same ordering.
	    for (unsigned int nn=0; nn<new_elem->n_nodes(); ++nn)
	      {
		// Get the correct node pointer, based on the id()
		Node* new_node = boundary_mesh.node_ptr(new_elem->node(nn));
		
		// sanity check: be sure that the new Nodes global id really matches
		libmesh_assert (new_node->id() == new_elem->node(nn));

		// Assign the new node pointer
		new_elem->set_node(nn) = new_node;
	      }
	  }
    } // end loop over active elements

  
  
  // When desired, copy the MeshData
  // to the boundary_mesh
  if ((boundary_mesh_data != NULL) && (this_mesh_data != NULL))
    boundary_mesh_data->assign(*this_mesh_data);

  // Don't repartition this mesh; we're using the processor_id values
  // as a hack to display bcids for now.
  boundary_mesh.partitioner().reset(NULL);

  // Trim any un-used nodes from the Mesh
  boundary_mesh.prepare_for_use(/*skip_renumber =*/ false);
}


			 

void BoundaryInfo::sync (const std::set<short int> &requested_boundary_ids,
			 BoundaryMesh& boundary_mesh)
{
  // Re-create the boundary mesh.
  boundary_mesh.clear();
    
  boundary_mesh.set_n_partitions() = _mesh.n_partitions();

  // Make individual copies of all the nodes in the current mesh
  // and add them to the boundary mesh.  Yes, this is overkill because
  // all of the current mesh nodes will not end up in the the boundary
  // mesh.  These nodes can be trimmed later via a call to prepare_for_use().
  {
    libmesh_assert (boundary_mesh.n_nodes() == 0);
    boundary_mesh.reserve_nodes(_mesh.n_nodes());
    
    MeshBase::const_node_iterator it  = _mesh.nodes_begin();
    MeshBase::const_node_iterator end = _mesh.nodes_end();
    
    for(; it != end; ++it)
      {
	const Node* node = *it;
	boundary_mesh.add_point(*node); // calls Node::build(Point, id)
      }
  }

  // Add additional sides that aren't flagged with boundary conditions
  MeshBase::const_element_iterator       el     = _mesh.active_elements_begin();
  const MeshBase::const_element_iterator end_el = _mesh.active_elements_end(); 

  for ( ; el != end_el; ++el)
    {
      const Elem* elem = *el;
      
      for (unsigned int s=0; s<elem->n_sides(); s++)
	if (elem->neighbor(s) == NULL) // on the boundary
	  {
	    // Get the top-level parent for this element
	    const Elem* top_parent = elem->top_parent();

	    // A convenient typedef
	    typedef
	      std::multimap<const Elem*, std::pair<unsigned short int, short int> >::const_iterator
	      Iter;
	      
	    // Find all the bcs asociated with top_parent
	    std::pair<Iter, Iter> pos = _boundary_side_id.equal_range(top_parent);

	    // look for a bcid which is (i) in the user-requested set, and
	    // (ii) matches the current side #s
	    while (pos.first != pos.second)
	      {
		// if this side is flagged with a boundary condition
		// and the user wants this id
		if ((pos.first->second.first == s) &&
		    (requested_boundary_ids.count(pos.first->second.second)))
		  {
		    // Build the side - do not use a "proxy" element here:
		    // This will be going into the BoundaryMesh and needs to
		    // stand on its own.
		    AutoPtr<Elem> side (elem->build_side(s, false));
		    
		    // inherit processor_id and subdomain_id from parent
		    side->subdomain_id() = elem->subdomain_id();
		    side->processor_id() = elem->processor_id();

		    // Add the side
		    Elem* new_elem = boundary_mesh.add_elem(side.release());

		    // and set the parent
		    new_elem->set_parent (const_cast<Elem*>(elem));

		    // This side's Node pointers still point to the nodes of the  original mesh.
		    // We need to re-point them to the boundary mesh's nodes!  Since we copied *ALL* of
		    // the original mesh's nodes over, we should be guaranteed to have the same ordering.
		    for (unsigned int nn=0; nn<new_elem->n_nodes(); ++nn)
		      {
			// Get the correct node pointer, based on the id()
			Node* new_node = boundary_mesh.node_ptr(new_elem->node(nn));
			
			// sanity check: be sure that the new Nodes global id really matches
			libmesh_assert (new_node->id() == new_elem->node(nn));

			// Assign the new node pointer
			new_elem->set_node(nn) = new_node;
		      }

		    // go on to the next side
		    break;
		  }
		
		++pos.first;
	      } // end loop over bcs matching top_parent
	    
	  } // end if neighbor is NULL
    } // end loop over active elements

  // Don't repartition this mesh; but rather inherit the partitioning
  boundary_mesh.partitioner().reset(NULL);
  
  // Trim any un-used nodes from the Mesh
  boundary_mesh.prepare_for_use(/*skip_renumber =*/ false);

  // and finally distribute element partitioning to the nodes
  Partitioner::set_node_processor_ids(boundary_mesh);
}




void BoundaryInfo::add_node(const unsigned int node,
			    const short int id)
{
  this->add_node (_mesh.node_ptr(node), id);
}



void BoundaryInfo::add_node(const Node* node,
			    const short int id)
{
  if (id == invalid_id)
    {
      libMesh::err << "ERROR: You may not set a boundary ID of "
		    << invalid_id << std::endl
		    << " That is reserved for internal use.\n"
		    << std::endl;

      libmesh_error();
    }

  // A convenient typedef
  typedef std::multimap<const Node*, short int>::const_iterator Iter;
	      
  // Don't add the same ID twice
  std::pair<Iter, Iter> pos = _boundary_node_id.equal_range(node);

  for (;pos.first != pos.second; ++pos.first)
    if (pos.first->second == id)
      return;

  std::pair<const Node*, short int> kv (node, id);
  
  _boundary_node_id.insert(kv);
  _boundary_ids.insert(id);
}



void BoundaryInfo::add_side(const unsigned int e,
			    const unsigned short int side,
			    const short int id)
{
  this->add_side (_mesh.elem(e), side, id);
}



void BoundaryInfo::add_side(const Elem* elem,
			    const unsigned short int side,
			    const short int id)
{
  libmesh_assert (elem != NULL);

  // Only add BCs for level-0 elements.
  libmesh_assert (elem->level() == 0);
  
  if (id == invalid_id)
    {
      libMesh::err << "ERROR: You may not set a boundary ID of "
		    << invalid_id << std::endl
		    << " That is reserved for internal use.\n"
		    << std::endl;

      libmesh_error();
    }
  
  // A convenient typedef
  typedef std::multimap<const Elem*, std::pair<unsigned short int, short int> >::const_iterator Iter;
	      
  // Don't add the same ID twice
  std::pair<Iter, Iter> pos = _boundary_side_id.equal_range(elem);

  for (;pos.first != pos.second; ++pos.first)
    if (pos.first->second.first == side &&
        pos.first->second.second == id)
      return;

  std::pair<unsigned short int, short int> p(side,id);
  std::pair<const Elem*, std::pair<unsigned short int, short int> >
    kv (elem, p);
  
  _boundary_side_id.insert(kv);
  _boundary_ids.insert(id);
}



void BoundaryInfo::add_side(const Elem* elem,
			    const unsigned short int side,
			    const std::vector<short int>& ids)
{
  if (ids.empty())
    return;

  libmesh_assert (elem != NULL);

  // Only add BCs for level-0 elements.
  libmesh_assert (elem->level() == 0);
  
  // A convenient typedef
  typedef std::multimap<const Elem*, std::pair<unsigned short int, short int> >::const_iterator Iter;
	      
  // Don't add the same ID twice
  std::pair<Iter, Iter> pos = _boundary_side_id.equal_range(elem);

  for (unsigned int i=0; i!= ids.size(); ++i)
    {
      short int id=ids[i];

      if (id == invalid_id)
        {
          libMesh::err << "ERROR: You may not set a boundary ID of "
		        << invalid_id << std::endl
		        << " That is reserved for internal use.\n"
		        << std::endl;

          libmesh_error();
        }
  
      for (Iter p = pos.first;p != pos.second; ++p)
        if (p->second.first == side &&
            p->second.second == id)
          continue;

      std::pair<unsigned short int, short int> p(side,id);
      std::pair<const Elem*, std::pair<unsigned short int, short int> >
        kv (elem, p);
  
      _boundary_side_id.insert(kv);
      _boundary_ids.insert(id);
    }
}



std::vector<short int> BoundaryInfo::boundary_ids(const Node* node) const
{
  std::vector<short int> ids;
  
  // A convenient typedef
  typedef std::multimap<const Node*, short int>::const_iterator Iter;
	      
  // Don't add the same ID twice
  std::pair<Iter, Iter> pos = _boundary_node_id.equal_range(node);

  for (;pos.first != pos.second; ++pos.first)
    ids.push_back(pos.first->second);

  return ids;
}



short int BoundaryInfo::boundary_id(const Elem* const elem,
				    const unsigned short int side) const
{
  libmesh_assert (elem != NULL);

  // Only level-0 elements store BCs.  If this is not a level-0
  // element get its level-0 parent and infer the BCs.
  const Elem*  searched_elem = elem;
  if (elem->level() != 0)
  {
    if (elem->side(side).get() == NULL)
      searched_elem = elem->top_parent ();
    else
      while (searched_elem->parent() != NULL) {
        const Elem * parent = searched_elem->parent();
        if (parent->is_child_on_side(parent->which_child_am_i(searched_elem), side) == false)
          return invalid_id;
        searched_elem = parent;
      }
  }
  
  std::pair<std::multimap<const Elem*,
                          std::pair<unsigned short int, short int> >::const_iterator,
            std::multimap<const Elem*,
                          std::pair<unsigned short int, short int> >::const_iterator > 
    e = _boundary_side_id.equal_range(searched_elem);

  // elem not in the data structure
  if (e.first == e.second)
    return invalid_id;

  // elem is there, maybe multiple occurances
  while (e.first != e.second)
    {
      // if this is true we found the requested side
      // of the element and want to return the id
      if (e.first->second.first == side)
	return e.first->second.second;

      ++e.first;
    }

  // if we get here, we found elem in the data structure but not
  // the requested side, so return the default value
  return invalid_id;  
}



std::vector<short int> BoundaryInfo::boundary_ids (const Elem* const elem,
                                                   const unsigned short int side) const
{
  libmesh_assert (elem != NULL);

  std::vector<short int> ids;

  // Only level-0 elements store BCs.  If this is not a level-0
  // element get its level-0 parent and infer the BCs.
  const Elem*  searched_elem = elem;
  if (elem->level() != 0)
  {
    if (elem->side(side).get() == NULL)
      searched_elem = elem->top_parent ();
    else
      while (searched_elem->parent() != NULL) {
        const Elem * parent = searched_elem->parent();
        if (parent->is_child_on_side(parent->which_child_am_i(searched_elem), side) == false)
          return ids;
        searched_elem = parent;
      }
  }

  std::pair<std::multimap<const Elem*,
                          std::pair<unsigned short int, short int> >::const_iterator,
            std::multimap<const Elem*,
                          std::pair<unsigned short int, short int> >::const_iterator >
    e = _boundary_side_id.equal_range(searched_elem);

  // elem not in the data structure
  if (e.first == e.second)
    return ids;

  // elem is there, maybe multiple occurances
  while (e.first != e.second)
    {
      // if this is true we found the requested side of the element
      if (e.first->second.first == side)
        ids.push_back(e.first->second.second);

      ++e.first;
    }

  // if we get here, we found elem in the data structure but not
  // the requested side, so return the default value
  return ids;

}



std::vector<short int> BoundaryInfo::raw_boundary_ids (const Elem* const elem,
                                                       const unsigned short int side) const
{
  libmesh_assert (elem != NULL);

  std::vector<short int> ids;

  // Only level-0 elements store BCs.
  if (elem->parent())
    return ids;

  std::pair<std::multimap<const Elem*,
                          std::pair<unsigned short int, short int> >::const_iterator,
            std::multimap<const Elem*,
                          std::pair<unsigned short int, short int> >::const_iterator >
    e = _boundary_side_id.equal_range(elem);

  // Check any occurances
  while (e.first != e.second)
    {
      // if this is true we found the requested side of the element
      if (e.first->second.first == side)
        ids.push_back(e.first->second.second);

      ++e.first;
    }

  // if nothing got pushed back, we didn't find elem in the data
  // structure with the requested side, so return the default empty
  // vector 
  return ids;
}



void BoundaryInfo::remove_side (const Elem* elem,
                                const unsigned short int side)
{
  libmesh_assert (elem != NULL);

  // The user shouldn't be trying to remove only one child's boundary
  // id
  libmesh_assert (elem->level() == 0);
  
  std::pair<std::multimap<const Elem*,
                          std::pair<unsigned short int, short int> >::iterator,
            std::multimap<const Elem*,
                          std::pair<unsigned short int, short int> >::iterator > 
    e = _boundary_side_id.equal_range(elem);

  // elem may be there, maybe multiple occurances
  while (e.first != e.second)
    {
      // if this is true we found the requested side
      // of the element and want to erase the id
      if (e.first->second.first == side)
	{
	  // (postfix++ - increment the iterator before it's invalid)
          _boundary_side_id.erase(e.first++);
        }
      else
        ++e.first;
    }
}



void BoundaryInfo::remove_side (const Elem* elem,
                                const unsigned short int side,
                                const short int id)
{
  libmesh_assert (elem != NULL);

  // The user shouldn't be trying to remove only one child's boundary
  // id
  libmesh_assert (elem->level() == 0);
  
  std::pair<std::multimap<const Elem*,
                          std::pair<unsigned short int, short int> >::iterator,
            std::multimap<const Elem*,
                          std::pair<unsigned short int, short int> >::iterator > 
    e = _boundary_side_id.equal_range(elem);

  // elem may be there, maybe multiple occurances
  while (e.first != e.second)
    {
      // if this is true we found the requested side
      // of the element and want to erase the requested id
      if (e.first->second.first == side &&
          e.first->second.second == id)
	{
	  // (postfix++ - increment the iterator before it's invalid)
          _boundary_side_id.erase(e.first++);
        }
      else
        ++e.first;
    }
}



unsigned int BoundaryInfo::side_with_boundary_id(const Elem* const elem,
                                                 const unsigned short int boundary_id) const
{
  const Elem* searched_elem = elem;
  if (elem->level() != 0)
    searched_elem = elem->top_parent();

  std::pair<std::multimap<const Elem*,
                          std::pair<unsigned short int, short int> >::const_iterator,
            std::multimap<const Elem*,
                          std::pair<unsigned short int, short int> >::const_iterator > 
    e = _boundary_side_id.equal_range(searched_elem);

  // elem may have zero or multiple occurances
  for (; e.first != e.second; ++e.first)
    {
      // if this is true we found the requested boundary_id
      // of the element and want to return the side
      if (e.first->second.second == boundary_id)
        {
         unsigned int side = e.first->second.first;

         // If we're on this external boundary then we share this
         // external boundary id
         if (elem->neighbor(side) == NULL)
           return side;

         // If we're on an internal boundary then we need to be sure
         // it's the same internal boundary as our top_parent
         const Elem *p = elem;
         while (p != NULL)
           {
             const Elem *parent = p->parent();
             if (!parent->is_child_on_side(parent->which_child_am_i(p), side))
               break;
             p = parent;
           }
         // We're on that side of our top_parent; return it
         if (!p)
           return side;
       }
    }

  // if we get here, we found elem in the data structure but not
  // the requested boundary id, so return the default value
  return libMesh::invalid_uint;  
}

void BoundaryInfo::build_node_boundary_ids(std::vector<short int> &b_ids)
{
  b_ids.clear();

  std::multimap<const Node*, short int>::const_iterator pos
    = _boundary_node_id.begin();

  for (; pos != _boundary_node_id.end(); ++pos)
    {
      short int id = pos->second;
      
      if(std::find(b_ids.begin(),b_ids.end(),id) == b_ids.end())
        b_ids.push_back(id);
    }
}

void BoundaryInfo::build_side_boundary_ids(std::vector<short int> &b_ids)
{
  b_ids.clear();

  std::multimap<const Elem*, std::pair<unsigned short int, short int> >::const_iterator pos
    = _boundary_side_id.begin();

  for (; pos != _boundary_side_id.end(); ++pos)
    {
      short int id = pos->second.second;
      
      if(std::find(b_ids.begin(),b_ids.end(),id) == b_ids.end())
        b_ids.push_back(id);
    }
}

unsigned int BoundaryInfo::n_boundary_conds () const
{
  // in serial we know the number of bcs from the
  // size of the container
  if (_mesh.is_serial())
    return _boundary_side_id.size();

  // in parallel we need to sum the number of local bcs
  parallel_only();
  
  unsigned int nbcs=0;

  std::multimap<const Elem*,
                std::pair<unsigned short int,
                          short int> >::const_iterator pos;

  for (pos=_boundary_side_id.begin(); pos != _boundary_side_id.end(); ++pos)
    if (pos->first->processor_id() == libMesh::processor_id())
      nbcs++;
  
  Parallel::sum (nbcs);
  
  return nbcs;
}



void BoundaryInfo::build_node_list (std::vector<unsigned int>& nl,
				    std::vector<short int>&    il) const
{
  // Reserve the size, then use push_back
  nl.reserve (_boundary_node_id.size());
  il.reserve (_boundary_node_id.size());
  
  std::multimap<const Node*, short int>::const_iterator pos
    = _boundary_node_id.begin();

  for (; pos != _boundary_node_id.end(); ++pos)
    {
      nl.push_back (pos->first->id());
      il.push_back (pos->second);
    }
}


void
BoundaryInfo::build_node_list_from_side_list()
{
  std::multimap<const Elem*,
                std::pair<unsigned short int,
                          short int> >::const_iterator pos;
  
  //Loop over the side list
  for (pos=_boundary_side_id.begin(); pos != _boundary_side_id.end(); ++pos)
  {
    //Need to loop over the sides of any possible children
    std::vector< const Elem * > family;
#ifdef LIBMESH_ENABLE_AMR
    pos->first->active_family_tree_by_side (family, pos->second.first);
#else
    family.push_back(pos->first);
#endif

    for(unsigned int elem_it=0; elem_it < family.size(); elem_it++)
    {
      const Elem * cur_elem = family[elem_it];
      
      AutoPtr<Elem> side = cur_elem->build_side(pos->second.first);

      //Add each node node on the side with the side's boundary id
      for(unsigned int i=0; i<side->n_nodes(); i++)
      {
        Node * node = side->get_node(i);
        
        this->add_node(node, pos->second.second);
      }
    }
  }
}

void BoundaryInfo::build_side_list (std::vector<unsigned int>&       el,
				    std::vector<unsigned short int>& sl,
				    std::vector<short int>&          il) const
{
  // Reserve the size, then use push_back
  el.reserve (_boundary_side_id.size());
  sl.reserve (_boundary_side_id.size());
  il.reserve (_boundary_side_id.size());

  std::multimap<const Elem*,
                std::pair<unsigned short int,
                          short int> >::const_iterator pos;

  for (pos=_boundary_side_id.begin(); pos != _boundary_side_id.end();
       ++pos)
    {
      el.push_back (pos->first->id());
      sl.push_back (pos->second.first);
      il.push_back (pos->second.second);
    }
}



void BoundaryInfo::print_info(std::ostream& out) const
{
  // Print out the nodal BCs
  if (!_boundary_node_id.empty())
    {
      out << "Nodal Boundary conditions:" << std::endl
	  << "--------------------------" << std::endl
	  << "  (Node No., ID)               " << std::endl;

//       std::for_each(_boundary_node_id.begin(),
// 		    _boundary_node_id.end(),
// 		    PrintNodeInfo());

      std::multimap<const Node*, short int>::const_iterator it        = _boundary_node_id.begin();
      const std::multimap<const Node*, short int>::const_iterator end = _boundary_node_id.end();

      for (; it != end; ++it)
	out << "  (" << (*it).first->id()
	    << ", "  << (*it).second
	    << ")"  << std::endl;
    }
  
  // Print out the element BCs
  if (!_boundary_side_id.empty())
    {
      out << std::endl
	  << "Side Boundary conditions:" << std::endl
	  << "-------------------------" << std::endl
	  << "  (Elem No., Side No., ID)      " << std::endl;

//       std::for_each(_boundary_side_id.begin(),
// 		    _boundary_side_id.end(),
//   		    PrintSideInfo());

      std::multimap<const Elem*,
	std::pair<unsigned short int, short int> >::const_iterator it = _boundary_side_id.begin();
      const std::multimap<const Elem*,
	std::pair<unsigned short int, short int> >::const_iterator end = _boundary_side_id.end();

      for (; it != end; ++it)
        out << "  (" << (*it).first->id()
	    << ", "  << (*it).second.first
	    << ", "  << (*it).second.second 
	    << ")"   << std::endl;
    }
}



void BoundaryInfo::print_summary(std::ostream& out) const
{
  // Print out the nodal BCs
  if (!_boundary_node_id.empty())
    {
      out << "Nodal Boundary conditions:" << std::endl
	  << "--------------------------" << std::endl
	  << "  (ID, number of nodes)   " << std::endl;

      std::map<short int, unsigned int> ID_counts;

      std::multimap<const Node*, short int>::const_iterator it        = _boundary_node_id.begin();
      const std::multimap<const Node*, short int>::const_iterator end = _boundary_node_id.end();

      for (; it != end; ++it)
        ID_counts[(*it).second]++;

      std::map<short int, unsigned int>::const_iterator ID_it        = ID_counts.begin();
      const std::map<short int, unsigned int>::const_iterator ID_end = ID_counts.end();

      for (; ID_it != ID_end; ++ID_it)
	out << "  (" << (*ID_it).first
	    << ", "  << (*ID_it).second
	    << ")"  << std::endl;
    }
  
  // Print out the element BCs
  if (!_boundary_side_id.empty())
    {
      out << std::endl
	  << "Side Boundary conditions:" << std::endl
	  << "-------------------------" << std::endl
	  << "  (ID, number of sides)   " << std::endl;

      std::map<short int, unsigned int> ID_counts;

      std::multimap<const Elem*,
	std::pair<unsigned short int, short int> >::const_iterator it = _boundary_side_id.begin();
      const std::multimap<const Elem*,
	std::pair<unsigned short int, short int> >::const_iterator end = _boundary_side_id.end();

      for (; it != end; ++it)
        ID_counts[(*it).second.second]++;

      std::map<short int, unsigned int>::const_iterator ID_it        = ID_counts.begin();
      const std::map<short int, unsigned int>::const_iterator ID_end = ID_counts.end();

      for (; ID_it != ID_end; ++ID_it)
	out << "  (" << (*ID_it).first
	    << ", "  << (*ID_it).second
	    << ")"  << std::endl;
    }
}

} // namespace libMesh

// $Id: boundary_info.C,v 1.21 2003-05-15 23:34:35 benkirk Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2002  Benjamin S. Kirk, John W. Peterson
  
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
#include <algorithm>
#include <fstream>

// Local includes
#include "mesh_config.h"
#include "boundary_info.h"
#include "mesh_base.h"
#include "boundary_mesh.h"



//------------------------------------------------------
// BoundaryInfo static member initializations
const short int BoundaryInfo::invalid_id = -1234;



//------------------------------------------------------
// BoundaryInfo functions
BoundaryInfo::BoundaryInfo(const MeshBase& m) :
  mesh           (m)
{
}



BoundaryInfo::~BoundaryInfo()
{
  clear();
}



void BoundaryInfo::clear()
{
  boundary_node_id.clear();
  boundary_side_id.clear();
  boundary_ids.clear();
}



void BoundaryInfo::sync(BoundaryMesh& boundary_mesh)
{
  boundary_mesh.clear();

  /**
   * Re-create the boundary mesh.
   */
  
  std::map<unsigned int, unsigned int> new_node_numbers;
  
  boundary_mesh.set_n_subdomains() = n_boundary_ids();


  // Add sides to the structure.
  std::map<short int, unsigned int> id_map;

  // Original Code
  //     unsigned int cnt = 0;
  //     for (std::set<short int>::iterator pos = boundary_ids.begin();
  // 	 pos != boundary_ids.end(); ++pos)
  //       id_map[*pos] = cnt++;
  
  //     id_map[invalid_id] = cnt;

    
  // New code 
  std::for_each(boundary_ids.begin(),
		boundary_ids.end(),
		Fill(id_map));
    
    

  boundary_mesh.set_n_subdomains() = id_map.size();
  boundary_mesh.set_n_processors() = id_map.size();

  // Add additional sides that aren't flagged with boundary conditions
  const_active_elem_iterator       el     (mesh.elements_begin());
  const const_active_elem_iterator end_el (mesh.elements_end());
  
  for ( ; el != end_el; ++el)
    {
      const Elem* elem = *el;
      
      for (unsigned int s=0; s<elem->n_sides(); s++)
	if (elem->neighbor(s) == NULL) // on the boundary
	  {
	    // Build the side
	    AutoPtr<Elem> side (elem->build_side(s));
	    
	    // The side lives on the same processor as the parent
	    //side->set_processor_id() = elem->processor_id();
	    
	    // Find the right id number for that side
	    std::pair<std::multimap<const Elem*,
		      std::pair<unsigned short int, short int> >::iterator,
		      std::multimap<const Elem*,
		      std::pair<unsigned short int, short int> >::iterator > 
	      pos = boundary_side_id.equal_range(elem);

	    while (pos.first != pos.second)
	      {
		if (pos.first->second.first == s) // already flagged with a boundary condition
		  {
		    side->set_subdomain_id() =
		      id_map[pos.first->second.second];
		    
		    side->set_processor_id() =
		      side->subdomain_id();
		    break;
		  }
		
		++pos.first;
	      }

	    // either the element wasn't found or side s
	    // doesn't have a booundary condition
	    if (pos.first == pos.second)
	      {
		side->set_subdomain_id() = id_map[invalid_id];
	      }
	    
	    // Add the side
	    boundary_mesh.add_elem(side.release());
	  }
    }

  // Copy over the nodes
  boundary_mesh._nodes = mesh._nodes;
}



void BoundaryInfo::add_node(const unsigned int node,
			    const short int id)
{
  add_node (mesh.node_ptr(node), id);
}



void BoundaryInfo::add_node(const Node* node,
			    const short int id)
{
  if (id == invalid_id)
    {
      std::cerr << "ERROR: You may not set a boundary ID of "
		<< invalid_id << std::endl
		<< " That is reserved for internal use.\n"
		<< std::endl;

      error();
    }
  
  boundary_node_id[node] = id;
  boundary_ids.insert(id);
}



void BoundaryInfo::add_side(const unsigned int e,
			    const unsigned short int side,
			    const short int id)
{
  add_side (mesh.elem(e), side, id);
}



void BoundaryInfo::add_side(const Elem* elem,
			    const unsigned short int side,
			    const short int id)
{
  if (id == invalid_id)
    {
      std::cerr << "ERROR: You may not set a boundary ID of "
		<< invalid_id << std::endl
		<< " That is reserved for internal use.\n"
		<< std::endl;

      error();
    }
  
  std::pair<unsigned short int, short int> p(side,id);
  std::pair<const Elem*, std::pair<unsigned short int, short int> >
    kv (elem, p);
  
  boundary_side_id.insert(kv);
  boundary_ids.insert(id);

  // Possilby add the nodes of the side,
  // if they aren't already there. MGF meshes
  // seem to cause some trouble here, so don't
  // do this if the library is configured with
  // --enable-mgf-workaround
#ifndef ENABLE_MGF_WORKAROUND 
  {
    assert (side < elem->n_sides());
    
    AutoPtr<Elem> side_elem(elem->build_side(side));

    for (unsigned int n=0; n<side_elem->n_nodes(); n++)
      if (boundary_id(side_elem->get_node(n)) == invalid_id)
 	add_node(side_elem->get_node(n), id);
  }
#endif
  
}



short int BoundaryInfo::boundary_id(const Node* node) const
{ 
  std::map<const Node*, short int>::const_iterator n
    = boundary_node_id.find(node);

  // node not in the data structure
  if (n == boundary_node_id.end())
    return invalid_id;

  return n->second;
}



short int BoundaryInfo::boundary_id(const Elem* elem,
				    const unsigned short int side) const
{ 
  std::pair<std::multimap<const Elem*,
                          std::pair<unsigned short int, short int> >::const_iterator,
            std::multimap<const Elem*,
                          std::pair<unsigned short int, short int> >::const_iterator > 
    e = boundary_side_id.equal_range(elem);

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



void BoundaryInfo::build_node_list (std::vector<unsigned int>& nl,
				    std::vector<short int>&    il) const
{
  // Reserve the size, then use push_back
  nl.reserve (boundary_node_id.size());
  il.reserve (boundary_node_id.size());
  
  std::map<const Node*, short int>::const_iterator pos;

  for (pos=boundary_node_id.begin(); pos != boundary_node_id.end();
       ++pos)
    {
      nl.push_back (pos->first->id());
      il.push_back (pos->second);
    }
}



void BoundaryInfo::build_side_list (std::vector<unsigned int>&       el,
				    std::vector<unsigned short int>& sl,
				    std::vector<short int>&          il) const
{
  // Reserve the size, then use push_back
  el.reserve (boundary_side_id.size());
  sl.reserve (boundary_side_id.size());
  il.reserve (boundary_side_id.size());

  std::multimap<const Elem*,
                std::pair<unsigned short int,
                          short int> >::const_iterator pos;

  for (pos=boundary_side_id.begin(); pos != boundary_side_id.end();
       ++pos)
    {
      el.push_back (pos->first->id());
      sl.push_back (pos->second.first);
      il.push_back (pos->second.second);
    }
}



void BoundaryInfo::print_info() const
{
  // Print out the nodal BCs
  if (!boundary_node_id.empty())
    {
      std::cout << "Nodal Boundary conditions:" << std::endl
		<< "--------------------------" << std::endl
		<< "  (Node No., ID)               " << std::endl;

      std::for_each(boundary_node_id.begin(),
		    boundary_node_id.end(),
		    PrintNodeInfo());
    }

  // Print out the element BCs
  if (!boundary_side_id.empty())
    {
      std::cout << "Side Boundary conditions:" << std::endl
		<< "-------------------------" << std::endl
		<< "  (Elem No., Side No., ID)      " << std::endl;

      std::for_each(boundary_side_id.begin(),
		    boundary_side_id.end(),
  		    PrintSideInfo()); 
    }
}

// $Id: boundary_info.C,v 1.48 2007-02-12 20:29:39 jwpeterson Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2005  Benjamin S. Kirk, John W. Peterson
  
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

//------------------------------------------------------
// BoundaryInfo static member initializations
const short int BoundaryInfo::invalid_id = -1234;



//------------------------------------------------------
// BoundaryInfo functions
BoundaryInfo::BoundaryInfo(const MeshBase& m) :
  _mesh (m)
{
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



void BoundaryInfo::sync(BoundaryMesh& boundary_mesh,
			MeshData*     boundary_mesh_data,
			MeshData*     this_mesh_data)
{
  boundary_mesh.clear();

  /**
   * Re-create the boundary mesh.
   */
  
  std::map<unsigned int, unsigned int> new_node_numbers;
  
  boundary_mesh.set_n_subdomains() = this->n_boundary_ids();


  // Add sides to the structure.
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
    
    

  boundary_mesh.set_n_subdomains() = id_map.size();


  // Make individual copies of all the nodes in the current mesh
  // and add them to the boundary mesh.  Yes, this is overkill because
  // all of the current mesh nodes will not end up in the the boundary
  // mesh.  These nodes can be trimmed later via a call to prepare_for_use().
  {
    assert (boundary_mesh.n_nodes() == 0);
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
		    
		    side->processor_id() =
		      side->subdomain_id();
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
		assert (new_node->id() == new_elem->node(nn));

		// Assign the new node pointer
		new_elem->set_node(nn) = new_node;
	      }
	  }
    } // end loop over active elements

  
  
  // When desired, copy the MeshData
  // to the boundary_mesh
  if ((boundary_mesh_data != NULL) && (this_mesh_data != NULL))
    boundary_mesh_data->assign(*this_mesh_data);

  // Trim any un-used nodes from the Mesh
  // boundary_mesh.prepare_for_use();
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
      std::cerr << "ERROR: You may not set a boundary ID of "
		<< invalid_id << std::endl
		<< " That is reserved for internal use.\n"
		<< std::endl;

      error();
    }
  
  _boundary_node_id[node] = id;
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
  assert (elem != NULL);

  // Only add BCs for level-0 elements.
  assert (elem->level() == 0);
  
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
  
  _boundary_side_id.insert(kv);
  _boundary_ids.insert(id);

  // Possilby add the nodes of the side,
  // if they aren't already there.
  {
    // MGF meshes seem to cause some trouble here, so don't
    // do this if the library is initialized with
    // --enable-mgf-workaround as an argument
    // command line can't change at run-time, so use a static here.
    static bool mgf_workaround =
      libMesh::on_command_line("--enable-mgf-workaround");

    if (!mgf_workaround)
      {
	assert (side < elem->n_sides());
	
	AutoPtr<Elem> side_elem(elem->build_side(side));
	
	for (unsigned int n=0; n<side_elem->n_nodes(); n++)
	  if (this->boundary_id(side_elem->get_node(n)) == invalid_id)
	    this->add_node(side_elem->get_node(n), id);
      }
  }  
}



short int BoundaryInfo::boundary_id(const Node* node) const
{ 
  std::map<const Node*, short int>::const_iterator
    n = _boundary_node_id.find(node);

  // node not in the data structure
  if (n == _boundary_node_id.end())
    return invalid_id;

  return n->second;
}



short int BoundaryInfo::boundary_id(const Elem* const elem,
				    const unsigned short int side) const
{
  assert (elem != NULL);

  // Only level-0 elements store BCs.  If this is not a level-0
  // element get its level-0 parent and infer the BCs.
  const Elem*  searched_elem = elem;
  if (elem->level() != 0)
    searched_elem = elem->top_parent ();
  
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



void BoundaryInfo::build_node_list (std::vector<unsigned int>& nl,
				    std::vector<short int>&    il) const
{
  // Reserve the size, then use push_back
  nl.reserve (_boundary_node_id.size());
  il.reserve (_boundary_node_id.size());
  
  std::map<const Node*, short int>::const_iterator pos
    = _boundary_node_id.begin();

  for (; pos != _boundary_node_id.end(); ++pos)
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



void BoundaryInfo::print_info() const
{
  // Print out the nodal BCs
  if (!_boundary_node_id.empty())
    {
      std::cout << "Nodal Boundary conditions:" << std::endl
		<< "--------------------------" << std::endl
		<< "  (Node No., ID)               " << std::endl;

//       std::for_each(_boundary_node_id.begin(),
// 		    _boundary_node_id.end(),
// 		    PrintNodeInfo());

      std::map<const Node*, short int>::const_iterator it        = _boundary_node_id.begin();
      const std::map<const Node*, short int>::const_iterator end = _boundary_node_id.end();

      for (; it != end; ++it)
	std::cout << "  (" << (*it).first->id()
		  << ", "  << (*it).second
		  << ")"  << std::endl;
    }
  
  // Print out the element BCs
  if (!_boundary_side_id.empty())
    {
      std::cout << std::endl
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
       std::cout << "  (" << (*it).first->id()
		 << ", "  << (*it).second.first
		 << ", "  << (*it).second.second 
		 << ")"   << std::endl;
    }
}

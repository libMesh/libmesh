// $Id: boundary_info.C,v 1.12 2003-02-13 22:56:12 benkirk Exp $

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
#include "mesh.h"


//------------------------------------------------------
// BoundaryInfo static member initializations
const short int BoundaryInfo::invalid_id = -1234;


//------------------------------------------------------
// BoundaryInfo functions
BoundaryInfo::BoundaryInfo(unsigned int d,
			   const Mesh& m) :
  boundary_mesh(d-1, m.processor_id()),
  dim(d),
  mesh(m)
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
  boundary_values.clear();
  
  node_list.clear();
  elem_list.clear();
  side_list.clear();
  node_id_list.clear();
  elem_id_list.clear();

  boundary_mesh.clear();
}



void BoundaryInfo::sync()
{
  error();
  
  //   boundary_mesh.clear();
  
  //   /**
  //    * At this point we have a list of elements stored in
  //    * boundary_mesh, but no nodes.  Furthermore, the connectivity
  //    * for the stored elements is in terms of the _global_ node
  //    * numbers.  In this routine we will renumber that connectivity
  //    * and create the necessary nodes in the boundary_mesh.
  //    */

  //   std::map<unsigned int, unsigned int> new_node_numbers;
  //   unsigned int next_node_number=0;
  
  //   boundary_mesh.set_n_subdomains() = n_boundary_ids();


  //   // Add sides to the structure.
  //   {
  //     std::map<short int, unsigned int> id_map;

  //     // Original Code
  //     //     unsigned int cnt = 0;
  //     //     for (std::set<short int>::iterator pos = boundary_ids.begin();
  //     // 	 pos != boundary_ids.end(); ++pos)
  //     //       id_map[*pos] = cnt++;
    
  //     //     id_map[invalid_id] = cnt;

    
  //     // New code 
  //     std::for_each(boundary_ids.begin(),
  //      		  boundary_ids.end(),
  //      		  Fill(id_map));
    
    

  //     boundary_mesh.set_n_subdomains() = id_map.size();

  //     // Add additional sides that aren't flagged with boundary conditions
  //     for (unsigned int e=0; e<mesh.n_elem(); e++)
  //       if (mesh.elem(e)->active())
  // 	for (unsigned int s=0; s<mesh.elem(e)->n_sides(); s++)
  // 	  if (mesh.elem(e)->neighbor(s) == NULL) // on the boundary
  // 	    {
  // 	      // Add the side
  // 	      Elem* elem = &mesh.elem(e)->side(s);
  // 	      boundary_mesh.add_elem(&elem);

  // 	      // The side lives on the same processor as the parent
  // 	      elem->processor_id() = mesh.elem(e)->processor_id();

  // 	      // Find the right id number for that side
  // 	      std::pair<std::multimap<const Elem*,
  // 	                              std::pair<unsigned short int, short int> >::iterator,
  // 	                std::multimap<const Elem*,
  //                                       std::pair<unsigned short int, short int> >::iterator > 
  // 		pos = boundary_side_id.equal_range(mesh.elem(e));

  // 	      while (pos.first != pos.second)
  // 		{
  // 		  if (pos.first->second.first == s) // already flagged with a boundary condition
  // 		    {
  // 		      elem->subdomain_id() =
  // 			id_map[pos.first->second.second];
  // 		      break;
  // 		    }
		  
  // 		  ++pos.first;
  // 		}

  // 	      // either the element wasn't found or side s
  // 	      // doesn't have a booundary condition
  // 	      if (pos.first == pos.second)
  // 		{
  // 		  elem->subdomain_id() = id_map[invalid_id];
  // 		}
  // 	    }
  //   }


  
      
  //   for (unsigned int e=0; e<boundary_mesh.n_elem(); e++)
  //     for (unsigned int n=0; n<boundary_mesh.elem(e)->n_nodes(); n++)
  //       {
  // 	const unsigned int node_number = boundary_mesh.elem(e)->node(n);
	
  // 	std::map<unsigned int, unsigned int>::iterator
  // 	  pos = new_node_numbers.find(node_number);

  // 	if (pos == new_node_numbers.end())
  // 	  {
  // 	    new_node_numbers[node_number] =
  // 	      next_node_number++;
	    
  // 	    boundary_mesh.add_point(mesh.point(node_number));
  // 	  }
  //       }

  

  
  //   for (unsigned int e=0; e<boundary_mesh.n_elem(); e++)
  //     for (unsigned int n=0; n<boundary_mesh.elem(e)->n_nodes(); n++)
  //       {
  // 	const unsigned int old_number = boundary_mesh.elem(e)->node(n); 
	    
  // 	boundary_mesh.elem(e)->node(n) = new_node_numbers[old_number];
  //       }
}




void BoundaryInfo::add_node(const unsigned int node,
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

  node_list.push_back(node);
  node_id_list.push_back(id);
}




void BoundaryInfo::add_side(const unsigned int e,
			    const unsigned short int side,
			    const short int id)
{
  elem_list.push_back(e);
  side_list.push_back(side);
  elem_id_list.push_back(id);

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
      if (boundary_id(side_elem->node(n)) == invalid_id)
 	add_node(side_elem->node(n), id);
  }
#endif
  
}




void BoundaryInfo::read_shanee_boundary (const std::string& name)
{
  std::ifstream in(name.c_str());

  read_shanee_boundary(in);

  return;
}




short int BoundaryInfo::boundary_id(const unsigned int node) const
{ 
  std::map<unsigned int, short int>::const_iterator n
    = boundary_node_id.find(node);

  // node not in the data structure
  if (n == boundary_node_id.end())
    return invalid_id;

  return n->second;
}



short int BoundaryInfo::boundary_id(const unsigned int e,
				    const unsigned short int side) const
{
  return boundary_id(mesh.elem(e), side);
}



short int BoundaryInfo::boundary_id(const Elem* elem,
				    const unsigned short int side) const
{ 
  std::pair<std::multimap<const Elem*,
    std::pair<unsigned short int, short int> >::const_iterator,
    std::multimap<const Elem*,
    std::pair<unsigned short int, short int> >::const_iterator > 
    e=boundary_side_id.equal_range(elem);

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




void BoundaryInfo::read_shanee_boundary(std::istream& in)
{
  assert (in);
  assert (dim == 2);

  unsigned int n_elem, n_nds, n0, n1, id;
  
  in >> n_elem
     >> n_nds;

  for (unsigned int elem=0; elem<n_elem; elem++)
    {
      in >> n0 >> n1 >> id;
      
      //      Edge2 edge(n0,n1);

      //      add_edge(edge, static_cast<short int>(id));
      add_node(n0, static_cast<short int>(id));
      add_node(n1, static_cast<short int>(id));
      
    }
  
  for (unsigned int node=0; node<n_nds; node++)
    {
      in >> n0 >> id;
      
      add_node(n0, static_cast<short int>(id));      
    }    
}



void BoundaryInfo::add_boundary_values(const unsigned int node,
				       const std::vector<Real> values,
				       const short int id)
{
  add_node(node, id);
  boundary_values.push_back(std::make_pair(node, values));
}


std::vector<Real> BoundaryInfo::get_boundary_values(const unsigned int node) const
{
  std::vector<std::pair<unsigned int,
              std::vector<Real> > >::const_iterator pos;
  
  for (pos=boundary_values.begin(); pos!=boundary_values.end(); ++pos)
    {
      if (pos->first == node)
	{
	  return pos->second;
	}
    }

  std::cerr << "ERROR: No boundary values are specified for Node: "
	    << node << std::endl;

  error();
  std::vector<Real> v;
  return v;
}


void BoundaryInfo::print_info() const
{
  // Print out the nodal BCs
  if (!boundary_node_id.empty())
    {
      std::cout << "Nodal Boundary conditions:" << std::endl
		<< "--------------------------" << std::endl
		<< "  (Node No., ID)               " << std::endl;

      
      // Original Code
      // for (std::map<unsigned int, short int>::const_iterator
      // 	     it = boundary_node_id.begin(); it != boundary_node_id.end();
      // 	   ++it)
      // 	std::cout << "  (" << it->first
      // 		  << ", "  << it->second
      //

      // New code
      std::for_each(boundary_node_id.begin(),
		    boundary_node_id.end(),
		    PrintNodeInfo());
    }

  // Print out the element BCs
  if (!boundary_side_id.empty())
    {
      // This map must remain local to this scope.
      std::map<const Elem*, unsigned int> elem_star_to_num;

      for (unsigned int e=0; e<mesh.n_elem(); e++)
	elem_star_to_num[mesh.elem(e)] = e;      


      std::cout << "Side Boundary conditions:" << std::endl
		<< "-------------------------" << std::endl
		<< "  (Elem No., Side No., ID)      " << std::endl;

      // Original Code
      // for (std::map<const Elem*, std::pair<unsigned short int, short int> >::const_iterator
      // 	     it = boundary_side_id.begin(); it != boundary_side_id.end();
      // 	   ++it)
      // 	std::cout << "  (" << elem_star_to_num[it->first]
      // 		  << ", "  << it->second.first
      // 		  << ", "  << it->second.second 
      // 		  << ")"   << std::endl;

      // New code
      std::for_each(boundary_side_id.begin(),
		    boundary_side_id.end(),
  		    PrintSideInfo(elem_star_to_num));

      
    }
  
}




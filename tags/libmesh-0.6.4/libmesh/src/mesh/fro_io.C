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
#include <fstream>
#include <deque>
#include <map>

// Local includes
#include "libmesh_config.h"
#include "fro_io.h"
#include "mesh_base.h"
#include "boundary_info.h"
#include "elem.h"



// ------------------------------------------------------------
// FroIO  members
void FroIO::write (const std::string& fname)
{
  if (libMesh::processor_id() == 0)
    {
      // Open the output file stream
      std::ofstream out (fname.c_str());
      libmesh_assert (out.good());

      // Make sure it opened correctly
      if (!out.good())
        libmesh_file_error(fname.c_str());

      // Get a reference to the mesh
      const MeshBase& mesh = MeshOutput<MeshBase>::mesh();

      // Write the header
      out << mesh.n_elem()  << " "
	  << mesh.n_nodes() << " "
	  << "0 0 "
	  << mesh.boundary_info->n_boundary_ids()  << " 1\n";

      // Write the nodes -- 1-based!
      for (unsigned int n=0; n<mesh.n_nodes(); n++)
	out << n+1 << " \t"
	    << std::scientific
	    << std::setprecision(12)
	    << mesh.point(n)(0) << " \t"
	    << mesh.point(n)(1) << " \t"
	    << 0. << '\n';

      // Write the elements -- 1-based!
      MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
      const MeshBase::const_element_iterator end = mesh.active_elements_end();
      unsigned int e=0;
      
      for ( ; it != end; ++it)
	{
	  // .fro likes TRI3's
	  if ((*it)->type() != TRI3)
	    {
	      std::cerr << "ERROR:  .fro format only valid for triangles!\n"
			<< "  writing of " << fname << " aborted.\n"
			<< std::endl;
	      libmesh_error();
	    }
	  
	  out << ++e << " \t";

 	  for (unsigned int n=0; n<(*it)->n_nodes(); n++)
 	    out << (*it)->node(n)+1 << " \t";

// 	  // LHS -> RHS Mapping, for inverted triangles
// 	  out << (*it)->node(0)+1 << " \t";
// 	  out << (*it)->node(2)+1 << " \t";
// 	  out << (*it)->node(1)+1 << " \t";

	  out << "1\n";
	}

      // Write BCs.
      {
	const std::set<short int>& bc_ids =
	  mesh.boundary_info->get_boundary_ids();

 	std::vector<unsigned int>       el;
 	std::vector<unsigned short int> sl;
 	std::vector<short int>          il;
	
 	mesh.boundary_info->build_side_list (el, sl, il);


	// Map the boundary ids into [1,n_bc_ids],
	// treat them one at a time.
	short int bc_id=0;
	for (std::set<short int>::const_iterator id = bc_ids.begin();
	     id != bc_ids.end(); ++id)
	  {
	    std::deque<unsigned int> node_list;
	    
	    std::map<unsigned int, unsigned int>
	      forward_edges, backward_edges;
	    
	    // Get all sides on this element with the relevant BC id.
	    for (unsigned int e=0; e<el.size(); e++)
	      if (il[e] == *id)
		{
		  // need to build up node_list as a sorted array of edge nodes...
		  // for the following:
		  // a---b---c---d---e
		  // node_list [ a b c d e];
		  //
		  // the issue is just how to get this out of the elem/side based data structure.
		  // the approach is to build up 'chain links' like this:
		  // a---b b---c c---d d---e
		  // and piece them together.
		  //
		  // so, for an arbitray edge n0---n1, we build the
		  // "forward_edges"  map n0-->n1
		  // "backward_edges" map n1-->n0
		  // and then start with one chain link, and add on...
		  //
		  AutoPtr<Elem> side = mesh.elem(el[e])->build_side(sl[e]);

		  const unsigned int
		    n0 = side->node(0),
		    n1 = side->node(1);

		  // insert into forward-edge set
		  forward_edges.insert (std::make_pair(n0, n1));

		  // insert into backward-edge set
		  backward_edges.insert (std::make_pair(n1, n0));

		  // go ahead and add one edge to the list -- this will give us the beginning of a
		  // chain to work from!
		  if (node_list.empty())
		    {
		      node_list.push_front(n0);
		      node_list.push_back (n1);
		    }
		}

	    // we now have the node_list with one edge, the forward_edges, and the backward_edges
	    // the node_list will be filled when (node_list.size() == (n_edges+1))
	    // until that is the case simply add on to the beginning and end of the node_list,
	    // building up a chain of ordered nodes...
	    const unsigned int n_edges = forward_edges.size();
	    
	    while (node_list.size() != (n_edges+1))
	      {
		const unsigned int
		  front_node = node_list.front(),
		  back_node  = node_list.back();

		// look for front_pair in the backward_edges list
		{
		  std::map<unsigned int, unsigned int>::iterator
		    pos = backward_edges.find(front_node);

		  if (pos != backward_edges.end())
		    {
		      node_list.push_front(pos->second);
		       
		      backward_edges.erase(pos);
		    }
		}

		// look for back_pair in the forward_edges list
		{
		  std::map<unsigned int, unsigned int>::iterator
		    pos = forward_edges.find(back_node);

		  if (pos != forward_edges.end())
		    {
		      node_list.push_back(pos->second);
		       
		      forward_edges.erase(pos);
		    }
		}

// 		std::cout << "node_list.size()=" << node_list.size()
// 			  << ", n_edges+1=" << n_edges+1 << std::endl;		  
	      }
	    	    

	    out << ++bc_id << " " << node_list.size() << '\n';

	    std::deque<unsigned int>::iterator pos = node_list.begin();
	    for ( ; pos != node_list.end(); ++pos)
		out << *pos+1 << " \t0\n";
	  }	
      }      
    }
}

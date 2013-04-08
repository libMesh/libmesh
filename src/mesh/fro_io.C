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



// C++ includes
#include <fstream>
#include <iomanip>
#include <iostream>
#include <deque>
#include <map>

// Local includes
#include "libmesh/libmesh_config.h"
#include "libmesh/fro_io.h"
#include "libmesh/mesh_base.h"
#include "libmesh/boundary_info.h"
#include "libmesh/elem.h"

namespace libMesh
{



// ------------------------------------------------------------
// FroIO  members
void FroIO::write (const std::string& fname)
{
  // We may need to gather a ParallelMesh to output it, making that
  // const qualifier in our constructor a dirty lie
  MeshSerializer serialize(const_cast<MeshBase&>(this->mesh()), !_is_parallel_format);

  if (this->mesh().processor_id() == 0)
    {
      // Open the output file stream
      std::ofstream out_stream (fname.c_str());
      libmesh_assert (out_stream.good());

      // Make sure it opened correctly
      if (!out_stream.good())
        libmesh_file_error(fname.c_str());

      // Get a reference to the mesh
      const MeshBase& the_mesh = MeshOutput<MeshBase>::mesh();

      // Write the header
      out_stream << the_mesh.n_elem()  << " "
	  << the_mesh.n_nodes() << " "
	  << "0 0 "
	  << the_mesh.boundary_info->n_boundary_ids()  << " 1\n";

      // Write the nodes -- 1-based!
      for (unsigned int n=0; n<the_mesh.n_nodes(); n++)
	out_stream << n+1 << " \t"
	    << std::scientific
	    << std::setprecision(12)
	    << the_mesh.point(n)(0) << " \t"
	    << the_mesh.point(n)(1) << " \t"
	    << 0. << '\n';

      // Write the elements -- 1-based!
      MeshBase::const_element_iterator       it  = the_mesh.active_elements_begin();
      const MeshBase::const_element_iterator end = the_mesh.active_elements_end();

      for (unsigned int e=0 ; it != end; ++it)
	{
	  // .fro likes TRI3's
	  if ((*it)->type() != TRI3)
	    {
	      libMesh::err << "ERROR:  .fro format only valid for triangles!\n"
			    << "  writing of " << fname << " aborted.\n"
			    << std::endl;
	      libmesh_error();
	    }

	  out_stream << ++e << " \t";

 	  for (unsigned int n=0; n<(*it)->n_nodes(); n++)
            out_stream << (*it)->node(n)+1 << " \t";

// 	  // LHS -> RHS Mapping, for inverted triangles
// 	  out_stream << (*it)->node(0)+1 << " \t";
// 	  out_stream << (*it)->node(2)+1 << " \t";
// 	  out_stream << (*it)->node(1)+1 << " \t";

	  out_stream << "1\n";
	}

      // Write BCs.
      {
	const std::set<boundary_id_type>& bc_ids =
	  the_mesh.boundary_info->get_boundary_ids();

 	std::vector<dof_id_type>        el;
 	std::vector<unsigned short int> sl;
 	std::vector<boundary_id_type>   il;

        the_mesh.boundary_info->build_side_list (el, sl, il);


	// Map the boundary ids into [1,n_bc_ids],
	// treat them one at a time.
	boundary_id_type bc_id=0;
	for (std::set<boundary_id_type>::const_iterator id = bc_ids.begin();
	     id != bc_ids.end(); ++id)
	  {
	    std::deque<dof_id_type> node_list;

	    std::map<dof_id_type, dof_id_type>
	      forward_edges, backward_edges;

	    // Get all sides on this element with the relevant BC id.
	    for (std::size_t e=0; e<el.size(); e++)
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
		  AutoPtr<Elem> side = the_mesh.elem(el[e])->build_side(sl[e]);

		  const dof_id_type
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
	    const std::size_t n_edges = forward_edges.size();

	    while (node_list.size() != (n_edges+1))
	      {
		const dof_id_type
		  front_node = node_list.front(),
		  back_node  = node_list.back();

		// look for front_pair in the backward_edges list
		{
		  std::map<dof_id_type, dof_id_type>::iterator
		    pos = backward_edges.find(front_node);

		  if (pos != backward_edges.end())
		    {
		      node_list.push_front(pos->second);

		      backward_edges.erase(pos);
		    }
		}

		// look for back_pair in the forward_edges list
		{
		  std::map<dof_id_type, dof_id_type>::iterator
		    pos = forward_edges.find(back_node);

		  if (pos != forward_edges.end())
		    {
		      node_list.push_back(pos->second);

		      forward_edges.erase(pos);
		    }
		}

// 		libMesh::out << "node_list.size()=" << node_list.size()
// 			      << ", n_edges+1=" << n_edges+1 << std::endl;
	      }


	    out_stream << ++bc_id << " " << node_list.size() << '\n';

	    std::deque<dof_id_type>::iterator pos = node_list.begin();
	    for ( ; pos != node_list.end(); ++pos)
		out_stream << *pos+1 << " \t0\n";
	  }
      }
    }
}

} // namespace libMesh

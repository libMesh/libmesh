// $Id: mesh_smoother_laplace.C,v 1.5 2003-07-15 12:40:12 benkirk Exp $

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


// Local includes
#include "mesh_smoother_laplace.h"
#include <algorithm> // for std::copy, std::sort



// Member functions for the Laplace smoother
void LaplaceMeshSmoother::smooth(unsigned int n_iterations)
{
  if (!_initialized)
    this->init();
  
  // Don't smooth the nodes on the boundary...
  // this would change the mesh geometry which
  // is probably not something we want!
  std::vector<bool> on_boundary;
  _mesh.find_boundary_nodes(on_boundary);

  for (unsigned int n=0; n<n_iterations; n++)
    {
      for (unsigned int i=0; i<_mesh.n_nodes(); ++i)
	{
	  if (!on_boundary[i])
	    {
	      // Smooth!
	      Real avg_x = 0.;
	      Real avg_y = 0.;
	      Real avg_z = 0.;
	      for (unsigned int j=0; j<_graph[i].size(); ++j)
		{
		  // Get a reference to the current node in
		  // the graph
		  const Node& node = _mesh.node(_graph[i][j]);

		  avg_x += node(0);
		  avg_y += node(1);
		  avg_z += node(2);
		}

	      assert (_graph[i].size() != 0);
	      avg_x /= _graph[i].size();
	      avg_y /= _graph[i].size();
	      avg_z /= _graph[i].size();

	      // Update the location of node(i)
	      Node& node = _mesh.node(i);
	      node(0) = avg_x;
	      node(1) = avg_y;
	      node(2) = avg_z;
	    }
	}
    }
}





void LaplaceMeshSmoother::init()
{
  switch (_mesh.mesh_dimension())
    {
      


      //TODO:[BSK] Fix this to work for refined meshes...  I think
      // the implementation was done quickly for Damien, who did not have
      // refined grids.  Fix it here and in the original Mesh member.
      
    case 2: // Stolen directly from build_L_graph in mesh_base.C
      {
	// Initialize space in the graph.  It is n_nodes
	// long and each node is assumed to be connected to
	// approximately 4 neighbors.
	_graph.resize(_mesh.n_nodes());
	for (unsigned int i=0; i<_mesh.n_nodes(); ++i)
	  _graph[i].reserve(4);
	
	const_active_elem_iterator       el (_mesh.elements_begin());
	const const_active_elem_iterator end(_mesh.elements_end());
	
	for (; el != end; ++el)
	  {
	    // Shortcut notation for simplicity
	    const Elem* elem = *el;
	    
	    for (unsigned int s=0; s<elem->n_neighbors(); s++)
	      {
		// Only operate on sides which are on the
		// boundary or for which the current element's
		// id is greater than its neighbor's. 
		if ((elem->neighbor(s) == NULL) ||
		    (elem->id() > elem->neighbor(s)->id()))
		  {
		    AutoPtr<Elem> side(elem->build_side(s));
		  
		    _graph[side->node(0)].push_back(side->node(1));
		    _graph[side->node(1)].push_back(side->node(0));
		}
	      }
	  }
	_initialized = true;
	break;
      }

    case 3: // Stolen blatantly from build_L_graph in mesh_base.C
      {
	// Initialize space in the graph.  In 3D, I've assumed
	// that each node was connected to approximately 3 neighbors.
	_graph.resize(_mesh.n_nodes());
	for (unsigned int i=0; i<_mesh.n_nodes(); ++i)
	  _graph[i].reserve(8);
	
	const_active_elem_iterator       el (_mesh.elements_begin());
	const const_active_elem_iterator end(_mesh.elements_end());

	for (; el != end; ++el)
	  {
	    // Shortcut notation for simplicity
	    const Elem* elem = *el;
	    
	    for (unsigned int f=0; f<elem->n_neighbors(); f++) // Loop over faces
	      if ((elem->neighbor(f) == NULL) ||
		  (elem->id() > elem->neighbor(f)->id()))
		{
		  AutoPtr<Elem> face(elem->build_side(f));
		
		  for (unsigned int s=0; s<face->n_neighbors(); s++) // Loop over face's edges
		    {
		      AutoPtr<Elem> side(face->build_side(s));
		    
		      // At this point, we just insert the node numbers
		      // again.  At the end we'll call sort and unique
		      // to make sure there are no duplicates
		      _graph[side->node(0)].push_back(side->node(1));
		      _graph[side->node(1)].push_back(side->node(0));
		    }
		}
	  }

	// Now call sort and unique to remove duplicate entries.
	for (unsigned int i=0; i<_mesh.n_nodes(); ++i)
	  {
	    std::sort  (_graph[i].begin(), _graph[i].end());
	    _graph[i].erase(std::unique(_graph[i].begin(), _graph[i].end()), _graph[i].end());
	  }
	
	_initialized = true;
	break;
      }

    default:
      {
	std::cerr << "At this time it is not possible "
		  << "to smooth a dimension "
		  << _mesh.mesh_dimension()
		  << "mesh.  Aborting..."
		  << std::endl;
	error();
      }
      
    }
}




void LaplaceMeshSmoother::print_graph()
{
  for (unsigned int i=0; i<_graph.size(); ++i)
    {
      std::cout << i << ": ";
      std::copy(_graph[i].begin(),
		_graph[i].end(),
		std::ostream_iterator<unsigned int>(std::cout, " "));
      std::cout << std::endl;
    }
}

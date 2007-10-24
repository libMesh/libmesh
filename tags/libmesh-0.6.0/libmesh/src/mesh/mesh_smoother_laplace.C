// $Id: mesh_smoother_laplace.C,v 1.17 2007-02-15 17:04:49 jwpeterson Exp $

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
#include <algorithm> // for std::copy, std::sort


// Local includes
#include "mesh_smoother_laplace.h"
#include "mesh_tools.h"
#include "elem.h"
#include "mesh.h"

// Member functions for the Laplace smoother
void LaplaceMeshSmoother::smooth(unsigned int n_iterations)
{
  if (!_initialized)
    this->init();
  
  // Don't smooth the nodes on the boundary...
  // this would change the mesh geometry which
  // is probably not something we want!
  std::vector<bool> on_boundary;
  MeshTools::find_boundary_nodes(_mesh, on_boundary);

  // We can only update the nodes after all new positions were
  // determined. We store the new positions here
  std::vector<Point> new_positions;

  for (unsigned int n=0; n<n_iterations; n++)
    {
      new_positions.resize(_mesh.n_nodes());
      
      for (MeshBase::node_iterator it  = _mesh.nodes_begin();
	   it != _mesh.nodes_end();
	   ++it) 
	{
	  Node* node = *it;
          // leave the boundary intact
          // Only relocate the nodes which are vertices of an element
          // All other entries of _graph (the secondary nodes) are empty
	  if (!on_boundary[node->id()] && (_graph[node->id()].size() > 0) )
	    {
              Point avg_position(0.,0.,0.);
	      for (unsigned int j=0; j<_graph[node->id()].size(); ++j)
                avg_position.add(_mesh.node(_graph[node->id()][j]));
              new_positions[node->id()] = avg_position /
                static_cast<Real>(_graph[node->id()].size());
	    }
	}
      
      // now update the node positions
      for (MeshBase::node_iterator it  = _mesh.nodes_begin();
	   it != _mesh.nodes_end();
	   ++it)
	{
	  Node* node = *it;
	  
	  if (!on_boundary[node->id()] && (_graph[node->id()].size() > 0) )
	    {
	      // Should call Point::op=
	      _mesh.node(node->id()) = new_positions[node->id()];
	    }
	}
    }
  
  // finally adjust the second order nodes (those located between vertices)
  // these nodes will be located between their adjacent nodes
  // do this element-wise
  MeshBase::element_iterator       el  = _mesh.active_elements_begin();
  const MeshBase::element_iterator end = _mesh.active_elements_end(); 
	
  for (; el != end; ++el)
    {
      // Constant handle for the element
      const Elem* elem = *el;

      // get the second order nodes (son)
      // their element indices start at n_vertices and go to n_nodes
      const unsigned int son_begin = elem->n_vertices();
      const unsigned int son_end   = elem->n_nodes();
      
      // loop over all second order nodes (son)
      for (unsigned int son=son_begin; son<son_end; son++)
        {
	  // Don't smooth second-order nodes which are on the boundary
	  if (!on_boundary[elem->node(son)])
	    {
	      const unsigned int n_adjacent_vertices =
		elem->n_second_order_adjacent_vertices(son);

	      // calculate the new position which is the average of the
	      // position of the adjacent vertices
	      Point avg_position(0,0,0);
	      for (unsigned int v=0; v<n_adjacent_vertices; v++)
		avg_position +=
		  _mesh.point( elem->node( elem->second_order_adjacent_vertex(son,v) ) );

	      _mesh.node(elem->node(son)) = avg_position / n_adjacent_vertices;
	    }
        }
    }
}


void LaplaceMeshSmoother::init()
{


  switch (_mesh.mesh_dimension())
    {
      
      // TODO:[BSK] Fix this to work for refined meshes...  I think
      // the implementation was done quickly for Damien, who did not have
      // refined grids.  Fix it here and in the original Mesh member.
      
    case 2: // Stolen directly from build_L_graph in mesh_base.C
      {
	// Initialize space in the graph.  It is n_nodes
	// long and each node is assumed to be connected to
	// approximately 4 neighbors.
	_graph.resize(_mesh.n_nodes());
// 	for (unsigned int i=0; i<_mesh.n_nodes(); ++i)
// 	  _graph[i].reserve(4);
	
	MeshBase::element_iterator       el  = _mesh.active_elements_begin();
	const MeshBase::element_iterator end = _mesh.active_elements_end(); 
	
	for (; el != end; ++el)
	  {
	    // Constant handle for the element
	    const Elem* elem = *el;
	    
	    for (unsigned int s=0; s<elem->n_neighbors(); s++)
	      {
		// Only operate on sides which are on the
		// boundary or for which the current element's
		// id is greater than its neighbor's.
                // Sides get only built once.
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
// 	for (unsigned int i=0; i<_mesh.n_nodes(); ++i)
// 	  _graph[i].reserve(8);
	
	MeshBase::element_iterator       el  = _mesh.active_elements_begin();
	const MeshBase::element_iterator end = _mesh.active_elements_end(); 

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




void LaplaceMeshSmoother::print_graph() const
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

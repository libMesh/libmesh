// $Id: mesh_modification.C,v 1.3 2004-01-03 15:37:43 benkirk Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2004  Benjamin S. Kirk, John W. Peterson
  
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
#include <map>


// Local includes
#include "mesh.h"
#include "libmesh.h"
#include "face_inf_quad4.h"
#include "face_inf_quad6.h"
#include "cell_inf_prism6.h"
#include "cell_inf_prism12.h"
#include "cell_inf_hex8.h"
#include "cell_inf_hex16.h"
#include "cell_inf_hex18.h"
#include "libmesh_logging.h"


// ------------------------------------------------------------
// Mesh class member functions for mesh modification
void Mesh::all_second_order (const bool full_ordered)
{
  /*
   * when the mesh is not prepared,
   * at least renumber the nodes and 
   * elements, so that the node ids
   * are correct
   */
  if (!this->_is_prepared)
    this->renumber_nodes_and_elements ();

  // does this work also in parallel?
  assert (this->n_processors() == 1);

  START_LOG("all_second_order()", "MeshBase");

  /*
   * this map helps in identifying second order
   * nodes.  Namely, a second-order node:
   * - edge node
   * - face node
   * - bubble node
   * is uniquely defined through a set of adjacent
   * vertices.  This set of adjacent vertices is
   * used to identify already added higher-order
   * nodes.  We are safe to use node id's since we
   * make sure that these are correctly numbered.
   */
  std::map<std::vector<unsigned int>, Node*> adj_vertices_to_so_nodes;


  /*
   * for speed-up of the \p add_point() method, we
   * can reserve memory.  Guess the number of additional
   * nodes for different dimensions
   */
  switch (this->mesh_dimension())
  {
    case 1:
      /*
       * in 1D, there can only be order-increase from Edge2
       * to Edge3.  Something like 1/2 of n_nodes() have
       * to be added
       */
      this->_nodes.reserve(static_cast<unsigned int>(1.5*this->_nodes.size()));
      break;

    case 2:
      /*
       * in 2D, either refine from Tri3 to Tri6 (double the nodes)
       * or from Quad4 to Quad8 (again, double) or Quad9 (2.25 that much)
       */
      this->_nodes.reserve(static_cast<unsigned int>(2*this->_nodes.size()));
      break;


    case 3:
      /*
       * in 3D, either refine from Tet4 to Tet10 (factor = 2.5) up to
       * Hex8 to Hex27 (something  > 3).  Since in 3D there _are_ already
       * quite some nodes, and since we do not want to overburden the memory by
       * a too conservative guess, use the lower bound
       */
      this->_nodes.reserve(static_cast<unsigned int>(2.5*this->_nodes.size()));
      break;
	
    default:
      // Hm?
      error();
  }



  /*
   * form a vector that will hold the node id's of
   * the vertices that are adjacent to the son-th
   * second-order node.  Pull this outside of the
   * loop so that silly compilers don't repeatedly
   * create and destroy the vector.
   */
  std::vector<unsigned int> adjacent_vertices_ids;


  /**
   * Loop over the low-ordered elements in the _elements vector.
   * First make sure they _are_ indeed low-order, and then replace
   * them with an equivalent second-order element.  Don't
   * forget to delete the low-order element, or else it will leak!
   */
  for (unsigned int e=0; e<_elements.size(); e++)
    {
      // the linear-order element
      Elem* lo_elem = _elements[e];

      assert (lo_elem != NULL);

      // make sure it is linear order
      if (lo_elem->default_order() != FIRST)
        {	  
	  std::cerr << "ERROR: This is not a linear element: type=" 
		    << lo_elem->type() << std::endl;
	  error();
	}

      // this does _not_ work for refined elements
      assert (lo_elem->level () == 0);

      /*
       * build the second-order equivalent, add to
       * the new_elements list.  Note that this here
       * is the only point where \p full_ordered
       * is necessary.  The remaining code works well
       * for either type of seconrd-order equivalent, e.g.
       * Hex20 or Hex27, as equivalents for Hex8
       */
      Elem* so_elem = Elem::build ( Elem::second_order_equivalent_type(lo_elem->type(), 
								       full_ordered) );
      assert (lo_elem->n_vertices() == so_elem->n_vertices());


      /*
       * By definition the vertices of the linear and
       * second order element are identically numbered.
       * transfer these.
       */
      for (unsigned int v=0; v < lo_elem->n_vertices(); v++)
	so_elem->set_node(v) = lo_elem->get_node(v);

      /*
       * Now handle the additional mid-side nodes.  This
       * is simply handled through a map that remembers
       * the already-added nodes.  This map maps the global
       * ids of the vertices (that uniquely define this 
       * higher-order node) to the new node. 
       * Notation: son = second-order node
       */
      const unsigned int son_begin = so_elem->n_vertices();
      const unsigned int son_end   = so_elem->n_nodes();
      

      for (unsigned int son=son_begin; son<son_end; son++)
        {
	  const unsigned int n_adjacent_vertices =
	    so_elem->n_second_order_adjacent_vertices(son);

	  adjacent_vertices_ids.resize(n_adjacent_vertices);
	  
	  for (unsigned int v=0; v<n_adjacent_vertices; v++)
	    adjacent_vertices_ids[v] =
	      so_elem->node( so_elem->second_order_adjacent_vertex(son,v) );

	  /*
	   * \p adjacent_vertices_ids is now in order of the current
	   * side.  sort it, so that comparisons  with the 
	   * \p adjacent_vertices_ids created through other elements' 
	   * sides can match
	   */
	  std::sort(adjacent_vertices_ids.begin(),
		    adjacent_vertices_ids.end());


	  // does this set of vertices already has a mid-node added?
	  std::pair<std::map<std::vector<unsigned int>, Node*>::iterator,
                    std::map<std::vector<unsigned int>, Node*>::iterator>	    
	    pos = adj_vertices_to_so_nodes.equal_range (adjacent_vertices_ids);

	  // no, not added yet
	  if (pos.first == pos.second)
	    {
	      /*
	       * for this set of vertices, there is no 
	       * second_order node yet.  Add it.
	       *
	       * compute the location of the new node as
	       * the average over the adjacent vertices.
	       */
	      Point new_location = this->point(adjacent_vertices_ids[0]);
	      for (unsigned int v=1; v<n_adjacent_vertices; v++)
		new_location += this->point(adjacent_vertices_ids[v]);

	      new_location /= static_cast<Real>(n_adjacent_vertices);
	      
	      // add the new point to the mesh
	      Node* so_node = this->add_point (new_location);

	      /* 
	       * insert the new node with its defining vertex
	       * set into the map, and relocate pos to this
	       * new entry, so that the so_elem can use
	       * \p pos for inserting the node
	       */
	      adj_vertices_to_so_nodes.insert(pos.first,
					      std::make_pair(adjacent_vertices_ids,
							     so_node));

	      so_elem->set_node(son) = so_node;
	    }
	  // yes, already added.
	  else
	    {
	      assert (pos.first->second != NULL);
	      
	      so_elem->set_node(son) = pos.first->second;
	    }
	}


      /**
       * If the linear element had any boundary conditions they
       * should be transfered to the second-order element, and then
       * removed from the BoundaryInfo data structure.
       */
      {
	assert (lo_elem->n_sides() == so_elem->n_sides());
	
	for (unsigned int s=0; s<lo_elem->n_sides(); s++)
	  {
	    const short int boundary_id =
	      this->boundary_info.boundary_id (lo_elem, s);
	    
	    if (boundary_id != this->boundary_info.invalid_id)
	      this->boundary_info.add_side (so_elem, s, boundary_id);
	  }
	
	/**
	 * We have taken any boundary conditions the low-ordered
	 * element may have had.  Since we are about to delete
	 * the low-ordered element, we should first un-associate
	 * any boundary conditions it has.
	 */
	this->boundary_info.remove (lo_elem);
      }

      
      /*
       * The new second-order element is ready.
       * Delete the linear element and replace it with
       * the second-order element.
       */
      delete lo_elem;
      _elements[e] = so_elem;
    }


  // we can clear the map
  adj_vertices_to_so_nodes.clear();


  STOP_LOG("all_second_order()", "MeshBase");

  // renumber nodes, elements etc
  this->prepare_for_use();
}

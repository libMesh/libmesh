// $Id: mesh_base_modification.C,v 1.1 2003-08-07 19:25:31 ddreyer Exp $

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
#include <map>


// Local includes
#include "mesh_base.h"
#include "libmesh.h"
#include "elem.h"
#include "mesh_logging.h"


// ------------------------------------------------------------
// Mesh class member functions for mesh modification
void MeshBase::all_second_order ()
{
  /*
   * when the mesh is not prepared,
   * at least renumber the nodes and 
   * elements, so that the node ids
   * are correct
   */
  if (!this->_is_prepared)
      renumber_nodes_and_elements ();

  // does this work also in parallel?
  assert (this->n_processors() == 1);

  START_LOG("all_second_order()", "MeshBase");

  /*
   * the vector holding the new second-order
   * replacement elements
   */
  std::vector<Elem*> new_elements;
  new_elements.reserve(this->n_elem());

  /*
   * this map helps in identifying second order
   * nodes.  Namely, a second-order node:
   * - edge node
   * - face node
   * - bubble node
   * is uniquely defined through a set of adjacent
   * vertices.  This set of adjacent vertices is
   * used to identify already added higher-order
   * nodes
   */
//  std::map<std::vector<Node*>, Node*> adj_vertices_to_so_nodes;
  // since we assert _is_prepared, the nodes are correctly numbered
  // use these global id's, instead of Node* here
  std::map<std::vector<unsigned int>, Node*> adj_vertices_to_so_nodes;


  /*
   * we will have to add nodes, so remember
   * the next free node id
   */
  unsigned int next_free_node_id = this->n_nodes();



  /*
   * iterate over all elements contained in the 
   * mesh
   */
  elem_iterator       old_elements_it  (this->elements_begin());
  const elem_iterator old_elements_end (this->elements_end());

  for (; old_elements_it != old_elements_end; ++old_elements_it)
    {
      // the linear-order element
      Elem* lo_elem = *old_elements_it;

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
       * the new_elements list
       */
      Elem* so_elem = Elem::build (lo_elem->second_order_equivalent_type());
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
	  const unsigned int n_adjacent_vertices = so_elem->n_second_order_adjacent_vertices(son);

	  /*
	   * form a vector that will hold the node id's of
	   * the vertices that are adjacent to the son-th
	   * second-order node
	   */
	  std::vector<unsigned int> adjacent_vertices_ids;
	  adjacent_vertices_ids.resize(n_adjacent_vertices);
	  for (unsigned int v=0; v<n_adjacent_vertices; v++)
	      adjacent_vertices_ids[v] = so_elem->node( so_elem->second_order_adjacent_vertex(son,v) );


	  // does this set of vertices already has a mid-node added?
	  std::map<std::vector<unsigned int>, Node*>::const_iterator pos =  
	      adj_vertices_to_so_nodes.find(adjacent_vertices_ids);

	  if (pos == adj_vertices_to_so_nodes.end())
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

	      // build the node
	      Node* so_node = Node::build(new_location,
					  next_free_node_id++);

	      /* 
	       * insert the new node with its defining vertex
	       * set into the map, and relocate pos to this
	       * new entry, so that the so_elem can use
	       * \p pos for inserting the node
	       */
	      adj_vertices_to_so_nodes.insert(std::make_pair(adjacent_vertices_ids, so_node));

	      so_elem->set_node(son) = so_node;
	    }
	  else
	      so_elem->set_node(son) = pos->second;

	}


      /*
       * The new second-order element is ready.
       * Add it to the new_elements vector
       */
      new_elements.push_back(so_elem);

    }



  /*
   * Now, the \p _elements vector has to be replaced
   * by the \p new_elements vector.  Delete the
   * old element, then put the new element in place.
   */
  {
    for (unsigned int e=0; e<_elements.size(); e++)
      {
	assert (_elements[e] != NULL);
	delete _elements[e];
	_elements[e] = new_elements[e];
      }

    // now we can safely clear our local \p new_elements vector
    new_elements.clear();
  }


  /*
   * The nodes are not that easily updated as the 
   * elements.  We have to keep the old nodes,
   * @e and add the new nodes, for which pointers
   * are currently only stored in the map
   * \p adj_vertices_to_so_nodes
   */
  {
    // first reserve memory for @e all nodes
    const unsigned int n_all_nodes = this->n_nodes() + adj_vertices_to_so_nodes.size();
    _nodes.reserve(n_all_nodes);

    std::map<std::vector<unsigned int>, Node*>::const_iterator new_nodes_it        = adj_vertices_to_so_nodes.begin();
    const std::map<std::vector<unsigned int>, Node*>::const_iterator new_nodes_end = adj_vertices_to_so_nodes.end();

    for (; new_nodes_it != new_nodes_end; ++new_nodes_it)
	_nodes.push_back(new_nodes_it->second);


    // now we can clear the map
    adj_vertices_to_so_nodes.clear();
  }


  STOP_LOG("all_second_order()", "MeshBase");

  // renumber nodes, elements etc
  this->prepare_for_use();
}


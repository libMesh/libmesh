// $Id: tree.h,v 1.6 2004-07-26 15:01:30 jwpeterson Exp $

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



#ifndef __tree_h__
#define __tree_h__

// C++ includes

// Local includes
#include "tree_node.h"
#include "tree_base.h"
#include "mesh_base.h"

/**
 * This class defines a tree that may be used for fast point
 * location in space.
 *
 * @author Benjamin S. Kirk, 2002
 */

// ------------------------------------------------------------
// Tree class definition
template <unsigned int N>
class Tree : public TreeBase
{
public:
  
  /**
   * \p enum defining how to build the tree.  \p NODES will populate
   * the tree with nodes and then replace the nodes with element
   * connectivity, \p ELEMENTS will populate the tree with the elements
   * directly.
   */
  enum BuildType {NODES=0,
		  ELEMENTS,
		  INVALID_BUILD_TYPE };
  
  /**
   * Constructor.
   */
  Tree (const MeshBase& m, const unsigned int level, BuildType bt=NODES);

  /**
   * Copy-constructor.
   */
  Tree (const Tree<N>& other_tree);

  /**
   * Destructor.
   */
  ~Tree() {}
  
  /**
   * Prints the nodes.
   */
  void print_nodes() const
  { std::cout << "Printing nodes...\n"; root.print_nodes(); }

  /**
   * Prints the nodes.
   */
  void print_elements() const
  { std::cout << "Printing elements...\n"; root.print_elements(); }
  
  /**
   * @returns the number of active bins.
   */
  unsigned int n_active_bins() const { return root.n_active_bins(); }

  /**
   * @returns a pointer to the element containing point p.
   */
  const Elem* find_element(const Point& p) const;

  /**
   * @returns a pointer to the element containing point p.
   */
  const Elem* operator() (const Point& p) const;
  
  
private:

  
  /**
   * The tree root.
   */
  TreeNode<N> root;

  /**
   * How the tree is built.
   */
  const BuildType build_type;
  
};



/**
 * For convenience we define QuadTrees and OctTrees 
 * explicitly.
 */
namespace Trees
{
  /**
   * A QuadTree is a tree appropriate
   * for 2D meshes.
   */
  typedef Tree<4> QuadTree;

  /**
   * An OctTree is a tree appropriate
   * for 3D meshes.
   */
  typedef Tree<8> OctTree;
}



// ------------------------------------------------------------
// Tree class inline methods

// constructor
template <unsigned int N>
inline
Tree<N>::Tree (const MeshBase& m,
	       const unsigned int level,
	       const BuildType bt) :
  TreeBase(m),
  root(m,level),
  build_type(bt)
{
  // Set the root node bounding box equal to the bounding
  // box for the entire domain.
  root.set_bounding_box (mesh.bounding_box());


  if (build_type == NODES)
    {
      // Add all the nodes to the root node.  It will 
      // automagically build the tree for us.
      const_node_iterator       it (mesh.const_nodes_begin());
      const const_node_iterator end(mesh.const_nodes_end());
      
      for (; it != end; ++it)
	root.insert (*it);
      
      // Now the tree contains the nodes.
      // However, we want element pointers, so here we
      // convert between the two.
      std::vector<std::vector<const Elem*> > nodes_to_elem;
      
      mesh.build_nodes_to_elem_map (nodes_to_elem);      
      root.transform_nodes_to_elements (nodes_to_elem);
    }

  else if (build_type == ELEMENTS)
    {
      // Add all the elements to the root node.  It will
      // automagically build the tree for us.
      const_active_elem_iterator       it (mesh.const_elements_begin());
      const const_active_elem_iterator end(mesh.const_elements_end());

      for (; it != end; ++it)
	root.insert (*it);
    }
}



// copy-constructor
template <unsigned int N>
inline
Tree<N>::Tree (const Tree<N>& other_tree) :
  TreeBase   (other_tree),
  root       (other_tree.root),
  build_type (other_tree.build_type)
{
  error();
}



template <unsigned int N>
inline
const Elem* Tree<N>::operator() (const Point& p) const
{
  return this->find_element(p);
}



#endif

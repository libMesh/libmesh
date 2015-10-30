// $Id: tree_node.h,v 1.2 2004-01-03 15:37:42 benkirk Exp $

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



#ifndef __tree_node_h__
#define __tree_node_h__

// C++ includes
#include <set>

// Local includes
#include "libmesh_common.h"
#include "mesh_base.h"



/**
 * This class defines a node on a tree.  A tree node
 * contains a pointer to its parent (NULL if the node is
 * the root) and pointers to its children (NULL if the 
 * node is active.
 */
template <unsigned int N>
class TreeNode
{  
public:
  

  /**
   * Dummy constructor.
   */
  //TreeNode () {}


  /**
   * Constructor.  Takes a pointer to this node's 
   * parent.  The pointer should only be NULL 
   * for the top-level (root) node.
   */
  TreeNode (const MeshBase& m, 
	    const unsigned int l,	    
	    const TreeNode<N> *p = NULL);

  /**
   * Destructor.  Deletes all children, if any.  Thus
   * to delete a tree it is sufficient to explicitly 
   * delete the root node. 
   */
  ~TreeNode ();

  /**
   * @returns true if this node is the root node, false
   * otherwise.
   */
  bool is_root() const { return (parent == NULL); }

  /**
   * @returns true if this node is active (i.e. has no 
   * children), false otherwise.
   */
  bool active() const { return children.empty(); }

  /**
   * Inserts node number n into the TreeNode.
   */
  void insert (const unsigned int n);

  /**
   * Refine the tree node into N children if it contains
   * more than tol nodes.
   */
  void refine (); 

  /**
   * Sets the bounding box;
   */
  void set_bounding_box (const std::pair<Point, Point>& bbox);

  /**
   * @returns true if this TreeNode (or its children) contain node n,
   * false otherwise.
   */
  bool bounds_node (const unsigned int n) const { assert (n < mesh.n_nodes()); 
  return bounds_point(mesh.point(n)); }

  /**
   * @returns true if this TreeNode (or its children) contain point p,
   * false otherwise.
   */
  bool bounds_point (const Point &p) const;

  /**
   * @returns the level of the node.
   */
  unsigned int level () const; 
    
  /**
   * Prints the contents of the node_numbers vector if we
   * are active.
   */
  void print_nodes() const;
    
  /**
   * Prints the contents of the elements set if we
   * are active.
   */
  void print_elements() const;

  /**
   * Transforms node numbers to element pointers.
   */
  void transform_nodes_to_elements (std::vector<std::vector<unsigned int> >& 
				    nodes_to_elem);

  /**
   * @returns 1 if the node is active, otherwise
   * recursively calls this on the children.
   */
  unsigned int n_active_bins() const;

  /**
   * @returns an element containing point p.
   */
  Elem* find_element(const Point& p) const;
  
private:
  
  /**
   *
   */
  Elem* find_element_in_children(const Point& p) const;

  /**
   * Constructs the bounding box for child c.
   */
  std::pair<Point, Point> create_bounding_box (const unsigned int c) const;

  /**
   * Reference to the mesh.
   */
  const MeshBase& mesh;

  /**
   * The maximum number of things we should store before 
   * refining ourself.
   */
  const unsigned int max_level;

  /**
   * Pointer to this node's parent.
   */
  const TreeNode<N> *parent;

  /**
   * Pointers to our children.  This vector
   * is empty if the node is active.
   */
  std::vector<TreeNode<N>* > children;

  /**
   * The Cartesian bounding box for the node. 
   * The minimum point is stored as bounding_box.first,
   * the maximum point is stored as bounding_box.second.
   */
  std::pair<Point, Point> bounding_box;

  /**
   * Pointers to the elements in this tree node.
   */
  std::vector<Elem*> elements;

  /**
   * The node numbers contained in this portion of the tree.
   */
  std::vector<unsigned int> node_numbers;

};





// ------------------------------------------------------------
// TreeNode class inline methods

// constructor
template <unsigned int N>
inline
TreeNode<N>::TreeNode (const MeshBase& m, 
		       const unsigned int l,
		       const TreeNode<N>* p) :
  mesh(m),
  max_level(l),
  parent(p)
{
  // Assert our children are empty, thus we are active.
  assert (children.empty());
}



// destructor
template <unsigned int N>
inline
TreeNode<N>::~TreeNode ()
{
  // When we are destructed we must delete all of our
  // children.  They will this delete their children,
  // All the way down the line...
  for (unsigned int c=0; c<children.size(); c++)
    delete children[c];

  children.clear();
  elements.clear();
}



// level
template <unsigned int N>
inline
unsigned int TreeNode<N>::level () const 
{ 
  if (parent == NULL)
    return 0;
  else
    return parent->level()+1;

  error();
  return libMesh::invalid_uint;
}



#endif

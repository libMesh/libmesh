// The libMesh Finite Element Library.
// Copyright (C) 2002-2022 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_TREE_NODE_H
#define LIBMESH_TREE_NODE_H

// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/bounding_box.h"
#include "libmesh/point.h"

// C++ includes
#include <cstddef>
#include <set>
#include <unordered_map>
#include <vector>

namespace libMesh
{

// Forward Declarations
class MeshBase;
class Node;
class Elem;

/**
 * This class defines a node on a tree.  A tree node
 * contains a pointer to its parent (nullptr if the node is
 * the root) and pointers to its children (nullptr if the
 * node is active.
 *
 * \author Daniel Dreyer
 * \date 2003
 * \brief Base class for different Tree types.
 */
template <unsigned int N>
class TreeNode
{
public:
  /**
   * Constructor.  Takes a pointer to this node's
   * parent.  The pointer should only be nullptr
   * for the top-level (root) node.
   */
  TreeNode (const MeshBase & m,
            unsigned int tbs,
            const TreeNode<N> * p = nullptr);

  /**
   * Destructor.  Deletes all children, if any.  Thus
   * to delete a tree it is sufficient to explicitly
   * delete the root node.
   */
  ~TreeNode ();

  /**
   * \returns \p true if this node is the root node, false
   * otherwise.
   */
  bool is_root() const { return (parent == nullptr); }

  /**
   * \returns \p true if this node is active (i.e. has no
   * children), false otherwise.
   */
  bool active() const { return children.empty(); }

  /**
   * Tries to insert \p Node \p nd into the TreeNode.
   * \returns \p true iff \p nd is inserted into the TreeNode or one of
   * its children.
   */
  bool insert (const Node * nd);

  /**
   * Inserts \p Elem \p el into the TreeNode.
   * \returns \p true iff \p el is inserted into the TreeNode or one of
   * its children.
   */
  bool insert (const Elem * nd);

  /**
   * Refine the tree node into N children if it contains
   * more than tol nodes.
   */
  void refine ();

  /**
   * Sets the bounding box;
   */
  void set_bounding_box (const std::pair<Point, Point> & bbox);

  /**
   * \returns \p true if this TreeNode (or its children) contain node n
   * (within relative tolerance), false otherwise.
   */
  bool bounds_node (const Node * nd,
                    Real relative_tol = 0) const;

  /**
   * \returns \p true if this TreeNode (or its children) contain point p
   * (within relative tolerance), false otherwise.
   */
  bool bounds_point (const Point & p,
                     Real relative_tol = 0) const;

  /**
   * \returns The level of the node.
   */
  unsigned int level () const;

  /**
   * Prints the contents of the node_numbers vector if we
   * are active.
   */
  void print_nodes(std::ostream & out_stream=libMesh::out) const;

  /**
   * Prints the contents of the elements set if we
   * are active.
   */
  void print_elements(std::ostream & out_stream=libMesh::out) const;

  /**
   * Transforms node numbers to element pointers.
   */
  void transform_nodes_to_elements (std::vector<std::vector<const Elem *>> & nodes_to_elem);

  /**
   * Transforms node numbers to element pointers.
   */
  void transform_nodes_to_elements (std::unordered_map<dof_id_type, std::vector<const Elem *>> & nodes_to_elem);

  /**
   * \returns The number of active bins below
   * (including) this element.
   */
  unsigned int n_active_bins() const;

  /**
   * \returns An element containing point p,
   * optionally restricted to a set of allowed subdomains.
   */
  const Elem * find_element (const Point & p,
                             const std::set<subdomain_id_type> * allowed_subdomains = nullptr,
                             Real relative_tol = TOLERANCE) const;

  /**
   * Fills \p candidate_elements with any elements containing the
   * specified point \p p,
   * optionally restricted to a set of allowed subdomains,
   * optionally using a non-default relative tolerance for searches.
   */
  void find_elements (const Point & p,
                      std::set<const Elem *> & candidate_elements,
                      const std::set<subdomain_id_type> * allowed_subdomains = nullptr,
                      Real relative_tol = TOLERANCE) const;

private:
  /**
   * Look for point \p p in our children,
   * optionally restricted to a set of allowed subdomains.
   */
  const Elem * find_element_in_children (const Point & p,
                                         const std::set<subdomain_id_type> * allowed_subdomains,
                                         Real relative_tol) const;

  /**
   * Look for points in our children,
   * optionally restricted to a set of allowed subdomains.
   */
  void find_elements_in_children (const Point & p,
                                  std::set<const Elem *> & candidate_elements,
                                  const std::set<subdomain_id_type> * allowed_subdomains,
                                  Real relative_tol) const;

  /**
   * Constructs the bounding box for child \p c.
   */
  BoundingBox create_bounding_box (unsigned int c) const;

  /**
   * Reference to the mesh.
   */
  const MeshBase & mesh;

  /**
   * Pointer to this node's parent.
   */
  const TreeNode<N> * parent;

  /**
   * Pointers to our children.  This vector
   * is empty if the node is active.
   */
  std::vector<TreeNode<N> * > children;

  /**
   * The Cartesian bounding box for the node.
   */
  BoundingBox bounding_box;

  /**
   * Pointers to the elements in this tree node.
   */
  std::vector<const Elem *> elements;

  /**
   * The node numbers contained in this portion of the tree.
   */
  std::vector<const Node *> nodes;

  /**
   * The maximum number of things we should store before
   * refining ourself.
   */
  const unsigned int tgt_bin_size;

  /**
   * This specifies the refinement level beyond which we will
   * scale up the target bin size in child TreeNodes. We set
   * the default to be 10, which should be large enough such
   * that in most cases the target bin size does not need to
   * be increased.
   */
  unsigned int target_bin_size_increase_level;

  /**
   * Does this node contain any infinite elements.
   */
  bool contains_ifems;
};





// ------------------------------------------------------------
// TreeNode class inline methods
template <unsigned int N>
inline
TreeNode<N>::TreeNode (const MeshBase & m,
                       unsigned int tbs,
                       const TreeNode<N> * p) :
  mesh           (m),
  parent         (p),
  tgt_bin_size   (tbs),
  target_bin_size_increase_level(10),
  contains_ifems (false)
{
  // libmesh_assert our children are empty, thus we are active.
  libmesh_assert (children.empty());
  libmesh_assert (this->active());

  // Reserve space for the nodes & elements
  nodes.reserve    (tgt_bin_size);
  elements.reserve (tgt_bin_size);
}



template <unsigned int N>
inline
TreeNode<N>::~TreeNode ()
{
  // When we are destructed we must delete all of our
  // children.  They will thus delete their children,
  // All the way down the line...
  for (auto c : children)
    delete c;
}



template <unsigned int N>
inline
unsigned int TreeNode<N>::level () const
{
  if (parent != nullptr)
    return parent->level()+1;

  // if we have no parent, we are a level-0 box
  return 0;
}


} // namespace libMesh


#endif // LIBMESH_TREE_NODE_H

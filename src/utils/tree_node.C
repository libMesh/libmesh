// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include <set>

// Local includes
#include "libmesh/libmesh_config.h"
#include "libmesh/tree_node.h"
#include "libmesh/mesh_base.h"
#include "libmesh/elem.h"

namespace libMesh
{

// ------------------------------------------------------------
// TreeNode class methods
template <unsigned int N>
bool TreeNode<N>::insert (const Node * nd)
{
  libmesh_assert(nd);
  libmesh_assert_less (nd->id(), mesh.n_nodes());

  // Return if we don't bound the node
  if (!this->bounds_node(nd))
    return false;

  // Add the node to ourself if we are active
  if (this->active())
    {
      nodes.push_back (nd);

      // Refine ourself if we reach the target bin size for a TreeNode.
      if (nodes.size() == tgt_bin_size)
        this->refine();

      return true;
    }

  // If we are not active simply pass the node along to
  // our children
  libmesh_assert_equal_to (children.size(), N);

  bool was_inserted = false;
  for (unsigned int c=0; c<N; c++)
    if (children[c]->insert (nd))
      was_inserted = true;
  return was_inserted;
}



template <unsigned int N>
bool TreeNode<N>::insert (const Elem * elem)
{
  libmesh_assert(elem);

  // We first want to find the corners of the cuboid surrounding the cell.
  const BoundingBox bbox = elem->loose_bounding_box();

  // Next, find out whether this cuboid has got non-empty intersection
  // with the bounding box of the current tree node.
  //
  // If not, we should not care about this element.
  if (!this->bounding_box.intersect(bbox))
    return false;

  // Only add the element if we are active
  if (this->active())
    {
      elements.push_back (elem);

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

      // flag indicating this node contains
      // infinite elements
      if (elem->infinite())
        this->contains_ifems = true;

#endif

      unsigned int element_count = elements.size();
      if(!mesh.get_count_lower_dim_elems_in_point_locator())
      {
        const std::set<unsigned char> & elem_dimensions = mesh.elem_dimensions();
        if(elem_dimensions.size() > 1)
        {
          element_count = 0;
          unsigned char highest_dim_elem = *elem_dimensions.rbegin();
          for (std::size_t i=0; i<elements.size(); i++)
          {
            if(elements[i]->dim() == highest_dim_elem)
            {
              element_count++;
            }
          }
        }
      }

      // Refine ourself if we reach the target bin size for a TreeNode.
      if (element_count == tgt_bin_size)
        this->refine();

      return true;
    }

  // If we are not active simply pass the element along to
  // our children
  libmesh_assert_equal_to (children.size(), N);

  bool was_inserted = false;
  for (unsigned int c=0; c<N; c++)
    if (children[c]->insert (elem))
      was_inserted = true;
  return was_inserted;
}



template <unsigned int N>
void TreeNode<N>::refine ()
{
  // Huh?  better be active...
  libmesh_assert (this->active());
  libmesh_assert (children.empty());

  // A TreeNode<N> has by definition N children
  children.resize(N);

  // Scale up the target bin size in child TreeNodes if we have reached
  // the maximum number of refinement levels.
  unsigned int new_target_bin_size = tgt_bin_size;
  if(level() >= target_bin_size_increase_level)
    {
      new_target_bin_size *= 2;
    }

  for (unsigned int c=0; c<N; c++)
    {
      // Create the child and set its bounding box.
      children[c] = new TreeNode<N> (mesh, new_target_bin_size, this);
      children[c]->set_bounding_box(this->create_bounding_box(c));

      // Pass off our nodes to our children
      for (std::size_t n=0; n<nodes.size(); n++)
        children[c]->insert(nodes[n]);

      // Pass off our elements to our children
      for (std::size_t e=0; e<elements.size(); e++)
        children[c]->insert(elements[e]);
    }

  // We don't need to store nodes or elements any more, they have been
  // added to the children.  Use the "swap trick" to actually reduce
  // the capacity of these vectors.
  std::vector<const Node *>().swap(nodes);
  std::vector<const Elem *>().swap(elements);

  libmesh_assert_equal_to (nodes.capacity(), 0);
  libmesh_assert_equal_to (elements.capacity(), 0);
}



template <unsigned int N>
void TreeNode<N>::set_bounding_box (const std::pair<Point, Point> & bbox)
{
  bounding_box = bbox;
}



template <unsigned int N>
bool TreeNode<N>::bounds_node (const Node * nd,
                               Real relative_tol) const
{
  libmesh_assert(nd);
  return bounds_point(*nd, relative_tol);
}



template <unsigned int N>
bool TreeNode<N>::bounds_point (const Point & p,
                                Real relative_tol) const
{
  const Point & min = bounding_box.first;
  const Point & max = bounding_box.second;

  const Real tol = (max - min).norm() * relative_tol;

  if ((p(0) >= min(0) - tol)
      && (p(0) <= max(0) + tol)
#if LIBMESH_DIM > 1
      && (p(1) >= min(1) - tol)
      && (p(1) <= max(1) + tol)
#endif
#if LIBMESH_DIM > 2
      && (p(2) >= min(2) - tol)
      && (p(2) <= max(2) + tol)
#endif
      )
    return true;

  return false;
}



template <unsigned int N>
BoundingBox
TreeNode<N>::create_bounding_box (unsigned int c) const
{
  switch (N)
    {
      // How to refine an OctTree Node
    case 8:
      {
        const Real xmin = bounding_box.first(0);
        const Real ymin = bounding_box.first(1);
        const Real zmin = bounding_box.first(2);

        const Real xmax = bounding_box.second(0);
        const Real ymax = bounding_box.second(1);
        const Real zmax = bounding_box.second(2);

        const Real xc = .5*(xmin + xmax);
        const Real yc = .5*(ymin + ymax);
        const Real zc = .5*(zmin + zmax);

        switch (c)
          {
          case 0:
            return BoundingBox (Point(xmin, ymin, zmin),
                                Point(xc,   yc,   zc));
          case 1:
            return BoundingBox (Point(xc,   ymin, zmin),
                                Point(xmax, yc,   zc));
          case 2:
            return BoundingBox (Point(xmin, yc,   zmin),
                                Point(xc,   ymax, zc));
          case 3:
            return BoundingBox (Point(xc,   yc,   zmin),
                                Point(xmax, ymax, zc));
          case 4:
            return BoundingBox (Point(xmin, ymin, zc),
                                Point(xc,   yc,   zmax));
          case 5:
            return BoundingBox (Point(xc,   ymin, zc),
                                Point(xmax, yc,   zmax));
          case 6:
            return BoundingBox (Point(xmin, yc,   zc),
                                Point(xc,   ymax, zmax));
          case 7:
            return BoundingBox (Point(xc,   yc,   zc),
                                Point(xmax, ymax, zmax));
          default:
            libmesh_error_msg("c >= N! : " << c);
          }

        break;
      } // case 8

      // How to refine an QuadTree Node
    case 4:
      {
        const Real xmin = bounding_box.first(0);
        const Real ymin = bounding_box.first(1);

        const Real xmax = bounding_box.second(0);
        const Real ymax = bounding_box.second(1);

        const Real xc = .5*(xmin + xmax);
        const Real yc = .5*(ymin + ymax);

        switch (c)
          {
          case 0:
            return BoundingBox (Point(xmin, ymin),
                                Point(xc,   yc));
          case 1:
            return BoundingBox (Point(xc,   ymin),
                                Point(xmax, yc));
          case 2:
            return BoundingBox (Point(xmin, yc),
                                Point(xc,   ymax));
          case 3:
            return BoundingBox (Point(xc,   yc),
                                Point(xmax, ymax));
          default:
            libmesh_error_msg("c >= N!");
          }

        break;
      } // case 4

      // How to refine a BinaryTree Node
    case 2:
      {
        const Real xmin = bounding_box.first(0);

        const Real xmax = bounding_box.second(0);

        const Real xc = .5*(xmin + xmax);

        switch (c)
          {
          case 0:
            return BoundingBox (Point(xmin),
                                Point(xc));
          case 1:
            return BoundingBox (Point(xc),
                                Point(xmax));
          default:
            libmesh_error_msg("c >= N!");
          }

        break;
      } // case 2

    default:
      libmesh_error_msg("Only implemented for Octrees, QuadTrees, and Binary Trees!");
    }

  libmesh_error_msg("We'll never get here!");
  Point min, max;
  return BoundingBox (min, max);
}



template <unsigned int N>
void TreeNode<N>::print_nodes(std::ostream & out_stream) const
{
  if (this->active())
    {
      out_stream << "TreeNode Level: " << this->level() << std::endl;

      for (std::size_t n=0; n<nodes.size(); n++)
        out_stream << " " << nodes[n]->id();

      out_stream << std::endl << std::endl;
    }
  else
    {
      for (std::size_t child=0; child<children.size(); child++)
        children[child]->print_nodes();
    }
}



template <unsigned int N>
void TreeNode<N>::print_elements(std::ostream & out_stream) const
{
  if (this->active())
    {
      out_stream << "TreeNode Level: " << this->level() << std::endl;

      for (std::vector<const Elem *>::const_iterator pos=elements.begin();
           pos != elements.end(); ++pos)
        out_stream << " " << *pos;

      out_stream << std::endl << std::endl;
    }
  else
    {
      for (std::size_t child=0; child<children.size(); child++)
        children[child]->print_elements();
    }
}



template <unsigned int N>
void TreeNode<N>::transform_nodes_to_elements (std::vector<std::vector<const Elem *> > & nodes_to_elem)
{
  if (this->active())
    {
      elements.clear();

      // Temporarily use a set. Since multiple nodes
      // will likely map to the same element we use a
      // set to eliminate the duplication.
      std::set<const Elem *> elements_set;

      for (std::size_t n=0; n<nodes.size(); n++)
        {
          // the actual global node number we are replacing
          // with the connected elements
          const dof_id_type node_number = nodes[n]->id();

          libmesh_assert_less (node_number, mesh.n_nodes());
          libmesh_assert_less (node_number, nodes_to_elem.size());

          for (std::size_t e=0; e<nodes_to_elem[node_number].size(); e++)
            elements_set.insert(nodes_to_elem[node_number][e]);
        }

      // Done with the nodes.
      std::vector<const Node *>().swap(nodes);

      // Now the set is built.  We can copy this to the
      // vector.  Note that the resulting vector will
      // already be sorted, and will require less memory
      // than the set.
      elements.reserve(elements_set.size());

      for (std::set<const Elem *>::iterator pos=elements_set.begin();
           pos != elements_set.end(); ++pos)
        {
          elements.push_back(*pos);

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

          // flag indicating this node contains
          // infinite elements
          if ((*pos)->infinite())
            this->contains_ifems = true;

#endif
        }
    }
  else
    {
      for (std::size_t child=0; child<children.size(); child++)
        children[child]->transform_nodes_to_elements (nodes_to_elem);
    }

}



template <unsigned int N>
unsigned int TreeNode<N>::n_active_bins() const
{
  if (this->active())
    return 1;

  else
    {
      unsigned int sum=0;

      for (std::size_t c=0; c<children.size(); c++)
        sum += children[c]->n_active_bins();

      return sum;
    }
}



template <unsigned int N>
const Elem *
TreeNode<N>::find_element (const Point & p,
                           const std::set<subdomain_id_type> * allowed_subdomains,
                           Real relative_tol) const
{
  if (this->active())
    {
      // Only check our children if the point is in our bounding box
      // or if the node contains infinite elements
      if (this->bounds_point(p, relative_tol) || this->contains_ifems)
        // Search the active elements in the active TreeNode.
        for (std::vector<const Elem *>::const_iterator pos=elements.begin();
             pos != elements.end(); ++pos)
          if (!allowed_subdomains || allowed_subdomains->count((*pos)->subdomain_id()))
            if ((*pos)->active() && (*pos)->contains_point(p, relative_tol))
              return *pos;

      // The point was not found in any element
      return libmesh_nullptr;
    }
  else
    return this->find_element_in_children(p,allowed_subdomains,
                                          relative_tol);

  libmesh_error_msg("We'll never get here!");
  return libmesh_nullptr;
}




template <unsigned int N>
const Elem * TreeNode<N>::find_element_in_children (const Point & p,
                                                    const std::set<subdomain_id_type> * allowed_subdomains,
                                                    Real relative_tol) const
{
  libmesh_assert (!this->active());

  searched_child.assign(children.size(), false);

  // First only look in the children whose bounding box
  // contain the point p.
  for (std::size_t c=0; c<children.size(); c++)
    if (children[c]->bounds_point(p, relative_tol))
      {
        const Elem * e =
          children[c]->find_element(p,allowed_subdomains,
                                    relative_tol);

        if (e != libmesh_nullptr)
          return e;

        // If we get here then a child that bounds the
        // point does not have any elements that contain
        // the point.  So, we will search all our children.
        // However, we have already searched child c so there
        // is no use searching her again.
        searched_child[c] = true;
      }


  // If we get here then our child whose bounding box
  // was searched and did not find any elements containing
  // the point p.  So, let's look at the other children
  // but exclude the one we have already searched.
  for (std::size_t c=0; c<children.size(); c++)
    if (!searched_child[c])
      {
        const Elem * e =
          children[c]->find_element(p,allowed_subdomains,
                                    relative_tol);

        if (e != libmesh_nullptr)
          return e;
      }

  // If we get here we have searched all our children.
  // Since this process was started at the root node then
  // we have searched all the elements in the tree without
  // success.  So, we should return NULL since at this point
  // _no_ elements in the tree claim to contain point p.

  return libmesh_nullptr;
}



// ------------------------------------------------------------
// Explicit Instantiations
template class TreeNode<2>;
template class TreeNode<4>;
template class TreeNode<8>;

} // namespace libMesh

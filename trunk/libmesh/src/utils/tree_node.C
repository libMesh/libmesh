// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2007  Benjamin S. Kirk, John W. Peterson
  
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
#include "libmesh_config.h"
#include "tree_node.h"
#include "mesh_base.h"
#include "elem.h"

// ------------------------------------------------------------
// TreeNode class methods
template <unsigned int N>
void TreeNode<N>::insert (const Node* nd)
{
  assert (nd != NULL);
  assert (nd->id() < mesh.n_nodes());

  // Return if we don't bound the node
  if (!this->bounds_node(nd))
    return;
  
  // Only add the node if we are active
  if (this->active())
    {
      nodes.push_back (nd);

      // Refine ourself if we reach the target bin size for a TreeNode.
      if (nodes.size() == tgt_bin_size)
	this->refine();
    }
  
  // If we are not active simply pass the node along to
  // our children
  else
    {
      assert (children.size() == N);
      
      for (unsigned int c=0; c<N; c++)
	children[c]->insert (nd);
    }
}



template <unsigned int N>
void TreeNode<N>::insert (const Elem* elem)
{
  assert (elem != NULL);

  /* We first want to find the corners of the cuboid surrounding the
     cell.  */
  Point minCoord = elem->point(0);
  Point maxCoord = minCoord;
  unsigned int dim = elem->dim();
  for(unsigned int i=elem->n_nodes()-1; i>0; i--)
    {
      Point p = elem->point(i);
      for(unsigned int d=0; d<dim; d++)
	{
	  if(minCoord(d)>p(d)) minCoord(d) = p(d);
	  if(maxCoord(d)<p(d)) maxCoord(d) = p(d);
	}
    }

  /* Next, find out whether this cuboid has got non-empty intersection
     with the bounding box of the current tree node.  */
  bool intersects = true;
  for(unsigned int d=0; d<dim; d++)
    {
      if(maxCoord(d)<this->bounding_box.first(d) ||
	 minCoord(d)>this->bounding_box.second(d))
	{
	  intersects = false;
	}
    }

  /* If not, we should not care about this element.  */
  if(!intersects)
    {
      return;
    }

  // Only add the element if we are active
  if (this->active())
    {
      elements.push_back (elem);

#ifdef ENABLE_INFINITE_ELEMENTS

      // flag indicating this node contains
      // infinite elements
      if (elem->infinite())	
	this->contains_ifems = true;
      
#endif
      
      // Refine ourself if we reach the target bin size for a TreeNode.
      if (elements.size() == tgt_bin_size)
	this->refine();
    }
  
  // If we are not active simply pass the element along to
  // our children
  else
    {
      assert (children.size() == N);
      
      for (unsigned int c=0; c<N; c++)
	children[c]->insert (elem);
    }
}



template <unsigned int N>
void TreeNode<N>::refine ()
{
  // Huh?  better be active...
  assert (this->active());
  assert (children.empty());
  
  // A TreeNode<N> has by definition N children
  children.resize(N);

  for (unsigned int c=0; c<N; c++)
    {
      // Create the child and set its bounding box.
      children[c] = new TreeNode<N> (mesh, tgt_bin_size, this);
      children[c]->set_bounding_box(this->create_bounding_box(c));

      // Pass off our nodes to our children
      for (unsigned int n=0; n<nodes.size(); n++)
	children[c]->insert(nodes[n]);

      // Pass off our elements to our children
      for (unsigned int e=0; e<elements.size(); e++)
	children[c]->insert(elements[e]);
    }

  // We don't need to store nodes or elements any more,
  // they have been added to the children.
  // Note that we cannot use std::vector<>::clear() here
  // since that in general does not reduce capacity!!
  // That would be a bad, bad thing.
  std::vector<const Node*>().swap(nodes);
  std::vector<const Elem*>().swap(elements);

  assert (nodes.capacity()    == 0);
  assert (elements.capacity() == 0);
}



template <unsigned int N>
void TreeNode<N>::set_bounding_box (const std::pair<Point, Point>& bbox)
{
  bounding_box = bbox;
}



template <unsigned int N>
bool TreeNode<N>::bounds_point (const Point& p) const
{
  const Point& min = bounding_box.first;
  const Point& max = bounding_box.second;


  if ((p(0) >= min(0)) &&
      (p(1) >= min(1)) &&
      (p(2) >= min(2)) &&
      
      (p(0) <= max(0)) &&
      (p(1) <= max(1)) &&
      (p(2) <= max(2)))
    return true;   

  return false;
}



template <unsigned int N>
std::pair<Point, Point> 
TreeNode<N>::create_bounding_box (const unsigned int c) const
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

	const Real xc = xmin + .5*(xmax - xmin);
	const Real yc = ymin + .5*(ymax - ymin);
	const Real zc = zmin + .5*(zmax - zmin);

	
	switch (c)
	  {

	  case 0:
	    {
	      const Point min(xmin, ymin, zmin);
	      const Point max(xc,   yc,   zc);
	      return std::make_pair (min, max);
	      break;
	    }

	  case 1:
	    {
	      const Point min(xc,   ymin, zmin);
	      const Point max(xmax, yc,   zc);
	      return std::make_pair (min, max);
	      break;
	    }

	  case 2:
	    {
	      const Point min(xmin, yc,   zmin);
	      const Point max(xc,   ymax, zc);
	      return std::make_pair (min, max);
	      break;
	    }

	  case 3:
	    {
	      const Point min(xc,   yc,   zmin);
	      const Point max(xmax, ymax, zc);
	      return std::make_pair (min, max);
	      break;
	    }

	  case 4:
	    {
	      const Point min(xmin, ymin, zc);
	      const Point max(xc,   yc,   zmax);
	      return std::make_pair (min, max);
	      break;
	    }

	  case 5:
	    {
	      const Point min(xc,   ymin, zc);
	      const Point max(xmax, yc,   zmax);
	      return std::make_pair (min, max);
	      break;
	    }

	  case 6:
	    {
	      const Point min(xmin, yc,   zc);
	      const Point max(xc,   ymax, zmax);
	      return std::make_pair (min, max);
	      break;
	    }

	  case 7:
	    {
	      const Point min(xc,   yc,   zc);
	      const Point max(xmax, ymax, zmax);
	      return std::make_pair (min, max);
	      break;
	    }
    
	  default:
	    std::cerr << "c >= N! : " << c
		      << std::endl;
	    libmesh_error();
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

	const Real xc = xmin + .5*(xmax - xmin);
	const Real yc = ymin + .5*(ymax - ymin);

	switch (c)
	  {
	  case 0:
	    {
	      const Point min(xmin, ymin);
	      const Point max(xc,   yc);
	      return std::make_pair (min, max);
	      break;
	    }

	  case 1:
	    {
	      const Point min(xc,   ymin);
	      const Point max(xmax, yc);
	      return std::make_pair (min, max);
	      break;
	    }

	  case 2:
	    {
	      const Point min(xmin, yc);
	      const Point max(xc,   ymax);
	      return std::make_pair (min, max);
	      break;
	    }

	  case 3:
	    {
	      const Point min(xc,   yc);
	      const Point max(xmax, ymax);
	      return std::make_pair (min, max);
	      break;
	    }
	    
	  default:
	    std::cerr << "c >= N!" << std::endl;
	    libmesh_error();
	    
	  }

	break;
      } // case 4

    default:
      std::cerr << "Only implemented for Octrees and QuadTrees!" << std::endl;
      libmesh_error();

    }

  // How did we get here?
  libmesh_error();

  Point min, max;  
  return std::make_pair (min, max);
}



template <unsigned int N>
void TreeNode<N>::print_nodes() const
{
  if (this->active())
    {
      std::cout << "TreeNode Level: " << this->level() << std::endl;
      
      for (unsigned int n=0; n<nodes.size(); n++)
	std::cout << " " << nodes[n]->id();
      
      std::cout << std::endl << std::endl;
	
    }
  else
    {
      for (unsigned int child=0; child<children.size(); child++)
	children[child]->print_nodes();
    }
}



template <unsigned int N>
void TreeNode<N>::print_elements() const
{
  if (this->active())
    {
      std::cout << "TreeNode Level: " << this->level() << std::endl;
     
      for (std::vector<const Elem*>::const_iterator pos=elements.begin();
	   pos != elements.end(); ++pos)
	std::cout << " " << *pos;

      std::cout << std::endl << std::endl;	
    }
  else
    {
      for (unsigned int child=0; child<children.size(); child++)
	children[child]->print_elements();
    }
}



template <unsigned int N>
void TreeNode<N>::transform_nodes_to_elements (std::vector<std::vector<const Elem*> >& nodes_to_elem)
{
   if (this->active())
    {
      elements.clear();

      // Temporarily use a set. Since multiple nodes
      // will likely map to the same element we use a
      // set to eliminate the duplication.
      std::set<const Elem*> elements_set;
      
      for (unsigned int n=0; n<nodes.size(); n++)
	{
	  // the actual global node number we are replacing
	  // with the connected elements
	  const unsigned int node_number = nodes[n]->id();

	  assert (node_number < mesh.n_nodes());
	  assert (node_number < nodes_to_elem.size());
	  
	  for (unsigned int e=0; e<nodes_to_elem[node_number].size(); e++)
	    elements_set.insert(nodes_to_elem[node_number][e]);
	}

      // Done with the nodes.
      std::vector<const Node*>().swap(nodes);

      // Now the set is built.  We can copy this to the
      // vector.  Note that the resulting vector will 
      // already be sorted, and will require less memory
      // than the set.
      elements.reserve(elements_set.size());

      for (std::set<const Elem*>::iterator pos=elements_set.begin();
	   pos != elements_set.end(); ++pos)
	{
	  elements.push_back(*pos);
	  
#ifdef ENABLE_INFINITE_ELEMENTS

	  // flag indicating this node contains
	  // infinite elements
	  if ((*pos)->infinite())
	    this->contains_ifems = true;
	  
#endif
	}
    }
  else
    {
      for (unsigned int child=0; child<children.size(); child++)
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

      for (unsigned int c=0; c<children.size(); c++)
	sum += children[c]->n_active_bins();

      return sum;
    }
}



template <unsigned int N>
const Elem* TreeNode<N>::find_element(const Point& p) const
{
  if (this->active())
    {
      // Only check our children if the point is in our bounding box
      // or if the node contains infinite elements
      if (this->bounds_point(p) || this->contains_ifems)
	// Search the active elements in the active TreeNode.
	for (std::vector<const Elem*>::const_iterator pos=elements.begin();
	     pos != elements.end(); ++pos)
	  if ((*pos)->active())
	    if ((*pos)->contains_point(p))
	      return *pos;
      
      // The point was not found in any element
      return NULL;	    
    }
  else
    return this->find_element_in_children(p);
  
    

  // Should never get here.  See if-else structure
  // above with return statements that must get executed.
  libmesh_error();
  
  return NULL;
}




template <unsigned int N>
const Elem* TreeNode<N>::find_element_in_children(const Point& p) const
{
  assert (!this->active());

  unsigned int excluded_child = libMesh::invalid_uint;
  
  // First only look in the children whose bounding box
  // contain the point p.  Note that only one child will
  // bound the point since the bounding boxes are not
  // overlapping
  for (unsigned int c=0; c<children.size(); c++)
    if (children[c]->bounds_point(p))
      {
	if (children[c]->active())
	  {
	    const Elem* e = children[c]->find_element(p);
	    
	    if (e != NULL)
	      return e;	  
	  }
	else
	  {
	    const Elem* e = children[c]->find_element_in_children(p);

	    if (e != NULL)
	      return e;
	  }

	// If we get here than the child that bounds the
	// point does not have any elements that contain
	// the point.  So, we will search all our children.
	// However, we have already searched child c so there
	// is no use searching her again.
	excluded_child = c;
      }
	 

  // If we get here then our child whose bounding box
  // was searched and did not find any elements containing
  // the point p.  So, let's look at the other children
  // but exclude the one we have already searched.
  for (unsigned int c=0; c<children.size(); c++)
    if (c != excluded_child)    
      if (children[c]->active())
	{
	  const Elem* e = children[c]->find_element(p);
	  
	  if (e != NULL)
	    return e;	  
	}
      else
	{
	  const Elem* e = children[c]->find_element_in_children(p);
	  
	  if (e != NULL)
	    return e;
	}

  // If we get here we have searched all our children.
  // Since this process was started at the root node then
  // we have searched all the elements in the tree without
  // success.  So, we should return NULL since at this point
  // _no_ elements in the tree claim to contain point p.
  
  return NULL;     
}



// ------------------------------------------------------------
// Explicit Instantiations
template class TreeNode<4>;
template class TreeNode<8>;








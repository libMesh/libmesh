// $Id: tree.C,v 1.1.1.1 2003-01-10 16:17:48 libmesh Exp $

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
#include <set>
#include <algorithm>

// Local includes
#include "tree.h"
#include "elem.h"



// ------------------------------------------------------------
// TreeNode class methods
template <unsigned int N>
void TreeNode<N>::insert (const unsigned int n)
{
  assert (n < mesh.n_nodes());
  
  // Only consider the node for addition
  // if we are active
  if (active())
    {
      // Only add the node if we bound it
      if (bounds_node(n))
	node_numbers.push_back (n);
    }
    
  // If we are not active simply pass the node along to
  // our children
  else
    {
      for (unsigned int c=0; c<children.size(); c++)
	children[c]->insert (n);
    }

  // Refine ourself if we exceed the maximim 
  // number of nodes in a TreeNode.
  if (node_numbers.size() > max_level)
    {
      refine();
    }
};



template <unsigned int N>
void TreeNode<N>::refine ()
{
  // A TreeNode<N> has by definition N children
  children.resize(N);

  for (unsigned int c=0; c<children.size(); c++)
    {
      // Create the child and set its bounding box.
      children[c] = new TreeNode<N> (mesh, max_level, this);
      children[c]->set_bounding_box(create_bounding_box(c));

      // Pass off our nodes to our children
      for (unsigned int node=0; node<node_numbers.size(); node++)
	if (children[c]->bounds_node(node_numbers[node]))
	  children[c]->insert(node_numbers[node]);
    }

  // We don't need to store any nodes any more
  node_numbers.clear(); 
  elements.clear();
};



template <unsigned int N>
void TreeNode<N>::set_bounding_box (const std::pair<Point, Point>& bbox)
{
  bounding_box = bbox;
};



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
      (p(2) <= max(2))
      )
    return true;
   

  return false;
};



template <unsigned int N>
std::pair<Point, Point> 
TreeNode<N>::create_bounding_box (const unsigned int c) const
{
  switch (N)
    {

      // How to refine an OctTree Node
    case 8:
      {
	const real xmin = bounding_box.first(0);
	const real ymin = bounding_box.first(1);
	const real zmin = bounding_box.first(2);

	const real xmax = bounding_box.second(0);
	const real ymax = bounding_box.second(1);
	const real zmax = bounding_box.second(2);

	const real xc = xmin + .5*(xmax - xmin);
	const real yc = ymin + .5*(ymax - ymin);
	const real zc = zmin + .5*(zmax - zmin);

	
	switch (c)
	  {

	  case 0:
	    {
	      const Point min(xmin, ymin, zmin);
	      const Point max(xc,   yc,   zc);
	      return std::pair<Point, Point> (min, max);
	      break;
	    }

	  case 1:
	    {
	      const Point min(xc,   ymin, zmin);
	      const Point max(xmax, yc,   zc);
	      return std::pair<Point, Point> (min, max);
	      break;
	    }

	  case 2:
	    {
	      const Point min(xmin, yc,   zmin);
	      const Point max(xc,   ymax, zc);
	      return std::pair<Point, Point> (min, max);
	      break;
	    }

	  case 3:
	    {
	      const Point min(xc,   yc,   zmin);
	      const Point max(xmax, ymax, zc);
	      return std::pair<Point, Point> (min, max);
	      break;
	    }

	  case 4:
	    {
	      const Point min(xmin, ymin, zc);
	      const Point max(xc,   yc,   zmax);
	      return std::pair<Point, Point> (min, max);
	      break;
	    }

	  case 5:
	    {
	      const Point min(xc,   ymin, zc);
	      const Point max(xmax, yc,   zmax);
	      return std::pair<Point, Point> (min, max);
	      break;
	    }

	  case 6:
	    {
	      const Point min(xmin, yc,   zc);
	      const Point max(xc,   ymax, zmax);
	      return std::pair<Point, Point> (min, max);
	      break;
	    }

	  case 7:
	    {
	      const Point min(xc,   yc,   zc);
	      const Point max(xmax, ymax, zmax);
	      return std::pair<Point, Point> (min, max);
	      break;
	    }

    
	  default:
	    std::cerr << "c >= N!" << std::endl;
	    error();
	  };



	break;
      }

    default:
      std::cerr << "Only implemented for Octrees!" << std::endl;
      error();

    };

  // How did we get here?
  error();

  Point min, max;
  
  return std::pair<Point, Point> (min, max);
};



template <unsigned int N>
void TreeNode<N>::print_nodes() const
{
  if (active())
    {
      std::cout << "TreeNode Level: " << level() << std::endl;
      
      for (unsigned int node=0; node<node_numbers.size(); node++)
	std::cout << " " << node_numbers[node];
      
      std::cout << std::endl << std::endl;
	
    }
  else
    {
      for (unsigned int child=0; child<children.size(); child++)
	children[child]->print_nodes();
    };
};



template <unsigned int N>
void TreeNode<N>::print_elements() const
{
  if (active())
    {
      std::cout << "TreeNode Level: " << level() << std::endl;
     
      for (std::vector<Elem*>::const_iterator pos=elements.begin();
	   pos != elements.end(); ++pos)
	std::cout << " " << *pos;

      std::cout << std::endl << std::endl;	
    }
  else
    {
      for (unsigned int child=0; child<children.size(); child++)
	children[child]->print_elements();
    };
};



template <unsigned int N>
void TreeNode<N>::transform_nodes_to_elements (std::vector<std::vector<unsigned int> >& nodes_to_elem)
{
   if (active())
    {
      elements.clear();

      // Temporarily use a set. Since multiple nodes
      // will likely map to the same element we use a
      // set to eliminate the duplication.
      std::set<Elem*> elements_set;
      
      for (unsigned int node=0; node<node_numbers.size(); node++)
	{
	  // the actual global node number we are replacing
	  // with the connected elements
	  const unsigned int node_number = node_numbers[node];

	  assert (node_number < mesh.n_nodes());
	  assert (node_number < nodes_to_elem.size());
	  
	  for (unsigned int e=0; e<nodes_to_elem[node_number].size(); e++)
	    elements_set.insert(mesh.elem(nodes_to_elem[node_number][e]));
	};

      // Done with the node numbers.
      node_numbers.clear();

      // Now the set is built.  We can copy this to the
      // vector.  Note that the resulting vector will 
      // already be sorted, and will require less memory
      // than the set.
      elements.resize(elements_set.size());

      unsigned int cnt=0;		      

      for (std::set<Elem*>::iterator pos=elements_set.begin();
	   pos != elements_set.end(); ++pos)
	elements[cnt++] = *pos;
    }
  else
    {
      for (unsigned int child=0; child<children.size(); child++)
	children[child]->transform_nodes_to_elements (nodes_to_elem);
    };
 
};



template <unsigned int N>
unsigned int TreeNode<N>::n_active_bins() const
{
  if (active())
    return 1;
  else
    {
      unsigned int sum=0;

      for (unsigned int c=0; c<children.size(); c++)
	sum += children[c]->n_active_bins();

      return sum;
    }
};



template <unsigned int N>
Elem* TreeNode<N>::find_element(const Point& p) const
{
  if (active())
    {
      // Search the active elements in the active TreeNode.
      for (std::vector<Elem*>::const_iterator pos=elements.begin();
	   pos != elements.end(); ++pos)
	if ((*pos)->active())
	  if ((*pos)->contains_point(mesh, p))
	    return *pos;

      // The point was not found in any element
      return NULL;	    
    }
  else
    {
      return find_element_in_children(p);
    }
    

  // Should never get here.  See if-else structure
  // above with return statements that must get executed.
  error();
  
  return NULL;
};




template <unsigned int N>
Elem* TreeNode<N>::find_element_in_children(const Point& p) const
{
  assert (!active());

  unsigned int excluded_child = static_cast<unsigned int>(-1);
  
  // First only look in the children whose bounding box
  // contain the point p.  Note that only one child will
  // bound the point since the bounding boxes are not
  // overlapping
  for (unsigned int c=0; c<children.size(); c++)
    if (children[c]->bounds_point(p))
      {
	if (children[c]->active())
	  {
	    Elem* e = children[c]->find_element(p);
	    
	    if (e != NULL)
	      return e;	  
	  }
	else
	  {
	    Elem* e = children[c]->find_element_in_children(p);

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
	  Elem* e = children[c]->find_element(p);
	  
	  if (e != NULL)
	    return e;	  
	}
      else
	{
	  Elem* e = children[c]->find_element_in_children(p);
	  
	  if (e != NULL)
	    return e;
	}

  // If we get here we have searched all our children.
  // Since this process was started at the root node then
  // we have searched all the elements in the tree without
  // success.  So, we should return NULL since at this point
  // _no_ elements in the tree claim to contain point p.
  
  return NULL;     
};



// ------------------------------------------------------------
// Tree class method
template <unsigned int N>
Elem* Tree<N>::find_element(const Point& p) const
{
  return root.find_element(p);
};


// ------------------------------------------------------------
// Explicit Instantiations
template class TreeNode<4>;
template class TreeNode<8>;
template class Tree<4>;
template class Tree<8>;










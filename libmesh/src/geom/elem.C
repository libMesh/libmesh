// $Id: elem.C,v 1.33 2004-03-24 05:49:12 jwpeterson Exp $

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
#include <algorithm> // for std::sort
#include <iterator>  // for std::ostream_iterator

// Local includes
#include "elem.h"
#include "fe_type.h"
#include "fe_interface.h"
#include "edge_edge2.h"
#include "edge_edge3.h"
#include "edge_inf_edge2.h"
#include "face_tri3.h"
#include "face_tri6.h"
#include "face_quad4.h"
#include "face_quad8.h"
#include "face_quad9.h"
#include "face_inf_quad.h"
#include "face_inf_quad4.h"
#include "face_inf_quad6.h"
#include "cell_tet.h"
#include "cell_tet4.h"
#include "cell_tet10.h"
#include "cell_hex.h"
#include "cell_hex8.h"
#include "cell_hex20.h"
#include "cell_hex27.h"
#include "cell_inf_hex.h"
#include "cell_inf_hex8.h"
#include "cell_inf_hex16.h"
#include "cell_inf_hex18.h"
#include "cell_prism.h"
#include "cell_prism6.h"
#include "cell_prism15.h"
#include "cell_prism18.h"
#include "cell_inf_prism.h"
#include "cell_inf_prism6.h"
#include "cell_inf_prism12.h"
#include "cell_pyramid.h"
#include "cell_pyramid5.h"



// ------------------------------------------------------------
// Elem class member funcions
Elem* Elem::build(const ElemType type,
		  const Elem* p)
{
  switch (type)
    {
      // 1D elements 
    case EDGE2:
      {
	return new Edge2(p);
      }
    case EDGE3:
      {
	return new Edge3(p);
      }


      
      // 2D elements
    case TRI3:
      {
	return new Tri3(p);
      }
    case TRI6:
      {
	return new Tri6(p);
      }
    case QUAD4:
      {
	return new Quad4(p);
      }
    case QUAD8:
      {
	return new Quad8(p);
      }
    case QUAD9:
      {
	return new Quad9(p);
      }
   

      // 3D elements
    case TET4:
      {
	return new Tet4(p);
      }
    case TET10:
      {
	return new Tet10(p);
      }
    case HEX8:
      {
	return new Hex8(p);
      }
    case HEX20:
      {
	return new Hex20(p);
      }
    case HEX27:
      {
	return new Hex27(p);
      }
    case PRISM6:
      {
	return new Prism6(p);
      }
    case PRISM15:
      {
	return new Prism15(p);
      }
    case PRISM18:
      {
	return new Prism18(p);
      }
    case PYRAMID5:
      {
	return new Pyramid5(p);
      }



#ifdef ENABLE_INFINITE_ELEMENTS

      // 1D infinite elements
    case INFEDGE2:
      {
	return new InfEdge2(p);
      }


      // 2D infinite elements
    case INFQUAD4:
      {
	return new InfQuad4(p);
      }
    case INFQUAD6:
      {
	return new InfQuad6(p);
      }

   
    // 3D infinite elements
    case INFHEX8:
      {
	return new InfHex8(p);
      }
    case INFHEX16:
      {
	return new InfHex16(p);
      }
    case INFHEX18:
      {
	return new InfHex18(p);
      }
    case INFPRISM6:
      {
	return new InfPrism6(p);
      }
    case INFPRISM12:
      {
	return new InfPrism12(p);
      }

#endif
           
    default:
      {
	std::cerr << "ERROR: Undefined element type!." << std::endl;
	error();
      }
    }
    

  // This will never happen...  Look at the case-structure above.
  // error() will abort, so we won't get here...
  error();
  
  return NULL;
}



void Elem::write_tecplot_connectivity(std::ostream& out) const
{
  assert (!out.bad());
  assert (_nodes != NULL);
  
  for (unsigned int sc=0; sc <this->n_sub_elem(); sc++)
    {
      std::vector<unsigned int> conn =
	this->tecplot_connectivity(sc);
      
      std::copy(conn.begin(),
 		conn.end(),
 		std::ostream_iterator<unsigned int>(out, " "));
      
      out << std::endl;
    }
}



unsigned int Elem::key() const
{
  const unsigned int nv = this->n_vertices();
  
  assert (nv != 0);

  std::vector<unsigned int> vec (nv, 0);

  for (unsigned int v=0; v<nv; v++)
    vec[v] = this->node(v);

  
  std::sort(vec.begin(), vec.end());
  
  unsigned int n       = vec[0];
  const unsigned int m = vec[0];
  
  for (unsigned int i=1; i<nv; i++)
    n = n^vec[i];
  

  return n + m*m;  
} 



Point Elem::centroid() const
{
  Point cp;

  for (unsigned int n=0; n<this->n_vertices(); n++)
    cp.add (this->point(n));

  return (cp /= static_cast<Real>(this->n_vertices()));    
}



Real Elem::hmin() const
{
  Real h_min=1.e30;

  for (unsigned int n_outer=0; n_outer<this->n_vertices(); n_outer++)
    for (unsigned int n_inner=n_outer+1; n_inner<this->n_vertices(); n_inner++)
      {
	const Point diff = (this->point(n_outer) - this->point(n_inner));
	
	h_min = std::min(h_min,diff.size());
      }

  return h_min;
}



Real Elem::hmax() const
{
  Real h_max=0;

  for (unsigned int n_outer=0; n_outer<this->n_vertices(); n_outer++)
    for (unsigned int n_inner=n_outer+1; n_inner<this->n_vertices(); n_inner++)
      {
	const Point diff = (this->point(n_outer) - this->point(n_inner));
	
	h_max = std::max(h_max,diff.size());
      }

  return h_max;
}



Real Elem::length(const unsigned int n1, 
		  const unsigned int n2) const
{
  assert ( n1 < this->n_vertices() );
  assert ( n2 < this->n_vertices() );

  return (this->point(n1) - this->point(n2)).size();
}



bool Elem::operator == (const Elem& rhs) const
{
//   assert (n_nodes());
//   assert (rhs.n_nodes());

//   // Elements can only be equal if they
//   // contain the same number of nodes.
//   if (this->n_nodes() == rhs.n_nodes())
//     {
//       // Create a set that contains our global
//       // node numbers and those of our neighbor.
//       // If the set is the same size as the number
//       // of nodes in both elements then they must
//       // be connected to the same nodes.     
//       std::set<unsigned int> nodes_set;

//       for (unsigned int n=0; n<this->n_nodes(); n++)
//         {
//           nodes_set.insert(this->node(n));
//           nodes_set.insert(rhs.node(n));
//         }

//       // If this passes the elements are connected
//       // to the same global nodes
//       if (nodes_set.size() == this->n_nodes())
//         return true;
//     }

//   // If we get here it is because the elements either
//   // do not have the same number of nodes or they are
//   // connected to different nodes.  Either way they
//   // are not the same element
//   return false;
  
  // Useful typedefs
  typedef std::vector<unsigned int>::iterator iterator;

  
  // Elements can only be equal if they
  // contain the same number of nodes.
  // However, we will only test the vertices,
  // which is sufficient & cheaper
  if (this->n_nodes() == rhs.n_nodes())
    {
      // The number of nodes in the element
      const unsigned int nn = this->n_nodes();
      
      // Create a vector that contains our global
      // node numbers and those of our neighbor.
      // If the sorted, unique vector is the same size
      // as the number of nodes in both elements then
      // they must be connected to the same nodes.
      //
      // The vector will be no larger than 2*n_nodes(),
      // so we might as well reserve the space.
      std::vector<unsigned int> common_nodes;
      common_nodes.reserve (2*nn);

      // Add the global indices of the nodes
      for (unsigned int n=0; n<nn; n++)
	{
	  common_nodes.push_back (this->node(n));
	  common_nodes.push_back (rhs.node(n));
	}

      // Sort the vector and find out how long
      // the sorted vector is.
      std::sort (common_nodes.begin(), common_nodes.end());
      
      iterator new_end = std::unique (common_nodes.begin(),
				      common_nodes.end());
      
      const int new_size = std::distance (common_nodes.begin(),
					  new_end);
      
      // If this passes the elements are connected
      // to the same global vertex nodes
      if (new_size == static_cast<int>(nn))
	return true;
    }

  // If we get here it is because the elements either
  // do not have the same number of nodes or they are
  // connected to different nodes.  Either way they
  // are not the same element
  return false;
}



void Elem::write_ucd_connectivity(std::ostream &out) const
{
  assert (out);
  assert (_nodes != NULL);

  for (unsigned int i=0; i<this->n_nodes(); i++)
    out << this->node(i)+1 << "\t";

  out << std::endl;
}



Real Elem::quality (const ElemQuality q) const
{
  switch (q)
    {    
      /**
       * I don't know what to do for this metric. 
       */
    default:
      {
	here();

	std::cerr << "ERROR:  unknown quality metric: "
		  << q 
		  << std::endl
		  << "Cowardly returning 1."
		  << std::endl;

	return 1.;
      }
    }

    
    // Will never get here...
    error();
    return 0.;
}


#ifdef ENABLE_AMR

void Elem::family_tree (std::vector<const Elem*>& family,
			const bool reset) const
{
  // Clear the vector if the flag reset tells us to.
  if (reset)
    family.clear();

  // Add this element to the family tree.
  family.push_back(this);

  // Recurse into the elements children, if it has them.
  // Do not clear the vector any more.
  if (!this->active())
    for (unsigned int c=0; c<this->n_children(); c++)
      this->child(c)->family_tree (family, false);
}



void Elem::active_family_tree (std::vector<const Elem*>& active_family,
			       const bool reset) const
{
  // Clear the vector if the flag reset tells us to.
  if (reset)
    active_family.clear();

  // Add this element to the family tree if it is active
  if (this->active())
    active_family.push_back(this);

  // Otherwise recurse into the element's children.
  // Do not clear the vector any more.
  else 
    for (unsigned int c=0; c<this->n_children(); c++)
      this->child(c)->active_family_tree (active_family, false);

}

#endif // #ifdef ENABLE_AMR



bool Elem::contains_point (const Point& p) const
{
  // Declare a basic FEType.  Will ue a Lagrange
  // element by default.
  FEType fe_type(this->default_order());
  
  const Point mapped_point = FEInterface::inverse_map(this->dim(),
						      fe_type,
						      this,
						      p,
						      1.e-4,
						      false);

  return FEInterface::on_reference_element(mapped_point, this->type());
}



void Elem::nullify_neighbors ()
{
  // Tell any of my neighbors about my death...
  // Looks strange, huh?
  for (unsigned int n=0; n<this->n_neighbors(); n++)
    if (this->neighbor(n) != NULL)
      {
	Elem* neighbor = this->neighbor(n);

	// Note:  it is possible that I see the neighbor
	// (which is coarser than me)
	// but they don't see me, so avoid that case.
	if (neighbor->level() == this->level())
	  {	
	    const unsigned int w_n_a_i = neighbor->which_neighbor_am_i(this);
	    neighbor->set_neighbor(w_n_a_i, NULL);
	    this->set_neighbor(n, NULL);
	  }
      }
}




unsigned int Elem::n_second_order_adjacent_vertices (const unsigned int) const
{
  // for linear elements, always return 0
  return 0;
}



unsigned short int Elem::second_order_adjacent_vertex (const unsigned int,
						       const unsigned int) const
{
  // for linear elements, always return 0
  return 0;
}




ElemType Elem::second_order_equivalent_type (const ElemType et,
					     const bool full_ordered)
{ 
  /* for second-order elements, always return \p INVALID_ELEM
   * since second-order elements should not be converted
   * into something else.  Only linear elements should 
   * return something sensible here
   */
  switch (et)
    {
    case EDGE2:
      {
	// full_ordered not relevant
	return EDGE3;
      }

    case TRI3:
      {
	// full_ordered not relevant
	return TRI6;
      }

    case QUAD4:
      {
	if (full_ordered)
	  return QUAD9;
	else 
	  return QUAD8;
      }

    case TET4:
      {
	// full_ordered not relevant
	return TET10;
      }

    case HEX8:
      {
	// see below how this correlates with INFHEX8
	if (full_ordered)
	  return HEX27;
	else 
	  return HEX20;
      }

    case PRISM6:
      {
	if (full_ordered)
	  return PRISM18;
	else 
	  return PRISM15;
      }

    case PYRAMID5:
      {
	// error(); 
	return INVALID_ELEM;
      }



#ifdef ENABLE_INFINITE_ELEMENTS

    // infinite elements
    case INFEDGE2:
      {
	// error(); 
	return INVALID_ELEM;
      }

    case INFQUAD4:
      {
	// full_ordered not relevant
	return INFQUAD6;
      }

    case INFHEX8:
      {
	/*
	 * Note that this matches with \p Hex8:
	 * For full-ordered, \p InfHex18 and \p Hex27
	 * belong together, and for not full-ordered,
	 * \p InfHex16 and \p Hex20 belong together.
	 */
	if (full_ordered)
	  return INFHEX18;
	else 
	  return INFHEX16;
      }

    case INFPRISM6:
      {
	// full_ordered not relevant
	return INFPRISM12;
      }

#endif


    default:
      {
	// second-order element
	return INVALID_ELEM;
      }
    }
}




// $Id: elem.C,v 1.60 2006-10-21 18:28:18 roystgnr Exp $

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
#include <algorithm> // for std::sort
#include <iterator>  // for std::ostream_iterator

// Local includes
#include "elem.h"
#include "fe_type.h"
#include "fe_interface.h"
#include "edge_edge2.h"
#include "edge_edge3.h"
#include "edge_edge4.h"
#include "edge_inf_edge2.h"
#include "face_tri3.h"
#include "face_tri6.h"
#include "face_quad4.h"
#include "face_quad8.h"
#include "face_quad9.h"
#include "face_inf_quad4.h"
#include "face_inf_quad6.h"
#include "cell_tet4.h"
#include "cell_tet10.h"
#include "cell_hex8.h"
#include "cell_hex20.h"
#include "cell_hex27.h"
#include "cell_inf_hex8.h"
#include "cell_inf_hex16.h"
#include "cell_inf_hex18.h"
#include "cell_prism6.h"
#include "cell_prism15.h"
#include "cell_prism18.h"
#include "cell_inf_prism6.h"
#include "cell_inf_prism12.h"
#include "cell_pyramid5.h"
#include "fe_base.h"
#include "quadrature_gauss.h"

// Initialize static member variables
const unsigned int Elem::_bp1 = 65449;
const unsigned int Elem::_bp2 = 48661;

// ------------------------------------------------------------
// Elem class member funcions
AutoPtr<Elem> Elem::build(const ElemType type,
			  Elem* p)
{
  Elem* elem = NULL;
 
  switch (type)
    {
      // 1D elements 
    case EDGE2:
      {
	elem = new Edge2(p);
	break;
      }
    case EDGE3:
      {
	elem = new Edge3(p);
	break;
      }
    case EDGE4:
      {
        elem = new Edge4(p);
        break;
      }


      
      // 2D elements
    case TRI3:
      {
	elem = new Tri3(p);
	break;
      }
    case TRI6:
      {
	elem = new Tri6(p);
	break;
      }
    case QUAD4:
      {
	elem = new Quad4(p);
	break;
      }
    case QUAD8:
      {
	elem = new Quad8(p);
	break;
      }
    case QUAD9:
      {
	elem = new Quad9(p);
	break;
      }
   

      // 3D elements
    case TET4:
      {
	elem = new Tet4(p);
	break;
      }
    case TET10:
      {
	elem = new Tet10(p);
	break;
      }
    case HEX8:
      {
	elem = new Hex8(p);
	break;
      }
    case HEX20:
      {
	elem = new Hex20(p);
	break;
      }
    case HEX27:
      {
	elem = new Hex27(p);
	break;
      }
    case PRISM6:
      {
	elem = new Prism6(p);
	break;
      }
    case PRISM15:
      {
	elem = new Prism15(p);
	break;
      }
    case PRISM18:
      {
	elem = new Prism18(p);
	break;
      }
    case PYRAMID5:
      {
	elem = new Pyramid5(p);
	break;
      }



#ifdef ENABLE_INFINITE_ELEMENTS

      // 1D infinite elements
    case INFEDGE2:
      {
	elem = new InfEdge2(p);
	break;
      }


      // 2D infinite elements
    case INFQUAD4:
      {
	elem = new InfQuad4(p);
	break;
      }
    case INFQUAD6:
      {
	elem = new InfQuad6(p);
	break;
      }

   
    // 3D infinite elements
    case INFHEX8:
      {
	elem = new InfHex8(p);
	break;
      }
    case INFHEX16:
      {
	elem = new InfHex16(p);
	break;
      }
    case INFHEX18:
      {
	elem = new InfHex18(p);
	break;
      }
    case INFPRISM6:
      {
	elem = new InfPrism6(p);
	break;
      }
    case INFPRISM12:
      {
	elem = new InfPrism12(p);
	break;
      }

#endif
           
    default:
      {
	std::cerr << "ERROR: Undefined element type!." << std::endl;
	error();
      }
    }
  

  AutoPtr<Elem> ap(elem);
  return ap;
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



bool Elem::operator == (const DofObject& rhs) const
{

    // Cast rhs to an Elem*
    const Elem* rhs_elem = dynamic_cast<const Elem*>(&rhs);

    // If we cannot cast to an Elem*, rhs must be a Node
    if(rhs_elem == static_cast<const Elem*>(NULL))
        return false;

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
  if (this->n_nodes() == rhs_elem->n_nodes())
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
	  common_nodes.push_back (rhs_elem->node(n));
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





void Elem::write_connectivity (std::ostream& out,
			       const IOPackage iop) const
{
  assert (out.good());
  assert (_nodes != NULL);
  assert (iop != INVALID_IO_PACKAGE);

  switch (iop)
    {
    case TECPLOT:
      {
	// This connectivity vector will be used repeatedly instead
	// of being reconstructed inside the loop.
	std::vector<unsigned int> conn;
	for (unsigned int sc=0; sc <this->n_sub_elem(); sc++)
	  {
	    this->connectivity(sc, TECPLOT, conn);
	    
	    std::copy(conn.begin(),
		      conn.end(),
		      std::ostream_iterator<unsigned int>(out, " "));
	    
	    out << '\n';
	  }
	return;
      }

    case UCD:
      {
	for (unsigned int i=0; i<this->n_nodes(); i++)
	  out << this->node(i)+1 << "\t";
	
	out << '\n';
	return;
      }

    default:
      error();
    }

  error();
}


// void Elem::write_tecplot_connectivity(std::ostream& out) const
// {
//   assert (!out.bad());
//   assert (_nodes != NULL);

//   // This connectivity vector will be used repeatedly instead
//   // of being reconstructed inside the loop.
//   std::vector<unsigned int> conn;
//   for (unsigned int sc=0; sc <this->n_sub_elem(); sc++)
//     {
//       this->connectivity(sc, TECPLOT, conn);
      
//       std::copy(conn.begin(),
//  		conn.end(),
//  		std::ostream_iterator<unsigned int>(out, " "));
      
//       out << std::endl;
//     }
// }



// void Elem::write_ucd_connectivity(std::ostream &out) const
// {
//   assert (out);
//   assert (_nodes != NULL);

//   for (unsigned int i=0; i<this->n_nodes(); i++)
//     out << this->node(i)+1 << "\t";

//   out << std::endl;
// }



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



bool Elem::ancestor() const
{
#ifdef ENABLE_AMR

  if (this->active())
    return false;

if (!this->has_children())
    return false;
  if (this->child(0)->active())
    return true;

  return this->child(0)->ancestor();
#else
  return false;
#endif
}



void Elem::add_child (Elem* elem)
{
  if(_children == NULL)
  {
    _children = new Elem*[this->n_children()];
    
    for (unsigned int c=0; c<this->n_children(); c++)
      _children[c] = NULL;
  }

  for (unsigned int c=0; c<this->n_children(); c++)
  {
    if(_children[c] == NULL)
    {
      _children[c] = elem;
      return;
    }
  }

  std::cerr << "Error: Tried to add a child to an element with full children array"
            << std::endl;
  error();
}





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



void Elem::family_tree_by_neighbor (std::vector<const Elem*>& family,
                                    const Elem* neighbor,
			            const bool reset) const
{
  // Clear the vector if the flag reset tells us to.
  if (reset)
    family.clear();

  // This only makes sense if we're already a neighbor
  assert (this->is_neighbor(neighbor));

  // Add this element to the family tree.
  family.push_back(this);

  // Recurse into the elements children, if it has them.
  // Do not clear the vector any more.
  if (!this->active())
    for (unsigned int c=0; c<this->n_children(); c++)
      if (this->child(c)->is_neighbor(neighbor))
        this->child(c)->family_tree_by_neighbor (family, false);
}



void Elem::active_family_tree_by_neighbor (std::vector<const Elem*>& family,
                                           const Elem* neighbor,
			                   const bool reset) const
{
  // Clear the vector if the flag reset tells us to.
  if (reset)
    family.clear();

  // This only makes sense if we're already a neighbor
  assert (this->is_neighbor(neighbor));

  // Add an active element to the family tree.
  if (this->active())
    family.push_back(this);

  // Or recurse into an ancestor element's children.
  // Do not clear the vector any more.
  else
    for (unsigned int c=0; c<this->n_children(); c++)
      if (this->child(c)->is_neighbor(neighbor))
        this->child(c)->active_family_tree_by_neighbor (family, false);
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




ElemType Elem::first_order_equivalent_type (const ElemType et)
{ 
  switch (et)
    {
    case EDGE4:
    case EDGE3:
      return EDGE2;
    case TRI6:
      return TRI3;
    case QUAD8:
    case QUAD9:
      return QUAD4;
    case TET10:
      return TET4;
    case HEX27:
    case HEX20:
      return HEX8;
    case PRISM15:
    case PRISM18:
      return PRISM6;

#ifdef ENABLE_INFINITE_ELEMENTS

    case INFQUAD6:
      return INFQUAD4;
    case INFHEX16:
    case INFHEX18:
      return INFHEX8;
    case INFPRISM12:
      return INFPRISM6;

#endif

    default:
      // first-order or unknown element
      return INVALID_ELEM;
    }
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



Elem::side_iterator Elem::boundary_sides_begin()
{
  Predicates::BoundarySide<SideIter> bsp;
  return side_iterator(this->_first_side(), this->_last_side(), bsp);
}




Elem::side_iterator Elem::boundary_sides_end()
{
  Predicates::BoundarySide<SideIter> bsp;
  return side_iterator(this->_last_side(), this->_last_side(), bsp);
}


bool Elem::is_child_on_edge(const unsigned int c,
                            const unsigned int e) const
{
  assert (c < this->n_children());
  assert (e < this->n_edges());

  AutoPtr<Elem> my_edge = this->build_edge(e);
  AutoPtr<Elem> child_edge = this->build_edge(e);

  // We're assuming that an overlapping child edge has the same
  // number and orientation as its parent
  return (child_edge->node(0) == my_edge->node(0) ||
      child_edge->node(1) == my_edge->node(1));
}


bool Elem::is_child_on_side(const unsigned int c,
                            const unsigned int s) const
{
  assert (c < this->n_children());
  assert (s < this->n_sides());

  Elem *child = this->child(c);
  return ((child->neighbor(s) == NULL) // on boundary
      || (child->neighbor(s)->parent() != this));
}


unsigned int Elem::min_p_level_by_neighbor(const Elem* neighbor,
                                           unsigned int current_min) const
{
  assert(!this->subactive());
  assert(neighbor->active());

  // If we're an active element this is simple
  if (this->active())
    return std::min(current_min, this->p_level());

  assert(is_neighbor(neighbor));

  // The p_level() of an ancestor element is already the minimum
  // p_level() of its children - so if that's high enough, we don't
  // need to examine any children.
  if (current_min <= this->p_level())
    return current_min;

  unsigned int min_p_level = current_min;

  for (unsigned int c=0; c<this->n_children(); c++)
    {
      const Elem* const child = this->child(c);
      if (child->is_neighbor(neighbor))
        min_p_level =
	  child->min_p_level_by_neighbor(neighbor,
                                         min_p_level);
    }

  return min_p_level;
}


unsigned int Elem::min_new_p_level_by_neighbor(const Elem* neighbor,
                                               unsigned int current_min) const
{
  assert(!this->subactive());
  assert(neighbor->active());

  // If we're an active element this is simple
  if (this->active())
    {
      unsigned int new_p_level = this->p_level();
      if (this->p_refinement_flag() == Elem::REFINE)
        new_p_level += 1;
      if (this->p_refinement_flag() == Elem::COARSEN)
        {
          assert (new_p_level > 0);
          new_p_level -= 1;
        }
      return std::min(current_min, new_p_level);
    }

  assert(is_neighbor(neighbor));

  unsigned int min_p_level = current_min;

  for (unsigned int c=0; c<this->n_children(); c++)
    {
      const Elem* const child = this->child(c);
      if (child->is_neighbor(neighbor))
        min_p_level =
	  child->min_new_p_level_by_neighbor(neighbor,
                                             min_p_level);
    }

  return min_p_level;
}




Real Elem::volume () const
{
  // The default implementation builds a finite element of the correct
  // order and sums up the JxW contributions.  This can be expensive,
  // so the various element types can overload this method and compute
  // the volume more efficiently.
  FEType fe_type (this->default_order() , LAGRANGE);

  AutoPtr<FEBase> fe (FEBase::build(this->dim(),
				    fe_type));

   const std::vector<Real>& JxW = fe->get_JxW();
   
  // The default quadrature rule should integrate the mass matrix,
  // thus it should be plenty to compute the area
  QGauss qrule (this->dim(), fe_type.default_quadrature_order());

  fe->attach_quadrature_rule(&qrule);

  fe->reinit(this);

  Real vol=0.;
  for (unsigned int qp=0; qp<qrule.n_points(); ++qp)
    vol += JxW[qp];
  
  return vol;
  
}

// $Id: fe_monomial.C,v 1.15 2003-11-30 07:10:16 benkirk Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2002-2003  Benjamin S. Kirk, John W. Peterson
  
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



// Local includes
#include "fe.h"
#include "elem.h"




// ------------------------------------------------------------
// Monomials-specific implementations
template <unsigned int Dim, FEFamily T>
void FE<Dim,T>::nodal_soln(const Elem* elem,
			   const Order order,
			   const std::vector<Number>& elem_soln,
			   std::vector<Number>&       nodal_soln)
{
  const unsigned int n_nodes = elem->n_nodes();
  
  const ElemType type = elem->type();

  nodal_soln.resize(n_nodes);


  
  switch (order)
    {
      // Constant shape functions
    case CONSTANT:
      {
	assert (elem_soln.size() == 1);
	
	const Number val = elem_soln[0];
	
	for (unsigned int n=0; n<n_nodes; n++)
	  nodal_soln[n] = val;
	
	return;
      }


      // For other bases do interpolation at the nodes
      // explicitly.
    case FIRST:
    case SECOND:
    case THIRD:
    case FOURTH:
      {

	const unsigned int n_sf =
	  FE<Dim,T>::n_shape_functions(type, order);
	
	for (unsigned int n=0; n<n_nodes; n++)
	  {
	    const Point mapped_point = FE<Dim,T>::inverse_map(elem,
							      elem->point(n));

	    assert (elem_soln.size() == n_sf);

	    // Zero before summation
	    nodal_soln[n] = 0;

	    // u_i = Sum (alpha_i phi_i)
	    for (unsigned int i=0; i<n_sf; i++)
	      nodal_soln[n] += elem_soln[i]*FE<Dim,T>::shape(elem,
							     order,
							     i,
							     mapped_point);	    
	  }

	return;
      }
      
    default:
      {
	error();
      }
    }
}



template <unsigned int Dim, FEFamily T>
unsigned int FE<Dim,T>::n_dofs(const ElemType t, const Order o)
{
  switch (o)
    {

      // constant shape functions
      // no matter what shape there is only one DOF.
    case CONSTANT:
      return 1;


      // Discontinuous linear shape functions
      // expressed in the monomials.
    case FIRST:
      {
	switch (t)
	  {
	  case EDGE2:
	  case EDGE3:
	  case EDGE4:
	    return 2;

	  case TRI3:
	  case TRI6:
	  case QUAD4:
	  case QUAD8:
	  case QUAD9:
	    return 3;

	  case TET4:
	  case TET10:
	  case HEX8:
	  case HEX20:
	  case HEX27:
	  case PRISM6:
	  case PRISM15:
	  case PRISM18:
	  case PYRAMID5:
	    return 4;
	    
	  default:
	    {
#ifdef DEBUG
	      std::cerr << "ERROR: Bad ElemType = " << t
			<< " for " << o << "th order approximation!" 
			<< std::endl;
#endif
	      error();	    
	    }
	  }
      }


      // Discontinuous quadratic shape functions
      // expressed in the monomials.
    case SECOND:
      {
	switch (t)
	  {
	  case EDGE2:
	  case EDGE3:
	  case EDGE4:
	    return 3;

	  case TRI3:
	  case TRI6:
	  case QUAD4:
	  case QUAD8:
	  case QUAD9:
	    return 6;

	  case TET4:
	  case TET10:
	  case HEX8:
	  case HEX20:
	  case HEX27:
	  case PRISM6:
	  case PRISM15:
	  case PRISM18:
	  case PYRAMID5:
	    return 10;
	    
	  default:
	    {
#ifdef DEBUG
	      std::cerr << "ERROR: Bad ElemType = " << t
			<< " for " << o << "th order approximation!" 
			<< std::endl;
#endif
	      error();	    
	    }
	  }
      }


      // Discontinuous cubic shape functions
      // expressed in the monomials.
    case THIRD:
      {
	switch (t)
	  {
	  case EDGE2:
	  case EDGE3:
	  case EDGE4:
	    return 4;

	  case TRI3:
	  case TRI6:
	  case QUAD4:
	  case QUAD8:
	  case QUAD9:
	    return 10;

	  case TET4:
	  case TET10:
	  case HEX8:
	  case HEX20:
	  case HEX27:
	  case PRISM6:
	  case PRISM15:
	  case PRISM18:
	  case PYRAMID5:
	    return 20;
	    
	  default:
	    {
#ifdef DEBUG
	      std::cerr << "ERROR: Bad ElemType = " << t
			<< " for " << o << "th order approximation!" 
			<< std::endl;
#endif
	      error();	    
	    }
	  }
      }


      // Discontinuous quartic shape functions
      // expressed in the monomials.
    case FOURTH:
      {
	switch (t)
	  {
	  case EDGE2:
	  case EDGE3:
	    return 5;

	  case TRI3:
	  case TRI6:
	  case QUAD4:
	  case QUAD8:
	  case QUAD9:
	    return 15;

	  case TET4:
	  case TET10:
	  case HEX8:
	  case HEX20:
	  case HEX27:
	  case PRISM6:
	  case PRISM15:
	  case PRISM18:
	  case PYRAMID5:
	    return 35;
	    
	  default:
	    {
#ifdef DEBUG
	      std::cerr << "ERROR: Bad ElemType = " << t
			<< " for " << o << "th order approximation!" 
			<< std::endl;
#endif
	      error();	    
	    }
	  }
      }

      
    default:
      {
	error();
      }
    }
  
  error();
  
  return 0;
}



template <unsigned int Dim, FEFamily T>
unsigned int FE<Dim,T>::n_dofs_at_node(const ElemType,
				       const Order,
				       const unsigned int)
{
  // Monomials elements have no dofs at nodes
  // (just on the element)
  return 0;
}



template <unsigned int Dim, FEFamily T>
unsigned int FE<Dim,T>::n_dofs_per_elem(const ElemType t,
					const Order o)
{
  switch (o)
    {
      // constant shape functions always have 1 DOF per element
    case CONSTANT:
      return 1;

      
      // Discontinuous linear shape functions
      // expressed in the monomials.
    case FIRST:
      {
	switch (t)
	  {
	    // 1D linears have 2 DOFs per element
	  case EDGE2:
	  case EDGE3:
	  case EDGE4:
	    return 2;

	    // 2D linears have 3 DOFs per element
	  case TRI3:
	  case TRI6:
	  case QUAD4:
	  case QUAD8:
	  case QUAD9:
	    return 3;

	    // 3D linears have 4 DOFs per element
 	  case TET4:
	  case TET10:
	  case HEX8:
	  case HEX20:
	  case HEX27:
	  case PRISM6:
	  case PRISM15:
	  case PRISM18:
	  case PYRAMID5:
	    return 4;
	    
	  default:
	    {
#ifdef DEBUG
	      std::cerr << "ERROR: Bad ElemType = " << t
			<< " for " << o << "th order approximation!" 
			<< std::endl;
#endif
	      error();	    
	    }
	  }
      }


      // Discontinuous quadratic shape functions
      // expressed in the monomials.
    case SECOND:
      {
	switch (t)
	  {
	    // 1D quadratics have 3 DOFs per element
	  case EDGE2:
	  case EDGE3:
	  case EDGE4:
	    return 3;

	    // 2D quadratics have 6 DOFs per element
	  case TRI3:
	  case TRI6:
	  case QUAD4:
	  case QUAD8:
	  case QUAD9:
	    return 6;

	    // 3D quadratics have 10 DOFs per element
	  case TET4:
	  case TET10:
	  case HEX8:
	  case HEX20:
	  case HEX27:
	  case PRISM6:
	  case PRISM15:
	  case PRISM18:
	  case PYRAMID5:
	    return 10;
	    
	  default:
	    {
#ifdef DEBUG
	      std::cerr << "ERROR: Bad ElemType = " << t
			<< " for " << o << "th order approximation!" 
			<< std::endl;
#endif
	      error();	    
	    }
	  }
      }


      // Discontinuous cubic shape functions
      // expressed in the monomials.
    case THIRD:
      {
	switch (t)
	  {
	  case EDGE2:
	  case EDGE3:
	  case EDGE4:
	    return 4;

	  case TRI3:
	  case TRI6:
	  case QUAD4:
	  case QUAD8:
	  case QUAD9:
	    return 10;

	  case TET4:
	  case TET10:
	  case HEX8:
	  case HEX20:
	  case HEX27:
	  case PRISM6:
	  case PRISM15:
	  case PRISM18:
	  case PYRAMID5:
	    return 20;
	    
	  default:
	    {
#ifdef DEBUG
	      std::cerr << "ERROR: Bad ElemType = " << t
			<< " for " << o << "th order approximation!" 
			<< std::endl;
#endif
	      error();	    
	    }
	  }
      }


      // Discontinuous quartic shape functions
      // expressed in the monomials.
    case FOURTH:
      {
	switch (t)
	  {
	  case EDGE2:
	  case EDGE3:
	  case EDGE4:
	    return 5;

	  case TRI3:
	  case TRI6:
	  case QUAD4:
	  case QUAD8:
	  case QUAD9:
	    return 15;

	  case TET4:
	  case TET10:
	  case HEX8:
	  case HEX20:
	  case HEX27:
	  case PRISM6:
	  case PRISM15:
	  case PRISM18:
	  case PYRAMID5:
	    return 35;
	    
	  default:
	    {
#ifdef DEBUG
	      std::cerr << "ERROR: Bad ElemType = " << t
			<< " for " << o << "th order approximation!" 
			<< std::endl;
#endif
	      error();	    
	    }
	  }
      }


      
      // Otherwise no DOFS per element
    default:
      return 0;
    }
}



template <unsigned int Dim, FEFamily T>
void FE<Dim,T>::compute_constraints (std::map<unsigned int,
				            std::map<unsigned int,
				                     float> > &,
				     const unsigned int,
				     const unsigned int,
				     const FEType&,
				     const Elem*)
{
  // Monomials are discontinuous...  No constraints.
  return;
}



template <unsigned int Dim, FEFamily T>
bool FE<Dim,T>::shapes_need_reinit() const
{
  return false;
}



//--------------------------------------------------------------
// Explicit instantiations
template class FE<1,MONOMIAL>;
template class FE<2,MONOMIAL>;
template class FE<3,MONOMIAL>;

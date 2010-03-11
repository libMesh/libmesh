// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2008 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner
  
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
#include "fe_macro.h"
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

  const Order totalorder = static_cast<Order>(order+elem->p_level());
  
  switch (totalorder)
    {
      // Constant shape functions
    case CONSTANT:
      {
	libmesh_assert (elem_soln.size() == 1);
	
	const Number val = elem_soln[0];
	
	for (unsigned int n=0; n<n_nodes; n++)
	  nodal_soln[n] = val;
	
	return;
      }


      // For other bases do interpolation at the nodes
      // explicitly.
    default:
      {

	const unsigned int n_sf =
	  FE<Dim,T>::n_shape_functions(type, totalorder);
	
	for (unsigned int n=0; n<n_nodes; n++)
	  {
	    const Point mapped_point = FE<Dim,T>::inverse_map(elem,
							      elem->point(n));

	    libmesh_assert (elem_soln.size() == n_sf);

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
	  case NODEELEM:
	    return 1;

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
	      libMesh::err << "ERROR: Bad ElemType = " << t
			    << " for " << o << "th order approximation!" 
			    << std::endl;
#endif
	      libmesh_error();	    
	    }
	  }
      }


      // Discontinuous quadratic shape functions
      // expressed in the monomials.
    case SECOND:
      {
	switch (t)
	  {
	  case NODEELEM:
	    return 1;

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
	      libMesh::err << "ERROR: Bad ElemType = " << t
			    << " for " << o << "th order approximation!" 
			    << std::endl;
#endif
	      libmesh_error();	    
	    }
	  }
      }


      // Discontinuous cubic shape functions
      // expressed in the monomials.
    case THIRD:
      {
	switch (t)
	  {
	  case NODEELEM:
	    return 1;

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
	      libMesh::err << "ERROR: Bad ElemType = " << t
			    << " for " << o << "th order approximation!" 
			    << std::endl;
#endif
	      libmesh_error();	    
	    }
	  }
      }


      // Discontinuous quartic shape functions
      // expressed in the monomials.
    case FOURTH:
      {
	switch (t)
	  {
	  case NODEELEM:
	    return 1;

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
	      libMesh::err << "ERROR: Bad ElemType = " << t
			    << " for " << o << "th order approximation!" 
			    << std::endl;
#endif
	      libmesh_error();	    
	    }
	  }
      }

      
    default:
      {
        const unsigned int order = static_cast<unsigned int>(o);
	switch (t)
	  {
	  case NODEELEM:
	    return 1;

	  case EDGE2:
	  case EDGE3:
	    return (order+1);

	  case TRI3:
	  case TRI6:
	  case QUAD4:
	  case QUAD8:
	  case QUAD9:
	    return (order+1)*(order+2)/2;

	  case TET4:
	  case TET10:
	  case HEX8:
	  case HEX20:
	  case HEX27:
	  case PRISM6:
	  case PRISM15:
	  case PRISM18:
	  case PYRAMID5:
	    return (order+1)*(order+2)*(order+3)/6;
	    
	  default:
	    {
#ifdef DEBUG
	      libMesh::err << "ERROR: Bad ElemType = " << t
			    << " for " << o << "th order approximation!" 
			    << std::endl;
#endif
	      libmesh_error();	    
	    }
	  }
      }
    }
  
  libmesh_error();
  
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
	  case NODEELEM:
	    return 1;

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
	      libMesh::err << "ERROR: Bad ElemType = " << t
			    << " for " << o << "th order approximation!" 
			    << std::endl;
#endif
	      libmesh_error();	    
	    }
	  }
      }


      // Discontinuous quadratic shape functions
      // expressed in the monomials.
    case SECOND:
      {
	switch (t)
	  {
	  case NODEELEM:
	    return 1;

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
	      libMesh::err << "ERROR: Bad ElemType = " << t
			    << " for " << o << "th order approximation!" 
			    << std::endl;
#endif
	      libmesh_error();	    
	    }
	  }
      }


      // Discontinuous cubic shape functions
      // expressed in the monomials.
    case THIRD:
      {
	switch (t)
	  {
	  case NODEELEM:
	    return 1;

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
	      libMesh::err << "ERROR: Bad ElemType = " << t
			    << " for " << o << "th order approximation!" 
			    << std::endl;
#endif
	      libmesh_error();	    
	    }
	  }
      }


      // Discontinuous quartic shape functions
      // expressed in the monomials.
    case FOURTH:
      {
	switch (t)
	  {
	  case NODEELEM:
	    return 1;

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
	      libMesh::err << "ERROR: Bad ElemType = " << t
			    << " for " << o << "th order approximation!" 
			    << std::endl;
#endif
	      libmesh_error();	    
	    }
	  }
      }
      
    default:
      {
        const unsigned int order = static_cast<unsigned int>(o);
	switch (t)
	  {
	  case NODEELEM:
	    return 1;

	  case EDGE2:
	  case EDGE3:
	    return (order+1);

	  case TRI3:
	  case TRI6:
	  case QUAD4:
	  case QUAD8:
	  case QUAD9:
	    return (order+1)*(order+2)/2;

	  case TET4:
	  case TET10:
	  case HEX8:
	  case HEX20:
	  case HEX27:
	  case PRISM6:
	  case PRISM15:
	  case PRISM18:
	  case PYRAMID5:
	    return (order+1)*(order+2)*(order+3)/6;
	    
	  default:
	    {
#ifdef DEBUG
	      libMesh::err << "ERROR: Bad ElemType = " << t
			    << " for " << o << "th order approximation!" 
			    << std::endl;
#endif
	      libmesh_error();	    
	    }
	  }
      }

      return 0;
    }
}



template <unsigned int Dim, FEFamily T>
FEContinuity FE<Dim,T>::get_continuity() const
{
  return DISCONTINUOUS;
}



template <unsigned int Dim, FEFamily T>
bool FE<Dim,T>::is_hierarchic() const
{
  return true;
}



#ifdef LIBMESH_ENABLE_AMR
template <unsigned int Dim, FEFamily T>
void FE<Dim,T>::compute_constraints (DofConstraints &,
				     DofMap &,
				     const unsigned int,
				     const Elem*)
{
  // Monomials are discontinuous...  No constraints.
  return;
}
#endif // #ifdef LIBMESH_ENABLE_AMR



template <unsigned int Dim, FEFamily T>
bool FE<Dim,T>::shapes_need_reinit() const
{
  return false;
}



//--------------------------------------------------------------
// Explicit instantiation of member functions
INSTANTIATE_MBRF(0,MONOMIAL);
INSTANTIATE_MBRF(1,MONOMIAL);
INSTANTIATE_MBRF(2,MONOMIAL);
INSTANTIATE_MBRF(3,MONOMIAL);

#ifdef LIBMESH_ENABLE_AMR
template void FE<2,MONOMIAL>::compute_constraints(DofConstraints&, DofMap&, 
						  const unsigned int,
						  const Elem*);
template void FE<3,MONOMIAL>::compute_constraints(DofConstraints&, DofMap&, 
						  const unsigned int,
						  const Elem*);
#endif // #ifdef LIBMESH_ENABLE_AMR

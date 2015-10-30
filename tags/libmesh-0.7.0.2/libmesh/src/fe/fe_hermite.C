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
#include "elem.h"
#include "fe.h"
#include "fe_macro.h"

namespace libMesh
{




// ------------------------------------------------------------
// Hierarchic-specific implementations
template <unsigned int Dim, FEFamily T>
void FE<Dim,T>::nodal_soln(const Elem* elem,
			   const Order order,
			   const std::vector<Number>& elem_soln,
			   std::vector<Number>&       nodal_soln)
{
  const unsigned int n_nodes = elem->n_nodes();
  
  const ElemType type = elem->type();

  nodal_soln.resize(n_nodes);

  const Order totalorder = static_cast<Order>(order + elem->p_level());
  
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
}



template <unsigned int Dim, FEFamily T>
unsigned int FE<Dim,T>::n_dofs(const ElemType t, const Order o)
{
  libmesh_assert (o > 2);
  // Piecewise (bi/tri)cubic C1 Hermite splines
  switch (t)
    {
    case NODEELEM:
      return 1;
    case EDGE2:
      libmesh_assert (o < 4);
    case EDGE3:
      return (o+1);
	    
    case QUAD4:
    case QUAD8:
      libmesh_assert (o < 4);
    case QUAD9:
      return ((o+1)*(o+1));
	    
    case HEX8:
    case HEX20:
      libmesh_assert (o < 4);
    case HEX27:
      return ((o+1)*(o+1)*(o+1));
	    
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
  
  libmesh_error();  
  return 0;
}



template <unsigned int Dim, FEFamily T>
unsigned int FE<Dim,T>::n_dofs_at_node(const ElemType t,
				       const Order o,
				       const unsigned int n)
{
  libmesh_assert (o > 2);
  // Piecewise (bi/tri)cubic C1 Hermite splines
  switch (t)
    {
    case NODEELEM:
      return 1;
    case EDGE2:
    case EDGE3:
      {
        switch (n)
	{
	  case 0:
	  case 1:
	    return 2;
	  case 2:
//          Interior DoFs are carried on Elems
//	    return (o-3);
            return 0;

	  default:
	    libmesh_error();
	}
      }
	    
    case QUAD4:
      libmesh_assert (o < 4);
    case QUAD8:
    case QUAD9:
      {
	switch (n)
	  {
          // Vertices
	  case 0:
	  case 1:
	  case 2:
	  case 3:
	    return 4;
          // Edges
	  case 4:
	  case 5:
	  case 6:
	  case 7:
	    return (2*(o-3));
	  case 8:
//          Interior DoFs are carried on Elems
//	    return ((o-3)*(o-3));
            return 0;

	  default:
	    libmesh_error();
	  }
      }
	    
    case HEX8:
    case HEX20:
      libmesh_assert (o < 4);
    case HEX27:
      {
        switch (n)
	  {
          // Vertices
	  case 0:
	  case 1:
	  case 2:
	  case 3:
	  case 4:
	  case 5:
	  case 6:
	  case 7:
	    return 8;
          // Edges
	  case 8:
	  case 9:
	  case 10:
	  case 11:
	  case 12:
	  case 13:
	  case 14:
	  case 15:
	  case 16:
	  case 17:
	  case 18:
	  case 19:
	    return (4*(o-3));
          // Faces
	  case 20:
	  case 21:
	  case 22:
	  case 23:
	  case 24:
	  case 25:
	    return (2*(o-3)*(o-3));
	  case 26:
          // Interior DoFs are carried on Elems
//	    return ((o-3)*(o-3)*(o-3));
	    return 0;

	  default:
	    libmesh_error();
	  }
      }
	    
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
  
  libmesh_error();
  
  return 0;
}



template <unsigned int Dim, FEFamily T>
unsigned int FE<Dim,T>::n_dofs_per_elem(const ElemType t,
					const Order o)
{
  libmesh_assert (o > 2);

  switch (t)
    {
    case NODEELEM:
      return 0;
    case EDGE2:
    case EDGE3:
      return (o-3);
    case QUAD4:
      libmesh_assert (o < 4);
    case QUAD8:
    case QUAD9:
      return ((o-3)*(o-3));
    case HEX8:
      libmesh_assert (o < 4);
    case HEX20:
    case HEX27:
      return ((o-3)*(o-3)*(o-3));

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
 
  // Will never get here...
  libmesh_error();
  return 0;
}



template <unsigned int Dim, FEFamily T>
FEContinuity FE<Dim,T>::get_continuity() const
{
  return C_ONE;
}



template <unsigned int Dim, FEFamily T>
bool FE<Dim,T>::is_hierarchic() const
{
  return true;
}



#ifdef LIBMESH_ENABLE_AMR
template <unsigned int Dim, FEFamily T>
void FE<Dim,T>::compute_constraints (DofConstraints &constraints,
				     DofMap &dof_map,
				     const unsigned int variable_number,
				     const Elem* elem)
{
  compute_proj_constraints(constraints, dof_map, variable_number, elem);
}
#endif // #ifdef LIBMESH_ENABLE_AMR



template <unsigned int Dim, FEFamily T>
bool FE<Dim,T>::shapes_need_reinit() const
{
  return true;
}


//--------------------------------------------------------------
// Explicit instantiation of member functions
INSTANTIATE_MBRF(0,HERMITE);
INSTANTIATE_MBRF(1,HERMITE);
INSTANTIATE_MBRF(2,HERMITE);
INSTANTIATE_MBRF(3,HERMITE);

#ifdef LIBMESH_ENABLE_AMR
template void FE<2,HERMITE>::compute_constraints(DofConstraints&, DofMap&,
                                                 const unsigned int,
                                                 const Elem*);
template void FE<3,HERMITE>::compute_constraints(DofConstraints&, DofMap&,
                                                 const unsigned int,
                                                 const Elem*);
#endif // #ifdef LIBMESH_ENABLE_AMR

} // namespace libMesh

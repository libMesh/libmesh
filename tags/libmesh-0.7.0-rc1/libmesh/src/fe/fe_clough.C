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

  
  switch (totalorder)
    {
      // Piecewise cubic shape functions with linear flux edges
    case SECOND:
      // Piecewise cubic shape functions
    case THIRD:
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
      
    default:
      {
	libmesh_error();
      }
    }
}



template <unsigned int Dim, FEFamily T>
unsigned int FE<Dim,T>::n_dofs(const ElemType t, const Order o)
{
  switch (o)
    {
      // Piecewise cubic shape functions with linear flux edges
    case SECOND:
      {
	switch (t)
	  {
	  case TRI6:
	    return 9;
	    
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
      // Piecewise cubic Clough-Tocher element
    case THIRD:
      {
	switch (t)
	  {
	  case TRI6:
	    return 12;
	    
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
  switch (o)
    {
      // Piecewise cubic shape functions with linear flux edges
    case SECOND:
      {
	switch (t)
	  {
	    // The 2D Clough-Tocher defined on a 6-noded triangle
	  case TRI6:
	    {
	      switch (n)
		{
		case 0:
		case 1:
		case 2:
		  return 3;

		case 3:
		case 4:
		case 5:
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
      }
      // The third-order hierarchic shape functions
    case THIRD:
      {
	switch (t)
	  {
	    // The 2D Clough-Tocher defined on a 6-noded triangle
	  case TRI6:
	    {
	      switch (n)
		{
		case 0:
		case 1:
		case 2:
		  return 3;

		case 3:
		case 4:
		case 5:
		  return 1;

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
      }
    default:
      {
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
  switch (o)
    {
      // Piecewise cubic shape functions with linear flux edges
    case SECOND:
      // The third-order Clough-Tocher shape functions
    case THIRD:
      {
	switch (t)
	  {
	    // The 2D hierarchic defined on a 6-noded triangle
	  case TRI6:
	    return 0;

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
      // Otherwise no DOFS per element
    default:
      libmesh_error();	    
      return 0;
    }
}



template <unsigned int Dim, FEFamily T>
FEContinuity FE<Dim,T>::get_continuity() const
{
  return C_ONE;
}



template <unsigned int Dim, FEFamily T>
bool FE<Dim,T>::is_hierarchic() const
{
  return false;  // FIXME - this will be changed
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
INSTANTIATE_MBRF(0,CLOUGH);
INSTANTIATE_MBRF(1,CLOUGH);
INSTANTIATE_MBRF(2,CLOUGH);
INSTANTIATE_MBRF(3,CLOUGH);

#ifdef LIBMESH_ENABLE_AMR
template void FE<2,CLOUGH>::compute_constraints(DofConstraints&, DofMap&,
                                                const unsigned int,
                                                const Elem*);
template void FE<3,CLOUGH>::compute_constraints(DofConstraints&, DofMap&,
                                                const unsigned int,
                                                const Elem*);
#endif // #ifdef LIBMESH_ENABLE_AMR

} // namespace libMesh

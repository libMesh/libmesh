// $Id$

// The Next Great Finite Element Library.
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
#include "libmesh_config.h"
#ifdef LIBMESH_ENABLE_HIGHER_ORDER_SHAPES

#include "fe.h"
#include "fe_macro.h"
#include "elem.h"




// ------------------------------------------------------------
// Bernstein-specific implementations
// Steffen Petersen 2005
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
    case FIRST:
    case SECOND:
    case THIRD:
    case FOURTH:
    case FIFTH:
    case SIXTH:
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
	return;
      }
    }

  
  // How did we get here?
  libmesh_error();  
  return;
}



template <unsigned int Dim, FEFamily T>
unsigned int FE<Dim,T>::n_dofs(const ElemType t, const Order o)
{

  switch (t)
    {
    case EDGE2:
    case EDGE3:
      return (o+1);
    case QUAD4:
      libmesh_assert(o < 2);
    case QUAD8:
      {
	if (o == 1)
	  return 4;
	else if (o == 2)
	  return 8;
	else
	  libmesh_error();  
      }
    case QUAD9:
      return ((o+1)*(o+1));
    case HEX8:
      libmesh_assert(o < 2);
    case HEX20:
      {
	if (o == 1)
	  return 8;
	else if (o == 2)
	  return 20;
	else
	  libmesh_error();
      }
    case HEX27:
      return ((o+1)*(o+1)*(o+1));
    case TRI3:
      libmesh_assert (o<2);
    case TRI6:
      return ((o+1)*(o+2)/2);
    case TET4:
      libmesh_assert (o<2);
    case TET10:
      {
	libmesh_assert (o<3);
	return ((o+1)*(o+2)*(o+3)/6); 
      }
    default:
      libmesh_error();
    }

  libmesh_error();
  return 0;

  // old code
//   switch (o)
//     {
//       // Bernstein 1st-order polynomials.
//     case FIRST:
//       {
// 	switch (t)
// 	  {
// 	  case EDGE2:
// 	  case EDGE3:
// 	    return 2;

// 	  case TRI6:
// 	    return 3;

// 	  case QUAD8:
// 	  case QUAD9:
// 	    return 4;

// 	  case TET4:
// 	  case TET10:
// 	    return 4;

// 	  case HEX20:
// 	  case HEX27:
// 	    return 8;

// 	  default:
// 	    libmesh_error();
// 	  }
//       }


//       // Bernstein 2nd-order polynomials.
//     case SECOND:
//       {
// 	switch (t)
// 	  {
// 	  case EDGE2:
// 	  case EDGE3:
// 	    return 3;

// 	  case TRI6:
// 	    return 6;

// 	  case QUAD8:
// 	    return 8;

// 	  case QUAD9:
// 	    return 9;

// 	  case TET10:
// 	    return 10;
	    
// 	  case HEX20:
// 	    return 20;

// 	  case HEX27:
// 	    return 27;

// 	  default:
// 	    libmesh_error();
// 	  }
//       }


//       // Bernstein 3rd-order polynomials.
//     case THIRD:
//       {
// 	switch (t)
// 	  {
// 	  case EDGE2:
// 	  case EDGE3:
// 	    return 4;

// 	  case TRI6:
// 	    return 10;
	    
// 	  case QUAD8:
// 	  case QUAD9:
// 	    return 16;

// 	  case TET10:
// 	    return 20;

// 	  case HEX27:
// 	    return 64;

// 	  default:
// 	    libmesh_error();
// 	  }
//       }


//       // Bernstein 4th-order polynomials.
//     case FOURTH:
//       {
// 	switch (t)
// 	  {
// 	  case EDGE2:
// 	  case EDGE3:
// 	    return 5;

// 	  case TRI6:
// 	    return 15;
	    
// 	  case QUAD8:
// 	  case QUAD9:
// 	    return 25;

// 	  case TET10:
// 	    return 35;

// 	  case HEX27:
// 	    return 125;

// 	  default:
// 	    libmesh_error();
// 	  }
//       }


//       // Bernstein 5th-order polynomials.
//     case FIFTH:
//       {
// 	switch (t)
// 	  {
// 	  case EDGE2:
// 	  case EDGE3:
// 	    return 6;

// 	  case TRI6:
// 	    return 21;
	      
// 	  case QUAD8:
// 	  case QUAD9:
// 	    return 36;

// 	  case HEX27:
// 	    return 216;

// 	  default:
// 	    libmesh_error();
// 	  }
//       }

     
//       // Bernstein 6th-order polynomials.
//     case SIXTH:
//       {
// 	switch (t)
// 	  {
// 	  case EDGE2:
// 	  case EDGE3:
// 	    return 7;

// 	  case TRI6:
// 	    return 28;
	      
// 	  case QUAD8:
// 	  case QUAD9:
// 	    return 49;

// 	  case HEX27:
// 	    return 343;

// 	  default:
// 	    libmesh_error();
// 	  }
//       }
      
//     default:
//       {
// 	libmesh_error();
//       }
//     }
  
//   libmesh_error();  
//   return 0;
}



template <unsigned int Dim, FEFamily T>
unsigned int FE<Dim,T>::n_dofs_at_node(const ElemType t,
				       const Order o,
				       const unsigned int n)
{
  switch (t)
    {

    case EDGE2:
    case EDGE3:
      switch (n)
        {
	case 0:
	case 1:
	  return 1;
	case 2:
	  libmesh_assert (o>1);
	  return (o-1);
	default:
	  libmesh_error();		  
	}
    case TRI6:
      switch (n)
	{
	case 0:
	case 1:
	case 2:
	  return 1;

	case 3:
	case 4:
	case 5:
	  return (o-1);
        // Internal DoFs are associated with the elem, not its nodes
	default:
	  libmesh_error();
	}
    case QUAD8:
      libmesh_assert (n<8);
      libmesh_assert (o<3);
    case QUAD9:
      {
	switch (n)
	  {
	  case 0:
	  case 1:
	  case 2:
	  case 3:
	    return 1;
	    
	  case 4:
	  case 5:
	  case 6:
	  case 7:
	    return (o-1);
	    
	    // Internal DoFs are associated with the elem, not its nodes
	  case 8:
	    return 0;
	    
	  default:
	    libmesh_error();		  
	  }
      }
    case HEX8:
      libmesh_assert (n < 8);
      libmesh_assert (o < 2);
    case HEX20:
      libmesh_assert (n < 20);
      libmesh_assert (o < 3);
    case HEX27:
      switch (n)
	{
	case 0:
	case 1:
	case 2:
	case 3:
	case 4:
	case 5:
	case 6:
	case 7:
	  return 1;

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
	  return (o-1);
		  
	case 20:
	case 21:
	case 22:
	case 23:
	case 24:
	case 25:
	  return ((o-1)*(o-1));
	 
        // Internal DoFs are associated with the elem, not its nodes
	case 26:
	  return 0;
	  
	default:
	  libmesh_error();
        }
    case TET4:
      libmesh_assert(n<4);
      libmesh_assert(o<2);
    case TET10:
      libmesh_assert (o<3);
      libmesh_assert (n<10);
      switch (n)
	{
	case 0:
	case 1:
	case 2:
	case 3:
	  return 1;
	  
	case 4:
	case 5:
	case 6:
	case 7:
	case 8:
	case 9:
	  return (o-1);
	  
	default:
	  libmesh_error();
	}
      
    default:
      libmesh_error();
      
    }

  libmesh_error();
  return 0;

//   switch (o)
//     {
//       // The first-order Bernstein shape functions
//     case FIRST:
//       {
// 	switch (t)
// 	  {
// 	    // The 1D Bernstein defined on a three-noded edge
// 	  case EDGE2:
// 	  case EDGE3:
// 	    {
// 	      switch (n)
// 		{
// 		case 0:
// 		case 1:
// 		  return 1;

// 		case 2:
// 		  return 0;
		  
// 		default:
// 		  libmesh_error();		  
// 		}
// 	    }

	    
// 	    // The 2D Bernstein defined on a 6-noded triangle
// 	  case TRI6:
// 	    {
// 	      switch (n)
// 		{
// 		case 0:
// 		case 1:
// 		case 2:
// 		  return 1;

// 		case 3:
// 		case 4:
// 		case 5:
// 		case 6:
// 		  return 0;

// 		default:
// 		  libmesh_error();
// 		}
// 	    }
	    

// 	    // The 2D tensor-product Bernsteins defined on a
// 	    // nine-noded quadrilateral.
// 	  case QUAD8:
// 	  case QUAD9:
// 	    {
// 	      switch (n)
// 		{
// 		case 0:
// 		case 1:
// 		case 2:
// 		case 3:
// 		  return 1;

// 		case 4:
// 		case 5:
// 		case 6:
// 		case 7:
// 		  return 0;
		  
// 		case 8:
// 		  return 0;
		  
// 		default:
// 		  libmesh_error();
// 		}
// 	    }


// 	    // The 3D Bernsteins defined on a ten-noded tetrahedral.
// 	  case TET4:
// 	  case TET10:
// 	    {
// 	      switch (n)
// 		{
// 		case 0:
// 		case 1:
// 		case 2:
// 		case 3:
// 		  return 1;

// 		case 4:
// 		case 5:
// 		case 6:
// 		case 7:		 
// 		case 8:
// 		case 9:
// 		  return 0;
		  
// 		default:
// 		  libmesh_error();
// 		}
// 	    }


// 	    // The 3D Bernsteins defined on a hexahedral.
// 	  case HEX20:
// 	  case HEX27:
// 	    {
// 	      switch (n)
// 		{
// 		case 0:
// 		case 1:
// 		case 2:
// 		case 3:
// 		case 4:
// 		case 5:
// 		case 6:
// 		case 7:
// 		  return 1;

// 		case 8:
// 		case 9:
// 		case 10:
// 		case 11:
// 		case 12:
// 		case 13:
// 		case 14:
// 		case 15:
// 		case 16:
// 		case 17:
// 		case 18:
// 		case 19:
// 		  return 0;
		  
// 		case 20:
// 		case 21:
// 		case 22:
// 		case 23:
// 		case 24:
// 		case 25:
// 		  return 0;
		  
// 		case 26:
// 		  return 0;
// 		}
// 	    }
	    
// 	  default:
// 	    libmesh_error();
	    
// 	  }
//       }



//       // The second-order Bernstein shape functions
//     case SECOND:
//       {
// 	switch (t)
// 	  {
// 	    // The 1D Bernstein defined on a three-noded edge
// 	  case EDGE2:
// 	  case EDGE3:
// 	    {
// 	      libmesh_assert (n<3);
// 	      return 1;
// 	    }

	    
// 	    // The 2D Bernstein defined on a 6-noded triangle
// 	  case TRI6:
// 	    {
// 	      libmesh_assert (n<6);
// 	      return 1;
// 	    }

// 	    // The 8-noded quadrilateral
// 	  case QUAD8:
// 	    {
// 	      libmesh_assert (n<8);
// 	      return 1;
// 	    }
// 	    // The 2D tensor-product Bernsteins defined on a
// 	    // nine-noded quadrilateral.
// 	  case QUAD9:
// 	    {
// 	      libmesh_assert (n<9);
// 	      return 1;
// 	    }


// 	    // The 3D Bernsteins defined on a ten-noded tetrahedral.
// 	  case TET10:
// 	    {
// 	      libmesh_assert(n<10);
// 	      return 1;
// 	    }

// 	    // The 3D Bernsteins defined on a
// 	    // 20-noded hexahedral.
// 	  case HEX20:
// 	    {
// 	      libmesh_assert(n<20);
// 	      return 1;
// 	    }

// 	    // The 3D tensor-product Bernsteins defined on a
// 	    // twenty-seven noded hexahedral.
// 	  case HEX27:
// 	    {
// 	      libmesh_assert(n<27);
// 	      return 1;
// 	    }
	    
// 	  default:
// 	    libmesh_error();
	    
// 	  }
//       }



//       // The third-order Bernstein shape functions
//     case THIRD:
//       {
// 	switch (t)
// 	  {
// 	    // The 1D Bernstein defined on a three-noded edge
// 	  case EDGE2:
// 	  case EDGE3:
// 	    {
// 	      switch (n)
// 		{
// 		case 0:
// 		case 1:
// 		  return 1;

// 		case 2:
// 		  return 2;
		  
// 		default:
// 		  libmesh_error();		  
// 		}
// 	    }

	    
// 	    // The 2D Bernstein defined on a 6-noded triangle
// 	  case TRI6:
// 	    {
// 	      switch (n)
// 		{
// 		case 0:
// 		case 1:
// 		case 2:
// 		  return 1;

// 		case 3:
// 		case 4:
// 		case 5:
// 		  return 2;

// 		case 6:
// 	      return 1;

// 		default:
// 		  libmesh_error();
// 		}
// 	    }


// 	    // The 2D tensor-product Bernsteins defined on a
// 	    // nine-noded quadrilateral.
// 	  case QUAD8:
// 	  case QUAD9:
// 	    {
// 	      switch (n)
// 		{
// 		case 0:
// 		case 1:
// 		case 2:
// 		case 3:
// 		  return 1;

// 		case 4:
// 		case 5:
// 		case 6:
// 		case 7:
// 		  return 2;
		  
// 		case 8:
// 		  return 4;
		  
// 		default:
// 		  libmesh_error();
// 		}
// 	    }


// 	    // The 3D Bernsteins defined on a ten-noded tetrahedral.
// 	  case TET10:
// 	    {
// 	      switch (n)
// 		{
// 		case 0:
// 		case 1:
// 		case 2:
// 		case 3:
// 		  return 1;

// 		case 4:
// 		case 5:
// 		case 6:
// 		case 7:
// 		case 8:
// 		case 9:
// 		  return 2;
		  
// 		default:
// 		  libmesh_error();
// 		}
// 	    }

// 	    // The 3D tensor-product Bernsteins defined on a
// 	    // twenty-seven noded hexahedral.
// 	  case HEX27:
// 	    {
// 	      switch (n)
// 		{
// 		case 0:
// 		case 1:
// 		case 2:
// 		case 3:
// 		case 4:
// 		case 5:
// 		case 6:
// 		case 7:
// 		  return 1;

// 		case 8:
// 		case 9:
// 		case 10:
// 		case 11:
// 		case 12:
// 		case 13:
// 		case 14:
// 		case 15:
// 		case 16:
// 		case 17:
// 		case 18:
// 		case 19:
// 		  return 2;
		  
// 		case 20:
// 		case 21:
// 		case 22:
// 		case 23:
// 		case 24:
// 		case 25:
// 		  return 4;
		  
// 		case 26:
// 		  return 8;
// 		}
// 	    }
	    
// 	  default:
// 	    libmesh_error();
	    
// 	  }
//       }



//       // The fourth-order Bernstein shape functions
//     case FOURTH:
//       {
// 	switch (t)
// 	  {
// 	    // The 1D Bernstein defined on a three-noded edge
// 	  case EDGE2:
// 	  case EDGE3:
// 	    {
// 	      switch (n)
// 		{
// 		case 0:
// 		case 1:
// 		  return 1;

// 		case 2:
// 		  return 3;
		  
// 		default:
// 		  libmesh_error();		  
// 		}
// 	    }

	    
// 	    // The 2D Bernstein defined on a 6-noded triangle
// 	  case TRI6:
// 	    {
// 	      switch (n)
// 		{
// 		case 0:
// 		case 1:
// 		case 2:
// 		  return 1;

// 		case 3:
// 		case 4:
// 		case 5:
// 		  return 3;

// 		case 6:
// 		  return 3;

// 		default:
// 		  libmesh_error();
// 		}
// 	    }


// 	    // The 2D tensor-product Bernsteins defined on a
// 	    // nine-noded quadrilateral.
// 	  case QUAD8:
// 	  case QUAD9:
// 	    {
// 	      switch (n)
// 		{
// 		case 0:
// 		case 1:
// 		case 2:
// 		case 3:
// 		  return 1;

// 		case 4:
// 		case 5:
// 		case 6:
// 		case 7:
// 		  return 3;
		  
// 		case 8:
// 		  return 9;
		  
// 		default:
// 		  libmesh_error();
// 		}
// 	    }

// 	    // The 3D Bernsteins defined on a ten-noded tetrahedral.
// 	  case TET10:
// 	    {
// 	      switch (n)
// 		{
// 		case 0:
// 		case 1:
// 		case 2:
// 		case 3:
// 		  return 1;

// 		case 4:
// 		case 5:
// 		case 6:
// 		case 7:
// 		case 8:
//         case 9:
// 		  return 3;
		  
// 		default:
// 		  libmesh_error();
// 		}
// 	    }

// 	    // The 3D tensor-product Bernsteins defined on a
// 	    // twenty-seven noded hexahedral.
// 	  case HEX27:
// 	    {
// 	      switch (n)
// 		{
// 		case 0:
// 		case 1:
// 		case 2:
// 		case 3:
// 		case 4:
// 		case 5:
// 		case 6:
// 		case 7:
// 		  return 1;

// 		case 8:
// 		case 9:
// 		case 10:
// 		case 11:
// 		case 12:
// 		case 13:
// 		case 14:
// 		case 15:
// 		case 16:
// 		case 17:
// 		case 18:
// 		case 19:
// 		  return 3;
		  
// 		case 20:
// 		case 21:
// 		case 22:
// 		case 23:
// 		case 24:
// 		case 25:
// 		  return 9;
		  
// 		case 26:
// 		  return 27;
// 		}
// 	    }
	    
// 	  default:
// 	    libmesh_error();
	    
// 	  }
//       }



//       // The fifth-order Bernstein shape functions
//     case FIFTH:
//       {
// 	switch (t)
// 	  {
// 	    // The 1D  Bernstein defined on a three-noded edge
// 	  case EDGE2:
// 	  case EDGE3:
// 	    {
// 	      switch (n)
// 		{
// 		case 0:
// 		case 1:
// 		  return 1;

// 		case 2:
// 		  return 4;
		  
// 		default:
// 		  libmesh_error();		  
// 		}
// 	    }

	    
// 	    // The 2D Bernstein defined on a 6-noded triangle
// 	  case TRI6:
// 	    {
// 	      switch (n)
// 		{
// 		case 0:
// 		case 1:
// 		case 2:
// 		  return 1;

// 		case 3:
// 		case 4:
// 		case 5:
// 		  return 4;

// 		case 6:
// 		  return 6;

// 		default:
// 		  libmesh_error();
// 		}
// 	    }


// 	    // The 2D tensor-product Bernsteins defined on a
// 	    // nine-noded quadrilateral.
// 	  case QUAD8:
// 	  case QUAD9:
// 	    {
// 	      switch (n)
// 		{
// 		case 0:
// 		case 1:
// 		case 2:
// 		case 3:
// 		  return 1;

// 		case 4:
// 		case 5:
// 		case 6:
// 		case 7:
// 		  return 4;
		  
// 		case 8:
// 		  return 16;
		  
// 		default:
// 		  libmesh_error();
// 		}
// 	    }

// 	    // The 3D tensor-product Bernsteins defined on a
// 	    // twenty-seven noded hexahedral.
// 	  case HEX27:
// 	    {
// 	      switch (n)
// 		{
// 		case 0:
// 		case 1:
// 		case 2:
// 		case 3:
// 		case 4:
// 		case 5:
// 		case 6:
// 		case 7:
// 		  return 1;

// 		case 8:
// 		case 9:
// 		case 10:
// 		case 11:
// 		case 12:
// 		case 13:
// 		case 14:
// 		case 15:
// 		case 16:
// 		case 17:
// 		case 18:
// 		case 19:
// 		  return 4;
		  
// 		case 20:
// 		case 21:
// 		case 22:
// 		case 23:
// 		case 24:
// 		case 25:
// 		  return 16;
		  
// 		case 26:
// 		  return 64;
// 		}
// 	    }
	    
// 	  default:
// 	    libmesh_error();
	    
// 	  } // switch ElemType
//       } // case FIFTH


//       // The sixth-order Bernstein shape functions
//     case SIXTH:
//       {
// 	switch (t)
// 	  {
// 	    // The 1D Bernstein defined on a three-noded edge
// 	  case EDGE2:
// 	  case EDGE3:
// 	    {
// 	      switch (n)
// 		{
// 		case 0:
// 		case 1:
// 		  return 1;

// 		case 2:
// 		  return 5;
		  
// 		default:
// 		  libmesh_error();		  
// 		}
// 	    }

	    
// 	    // The 2D Bernstein defined on a 6-noded triangle
// 	  case TRI6:
// 	    {
// 	      switch (n)
// 		{
// 		case 0:
// 		case 1:
// 		case 2:
// 		  return 1;

// 		case 3:
// 		case 4:
// 		case 5:
// 		  return 5;

// 		case 6:
// 		  return 10;

// 		default:
// 		  libmesh_error();
// 		}
// 	    }


// 	    // The 2D tensor-product Bernsteins defined on a
// 	    // nine-noded quadrilateral.
// 	  case QUAD8:
// 	  case QUAD9:
// 	    {
// 	      switch (n)
// 		{
// 		case 0:
// 		case 1:
// 		case 2:
// 		case 3:
// 		  return 1;

// 		case 4:
// 		case 5:
// 		case 6:
// 		case 7:
// 		  return 5;
		  
// 		case 8:
// 		  return 25;
		  
// 		default:
// 		  libmesh_error();
// 		}
// 	    }

// 	    // The 3D tensor-product Bernsteins defined on a
// 	    // twenty-seven noded hexahedral.
// 	  case HEX27:
// 	    {
// 	      switch (n)
// 		{
// 		case 0:
// 		case 1:
// 		case 2:
// 		case 3:
// 		case 4:
// 		case 5:
// 		case 6:
// 		case 7:
// 		  return 1;

// 		case 8:
// 		case 9:
// 		case 10:
// 		case 11:
// 		case 12:
// 		case 13:
// 		case 14:
// 		case 15:
// 		case 16:
// 		case 17:
// 		case 18:
// 		case 19:
// 		  return 5;
		  
// 		case 20:
// 		case 21:
// 		case 22:
// 		case 23:
// 		case 24:
// 		case 25:
// 		  return 25;
		  
// 		case 26:
// 		  return 125;
// 		}
// 	    }
	    
// 	  default:
// 	    libmesh_error();
	    
// 	  } // switch ElemType
//       } // case SIXTH

      
//     default:
//       libmesh_error();
//       
//     }
  
//   libmesh_error();
//   return 0;

}



template <unsigned int Dim, FEFamily T>
unsigned int FE<Dim,T>::n_dofs_per_elem(const ElemType t,
					const Order o)
{
  switch (t)
    {
    case EDGE2:
    case EDGE3:
      return 0;
    case TRI3:
    case QUAD4:
      return 0;
    case TRI6:
      return ((o-1)*(o-2)/2);
    case QUAD8:
      if (o <= 2)
	return 0;
    case QUAD9:
      return ((o-1)*(o-1));
    case HEX8:
      libmesh_assert(o < 2);
    case HEX20:
      libmesh_assert(o < 3);
      return 0;
    case HEX27:
      return ((o-1)*(o-1)*(o-1));
    case TET4:
      libmesh_assert (o<2);
    case TET10:
      libmesh_assert (o<3);
	return 0;
    default:
      libmesh_error();
    }

  libmesh_error();
  return 0;

//   switch (o)
//     {
//       // The first-order Bernstein shape functions
//     case FIRST:
//       {
// 	switch (t)
// 	  {
// 	    // The 1D Bernstein defined on a two-noded edge
// 	  case EDGE2:
// 	    return 0;
	    
// 	    // The 1D Bernstein defined on a three-noded edge
// 	  case EDGE3:
// 	    return 0;
	    
// 	    // The 2D Bernstein defined on a 6-noded triangle
// 	  case TRI6:
// 	    return 0;

// 	    // The 2D tensor-product Bernsteins defined on a
// 	    // nine-noded quadrilateral.
// 	  case QUAD8:
// 	    return 0;
	    
// 	    // The 2D tensor-product Bernsteins defined on a
// 	    // nine-noded quadrilateral.
// 	  case QUAD9:
// 	    return 0;
	    
// 	    // The 3D Bernsteins defined on a four-noded hexahedral.
// 	  case TET4:
// 	    return 0;
	    
// 	    // The 3D Bernsteins defined on a ten-noded hexahedral.
// 	  case TET10:
// 	    return 0;
	    
// 	    // The 3D Bernsteins defined on a 20-noded hexahedral.
// 	  case HEX20:
// 	    return 0;

// 	    // The 3D tensor-product Bernsteins defined on a
// 	    // twenty-seven noded hexahedral.
// 	  case HEX27:
// 	    return 0;
	    
	    
// 	  default:
// 	    libmesh_error();
	    
// 	  }
//       }



//       // The second-order Bernstein shape functions
//     case SECOND:
//       {
// 	switch (t)
// 	  {
// 	    // The 1D Bernstein defined on a two-noded edge
// 	  case EDGE2:
// 	    return 1;
	    
// 	    // The 1D Bernstein defined on a three-noded edge
// 	  case EDGE3:
// 	    return 0;
	    
// 	    // The 2D Bernstein defined on a 6-noded triangle
// 	  case TRI6:
// 	    return 0;

// 	    // The 2D tensor-product Bernsteins defined on a
// 	    // eight-noded quadrilateral.
// 	  case QUAD8:
// 	    return 0;

// 	    // The 2D tensor-product Bernsteins defined on a
// 	    // nine-noded quadrilateral.
// 	  case QUAD9:
// 	    return 0;
	    	    
// 	    // The 3D Bernsteins defined on a ten-noded hexahedral.
// 	  case TET10:
// 	    return 0;
	    	    
// 	    // The 3D Bernsteins defined on a 20-noded hexahedral.
// 	  case HEX20:
// 	    return 0;

// 	    // The 3D tensor-product Bernsteins defined on a
// 	    // twenty-seven noded hexahedral.
// 	  case HEX27:
// 	    return 0;


	    
// 	  default:
// 	    libmesh_error();
	    
// 	  }
//       }



//       // The third-order Bernstein shape functions
//     case THIRD:
//       {
// 	switch (t)
// 	  {
// 	    // The 1D Bernstein defined on a two-noded edge
// 	  case EDGE2:
// 	    return 2;
	    
// 	    // The 1D Bernstein defined on a three-noded edge
// 	  case EDGE3:
// 	    return 0;
	    
// 	    // The 2D Bernstein defined on a 6-noded triangle
// 	  case TRI6:
// 	    return 1;
	    
// 	    // The 2D tensor-product Bernsteins defined on a
// 	    // eight-noded quadrilateral.
// 	  case QUAD8:
// 	    return 4;

// 	    // The 2D tensor-product Bernsteins defined on a
// 	    // nine-noded quadrilateral.
// 	  case QUAD9:
// 	    return 0;
	    	    
// 	    // The 3D Bernsteins defined on a ten-noded hexahedral.
// 	  case TET10:
// 	    return 4;
	    
// 	    // The 3D tensor-product Bernsteins defined on a
// 	    // twenty-seven noded hexahedral.
// 	  case HEX27:
// 	    return 0;


	    
// 	  default:
// 	    libmesh_error();
	    
// 	  }
//       }



//       // The fourth-order Bernstein shape functions
//     case FOURTH:
//       {
// 	switch (t)
// 	  {
// 	    // The 1D Bernstein defined on a two-noded edge
// 	  case EDGE2:
// 	    return 3;
	    
// 	    // The 1D Bernstein defined on a three-noded edge
// 	  case EDGE3:
// 	    return 0;
	    
// 	    // The 2D Bernstein defined on a 6-noded triangle
// 	  case TRI6:
// 	    return 3;
	    
// 	    // The 2D tensor-product Bernsteins defined on a
// 	    // eight-noded quadrilateral.
// 	  case QUAD8:
// 	    return 9;

// 	    // The 2D tensor-product Bernsteins defined on a
// 	    // nine-noded quadrilateral.
// 	  case QUAD9:
// 	    return 0;
	    
// 	    // The 3D Bernsteins defined on a ten-noded hexahedral.
// 	  case TET10:
// 	    return 13;
	               	    
// 	    // The 3D tensor-product Bernsteins defined on a
// 	    // twenty-seven noded hexahedral.
// 	  case HEX27:
// 	    return 0;


	    
// 	  default:
// 	    libmesh_error();
	    
// 	  }
//       }



//       // The fifth-order Bernstein shape functions
//     case FIFTH:
//       {
// 	switch (t)
// 	  {
// 	    // The 1D Bernstein defined on a two-noded edge
// 	  case EDGE2:
// 	    return 4;
	    
// 	    // The 1D Bernstein defined on a three-noded edge
// 	  case EDGE3:
// 	    return 0;
	    
// 	    // The 2D Bernstein defined on a 6-noded triangle
// 	  case TRI6:
// 	    return 6;
	    
// 	    // The 2D tensor-product Bernsteins defined on a
// 	    // eight-noded quadrilateral.
// 	  case QUAD8:
// 	    return 16;

// 	    // The 2D tensor-product Bernsteins defined on a
// 	    // nine-noded quadrilateral.
// 	  case QUAD9:
// 	    return 0;
	    
// 	    // The 3D tensor-product Bernsteins defined on a
// 	    // twenty-seven noded hexahedral.
// 	  case HEX27:
// 	    return 0;


	    
// 	  default:
// 	    libmesh_error();
	    
// 	  }
//       }



//       // The sixth-order Bernstein shape functions
//     case SIXTH:
//       {
// 	switch (t)
// 	  {
// 	    // The 1D Bernstein defined on a two-noded edge
// 	  case EDGE2:
// 	    return 5;
	    
// 	    // The 1D Bernstein defined on a three-noded edge
// 	  case EDGE3:
// 	    return 0;
	    
// 	    // The 2D Bernstein defined on a 6-noded triangle
// 	  case TRI6:
// 	    return 10;
	    
// 	    // The 2D tensor-product Bernsteins defined on a
// 	    // eight-noded quadrilateral.
// 	  case QUAD8:
// 	    return 25;

// 	    // The 2D tensor-product Bernsteins defined on a
// 	    // nine-noded quadrilateral.
// 	  case QUAD9:
// 	    return 0;
	    
// 	    // The 3D tensor-product Bernsteins defined on a
// 	    // twenty-seven noded hexahedral.
// 	  case HEX27:
// 	    return 0;


	    
// 	  default:
// 	    libmesh_error();
	    
// 	  }
//       }


      
//       // Otherwise no DOFS per element
//     default:
//       return 0;
//     }

}

template <unsigned int Dim, FEFamily T>
FEContinuity FE<Dim,T>::get_continuity() const
{
  return C_ZERO;
}


template <unsigned int Dim, FEFamily T>
bool FE<Dim,T>::is_hierarchic() const
{
  return false;
}


template <unsigned int Dim, FEFamily T>
bool FE<Dim,T>::shapes_need_reinit() const
{
  // reinit is only necessary for
  // approximation orders >= 3,
  // but we might reach that with p refinement
//  if(this->fe_type.order == FIRST ||
//     this->fe_type.order == SECOND)
//    return false;
//  else
    return true;
}




//--------------------------------------------------------------
// Explicit instantiations for member functions
INSTANTIATE_MBRF(1,BERNSTEIN);
INSTANTIATE_MBRF(2,BERNSTEIN);
INSTANTIATE_MBRF(3,BERNSTEIN);

#endif //LIBMESH_ENABLE_HIGHER_ORDER_SHAPES

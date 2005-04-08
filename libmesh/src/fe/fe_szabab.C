// $Id: fe_szabab.C,v 1.4 2005-04-08 09:31:02 spetersen Exp $

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



// Local includes
#include "libmesh_config.h"
#ifdef ENABLE_HIGHER_ORDER_SHAPES

// Local includes
#include "fe.h"
#include "elem.h"



// ------------------------------------------------------------
// Szabo-Babuska-specific implementations
// Steffen Petersen 2004
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
    case FIFTH:
    case SIXTH:
    case SEVENTH:
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
	return;
      }
    }

  
  // We should never get here?
  error();  
  return;
}



template <unsigned int Dim, FEFamily T>
unsigned int FE<Dim,T>::n_dofs(const ElemType t, const Order o)
{
  switch (o)
    {
      // Szabo-Babuska 1st-order polynomials.
    case FIRST:
      {
	switch (t)
	  {
	  case EDGE2:
	  case EDGE3:
	    return 2;

	  case TRI6:
	    return 3;

	  case QUAD8:
	  case QUAD9:
	    return 4;

	  default:
	    error();
	  }
      }


      // Szabo-Babuska 2nd-order polynomials.
    case SECOND:
      {
	switch (t)
	  {
	  case EDGE2:
	  case EDGE3:
	    return 3;

	  case TRI6:
	    return 6;

	  case QUAD8:
	  case QUAD9:
	    return 9;
	    
	  default:
	    error();
	  }
      }


      // Szabo-Babuska 3rd-order polynomials.
    case THIRD:
      {
	switch (t)
	  {
	  case EDGE2:
	  case EDGE3:
	    return 4;

	  case TRI6:
	    return 10;
	    
	  case QUAD8:
	  case QUAD9:
	    return 16;

	  default:
	    error();
	  }
      }


      // Szabo-Babuska 4th-order polynomials.
    case FOURTH:
      {
	switch (t)
	  {
	  case EDGE2:
	  case EDGE3:
	    return 5;

	  case TRI6:
	    return 15;
	    
	  case QUAD8:
	  case QUAD9:
	    return 25;

	  default:
	    error();
	  }
      }


      // Szabo-Babuska 5th-order polynomials.
    case FIFTH:
      {
	switch (t)
	  {
	  case EDGE2:
	  case EDGE3:
	    return 6;

	  case TRI6:
	    return 21;
	      
	  case QUAD8:
	  case QUAD9:
	    return 36;

	  default:
	    error();
	  }
      }

      
      // Szabo-Babuska 6th-order polynomials.
    case SIXTH:
      {
	switch (t)
	  {
	  case EDGE2:
	  case EDGE3:
	    return 7;

	  case TRI6:
	    return 28;
	      
	  case QUAD8:
	  case QUAD9:
	    return 49;

	  default:
	    error();
	  }
      }

      // Szabo-Babuska 7th-order polynomials.
    case SEVENTH:
      {
	switch (t)
	  {
	  case EDGE2:
	  case EDGE3:
	    return 8;

	  case TRI6:
	    return 36;
	      
	  case QUAD8:
	  case QUAD9:
	    return 64;

	  default:
	    error();
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
unsigned int FE<Dim,T>::n_dofs_at_node(const ElemType t,
				       const Order o,
				       const unsigned int n)
{
  switch (o)
    {
      // The first-order Szabo-Babuska shape functions
    case FIRST:
      {
	switch (t)
	  {
	    // The 1D Szabo-Babuska defined on a three-noded edge
	  case EDGE2:
	  case EDGE3:
	    {
	      switch (n)
		{
		case 0:
		case 1:
		  return 1;

		case 2:
		  return 0;
		  
		default:
		  error();		  
		}
	    }

	    
	    // The 2D Szabo-Babuska defined on a 6-noded triangle
	  case TRI6:
	    {
	      switch (n)
		{
		case 0:
		case 1:
		case 2:
		  return 1;

		case 3:
		case 4:
		case 5:
		  return 0;

		default:
		  error();
		}
	    }
	    

	    // The 2D tensor-product Szabo-Babuska defined on a
	    // nine-noded quadrilateral.
	  case QUAD8:
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
		  return 0;
		  
		case 8:
		  return 0;
		  
		default:
		  error();
		}
	    }


	  default:
	    error();
	    
	  }
      }



      // The second-order Szabo-Babuska shape functions
    case SECOND:
      {
	switch (t)
	  {
	    // The 1D Szabo-Babuska defined on a three-noded edge
	  case EDGE2:
	  case EDGE3:
	    {
	      switch (n)
		{
		case 0:
		case 1:
		  return 1;

		case 2:
		  return 1;
		  
		default:
		  error();		  
		}
	    }

	    
	    // The 2D Szabo-Babuska defined on a 6-noded triangle
	  case TRI6:
	    {
	      switch (n)
		{
		case 0:
		case 1:
		case 2:
		  return 1;

		case 3:
		case 4:
		case 5:
		  return 1;

		default:
		  error();
		}
	    }


	    // The 2D tensor-product Szabo-Babuska defined on a
	    // nine-noded quadrilateral.
	  case QUAD8:
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
		  return 1;
		  
		case 8:
		  return 1;
		  
		default:
		  error();
		}
	    }

	    
	  default:
	    error();
	    
	  }
      }



      // The third-order Szabo-Babuska shape functions
    case THIRD:
      {
	switch (t)
	  {
	    // The 1D Szabo-Babuska defined on a three-noded edge
	  case EDGE2:
	  case EDGE3:
	    {
	      switch (n)
		{
		case 0:
		case 1:
		  return 1;

		case 2:
		  return 2;
		  
		default:
		  error();		  
		}
	    }

	    
	    // The 2D Szabo-Babuska defined on a 6-noded triangle
	  case TRI6:
	    {
	      switch (n)
		{
		case 0:
		case 1:
		case 2:
		  return 1;

		case 3:
		case 4:
		case 5:
		  return 2;

		default:
		  error();
		}
	    }


	    // The 2D tensor-product Szabo-Babuska defined on a
	    // nine-noded quadrilateral.
	  case QUAD8:
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
		  return 2;
		  
		case 8:
		  return 4;
		  
		default:
		  error();
		}
	    }


	  default:
	    error();
	    
	  }
      }



      // The fourth-order Szabo-Babuska shape functions
    case FOURTH:
      {
	switch (t)
	  {
	    // The 1D Szabo-Babuska defined on a three-noded edge
	  case EDGE2:
	  case EDGE3:
	    {
	      switch (n)
		{
		case 0:
		case 1:
		  return 1;

		case 2:
		  return 3;
		  
		default:
		  error();		  
		}
	    }

	    
	    // The 2D Szabo-Babuska defined on a 6-noded triangle
	  case TRI6:
	    {
	      switch (n)
		{
		case 0:
		case 1:
		case 2:
		  return 1;

		case 3:
		case 4:
		case 5:
		  return 3;

		default:
		  error();
		}
	    }


	    // The 2D tensor-product Szabo-Babuska defined on a
	    // nine-noded quadrilateral.
	  case QUAD8:
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
		  return 3;
		  
		case 8:
		  return 9;
		  
		default:
		  error();
		}
	    }


	  default:
	    error();
	    
	  }
      }



      // The fifth-order Szabo-Babuska shape functions
    case FIFTH:
      {
	switch (t)
	  {
	    // The 1D Szabo-Babuska defined on a three-noded edge
	  case EDGE2:
	  case EDGE3:
	    {
	      switch (n)
		{
		case 0:
		case 1:
		  return 1;

		case 2:
		  return 4;
		  
		default:
		  error();		  
		}
	    }

	    
	    // The 2D Szabo-Babuska defined on a 6-noded triangle
	  case TRI6:
	    {
	      switch (n)
		{
		case 0:
		case 1:
		case 2:
		  return 1;

		case 3:
		case 4:
		case 5:
		  return 4;

		default:
		  error();
		}
	    }


	    // The 2D tensor-product Szabo-Babuska defined on a
	    // nine-noded quadrilateral.
	  case QUAD8:
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
		  return 4;
		  
		case 8:
		  return 16;
		  
		default:
		  error();
		}
	    }


	  default:
	    error();
	    
	  }
      }



      // The sixth-order Szabo-Babuska shape functions
    case SIXTH:
      {
	switch (t)
	  {
	    // The 1D Szabo-Babuska defined on a three-noded edge
	  case EDGE2:
	  case EDGE3:
	    {
	      switch (n)
		{
		case 0:
		case 1:
		  return 1;

		case 2:
		  return 5;
		  
		default:
		  error();		  
		}
	    }

	    
	    // The 2D Szabo-Babuska defined on a 6-noded triangle
	  case TRI6:
	    {
	      switch (n)
		{
		case 0:
		case 1:
		case 2:
		  return 1;

		case 3:
		case 4:
		case 5:
		  return 5;

		default:
		  error();
		}
	    }


	    // The 2D tensor-product Szabo-Babuska defined on a
	    // nine-noded quadrilateral.
	  case QUAD8:
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
		  return 5;
		  
		case 8:
		  return 25;
		  
		default:
		  error();
		}
	    }


	  default:
	    error();
	    
	  }
      }


      // The seventh-order Szabo-Babuska shape functions
    case SEVENTH:
      {
	switch (t)
	  {
	    // The 1D Szabo-Babuska defined on a three-noded edge
	  case EDGE2:
	  case EDGE3:
	    {
	      switch (n)
		{
		case 0:
		case 1:
		  return 1;

		case 2:
		  return 6;
		  
		default:
		  error();		  
		}
	    }

	    
	    // The 2D Szabo-Babuska defined on a 6-noded triangle
	  case TRI6:
	    {
	      switch (n)
		{
		case 0:
		case 1:
		case 2:
		  return 1;

		case 3:
		case 4:
		case 5:
		  return 6;

		default:
		  error();
		}
	    }


	    // The 2D tensor-product Szabo-Babuska defined on a
	    // nine-noded quadrilateral.
	  case QUAD8:
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
		  return 6;
		  
		case 8:
		  return 36;
		  
		default:
		  error();
		}
	    }


	  default:
	    error();
	    
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
unsigned int FE<Dim,T>::n_dofs_per_elem(const ElemType t,
					const Order o)
{
  switch (o)
    {
      // The first-order Szabo-Babuska shape functions
    case FIRST:
      {
	switch (t)
	  {
	    // The 1D Szabo-Babuska defined on a two-noded edge
	  case EDGE2:
	    return 0;
	    
	    // The 1D Szabo-Babuska defined on a three-noded edge
	  case EDGE3:
	    return 0;
	    
	    // The 2D Szabo-Babuska defined on a 6-noded triangle
	  case TRI6:
	    return 0;

	    // The 2D tensor-product Szabo-Babuska defined on a
	    // nine-noded quadrilateral.
	  case QUAD8:
	    return 0;
	    
	    // The 2D tensor-product Szabo-Babuska defined on a
	    // nine-noded quadrilateral.
	    
	  case QUAD9:
	    return 0;
	    
	    
	  default:
	    error();
	    
	  }
      }



      // The second-order Szabo-Babuska shape functions
    case SECOND:
      {
	switch (t)
	  {
	    // The 1D Szabo-Babuska defined on a two-noded edge
	  case EDGE2:
	    return 1;
	    
	    // The 1D Szabo-Babuska defined on a three-noded edge
	  case EDGE3:
	    return 0;
	    
	    // The 2D Szabo-Babuska defined on a 6-noded triangle
	  case TRI6:
	    return 0;

	    // The 2D tensor-product Szabo-Babuska defined on a
	    // eight-noded quadrilateral.
	  case QUAD8:
	    return 1;

	    // The 2D tensor-product Szabo-Babuska defined on a
	    // nine-noded quadrilateral.
	  case QUAD9:
	    return 0;

	    
	  default:
	    error();
	    
	  }
      }



      // The third-order Szabo-Babuska shape functions
    case THIRD:
      {
	switch (t)
	  {
	    // The 1D Szabo-Babuska defined on a two-noded edge
	  case EDGE2:
	    return 2;
	    
	    // The 1D Szabo-Babuska defined on a three-noded edge
	  case EDGE3:
	    return 0;
	    
	    // The 2D Szabo-Babuska defined on a 6-noded triangle
	  case TRI6:
	    return 1;

	    // The 2D tensor-product Szabo-Babuska defined on a
	    // eight-noded quadrilateral.
	  case QUAD8:
	    return 4;

	    // The 2D tensor-product Szabo-Babuska defined on a
	    // nine-noded quadrilateral.
	  case QUAD9:
	    return 0;

	    
	  default:
	    error();
	    
	  }
      }



      // The fourth-order Szabo-Babuska shape functions
    case FOURTH:
      {
	switch (t)
	  {
	    // The 1D Szabo-Babuska defined on a two-noded edge
	  case EDGE2:
	    return 3;
	    
	    // The 1D Szabo-Babuska defined on a three-noded edge
	  case EDGE3:
	    return 0;
	    
	    // The 2D Szabo-Babuska defined on a 6-noded triangle
	  case TRI6:
	    return 3;

	    // The 2D tensor-product Szabo-Babuska defined on a
	    // eight-noded quadrilateral.
	  case QUAD8:
	    return 9;

	    // The 2D tensor-product Szabo-Babuska defined on a
	    // nine-noded quadrilateral.
	  case QUAD9:
	    return 0;

	    
	  default:
	    error();
	    
	  }
      }



      // The fifth-order Szabo-Babuska shape functions
    case FIFTH:
      {
	switch (t)
	  {
	    // The 1D Szabo-Babuska defined on a two-noded edge
	  case EDGE2:
	    return 4;
	    
	    // The 1D Szabo-Babuska defined on a three-noded edge
	  case EDGE3:
	    return 0;
	    
	    // The 2D Szabo-Babuska defined on a 6-noded triangle
	  case TRI6:
	    return 6;

	    // The 2D tensor-product Szabo-Babuska defined on a
	    // eight-noded quadrilateral.
	  case QUAD8:
	    return 16;

	    // The 2D tensor-product Szabo-Babuska defined on a
	    // nine-noded quadrilateral.
	  case QUAD9:
	    return 0;

	    
	  default:
	    error();
	    
	  }
      }


  // The sixth-order Szabo-Babuska shape functions
    case SIXTH:
      {
	switch (t)
	  {
	    // The 1D Szabo-Babuska defined on a two-noded edge
	  case EDGE2:
	    return 5;
	    
	    // The 1D Szabo-Babuska defined on a three-noded edge
	  case EDGE3:
	    return 0;
	    
	    // The 2D Szabo-Babuska defined on a 6-noded triangle
	  case TRI6:
	    return 10;

	    // The 2D tensor-product Szabo-Babuska defined on a
	    // eight-noded quadrilateral.
	  case QUAD8:
	    return 25;

	    // The 2D tensor-product Szabo-Babuska defined on a
	    // nine-noded quadrilateral.
	  case QUAD9:
	    return 0;

	    
	  default:
	    error();
	    
	  }
      }


  // The seventh-order Szabo-Babuska shape functions
    case SEVENTH:
      {
	switch (t)
	  {
	    // The 1D Szabo-Babuska defined on a two-noded edge
	  case EDGE2:
	    return 6;
	    
	    // The 1D Szabo-Babuska defined on a three-noded edge
	  case EDGE3:
	    return 0;
	    
	    // The 2D Szabo-Babuska defined on a 6-noded triangle
	  case TRI6:
	    return 15;

	    // The 2D tensor-product Szabo-Babuska defined on a
	    // eight-noded quadrilateral.
	  case QUAD8:
	    return 36;

	    // The 2D tensor-product Szabo-Babuska defined on a
	    // nine-noded quadrilateral.
	  case QUAD9:
	    return 0;

	    
	  default:
	    error();
	    
	  }
      }

      
      // Otherwise no DOFS per element
    default:
      return 0;
    }
}



template <unsigned int Dim, FEFamily T>
bool FE<Dim,T>::shapes_need_reinit() const
{
  // reinit is only necessary for
  // approximation orders >= 3
  if(this->fe_type.order == FIRST ||
     this->fe_type.order == SECOND)
    return false;
  else
    return true;
}




//--------------------------------------------------------------
// Explicit instantiations
template class FE<1,SZABAB>;
template class FE<2,SZABAB>;
template class FE<3,SZABAB>;


#endif //ENABLE_HIGHER_ORDER_SHAPES

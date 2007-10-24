// $Id: fe_hierarchic_shape_3D.C,v 1.12 2005-02-22 22:17:36 jwpeterson Exp $

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

// Local includes
#include "fe.h"
#include "elem.h"




template <>
Real FE<3,HIERARCHIC>::shape(const ElemType,
			     const Order,
			     const unsigned int,
			     const Point&)
{
  std::cerr << "Hierarchic polynomials require the element type\n"
	    << "because edge and face orientation is needed."
	    << std::endl;
  
  error();
  return 0.;
}



template <>
Real FE<3,HIERARCHIC>::shape(const Elem* elem,
			     const Order order,
			     const unsigned int i,
			     const Point& p)
{
#if DIM == 3
  
  assert (elem != NULL);
  const ElemType type = elem->type();
  
  switch (order)
    {

      

      // 1st & 2nd-order Hierarchics.
    case FIRST:
    case SECOND:
      {
	switch (type)
	  {
	    // Hierarchic shape functions on the hexahedral.
	  case HEX27:
	    {
	      assert (i<27);
	
	      // Compute hex shape functions as a tensor-product
	      const Real xi   = p(0);
	      const Real eta  = p(1);
	      const Real zeta = p(2);
	
	      // The only way to make any sense of this
	      // is to look at the mgflo/mg2/mgf documentation
	      // and make the cut-out cube!
	      //                                0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26
	      static const unsigned int i0[] = {0, 1, 1, 0, 0, 1, 1, 0, 2, 1, 2, 0, 0, 1, 1, 0, 2, 1, 2, 0, 2, 2, 1, 2, 0, 2, 2};
	      static const unsigned int i1[] = {0, 0, 1, 1, 0, 0, 1, 1, 0, 2, 1, 2, 0, 0, 1, 1, 0, 2, 1, 2, 2, 0, 2, 1, 2, 2, 2};
	      static const unsigned int i2[] = {0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 2, 2, 2, 2, 1, 1, 1, 1, 0, 2, 2, 2, 2, 1, 2};
	
	      return (FE<1,HIERARCHIC>::shape(EDGE3, order, i0[i], xi)*
		      FE<1,HIERARCHIC>::shape(EDGE3, order, i1[i], eta)*
		      FE<1,HIERARCHIC>::shape(EDGE3, order, i2[i], zeta));
	    }

	    
	  default:
	    error();
	  }
	
      }

      

      // 3rd-order Hierarchics.
    case THIRD:
      {
	switch (type)
	  {
	    // Hierarchic shape functions on the hexahedral.
	  case HEX27:
	    {
	      assert (i<64);
	
	      // Compute hex shape functions as a tensor-product
	      const Real xi    = p(0);
	      const Real eta   = p(1);
	      const Real zeta  = p(2);
	      Real xi_mapped   = p(0);
	      Real eta_mapped  = p(1);
	      Real zeta_mapped = p(2);
	
	      // The only way to make any sense of this
	      // is to look at the mgflo/mg2/mgf documentation
	      // and make the cut-out cube!
	      //  Nodes                         0  1  2  3  4  5  6  7  8  8  9  9 10 10 11 11 12 12 13 13 14 14 15 15 16 16 17 17 18 18 19 19 20 20 20 20 21 21 21 21 22 22 22 22 23 23 23 23 24 24 24 24 25 25 25 25 26 26 26 26 26 26 26 26 
	      //  DOFS                          0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 18 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 60 62 63
	      static const unsigned int i0[] = {0, 1, 1, 0, 0, 1, 1, 0, 2, 3, 1, 1, 2, 3, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 2, 3, 1, 1, 2, 3, 0, 0, 2, 3, 2, 3, 2, 3, 2, 3, 1, 1, 1, 1, 2, 3, 2, 3, 0, 0, 0, 0, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3};
	      static const unsigned int i1[] = {0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 2, 3, 1, 1, 2, 3, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 2, 3, 1, 1, 2, 3, 2, 2, 3, 3, 0, 0, 0, 0, 2, 3, 2, 3, 1, 1, 1, 1, 2, 3, 2, 3, 2, 2, 3, 3, 2, 2, 3, 3, 2, 2, 3, 3};
	      static const unsigned int i2[] = {0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 2, 3, 2, 3, 2, 3, 2, 3, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 2, 2, 3, 3, 2, 2, 3, 3, 2, 2, 3, 3, 2, 2, 3, 3, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3};



	      // handle the edge orientation
	      {
		// Edge 0
		if ((i1[i] == 0) && (i2[i] == 0))
		  {
		    if (elem->node(0) != std::min(elem->node(0), elem->node(1)))
		      xi_mapped = -xi;
		  }
		// Edge 1
		else if ((i0[i] == 1) && (i2[i] == 0))
		  {
		    if (elem->node(1) != std::min(elem->node(1), elem->node(2)))
		      eta_mapped = -eta;
		  }
		// Edge 2
		else if ((i1[i] == 1) && (i2[i] == 0))
		  {
		    if (elem->node(3) != std::min(elem->node(3), elem->node(2)))
		      xi_mapped = -xi;
		  }
		// Edge 3
		else if ((i0[i] == 0) && (i2[i] == 0))
		  {
		    if (elem->node(0) != std::min(elem->node(0), elem->node(3)))
		      eta_mapped = -eta;
		  }
		// Edge 4
		else if ((i0[i] == 0) && (i1[i] == 0))
		  {
		    if (elem->node(0) != std::min(elem->node(0), elem->node(4)))
		      zeta_mapped = -zeta;
		  }		
		// Edge 5
		else if ((i0[i] == 1) && (i1[i] == 0))
		  {
		    if (elem->node(1) != std::min(elem->node(1), elem->node(5)))
		      zeta_mapped = -zeta;
		  }		
		// Edge 6
		else if ((i0[i] == 1) && (i1[i] == 1))
		  {
		    if (elem->node(2) != std::min(elem->node(2), elem->node(6)))
		      zeta_mapped = -zeta;
		  }
		// Edge 7
		else if ((i0[i] == 0) && (i1[i] == 1))
		  {
		    if (elem->node(3) != std::min(elem->node(3), elem->node(7)))
		      zeta_mapped = -zeta;
		  }		
		// Edge 8
		else if ((i1[i] == 0) && (i2[i] == 1))
		  {
		    if (elem->node(4) != std::min(elem->node(4), elem->node(5)))
		      xi_mapped = -xi;
		  }
		// Edge 9
		else if ((i0[i] == 1) && (i2[i] == 1))
		  {
		    if (elem->node(5) != std::min(elem->node(5), elem->node(6)))
		      eta_mapped = -eta;
		  }		
		// Edge 10
		else if ((i1[i] == 1) && (i2[i] == 1))
		  {
		    if (elem->node(7) != std::min(elem->node(7), elem->node(6)))
		      xi_mapped = -xi;
		  }
		// Edge 11
		else if ((i0[i] == 0) && (i2[i] == 1))
		  {
		    if (elem->node(4) != std::min(elem->node(4), elem->node(7)))
		      eta_mapped = -eta;
		  }
	      }


	      // handle the face orientation
	      {
		// Face 0
		if (     (i2[i] == 0) && (i0[i] >= 2) && (i1[i] >= 2))
		  {
		    const unsigned int min_node = std::min(elem->node(1),
							   std::min(elem->node(2),
								    std::min(elem->node(0),
									     elem->node(3))));
		    if (elem->node(0) == min_node)
		      if (elem->node(1) == std::min(elem->node(1), elem->node(3)))
			{
			  // Case 1
			  xi_mapped  = xi;
			  eta_mapped = eta;
			}
		      else
			{
			  // Case 2
			  xi_mapped  = eta;
			  eta_mapped = xi;
			}

		    else if (elem->node(3) == min_node)
		      if (elem->node(0) == std::min(elem->node(0), elem->node(2)))
			{
			  // Case 3
			  xi_mapped  = -eta;
			  eta_mapped = xi;
			}
		      else
			{
			  // Case 4
			  xi_mapped  = xi;
			  eta_mapped = -eta;
			}

		    else if (elem->node(2) == min_node)
		      if (elem->node(3) == std::min(elem->node(3), elem->node(1)))
			{
			  // Case 5
			  xi_mapped  = -xi;
			  eta_mapped = -eta;
			}
		      else
			{
			  // Case 6
			  xi_mapped  = -eta;
			  eta_mapped = -xi;
			}

		    else if (elem->node(1) == min_node)
		      if (elem->node(2) == std::min(elem->node(2), elem->node(0)))
			{
			  // Case 7
			  xi_mapped  = eta;
			  eta_mapped = -xi;
			}
		      else
			{
			  // Case 8
			  xi_mapped  = -xi;
			  eta_mapped = eta;
			}
		  }

		
		// Face 1
		else if ((i1[i] == 0) && (i0[i] >= 2) && (i2[i] >= 2))
		  {
		    const unsigned int min_node = std::min(elem->node(0),
							   std::min(elem->node(1),
								    std::min(elem->node(5),
									     elem->node(4))));
		    if (elem->node(0) == min_node)
		      if (elem->node(1) == std::min(elem->node(1), elem->node(4)))
			{
			  // Case 1
			  xi_mapped   = xi;
			  zeta_mapped = zeta;
			}
		      else
			{
			  // Case 2
			  xi_mapped   = zeta;
			  zeta_mapped = xi;
			}

		    else if (elem->node(1) == min_node)
		      if (elem->node(5) == std::min(elem->node(5), elem->node(0)))
			{
			  // Case 3
			  xi_mapped   = zeta;
			  zeta_mapped = -xi;
			}
		      else
			{
			  // Case 4
			  xi_mapped   = -xi;
			  zeta_mapped = zeta;
			}

		    else if (elem->node(5) == min_node)
		      if (elem->node(4) == std::min(elem->node(4), elem->node(1)))
			{
			  // Case 5
			  xi_mapped   = -xi;
			  zeta_mapped = -zeta;
			}
		      else
			{
			  // Case 6
			  xi_mapped   = -zeta;
			  zeta_mapped = -xi;
			}

		    else if (elem->node(4) == min_node)
		      if (elem->node(0) == std::min(elem->node(0), elem->node(5)))
			{
			  // Case 7
			  xi_mapped   = -xi;
			  zeta_mapped = zeta;
			}
		      else
			{
			  // Case 8
			  xi_mapped   = xi;
			  zeta_mapped = -zeta;
			}
		  }

		
		// Face 2
		else if ((i0[i] == 1) && (i1[i] >= 2) && (i2[i] >= 2))
		  {
		    const unsigned int min_node = std::min(elem->node(1),
							   std::min(elem->node(2),
								    std::min(elem->node(6),
									     elem->node(5))));
		    if (elem->node(1) == min_node)
		      if (elem->node(2) == std::min(elem->node(2), elem->node(5)))
			{
			  // Case 1
			  eta_mapped  = eta;
			  zeta_mapped = zeta;
			}
		      else
			{
			  // Case 2
			  eta_mapped  = zeta;
			  zeta_mapped = eta;
			}

		    else if (elem->node(2) == min_node)
		      if (elem->node(6) == std::min(elem->node(6), elem->node(1)))
			{
			  // Case 3
			  eta_mapped  = zeta;
			  zeta_mapped = -eta;
			}
		      else
			{
			  // Case 4
			  eta_mapped  = -eta;
			  zeta_mapped = zeta;
			}

		    else if (elem->node(6) == min_node)
		      if (elem->node(5) == std::min(elem->node(5), elem->node(2)))
			{
			  // Case 5
			  eta_mapped  = -eta;
			  zeta_mapped = -zeta;
			}
		      else
			{
			  // Case 6
			  eta_mapped  = -zeta;
			  zeta_mapped = -eta;
			}

		    else if (elem->node(5) == min_node)
		      if (elem->node(1) == std::min(elem->node(1), elem->node(6)))
			{
			  // Case 7
			  eta_mapped  = -zeta;
			  zeta_mapped = eta;
			}
		      else
			{
			  // Case 8
			  eta_mapped   = eta;
			  zeta_mapped = -zeta;
			}
		  }

		
		// Face 3
		else if ((i1[i] == 1) && (i0[i] >= 2) && (i2[i] >= 2))
		  {
 		    const unsigned int min_node = std::min(elem->node(2),
							   std::min(elem->node(3),
								    std::min(elem->node(7),
									     elem->node(6))));
		    if (elem->node(3) == min_node)
		      if (elem->node(2) == std::min(elem->node(2), elem->node(7)))
			{
			  // Case 1
			  xi_mapped   = xi;
			  zeta_mapped = zeta;
			}
		      else
			{
			  // Case 2
			  xi_mapped   = zeta;
			  zeta_mapped = xi;
			}

		    else if (elem->node(7) == min_node)
		      if (elem->node(3) == std::min(elem->node(3), elem->node(6)))
			{
			  // Case 3
			  xi_mapped   = -zeta;
			  zeta_mapped = xi;
			}
		      else
			{
			  // Case 4
			  xi_mapped   = xi;
			  zeta_mapped = -zeta;
			}

		    else if (elem->node(6) == min_node)
		      if (elem->node(7) == std::min(elem->node(7), elem->node(2)))
			{
			  // Case 5
			  xi_mapped   = -xi;
			  zeta_mapped = -zeta;
			}
		      else
			{
			  // Case 6
			  xi_mapped   = -zeta;
			  zeta_mapped = -xi;
			}

		    else if (elem->node(2) == min_node)
		      if (elem->node(6) == std::min(elem->node(3), elem->node(6)))
			{
			  // Case 7
			  xi_mapped   = zeta;
			  zeta_mapped = -xi;
			}
		      else
			{
			  // Case 8
			  xi_mapped   = -xi;
			  zeta_mapped = zeta;
			}
		  }

		
		// Face 4
		else if ((i0[i] == 0) && (i1[i] >= 2) && (i2[i] >= 2))
		  {
 		    const unsigned int min_node = std::min(elem->node(3),
							   std::min(elem->node(0),
								    std::min(elem->node(4),
									     elem->node(7))));
		    if (elem->node(0) == min_node)
		      if (elem->node(3) == std::min(elem->node(3), elem->node(4)))
			{
			  // Case 1
			  eta_mapped  = eta;
			  zeta_mapped = zeta;
			}
		      else
			{
			  // Case 2
			  eta_mapped  = zeta;
			  zeta_mapped = eta;
			}

		    else if (elem->node(4) == min_node)
		      if (elem->node(0) == std::min(elem->node(0), elem->node(7)))
			{
			  // Case 3
			  eta_mapped  = -zeta;
			  zeta_mapped = eta;
			}
		      else
			{
			  // Case 4
			  eta_mapped  = eta;
			  zeta_mapped = -zeta;
			}

		    else if (elem->node(7) == min_node)
		      if (elem->node(4) == std::min(elem->node(4), elem->node(3)))
			{
			  // Case 5
			  eta_mapped  = -eta;
			  zeta_mapped = -zeta;
			}
		      else
			{
			  // Case 6
			  eta_mapped  = -zeta;
			  zeta_mapped = -eta;
			}

		    else if (elem->node(3) == min_node)
		      if (elem->node(7) == std::min(elem->node(7), elem->node(0)))
			{
			  // Case 7
			  eta_mapped   = zeta;
			  zeta_mapped = -eta;
			}
		      else
			{
			  // Case 8
			  eta_mapped  = -eta;
			  zeta_mapped = zeta;
			}
		  }

		
		// Face 5
		else if ((i2[i] == 1) && (i0[i] >= 2) && (i1[i] >= 2))
		  {
 		    const unsigned int min_node = std::min(elem->node(4),
							   std::min(elem->node(5),
								    std::min(elem->node(6),
									     elem->node(7))));
		    if (elem->node(4) == min_node)
		      if (elem->node(5) == std::min(elem->node(5), elem->node(7)))
			{
			  // Case 1
			  xi_mapped  = xi;
			  eta_mapped = eta;
			}
		      else
			{
			  // Case 2
			  xi_mapped  = eta;
			  eta_mapped = xi;
			}

		    else if (elem->node(5) == min_node)
		      if (elem->node(6) == std::min(elem->node(6), elem->node(4)))
			{
			  // Case 3
			  xi_mapped  = eta;
			  eta_mapped = -xi;
			}
		      else
			{
			  // Case 4
			  xi_mapped  = -xi;
			  eta_mapped = eta;
			}

		    else if (elem->node(6) == min_node)
		      if (elem->node(7) == std::min(elem->node(7), elem->node(5)))
			{
			  // Case 5
			  xi_mapped  = -xi;
			  eta_mapped = -eta;
			}
		      else
			{
			  // Case 6
			  xi_mapped  = -eta;
			  eta_mapped = -xi;
			}

		    else if (elem->node(7) == min_node)
		      if (elem->node(4) == std::min(elem->node(4), elem->node(6)))
			{
			  // Case 7
			  xi_mapped  = -eta;
			  eta_mapped = xi;
			}
		      else
			{
			  // Case 8
			  xi_mapped  = xi;
			  eta_mapped = eta;
			}
		  }

		
	      }
		  
	      
	      return (FE<1,HIERARCHIC>::shape(EDGE3, order, i0[i], xi_mapped)*
		      FE<1,HIERARCHIC>::shape(EDGE3, order, i1[i], eta_mapped)*
		      FE<1,HIERARCHIC>::shape(EDGE3, order, i2[i], zeta_mapped));
	    }

	    
	  default:
	    error();
	  }	
      }

    default:
      error();
    }
  
#endif
      
  error();
  return 0.;
}




template <>
Real FE<3,HIERARCHIC>::shape_deriv(const ElemType,
				   const Order,
				   const unsigned int,
				   const unsigned int,
				   const Point& )
{
  std::cerr << "Hierarchic polynomials require the element type\n"
	    << "because edge and face orientation is needed."
	    << std::endl;
  error();
  
  return 0.;
}



template <>
Real FE<3,HIERARCHIC>::shape_deriv(const Elem* elem,
				   const Order order,
				   const unsigned int i,
				   const unsigned int j,
				   const Point& p)
{
#if DIM == 3
  assert (elem != NULL);
  const ElemType type = elem->type();

  assert (j < 3);
  
  switch (order)
    {

      

      // 1st & 2nd-order Hierarchics.
    case FIRST:
    case SECOND:
      {
	switch (type)
	  {
	    // Hierarchic shape functions on the hexahedral.
	  case HEX27:
	    {
	      assert (i<27);
	
	      // Compute hex shape functions as a tensor-product
	      const Real xi   = p(0);
	      const Real eta  = p(1);
	      const Real zeta = p(2);
	
	      // The only way to make any sense of this
	      // is to look at the mgflo/mg2/mgf documentation
	      // and make the cut-out cube!
	      //                                0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26
	      static const unsigned int i0[] = {0, 1, 1, 0, 0, 1, 1, 0, 2, 1, 2, 0, 0, 1, 1, 0, 2, 1, 2, 0, 2, 2, 1, 2, 0, 2, 2};
	      static const unsigned int i1[] = {0, 0, 1, 1, 0, 0, 1, 1, 0, 2, 1, 2, 0, 0, 1, 1, 0, 2, 1, 2, 2, 0, 2, 1, 2, 2, 2};
	      static const unsigned int i2[] = {0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 2, 2, 2, 2, 1, 1, 1, 1, 0, 2, 2, 2, 2, 1, 2};

	      switch (j)
		{
		  // d()/dxi
		case 0:
		  return (FE<1,HIERARCHIC>::shape_deriv(EDGE3, order, i0[i], 0, xi)*
			  FE<1,HIERARCHIC>::shape      (EDGE3, order, i1[i],    eta)*
			  FE<1,HIERARCHIC>::shape      (EDGE3, order, i2[i],    zeta));

		  // d()/deta
		case 1:
		  return (FE<1,HIERARCHIC>::shape      (EDGE3, order, i0[i],     xi)*
			  FE<1,HIERARCHIC>::shape_deriv(EDGE3, order, i1[i], 0, eta)*
			  FE<1,HIERARCHIC>::shape      (EDGE3, order, i2[i],    zeta));

		  // d()/dzeta
		case 2:
		  return (FE<1,HIERARCHIC>::shape      (EDGE3, order, i0[i],    xi)*
			  FE<1,HIERARCHIC>::shape      (EDGE3, order, i1[i],    eta)*
			  FE<1,HIERARCHIC>::shape_deriv(EDGE3, order, i2[i], 0, zeta));

		default:
		  error();
		}
	    }


	  default:
	    error();
	  }
      }

      

      // 3rd-order Hierarchics.
    case THIRD:
      {
	switch (type)
	  {
	    // Hierarchic shape functions on the hexahedral.
	  case HEX27:
	    {
	      assert (i<64);
	
	      // Compute hex shape functions as a tensor-product
	      const Real xi    = p(0);
	      const Real eta   = p(1);
	      const Real zeta  = p(2);
	      Real xi_mapped   = p(0);
	      Real eta_mapped  = p(1);
	      Real zeta_mapped = p(2);
	
	      // The only way to make any sense of this
	      // is to look at the mgflo/mg2/mgf documentation
	      // and make the cut-out cube!
	      //  Nodes                         0  1  2  3  4  5  6  7  8  8  9  9 10 10 11 11 12 12 13 13 14 14 15 15 16 16 17 17 18 18 19 19 20 20 20 20 21 21 21 21 22 22 22 22 23 23 23 23 24 24 24 24 25 25 25 25 26 26 26 26 26 26 26 26 
	      //  DOFS                          0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 18 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 60 62 63
	      static const unsigned int i0[] = {0, 1, 1, 0, 0, 1, 1, 0, 2, 3, 1, 1, 2, 3, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 2, 3, 1, 1, 2, 3, 0, 0, 2, 3, 2, 3, 2, 3, 2, 3, 1, 1, 1, 1, 2, 3, 2, 3, 0, 0, 0, 0, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3};
	      static const unsigned int i1[] = {0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 2, 3, 1, 1, 2, 3, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 2, 3, 1, 1, 2, 3, 2, 2, 3, 3, 0, 0, 0, 0, 2, 3, 2, 3, 1, 1, 1, 1, 2, 3, 2, 3, 2, 2, 3, 3, 2, 2, 3, 3, 2, 2, 3, 3};
	      static const unsigned int i2[] = {0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 2, 3, 2, 3, 2, 3, 2, 3, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 2, 2, 3, 3, 2, 2, 3, 3, 2, 2, 3, 3, 2, 2, 3, 3, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3};



	      // handle the edge orientation
	      {
		// Edge 0
		if ((i1[i] == 0) && (i2[i] == 0))
		  {
		    if (elem->node(0) != std::min(elem->node(0), elem->node(1)))
		      xi_mapped = -xi;
		  }
		// Edge 1
		else if ((i0[i] == 1) && (i2[i] == 0))
		  {
		    if (elem->node(1) != std::min(elem->node(1), elem->node(2)))
		      eta_mapped = -eta;
		  }
		// Edge 2
		else if ((i1[i] == 1) && (i2[i] == 0))
		  {
		    if (elem->node(3) != std::min(elem->node(3), elem->node(2)))
		      xi_mapped = -xi;
		  }
		// Edge 3
		else if ((i0[i] == 0) && (i2[i] == 0))
		  {
		    if (elem->node(0) != std::min(elem->node(0), elem->node(3)))
		      eta_mapped = -eta;
		  }
		// Edge 4
		else if ((i0[i] == 0) && (i1[i] == 0))
		  {
		    if (elem->node(0) != std::min(elem->node(0), elem->node(4)))
		      zeta_mapped = -zeta;
		  }		
		// Edge 5
		else if ((i0[i] == 1) && (i1[i] == 0))
		  {
		    if (elem->node(1) != std::min(elem->node(1), elem->node(5)))
		      zeta_mapped = -zeta;
		  }		
		// Edge 6
		else if ((i0[i] == 1) && (i1[i] == 1))
		  {
		    if (elem->node(2) != std::min(elem->node(2), elem->node(6)))
		      zeta_mapped = -zeta;
		  }
		// Edge 7
		else if ((i0[i] == 0) && (i1[i] == 1))
		  {
		    if (elem->node(3) != std::min(elem->node(3), elem->node(7)))
		      zeta_mapped = -zeta;
		  }		
		// Edge 8
		else if ((i1[i] == 0) && (i2[i] == 1))
		  {
		    if (elem->node(4) != std::min(elem->node(4), elem->node(5)))
		      xi_mapped = -xi;
		  }
		// Edge 9
		else if ((i0[i] == 1) && (i2[i] == 1))
		  {
		    if (elem->node(5) != std::min(elem->node(5), elem->node(6)))
		      eta_mapped = -eta;
		  }		
		// Edge 10
		else if ((i1[i] == 1) && (i2[i] == 1))
		  {
		    if (elem->node(7) != std::min(elem->node(7), elem->node(6)))
		      xi_mapped = -xi;
		  }
		// Edge 11
		else if ((i0[i] == 0) && (i2[i] == 1))
		  {
		    if (elem->node(4) != std::min(elem->node(4), elem->node(7)))
		      eta_mapped = -eta;
		  }
	      }


	      // handle the face orientation
	      {
		// Face 0
		if (     (i2[i] == 0) && (i0[i] >= 2) && (i1[i] >= 2))
		  {
		    const unsigned int min_node = std::min(elem->node(1),
							   std::min(elem->node(2),
								    std::min(elem->node(0),
									     elem->node(3))));
		    if (elem->node(0) == min_node)
		      if (elem->node(1) == std::min(elem->node(1), elem->node(3)))
			{
			  // Case 1
			  xi_mapped  = xi;
			  eta_mapped = eta;
			}
		      else
			{
			  // Case 2
			  xi_mapped  = eta;
			  eta_mapped = xi;
			}

		    else if (elem->node(3) == min_node)
		      if (elem->node(0) == std::min(elem->node(0), elem->node(2)))
			{
			  // Case 3
			  xi_mapped  = -eta;
			  eta_mapped = xi;
			}
		      else
			{
			  // Case 4
			  xi_mapped  = xi;
			  eta_mapped = -eta;
			}

		    else if (elem->node(2) == min_node)
		      if (elem->node(3) == std::min(elem->node(3), elem->node(1)))
			{
			  // Case 5
			  xi_mapped  = -xi;
			  eta_mapped = -eta;
			}
		      else
			{
			  // Case 6
			  xi_mapped  = -eta;
			  eta_mapped = -xi;
			}

		    else if (elem->node(1) == min_node)
		      if (elem->node(2) == std::min(elem->node(2), elem->node(0)))
			{
			  // Case 7
			  xi_mapped  = eta;
			  eta_mapped = -xi;
			}
		      else
			{
			  // Case 8
			  xi_mapped  = -xi;
			  eta_mapped = eta;
			}
		  }

		
		// Face 1
		else if ((i1[i] == 0) && (i0[i] >= 2) && (i2[i] >= 2))
		  {
		    const unsigned int min_node = std::min(elem->node(0),
							   std::min(elem->node(1),
								    std::min(elem->node(5),
									     elem->node(4))));
		    if (elem->node(0) == min_node)
		      if (elem->node(1) == std::min(elem->node(1), elem->node(4)))
			{
			  // Case 1
			  xi_mapped   = xi;
			  zeta_mapped = zeta;
			}
		      else
			{
			  // Case 2
			  xi_mapped   = zeta;
			  zeta_mapped = xi;
			}

		    else if (elem->node(1) == min_node)
		      if (elem->node(5) == std::min(elem->node(5), elem->node(0)))
			{
			  // Case 3
			  xi_mapped   = zeta;
			  zeta_mapped = -xi;
			}
		      else
			{
			  // Case 4
			  xi_mapped   = -xi;
			  zeta_mapped = zeta;
			}

		    else if (elem->node(5) == min_node)
		      if (elem->node(4) == std::min(elem->node(4), elem->node(1)))
			{
			  // Case 5
			  xi_mapped   = -xi;
			  zeta_mapped = -zeta;
			}
		      else
			{
			  // Case 6
			  xi_mapped   = -zeta;
			  zeta_mapped = -xi;
			}

		    else if (elem->node(4) == min_node)
		      if (elem->node(0) == std::min(elem->node(0), elem->node(5)))
			{
			  // Case 7
			  xi_mapped   = -xi;
			  zeta_mapped = zeta;
			}
		      else
			{
			  // Case 8
			  xi_mapped   = xi;
			  zeta_mapped = -zeta;
			}
		  }

		
		// Face 2
		else if ((i0[i] == 1) && (i1[i] >= 2) && (i2[i] >= 2))
		  {
		    const unsigned int min_node = std::min(elem->node(1),
							   std::min(elem->node(2),
								    std::min(elem->node(6),
									     elem->node(5))));
		    if (elem->node(1) == min_node)
		      if (elem->node(2) == std::min(elem->node(2), elem->node(5)))
			{
			  // Case 1
			  eta_mapped  = eta;
			  zeta_mapped = zeta;
			}
		      else
			{
			  // Case 2
			  eta_mapped  = zeta;
			  zeta_mapped = eta;
			}

		    else if (elem->node(2) == min_node)
		      if (elem->node(6) == std::min(elem->node(6), elem->node(1)))
			{
			  // Case 3
			  eta_mapped  = zeta;
			  zeta_mapped = -eta;
			}
		      else
			{
			  // Case 4
			  eta_mapped  = -eta;
			  zeta_mapped = zeta;
			}

		    else if (elem->node(6) == min_node)
		      if (elem->node(5) == std::min(elem->node(5), elem->node(2)))
			{
			  // Case 5
			  eta_mapped  = -eta;
			  zeta_mapped = -zeta;
			}
		      else
			{
			  // Case 6
			  eta_mapped  = -zeta;
			  zeta_mapped = -eta;
			}

		    else if (elem->node(5) == min_node)
		      if (elem->node(1) == std::min(elem->node(1), elem->node(6)))
			{
			  // Case 7
			  eta_mapped  = -zeta;
			  zeta_mapped = eta;
			}
		      else
			{
			  // Case 8
			  eta_mapped   = eta;
			  zeta_mapped = -zeta;
			}
		  }

		
		// Face 3
		else if ((i1[i] == 1) && (i0[i] >= 2) && (i2[i] >= 2))
		  {
		    const unsigned int min_node = std::min(elem->node(2),
							   std::min(elem->node(3),
								    std::min(elem->node(7),
									     elem->node(6))));
		    if (elem->node(3) == min_node)
		      if (elem->node(2) == std::min(elem->node(2), elem->node(7)))
			{
			  // Case 1
			  xi_mapped   = xi;
			  zeta_mapped = zeta;
			}
		      else
			{
			  // Case 2
			  xi_mapped   = zeta;
			  zeta_mapped = xi;
			}

		    else if (elem->node(7) == min_node)
		      if (elem->node(3) == std::min(elem->node(3), elem->node(6)))
			{
			  // Case 3
			  xi_mapped   = -zeta;
			  zeta_mapped = xi;
			}
		      else
			{
			  // Case 4
			  xi_mapped   = xi;
			  zeta_mapped = -zeta;
			}

		    else if (elem->node(6) == min_node)
		      if (elem->node(7) == std::min(elem->node(7), elem->node(2)))
			{
			  // Case 5
			  xi_mapped   = -xi;
			  zeta_mapped = -zeta;
			}
		      else
			{
			  // Case 6
			  xi_mapped   = -zeta;
			  zeta_mapped = -xi;
			}

		    else if (elem->node(2) == min_node)
		      if (elem->node(6) == std::min(elem->node(3), elem->node(6)))
			{
			  // Case 7
			  xi_mapped   = zeta;
			  zeta_mapped = -xi;
			}
		      else
			{
			  // Case 8
			  xi_mapped   = -xi;
			  zeta_mapped = zeta;
			}
		  }

		
		// Face 4
		else if ((i0[i] == 0) && (i1[i] >= 2) && (i2[i] >= 2))
		  {
		    const unsigned int min_node = std::min(elem->node(3),
							   std::min(elem->node(0),
								    std::min(elem->node(4),
									     elem->node(7))));
		    if (elem->node(0) == min_node)
		      if (elem->node(3) == std::min(elem->node(3), elem->node(4)))
			{
			  // Case 1
			  eta_mapped  = eta;
			  zeta_mapped = zeta;
			}
		      else
			{
			  // Case 2
			  eta_mapped  = zeta;
			  zeta_mapped = eta;
			}

		    else if (elem->node(4) == min_node)
		      if (elem->node(0) == std::min(elem->node(0), elem->node(7)))
			{
			  // Case 3
			  eta_mapped  = -zeta;
			  zeta_mapped = eta;
			}
		      else
			{
			  // Case 4
			  eta_mapped  = eta;
			  zeta_mapped = -zeta;
			}

		    else if (elem->node(7) == min_node)
		      if (elem->node(4) == std::min(elem->node(4), elem->node(3)))
			{
			  // Case 5
			  eta_mapped  = -eta;
			  zeta_mapped = -zeta;
			}
		      else
			{
			  // Case 6
			  eta_mapped  = -zeta;
			  zeta_mapped = -eta;
			}

		    else if (elem->node(3) == min_node)
		      if (elem->node(7) == std::min(elem->node(7), elem->node(0)))
			{
			  // Case 7
			  eta_mapped   = zeta;
			  zeta_mapped = -eta;
			}
		      else
			{
			  // Case 8
			  eta_mapped  = -eta;
			  zeta_mapped = zeta;
			}
		  }

		
		// Face 5
		else if ((i2[i] == 1) && (i0[i] >= 2) && (i1[i] >= 2))
		  {
		    const unsigned int min_node = std::min(elem->node(4),
							   std::min(elem->node(5),
								    std::min(elem->node(6),
									     elem->node(7))));
		    if (elem->node(4) == min_node)
		      if (elem->node(5) == std::min(elem->node(5), elem->node(7)))
			{
			  // Case 1
			  xi_mapped  = xi;
			  eta_mapped = eta;
			}
		      else
			{
			  // Case 2
			  xi_mapped  = eta;
			  eta_mapped = xi;
			}

		    else if (elem->node(5) == min_node)
		      if (elem->node(6) == std::min(elem->node(6), elem->node(4)))
			{
			  // Case 3
			  xi_mapped  = eta;
			  eta_mapped = -xi;
			}
		      else
			{
			  // Case 4
			  xi_mapped  = -xi;
			  eta_mapped = eta;
			}

		    else if (elem->node(6) == min_node)
		      if (elem->node(7) == std::min(elem->node(7), elem->node(5)))
			{
			  // Case 5
			  xi_mapped  = -xi;
			  eta_mapped = -eta;
			}
		      else
			{
			  // Case 6
			  xi_mapped  = -eta;
			  eta_mapped = -xi;
			}

		    else if (elem->node(7) == min_node)
		      if (elem->node(4) == std::min(elem->node(4), elem->node(6)))
			{
			  // Case 7
			  xi_mapped  = -eta;
			  eta_mapped = xi;
			}
		      else
			{
			  // Case 8
			  xi_mapped  = xi;
			  eta_mapped = eta;
			}
		  }

		
	      }
		  
	      

	      assert (j < 3);

	      switch (j)
		{
		  // d()/dxi
		case 0:
		  return (FE<1,HIERARCHIC>::shape_deriv(EDGE3, order, i0[i], 0, xi_mapped)*
			  FE<1,HIERARCHIC>::shape      (EDGE3, order, i1[i],    eta_mapped)*
			  FE<1,HIERARCHIC>::shape      (EDGE3, order, i2[i],    zeta_mapped));

		  // d()/deta
		case 1:
		  return (FE<1,HIERARCHIC>::shape      (EDGE3, order, i0[i],    xi_mapped)*
			  FE<1,HIERARCHIC>::shape_deriv(EDGE3, order, i1[i], 0, eta_mapped)*
			  FE<1,HIERARCHIC>::shape      (EDGE3, order, i2[i],    zeta_mapped));

		  // d()/dzeta
		case 2:
		  return (FE<1,HIERARCHIC>::shape      (EDGE3, order, i0[i],    xi_mapped)*
			  FE<1,HIERARCHIC>::shape      (EDGE3, order, i1[i],    eta_mapped)*
			  FE<1,HIERARCHIC>::shape_deriv(EDGE3, order, i2[i], 0, zeta_mapped));

		default:
		  error();
		}
	    }

	    
	  default:
	    error();
	  }
      }

    default:
      error();
    }

#endif
  
  error();
  return 0.;
}



template <>
Real FE<3,HIERARCHIC>::shape_second_deriv(const ElemType,
				          const Order,
				          const unsigned int,
				          const unsigned int,
				          const Point& )
{
  std::cerr << "Hierarchic polynomials require the element type\n"
	    << "because edge and face orientation is needed."
	    << std::endl;
  error();
  
  return 0.;
}



template <>
Real FE<3,HIERARCHIC>::shape_second_deriv(const Elem*,
				          const Order,
				          const unsigned int,
				          const unsigned int,
				          const Point&)
{
  static bool warning_given = false;

  if (!warning_given)
  std::cerr << "Second derivatives for 3D hierarchics "
	    << " are not yet implemented!"
	    << std::endl;

  warning_given = true;
  return 0.;
}

// $Id: fe_hierarchic_shape_2D.C,v 1.15 2005-01-13 22:10:14 roystgnr Exp $

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


// C++ inlcludes

// Local includes
#include "fe.h"
#include "elem.h"
#include "utility.h"




template <>
Real FE<2,HIERARCHIC>::shape(const ElemType,
			     const Order,
			     const unsigned int,
			     const Point&)
{
  std::cerr << "Hierarchic polynomials require the element type\n"
	    << "because edge orientation is needed."
	    << std::endl;
  
  error();
  return 0.;
}



template <>
Real FE<2,HIERARCHIC>::shape(const Elem* elem,
			     const Order order,
			     const unsigned int i,
			     const Point& p)
{
  assert (elem != NULL);

  // Declare that we are using our own special power function
  // from the Utility namespace.  This saves typing later.
  using Utility::pow;

  const ElemType type = elem->type();
  
  switch (order)
    {      
      // 1st & 2nd-order Hierarchics.
    case FIRST:
    case SECOND:
      {
	switch (type)
	  {
	    // Hierarchic shape functions on the triangle.
	  case TRI6:
	    {
	      const Real zeta1 = p(0);
	      const Real zeta2 = p(1);
	      const Real zeta0 = 1. - zeta1 - zeta2;
	      
	      assert (i<6);


	      switch (i)
		{
		case 0:
		  return zeta0;
		  
		case 1:
		  return zeta1;
		  
		case 2:
		  return zeta2;
		  
		case 3:
		  return (pow<2>(zeta1 - zeta0) - pow<2>(zeta0 + zeta1))/2.;
		  
		case 4:
		  return (pow<2>(zeta2 - zeta1) - pow<2>(zeta1 + zeta2))/2.;
		  
		case 5:
		  return (pow<2>(zeta0 - zeta2) - pow<2>(zeta2 + zeta0))/2.;

		default:
		  error();
		}
	    }
	    

	    // Hierarchic shape functions on the quadrilateral.
	  case QUAD8:
	  case QUAD9:
	    {
	      // Compute quad shape functions as a tensor-product
	      const Real xi  = p(0);
	      const Real eta = p(1);
	      
	      assert (i < 9);
	      
	      //                                0  1  2  3  4  5  6  7  8
	      static const unsigned int i0[] = {0, 1, 1, 0, 2, 1, 2, 0, 2};
	      static const unsigned int i1[] = {0, 0, 1, 1, 0, 2, 1, 2, 2};
	      	      
	      return (FE<1,HIERARCHIC>::shape(EDGE3, order, i0[i], xi)*
		      FE<1,HIERARCHIC>::shape(EDGE3, order, i1[i], eta));
	      
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
	    // Hierarchic shape functions on the triangle.
	  case TRI6:
	    {
	      const Real zeta1 = p(0);
	      const Real zeta2 = p(1);
	      const Real zeta0 = 1. - zeta1 - zeta2;
	      
	      assert (i<10);

	      // Get  factors to account for edge-flipping
	      Real f0 = 1;
	      Real f1 = 1;
	      Real f2 = 1;

	      if (elem->node(0) > elem->node(1))
		f0 = -1.;
	      
	      if (elem->node(1) > elem->node(2))
		f1 = -1.;
	      
	      if (elem->node(2) > elem->node(0))
		f2 = -1.;


	      
	      switch (i)
		{
		case 0:
		  return zeta0;
		  
		case 1:
		  return zeta1;
		  
		case 2:
		  return zeta2;
		  
		  // Shape functions for edge 0
		case 3:
		  return    (pow<2>(zeta1 - zeta0) -
			     pow<2>(zeta0 + zeta1))/2.;

		case 4:
		  return f0*(pow<3>(zeta1 - zeta0) -
			     (zeta1 - zeta0)*pow<2>(zeta0 + zeta1))/6.;
		  
		  // Shape functions for edge 1
		case 5:
		  return    (pow<2>(zeta2 - zeta1) -
			     pow<2>(zeta1 + zeta2))/2.;

		case 6:
		  return f1*(pow<3>(zeta2 - zeta1) -
			     (zeta2 - zeta1)*pow<2>(zeta1 + zeta2))/6.;
		  
		  // Shape functions for edge 2
		case 7:
		  return    (pow<2>(zeta0 - zeta2) -
			     pow<2>(zeta2 + zeta0))/2.;

		case 8:
		  return f2*(pow<3>(zeta0 - zeta2) -
			     (zeta0 - zeta2)*pow<2>(zeta2 + zeta0))/6.;

		  // interior shape functions
		case 9:
		  return    (zeta0*zeta1*zeta2);
		  
		default:
		  error();
		}
	    }
	    

	    // Hierarchic shape functions on the quadrilateral.
	  case QUAD8:
	  case QUAD9:
	    {
	      // Compute quad shape functions as a tensor-product
	      const Real xi  = p(0);
	      const Real eta = p(1);
	      
	      assert (i < 16);

	      //                                0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15
	      static const unsigned int i0[] = {0,  1,  1,  0,  2,  3,  1,  1,  2,  3,  0,  0,  2,  3,  2,  3};
	      static const unsigned int i1[] = {0,  0,  1,  1,  0,  0,  2,  3,  1,  1,  2,  3,  2,  2,  3,  3};
	      
	      
	      // Get the factors
	      Real f0 = 1.;
	      Real f1 = 1.;
	      Real f2 = 1.;
	      Real f3 = 1.;

	      if (elem->node(0) > elem->node(1))
		f0 = -1.;
	      
	      if (elem->node(1) > elem->node(2))
		f1 = -1.;
	      
	      if (elem->node(3) > elem->node(2))
		f2 = -1.;
	      
	      if (elem->node(0) > elem->node(3))
		f3 = -1.;
	      

	      Real f = 1.;

	      
	      if ((i0[i] == 3) &&
		  (i1[i] == 0))
		f = f0;
	      else if ((i0[i] == 3) &&
		       (i1[i] == 1))
		f = f2;
	      
	      if ((i1[i] == 3) &&
		  (i0[i] == 0))
		f = f3;
	      else if ((i1[i] == 3) &&
		       (i0[i] == 1))
		f = f1;
	      

	      return f*(FE<1,HIERARCHIC>::shape(EDGE3, order, i0[i], xi)*
			FE<1,HIERARCHIC>::shape(EDGE3, order, i1[i], eta));
	    }

	  default:
	    error();
	  }
      }
	   

      

      // 4th-order Hierarchics.
    case FOURTH:
      {
	switch (type)
	  {
	    // Hierarchic shape functions on the triangle.
	  case TRI6:
	    {
	      const Real zeta1 = p(0);
	      const Real zeta2 = p(1);
	      const Real zeta0 = 1. - zeta1 - zeta2;
	      
	      assert (i<15);

	      // Get  factors to account for edge-flipping
	      Real f0 = 1;
	      Real f1 = 1;
	      Real f2 = 1;

	      if (elem->node(0) > elem->node(1))
		f0 = -1.;
	      
	      if (elem->node(1) > elem->node(2))
		f1 = -1.;
	      
	      if (elem->node(2) > elem->node(0))
		f2 = -1.;


	      
	      switch (i)
		{
		case 0:
		  return zeta0;
		  
		case 1:
		  return zeta1;
		  
		case 2:
		  return zeta2;

		  // Shape functions for edge 0
		case 3:
		  return    (pow<2>(zeta1 - zeta0) -
			     pow<2>(zeta0 + zeta1))/2.;
		   
		case 4:
		  return f0*(pow<3>(zeta1 - zeta0) -
			     (zeta1 - zeta0)*pow<2>(zeta0 + zeta1))/6.;

		case 5:
		  return    (pow<4>(zeta1 - zeta0) -
			     pow<4>(zeta0 + zeta1))/24.;
		  
		  // Shape functions for edge 1
		case 6:
		  return    (pow<2>(zeta2 - zeta1) -
			     pow<2>(zeta1 + zeta2))/2.;

		case 7:
		  return f1*(pow<3>(zeta2 - zeta1) -
			     (zeta2 - zeta1)*pow<2>(zeta1 + zeta2))/6.;

		case 8:
		  return    (pow<4>(zeta2 - zeta1) -
			     pow<4>(zeta1 + zeta2))/24.;
		  
		  // Shape functions for edge 2
		case 9:
		  return    (pow<2>(zeta0 - zeta2) -
			     pow<2>(zeta2 + zeta0))/2.;

		case 10:
		  return f2*(pow<3>(zeta0 - zeta2) -
			     (zeta0 - zeta2)*pow<2>(zeta2 + zeta0))/6.;

		case 11:
		  return    (pow<4>(zeta0 - zeta2) -
			     pow<4>(zeta2 + zeta0))/24.;

		  // interior shape functions
		case 12:
		  return    (zeta0*zeta0*zeta1*zeta2);

		case 13:
		  return    (zeta0*zeta1*zeta1*zeta2);
		  
		case 14:
		  return    (zeta0*zeta1*zeta2*zeta2);
		  
		default:
		  error();
		}
	    }
	    

	    // Hierarchic shape functions on the quadrilateral.
	  case QUAD8:
	  case QUAD9:
	    {
	      // Compute quad shape functions as a tensor-product
	      const Real xi  = p(0);
	      const Real eta = p(1);
	      
	      assert (i < 25);

	      //                                0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24
	      static const unsigned int i0[] = {0, 1, 1, 0, 2, 3, 4, 1, 1, 1, 2, 3, 4, 0, 0, 0, 2, 3, 4, 2, 3, 4, 2, 3, 4};
	      static const unsigned int i1[] = {0, 0, 1, 1, 0, 0, 0, 2, 3, 4, 1, 1, 1, 2, 3, 4, 2, 2, 2, 3, 3, 3, 4, 4, 4};
	      
	      
	      // Get the factors
	      Real f0 = 1.;
	      Real f1 = 1.;
	      Real f2 = 1.;
	      Real f3 = 1.;

	      if (elem->node(0) > elem->node(1))
		f0 = -1.;
	      
	      if (elem->node(1) > elem->node(2))
		f1 = -1.;
	      
	      if (elem->node(3) > elem->node(2))
		f2 = -1.;
	      
	      if (elem->node(0) > elem->node(3))
		f3 = -1.;
	      

	      Real f = 1.;

	      
	      if ((i0[i] == 3) &&
		  (i1[i] == 0))
		f = f0;
	      else if ((i0[i] == 3) &&
		       (i1[i] == 1))
		f = f2;
	      
	      if ((i1[i] == 3) &&
		  (i0[i] == 0))
		f = f3;
	      else if ((i1[i] == 3) &&
		       (i0[i] == 1))
		f = f1;
	      

	      return f*(FE<1,HIERARCHIC>::shape(EDGE3, order, i0[i], xi)*
			FE<1,HIERARCHIC>::shape(EDGE3, order, i1[i], eta));
	    }

	  default:
	    error();
	  }
      }
	   

      

      // 5th-order Hierarchics.
    case FIFTH:
      {
	switch (type)
	  {
	    // Hierarchic shape functions on the triangle.
	  case TRI6:
	    {
	      const Real zeta1 = p(0);
	      const Real zeta2 = p(1);
	      const Real zeta0 = 1. - zeta1 - zeta2;
	      
	      assert (i<21);

	      // Get  factors to account for edge-flipping
	      Real f0 = 1;
	      Real f1 = 1;
	      Real f2 = 1;

	      if (elem->node(0) > elem->node(1))
		f0 = -1.;
	      
	      if (elem->node(1) > elem->node(2))
		f1 = -1.;
	      
	      if (elem->node(2) > elem->node(0))
		f2 = -1.;


	      
	      switch (i)
		{
		case 0:
		  return zeta0;
		  
		case 1:
		  return zeta1;
		  
		case 2:
		  return zeta2;

		  // Shape functions for edge 0
		case 3:
		  return    (pow<2>(zeta1 - zeta0) -
			     pow<2>(zeta0 + zeta1))/2.;
		   
		case 4:
		  return f0*(pow<3>(zeta1 - zeta0) -
			     (zeta1 - zeta0)*pow<2>(zeta0 + zeta1))/6.;

		case 5:
		  return    (pow<4>(zeta1 - zeta0) -
			     pow<4>(zeta0 + zeta1))/24.;

		case 6:
		  return f0*(pow<5>(zeta1 - zeta0) -
			     (zeta1 - zeta0)*pow<4>(zeta0 + zeta1))/120.;
		  		  
		  // Shape functions for edge 1
		case 7:
		  return    (pow<2>(zeta2 - zeta1) -
			     pow<2>(zeta1 + zeta2))/2.;

		case 8:
		  return f1*(pow<3>(zeta2 - zeta1) -
			     (zeta2 - zeta1)*pow<2>(zeta1 + zeta2))/6.;

		case 9:
		  return    (pow<4>(zeta2 - zeta1) -
			     pow<4>(zeta1 + zeta2))/24.;

		case 10:
		  return f1*(pow<5>(zeta2 - zeta1) -
			     (zeta2 - zeta1)*pow<4>(zeta1 + zeta2))/120.;
		  
		  // Shape functions for edge 2
		case 11:
		  return    (pow<2>(zeta0 - zeta2) -
			     pow<2>(zeta2 + zeta0))/2.;

		case 12:
		  return f2*(pow<3>(zeta0 - zeta2) -
			     (zeta0 - zeta2)*pow<2>(zeta2 + zeta0))/6.;

		case 13:
		  return    (pow<4>(zeta0 - zeta2) -
			     pow<4>(zeta2 + zeta0))/24.;

		case 14:
		  return f2*(pow<5>(zeta0 - zeta2) -
			     (zeta0 - zeta2)*pow<4>(zeta2 + zeta0))/120.;
		  
		  // interior shape functions
		case 15:
		  return    pow<3>(zeta0)*zeta1*zeta2;

		case 16:
		  return    zeta0*pow<3>(zeta1)*zeta2;
		  
		case 17:
		  return    zeta0*zeta1*pow<3>(zeta2);

		case 18:
		  return    pow<2>(zeta0*zeta1)*zeta2;
		  
		case 19:
		  return    zeta0*pow<2>(zeta1*zeta2);
		  
		case 20:
		  return    zeta1*pow<2>(zeta0*zeta2);
		  
		default:
		  error();
		}
	    }
	    

	    // Hierarchic shape functions on the quadrilateral.
	  case QUAD8:
	  case QUAD9:
	    {
	      // Compute quad shape functions as a tensor-product
	      const Real xi  = p(0);
	      const Real eta = p(1);
	      
	      assert (i < 36);

	      //                                0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35
	      static const unsigned int i0[] = {0, 1, 1, 0, 2, 3, 4, 5, 1, 1, 1, 1, 2, 3, 4, 5, 0, 0, 0, 0, 2, 3, 4, 5, 2, 3, 4, 5, 2, 3, 4, 5, 2, 3, 4, 5};
	      static const unsigned int i1[] = {0, 0, 1, 1, 0, 0, 0, 0, 2, 3, 4, 5, 1, 1, 1, 1, 2, 3, 4, 5, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5};
	      
	      
	      // Get the factors
	      Real f0 = 1.;
	      Real f1 = 1.;
	      Real f2 = 1.;
	      Real f3 = 1.;

	      if (elem->node(0) > elem->node(1))
		f0 = -1.;
	      
	      if (elem->node(1) > elem->node(2))
		f1 = -1.;
	      
	      if (elem->node(3) > elem->node(2))
		f2 = -1.;
	      
	      if (elem->node(0) > elem->node(3))
		f3 = -1.;
	      

	      Real f = 1.;

	      
	      if ( ((i0[i] == 3) || (i0[i] == 5)) &&
		   (i1[i] == 0))
		f = f0;
	      else if ( ((i0[i] == 3) || (i0[i] == 5))&&
			(i1[i] == 1))
		f = f2;
	      
	      if ( ((i1[i] == 3) || (i1[i] == 5)) &&
		   (i0[i] == 0))
		f = f3;
	      else if ( ((i1[i] == 3) || (i1[i] == 5)) &&
			(i0[i] == 1))
		f = f1;
	      

	      return f*(FE<1,HIERARCHIC>::shape(EDGE3, order, i0[i], xi)*
			FE<1,HIERARCHIC>::shape(EDGE3, order, i1[i], eta));	      
	    }

	  default:
	    error();
	  }
      }



      
      // by default throw an error
    default:
      std::cerr << "ERROR: Unsupported polynomial order!" << std::endl;
      error();
    }

  
  error();
  return 0.;
}



template <>
Real FE<2,HIERARCHIC>::shape_deriv(const ElemType,
				   const Order,			    
				   const unsigned int,
				   const unsigned int,
				   const Point&)
{
  std::cerr << "Hierarchic polynomials require the element type\n"
	    << "because edge orientation is needed."
	    << std::endl;

  error();
  return 0.;
}



template <>
Real FE<2,HIERARCHIC>::shape_deriv(const Elem* elem,
				   const Order order,
				   const unsigned int i,
				   const unsigned int j,
				   const Point& p)
{
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
	    // Hierarchic shape functions on the triangle.
	  case TRI6:
	    {
	      // I have been lazy here and am using finite differences
	      // to compute the derivatives!
	      const Real eps = 1.e-6;
	      
	      assert (i < 6);
	      assert (j < 2);
	      
	      switch (j)
		{
		  //  d()/dxi
		case 0:
		  {
		    const Point pp(p(0)+eps, p(1));
		    const Point pm(p(0)-eps, p(1));

		    return (FE<2,HIERARCHIC>::shape(elem, order, i, pp) -
			    FE<2,HIERARCHIC>::shape(elem, order, i, pm))/2./eps;
		  }

		  // d()/deta
		case 1:
		  {
		    const Point pp(p(0), p(1)+eps);
		    const Point pm(p(0), p(1)-eps);

		    return (FE<2,HIERARCHIC>::shape(elem, order, i, pp) -
			    FE<2,HIERARCHIC>::shape(elem, order, i, pm))/2./eps;
		  }
		  

		default:
		  error();
		}
	    }

	    

	    // Hierarchic shape functions on the quadrilateral.
	  case QUAD8:
	  case QUAD9:
	    {
	      // Compute quad shape functions as a tensor-product
	      const Real xi  = p(0);
	      const Real eta = p(1);
	      
	      assert (i < 9);
	      
	      //                                0  1  2  3  4  5  6  7  8
	      static const unsigned int i0[] = {0, 1, 1, 0, 2, 1, 2, 0, 2};
	      static const unsigned int i1[] = {0, 0, 1, 1, 0, 2, 1, 2, 2};

	      switch (j)
		{
		  // d()/dxi
		case 0:		      
		  return (FE<1,HIERARCHIC>::shape_deriv(EDGE3, order, i0[i], 0, xi)*
			  FE<1,HIERARCHIC>::shape      (EDGE3, order, i1[i],    eta));

		  // d()/deta
		case 1:		      
		  return (FE<1,HIERARCHIC>::shape      (EDGE3, order, i0[i],    xi)*
			  FE<1,HIERARCHIC>::shape_deriv(EDGE3, order, i1[i], 0, eta));

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
	    // Hierarchic shape functions on the triangle.
	  case TRI6:
	    {
	      // I have been lazy here and am using finite differences
	      // to compute the derivatives!
	      const Real eps = 1.e-6;
	      
	      assert (i < 10);
	      assert (j < 2);
	      
	      switch (j)
		{
		  //  d()/dxi
		case 0:
		  {
		    const Point pp(p(0)+eps, p(1));
		    const Point pm(p(0)-eps, p(1));

		    return (FE<2,HIERARCHIC>::shape(elem, order, i, pp) -
			    FE<2,HIERARCHIC>::shape(elem, order, i, pm))/2./eps;
		  }

		  // d()/deta
		case 1:
		  {
		    const Point pp(p(0), p(1)+eps);
		    const Point pm(p(0), p(1)-eps);

		    return (FE<2,HIERARCHIC>::shape(elem, order, i, pp) -
			    FE<2,HIERARCHIC>::shape(elem, order, i, pm))/2./eps;
		  }
		  

		default:
		  error();
		}
	    }

	    

	    // Hierarchic shape functions on the quadrilateral.
	  case QUAD8:
	  case QUAD9:
	    {
	      // Compute quad shape functions as a tensor-product
	      const Real xi  = p(0);
	      const Real eta = p(1);
	      
	      assert (i < 16);

	      //                                0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15
	      static const unsigned int i0[] = {0,  1,  1,  0,  2,  3,  1,  1,  2,  3,  0,  0,  2,  3,  2,  3};
	      static const unsigned int i1[] = {0,  0,  1,  1,  0,  0,  2,  3,  1,  1,  2,  3,  2,  2,  3,  3};
	      
	      
	      // Get the factors
	      Real f0 = 1.;
	      Real f1 = 1.;
	      Real f2 = 1.;
	      Real f3 = 1.;

	      if (elem->node(0) > elem->node(1))
		f0 = -1.;
	      
	      if (elem->node(1) > elem->node(2))
		f1 = -1.;
	      
	      if (elem->node(3) > elem->node(2))
		f2 = -1.;
	      
	      if (elem->node(0) > elem->node(3))
		f3 = -1.;
	      

	      Real f = 1.;

	      
	      if ((i0[i] == 3) &&
		  (i1[i] == 0))
		f = f0;
	      else if ((i0[i] == 3) &&
		       (i1[i] == 1))
		f = f2;
	      
	      if ((i1[i] == 3) &&
		  (i0[i] == 0))
		f = f3;
	      else if ((i1[i] == 3) &&
		       (i0[i] == 1))
		f = f1;

	      
	      switch (j)
		{
		  // d()/dxi
		case 0:		  		  
		  return f*(FE<1,HIERARCHIC>::shape_deriv(EDGE3, order, i0[i], 0, xi)*
			    FE<1,HIERARCHIC>::shape      (EDGE3, order, i1[i],    eta));
	      
		  // d()/deta
		case 1:		  		  
		  return f*(FE<1,HIERARCHIC>::shape      (EDGE3, order, i0[i],    xi)*
			    FE<1,HIERARCHIC>::shape_deriv(EDGE3, order, i1[i], 0, eta));

		default:
		  error();
		}
	    }

	  default:
	    error();
	  }
      }
	   

      

      // 4th-order Hierarchics.
    case FOURTH:
      {
	switch (type)
	  {
	    // Hierarchic shape functions on the triangle.
	  case TRI6:
	    {
	      // I have been lazy here and am using finite differences
	      // to compute the derivatives!
	      const Real eps = 1.e-6;
	      
	      assert (i < 15);
	      assert (j < 2);
	      
	      switch (j)
		{
		  //  d()/dxi
		case 0:
		  {
		    const Point pp(p(0)+eps, p(1));
		    const Point pm(p(0)-eps, p(1));

		    return (FE<2,HIERARCHIC>::shape(elem, order, i, pp) -
			    FE<2,HIERARCHIC>::shape(elem, order, i, pm))/2./eps;
		  }

		  // d()/deta
		case 1:
		  {
		    const Point pp(p(0), p(1)+eps);
		    const Point pm(p(0), p(1)-eps);

		    return (FE<2,HIERARCHIC>::shape(elem, order, i, pp) -
			    FE<2,HIERARCHIC>::shape(elem, order, i, pm))/2./eps;
		  }
		  

		default:
		  error();
		}
	    }

	    

	    // Hierarchic shape functions on the quadrilateral.
	  case QUAD8:
	  case QUAD9:
	    {
	      // Compute quad shape functions as a tensor-product
	      const Real xi  = p(0);
	      const Real eta = p(1);
	      
	      assert (i < 25);

	      //                                0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24
	      static const unsigned int i0[] = {0, 1, 1, 0, 2, 3, 4, 1, 1, 1, 2, 3, 4, 0, 0, 0, 2, 3, 4, 2, 3, 4, 2, 3, 4};
	      static const unsigned int i1[] = {0, 0, 1, 1, 0, 0, 0, 2, 3, 4, 1, 1, 1, 2, 3, 4, 2, 2, 2, 3, 3, 3, 4, 4, 4};
	      
	      
	      // Get the factors
	      Real f0 = 1.;
	      Real f1 = 1.;
	      Real f2 = 1.;
	      Real f3 = 1.;

	      if (elem->node(0) > elem->node(1))
		f0 = -1.;
	      
	      if (elem->node(1) > elem->node(2))
		f1 = -1.;
	      
	      if (elem->node(3) > elem->node(2))
		f2 = -1.;
	      
	      if (elem->node(0) > elem->node(3))
		f3 = -1.;
	      

	      Real f = 1.;

	      
	      if ((i0[i] == 3) &&
		  (i1[i] == 0))
		f = f0;
	      else if ((i0[i] == 3) &&
		       (i1[i] == 1))
		f = f2;
	      
	      if ((i1[i] == 3) &&
		  (i0[i] == 0))
		f = f3;
	      else if ((i1[i] == 3) &&
		       (i0[i] == 1))
		f = f1;	      

	      
	      switch (j)
		{
		  // d()/dxi
		case 0:		  		  
		  return f*(FE<1,HIERARCHIC>::shape_deriv(EDGE3, order, i0[i], 0, xi)*
			    FE<1,HIERARCHIC>::shape      (EDGE3, order, i1[i],    eta));
	      
		  // d()/deta
		case 1:		  		  
		  return f*(FE<1,HIERARCHIC>::shape      (EDGE3, order, i0[i],    xi)*
			    FE<1,HIERARCHIC>::shape_deriv(EDGE3, order, i1[i], 0, eta));

		default:
		  error();
		}
	    }

	  default:
	    error();
	  }
      }
	   

      

      // 5th-order Hierarchics.
    case FIFTH:
      {
	// Hierarchic shape functions on the quadrilateral.
	switch (type)
	  {
	    // Hierarchic shape functions on the triangle.
	  case TRI6:
	    {
	      // I have been lazy here and am using finite differences
	      // to compute the derivatives!
	      const Real eps = 1.e-6;
	      
	      assert (i < 21);
	      assert (j < 2);
	      
	      switch (j)
		{
		  //  d()/dxi
		case 0:
		  {
		    const Point pp(p(0)+eps, p(1));
		    const Point pm(p(0)-eps, p(1));

		    return (FE<2,HIERARCHIC>::shape(elem, order, i, pp) -
			    FE<2,HIERARCHIC>::shape(elem, order, i, pm))/2./eps;
		  }

		  // d()/deta
		case 1:
		  {
		    const Point pp(p(0), p(1)+eps);
		    const Point pm(p(0), p(1)-eps);

		    return (FE<2,HIERARCHIC>::shape(elem, order, i, pp) -
			    FE<2,HIERARCHIC>::shape(elem, order, i, pm))/2./eps;
		  }
		  

		default:
		  error();
		}
	    }

	    

	  case QUAD8:
	  case QUAD9:
	    {
	      // Compute quad shape functions as a tensor-product
	      const Real xi  = p(0);
	      const Real eta = p(1);
	      
	      assert (i < 36);

	      //                                0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35
	      static const unsigned int i0[] = {0, 1, 1, 0, 2, 3, 4, 5, 1, 1, 1, 1, 2, 3, 4, 5, 0, 0, 0, 0, 2, 3, 4, 5, 2, 3, 4, 5, 2, 3, 4, 5, 2, 3, 4, 5};
	      static const unsigned int i1[] = {0, 0, 1, 1, 0, 0, 0, 0, 2, 3, 4, 5, 1, 1, 1, 1, 2, 3, 4, 5, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5};
	      
	      
	      // Get the factors
	      Real f0 = 1.;
	      Real f1 = 1.;
	      Real f2 = 1.;
	      Real f3 = 1.;

	      if (elem->node(0) > elem->node(1))
		f0 = -1.;
	      
	      if (elem->node(1) > elem->node(2))
		f1 = -1.;
	      
	      if (elem->node(3) > elem->node(2))
		f2 = -1.;
	      
	      if (elem->node(0) > elem->node(3))
		f3 = -1.;
	      

	      Real f = 1.;

	      
	      if ( ((i0[i] == 3) || (i0[i] == 5)) &&
		   (i1[i] == 0))
		f = f0;
	      else if ( ((i0[i] == 3) || (i0[i] == 5))&&
			(i1[i] == 1))
		f = f2;
	      
	      if ( ((i1[i] == 3) || (i1[i] == 5)) &&
		   (i0[i] == 0))
		f = f3;
	      else if ( ((i1[i] == 3) || (i1[i] == 5)) &&
			(i0[i] == 1))
		f = f1;	      

	      
	      switch (j)
		{
		  // d()/dxi
		case 0:		  		  
		  return f*(FE<1,HIERARCHIC>::shape_deriv(EDGE3, order, i0[i], 0, xi)*
			    FE<1,HIERARCHIC>::shape      (EDGE3, order, i1[i],    eta));
	      
		  // d()/deta
		case 1:		  		  
		  return f*(FE<1,HIERARCHIC>::shape      (EDGE3, order, i0[i],    xi)*
			    FE<1,HIERARCHIC>::shape_deriv(EDGE3, order, i1[i], 0, eta));

		default:
		  error();
		}
	    }

	  default:
	    error();
	  }
      }



      
      // by default throw an error;call the orientation-independent shape functions
    default:
      std::cerr << "ERROR: Unsupported polynomial order!" << std::endl;
      error();
    }

  
  error();
  return 0.;
}



template <>
Real FE<2,HIERARCHIC>::shape_second_deriv(const ElemType,
				          const Order,			    
				          const unsigned int,
				          const unsigned int,
				          const Point&)
{
  std::cerr << "Hierarchic polynomials require the element type\n"
	    << "because edge orientation is needed."
	    << std::endl;

  error();
  return 0.;
}



template <>
Real FE<2,HIERARCHIC>::shape_second_deriv(const Elem* elem,
				          const Order order,
				          const unsigned int i,
				          const unsigned int j,
				          const Point& p)
{
  assert (elem != NULL);

  // I have been lazy here and am using finite differences
  // to compute the derivatives!
  const Real eps = 1.e-6;
  Point pp, pm;
  unsigned int prevj;
	      
  switch (j)
  {
    //  d^2()/dxi^2
    case 0:
      {
        pp = Point(p(0)+eps, p(1));
        pm = Point(p(0)-eps, p(1));
        prevj = 0;
	break;
      }
  
    // d^2()/dxideta
    case 1:
      {
        pp = Point(p(0), p(1)+eps);
        pm = Point(p(0), p(1)-eps);
        prevj = 0;
	break;
      }

    // d^2()/deta^2
    case 2:
      {
        pp = Point(p(0), p(1)+eps);
        pm = Point(p(0), p(1)-eps);
        prevj = 1;
	break;
      }
    default:
      error();
  }
  return (FE<2,HIERARCHIC>::shape_deriv(elem, order, i, prevj, pp) -
	  FE<2,HIERARCHIC>::shape_deriv(elem, order, i, prevj,
					pm))/2./eps;
}

// $Id: fe_lagrange_shape_2D.C,v 1.5 2003-01-24 17:24:41 jwpeterson Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2002  Benjamin S. Kirk, John W. Peterson
  
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
#include <math.h>


// Local includes
#include "fe.h"
#include "elem.h"




template <>
real FE<2,LAGRANGE>::shape(const ElemType type,
			   const Order order,
			   const unsigned int i,
			   const Point& p)
{
#if DIM > 1
  
  switch (order)
    {
      // linear Lagrange shape functions
    case FIRST:
      {
	switch (type)
	  {
	  case QUAD4:
	  case QUAD8:
	  case QUAD9:
	    {
	      // Compute quad shape functions as a tensor-product
	      const real xi  = p(0);
	      const real eta = p(1);
	      
	      assert (i<4);
	      
	      //                                0  1  2  3  
	      static const unsigned int i0[] = {0, 1, 1, 0};
	      static const unsigned int i1[] = {0, 0, 1, 1};
	      
	      return (FE<1,LAGRANGE>::shape(EDGE2, FIRST, i0[i], xi)*
		      FE<1,LAGRANGE>::shape(EDGE2, FIRST, i1[i], eta));
	    };

	  case TRI3:
	  case TRI6:
	    {
	      const real zeta1 = p(0);
	      const real zeta2 = p(1);
	      const real zeta0 = 1. - zeta1 - zeta2;

	      assert (i<3);
	      
	      switch(i)
		{
		case 0:
		  return zeta0;
		  
		case 1:
		  return zeta1;
		  
		case 2:
		  return zeta2;
		  
		default:
		  error();
		  
		};
	    };
	    
	  default:
	    {
	      std::cerr << "ERROR: Unsupported 2D element type!: " << type
			<< std::endl;
	      error();
	    };
	  };
      };
      

      // quadratic Lagrange shape functions
    case SECOND:
      {
	switch (type)
	  {
	  case QUAD8:
	    {
	      const real xi  = p(0);
	      const real eta = p(1);

	      assert (i<8);

	      switch (i)
		{
		case 0:
		  return .25*(1. - xi)*(1. - eta)*(-1. - xi - eta);
	    
		case 1:
		  return .25*(1. + xi)*(1. - eta)*(-1. + xi - eta);
	    
		case 2:
		  return .25*(1. + xi)*(1. + eta)*(-1. + xi + eta);
	    
		case 3:
		  return .25*(1. - xi)*(1. + eta)*(-1. - xi + eta);

		case 4:
		  return .5*(1. - xi*xi)*(1. - eta);
	      
		case 5:
		  return .5*(1. + xi)*(1. - eta*eta);

		case 6:
		  return .5*(1. - xi*xi)*(1. + eta);

		case 7:
		  return .5*(1. - xi)*(1. - eta*eta);

		default:
		  error();
		};
	    };
	    
	  case QUAD9:
	    {
	      // Compute quad shape functions as a tensor-product
	      const real xi  = p(0);
	      const real eta = p(1);
	      
	      assert (i<9);
	      
	      //                                0  1  2  3  4  5  6  7  8
	      static const unsigned int i0[] = {0, 1, 1, 0, 2, 1, 2, 0, 2};
	      static const unsigned int i1[] = {0, 0, 1, 1, 0, 2, 1, 2, 2};
		  
	      return (FE<1,LAGRANGE>::shape(EDGE3, SECOND, i0[i], xi)*
		      FE<1,LAGRANGE>::shape(EDGE3, SECOND, i1[i], eta));
	    };
	    
	  case TRI6:
	    {
	      const real zeta1 = p(0);
	      const real zeta2 = p(1);
	      const real zeta0 = 1. - zeta1 - zeta2;
	      
	      assert (i<6);
	      
	      switch(i)
		{
		case 0:
		  return 2.*zeta0*(zeta0-0.5);
		  
		case 1:
		  return 2.*zeta1*(zeta1-0.5);
		  
		case 2:
		  return 2.*zeta2*(zeta2-0.5);
		  
		case 3:
		  return 4.*zeta0*zeta1;
		  
		case 4:
		  return 4.*zeta1*zeta2;
		  
		case 5:
		  return 4.*zeta2*zeta0;
	    
		default:
		  error();
	    
		};
	    };
	    
	  default:
	    {
	      std::cerr << "ERROR: Unsupported 2D element type!: " << type
			<< std::endl;
	      error();
	    };
	  };
      };


      
      // unsupported order
    default:
      {
	std::cerr << "ERROR: Unsupported 2D FE order!: " << order
		  << std::endl;
	error();
      };
    };

  error();
  return 0.;

#endif
};



template <>
real FE<2,LAGRANGE>::shape(const Elem* elem,
			   const Order order,
			   const unsigned int i,
			   const Point& p)
{
  assert (elem != NULL);

  // call the orientation-independent shape functions
  return FE<2,LAGRANGE>::shape(elem->type(), order, i, p);
};



template <>
real FE<2,LAGRANGE>::shape_deriv(const ElemType type,
				 const Order order,
				 const unsigned int i,
				 const unsigned int j,
				 const Point& p)
{
#if DIM > 1

  
  assert (j<2);

  switch (order)
    {
      // linear Lagrange shape functions
    case FIRST:
      {
	switch (type)
	  {
	  case QUAD4:
	  case QUAD8:
	  case QUAD9:
	    {
	      // Compute quad shape functions as a tensor-product
	      const real xi  = p(0);
	      const real eta = p(1);

	      assert (i<4);
	
	      //                                0  1  2  3  
	      static const unsigned int i0[] = {0, 1, 1, 0};
	      static const unsigned int i1[] = {0, 0, 1, 1};
		  
	      switch (j)
		{
		  // d()/dxi
		case 0:
		  return (FE<1,LAGRANGE>::shape_deriv(EDGE2, FIRST, i0[i], 0, xi)*
			  FE<1,LAGRANGE>::shape      (EDGE2, FIRST, i1[i], eta));
		  
		  // d()/deta
		case 1:
		  return (FE<1,LAGRANGE>::shape      (EDGE2, FIRST, i0[i], xi)*
			  FE<1,LAGRANGE>::shape_deriv(EDGE2, FIRST, i1[i], 0, eta));
		  
		default:
		  error();
		};
	    };

	  case TRI3:
	  case TRI6:
	    {
	      assert (i<3);
	      
	      const real dzeta0dxi  = -1.;
	      const real dzeta1dxi  = 1.;
	      const real dzeta2dxi  = 0.;
	      
	      const real dzeta0deta = -1.;
	      const real dzeta1deta = 0.;
	      const real dzeta2deta = 1.;
	      
	      switch (j)
		{
		  // d()/dxi
		case 0:
		  {
		    switch(i)
		      {
		      case 0:
			return dzeta0dxi;
			
		      case 1:
			return dzeta1dxi;
			
		      case 2:
			return dzeta2dxi;
			
		      default:
			error();
		      };
		  }
		  // d()/deta
		case 1:
		  {
		    switch(i)
		      {
		      case 0:
			return dzeta0deta;
			
		      case 1:
			return dzeta1deta;
			
		      case 2:
			return dzeta2deta;
			
		      default:
			error();
			
		      };
		  }
		default:
		  error();
		};
	    };
	    
	  default:
	    {
	      std::cerr << "ERROR: Unsupported 2D element type!: " << type
			<< std::endl;
	      error();
	    };
	  };
      };
      

      // quadratic Lagrange shape functions
    case SECOND:
      {
	switch (type)
	  {
	  case QUAD8:
	    {
	      const real xi  = p(0);
	      const real eta = p(1);
	
	      assert (i<8);
	
	      switch (j)
		{
		  // d/dxi
		case 0:
		  switch (i)
		    {
		    case 0:
		      return .25*(1. - eta)*((1. - xi)*(-1.) +
					     (-1.)*(-1. - xi - eta));
		
		    case 1:
		      return .25*(1. - eta)*((1. + xi)*(1.) +
					     (1.)*(-1. + xi - eta));
		
		    case 2:
		      return .25*(1. + eta)*((1. + xi)*(1.) +
					     (1.)*(-1. + xi + eta));
		
		    case 3:
		      return .25*(1. + eta)*((1. - xi)*(-1.) +
					     (-1.)*(-1. - xi + eta));
		
		    case 4:
		      return .5*(-2.*xi)*(1. - eta);
		
		    case 5:
		      return .5*(1.)*(1. - eta*eta);
		
		    case 6:
		      return .5*(-2.*xi)*(1. + eta);
		
		    case 7:
		      return .5*(-1.)*(1. - eta*eta);
		
		    default:
		      error();
		    };
	    
		  // d/deta
		case 1:
		  switch (i)
		    {
		    case 0:
		      return .25*(1. - xi)*((1. - eta)*(-1.) +
					    (-1.)*(-1. - xi - eta));
		
		    case 1:
		      return .25*(1. + xi)*((1. - eta)*(-1.) +
					    (-1.)*(-1. + xi - eta));
		
		    case 2:
		      return .25*(1. + xi)*((1. + eta)*(1.) +
					    (1.)*(-1. + xi + eta));
		
		    case 3:
		      return .25*(1. - xi)*((1. + eta)*(1.) +
					    (1.)*(-1. - xi + eta));
		
		    case 4:
		      return .5*(1. - xi*xi)*(-1.);
		
		    case 5:
		      return .5*(1. + xi)*(-2.*eta);
		
		    case 6:
		      return .5*(1. - xi*xi)*(1.);
		
		    case 7:
		      return .5*(1. - xi)*(-2.*eta);
		
		    default:
		      error();
		    };

		default:
		  error();
		};
	    }
	    
	  case QUAD9:
	    {
	      // Compute quad shape functions as a tensor-product
	      const real xi  = p(0);
	      const real eta = p(1);
	      
	      assert (i<9);
	      
	      //                                0  1  2  3  4  5  6  7  8
	      static const unsigned int i0[] = {0, 1, 1, 0, 2, 1, 2, 0, 2};
	      static const unsigned int i1[] = {0, 0, 1, 1, 0, 2, 1, 2, 2};
	
	      switch (j)
		{
		  // d()/dxi
		case 0:
		  return (FE<1,LAGRANGE>::shape_deriv(EDGE3, SECOND, i0[i], 0, xi)*
			  FE<1,LAGRANGE>::shape      (EDGE3, SECOND, i1[i], eta));
		  
		  // d()/deta
		case 1:
		  return (FE<1,LAGRANGE>::shape      (EDGE3, SECOND, i0[i], xi)*
			  FE<1,LAGRANGE>::shape_deriv(EDGE3, SECOND, i1[i], 0, eta));
		  
		default:
		  error();
		};
	    };
	    
	  case TRI6:
	    {
	      assert (i<6);
		  
	      const real zeta1 = p(0);
	      const real zeta2 = p(1);
	      const real zeta0 = 1. - zeta1 - zeta2;
	
	      const real dzeta0dxi  = -1.;
	      const real dzeta1dxi  = 1.;
	      const real dzeta2dxi  = 0.;
	
	      const real dzeta0deta = -1.;
	      const real dzeta1deta = 0.;
	      const real dzeta2deta = 1.;
	
	      switch(j)
		{
		case 0:
		  {
		    switch(i)
		      {
		      case 0:
			return (4.*zeta0-1.)*dzeta0dxi;
		 
		      case 1:
			return (4.*zeta1-1.)*dzeta1dxi;
		 
		      case 2:
			return (4.*zeta2-1.)*dzeta2dxi;
		 
		      case 3:
			return 4.*zeta1*dzeta0dxi + 4.*zeta0*dzeta1dxi;
		 
		      case 4:
			return 4.*zeta2*dzeta1dxi + 4.*zeta1*dzeta2dxi;
		 
		      case 5:
			return 4.*zeta2*dzeta0dxi + 4*zeta0*dzeta2dxi;
		 
		      default:
			error();
		      };
		  }
		      
		case 1:
		  {
		    switch(i)
		      {
		      case 0:
			return (4.*zeta0-1.)*dzeta0deta;
		 
		      case 1:
			return (4.*zeta1-1.)*dzeta1deta;
		 
		      case 2:
			return (4.*zeta2-1.)*dzeta2deta;
		 
		      case 3:
			return 4.*zeta1*dzeta0deta + 4.*zeta0*dzeta1deta;
		 
		      case 4:
			return 4.*zeta2*dzeta1deta + 4.*zeta1*dzeta2deta;
		 
		      case 5:
			return 4.*zeta2*dzeta0deta + 4*zeta0*dzeta2deta;
		 
		      default:
			error();
		      };
		  }
		default:
		  error();
		};
	    };
	    
	  default:
	    {
	      std::cerr << "ERROR: Unsupported 2D element type!: " << type
			<< std::endl;
	      error();
	    };
	  };
      };

      
      
      // unsupported order
    default:
      {
	std::cerr << "ERROR: Unsupported 2D FE order!: " << order
		  << std::endl;
	error();
      }
    };


  error();
  return 0.;

#endif
};



template <>
real FE<2,LAGRANGE>::shape_deriv(const Elem* elem,
				 const Order order,
				 const unsigned int i,
				 const unsigned int j,
				 const Point& p)
{
  assert (elem != NULL);


  // call the orientation-independent shape functions
  return FE<2,LAGRANGE>::shape_deriv(elem->type(), order, i, j, p);
};


// $Id: fe_szabab_shape_2D.C,v 1.5 2004-02-10 13:28:07 benkirk Exp $

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



// C++ includes
#include <math.h>


// Local includes
#include "libmesh_config.h"

#ifdef ENABLE_HIGHER_ORDER_SHAPES

#include "fe.h"
#include "elem.h"
#include "utility.h"


// Anonymous namespace to hold static sqrt values
namespace
{
  static const Real sqrt2  = sqrt(2.);
  static const Real sqrt6  = sqrt(6.);
  static const Real sqrt10 = sqrt(10.);
  static const Real sqrt14 = sqrt(14.);
  static const Real sqrt22 = sqrt(22.);
  static const Real sqrt26 = sqrt(26.);
}


template <>
Real FE<2,SZABAB>::shape(const ElemType,
			 const Order,
			 const unsigned int,
			 const Point&)
{
  std::cerr << "Szabo-Babuska polynomials require the element type\n"
	    << "because edge orientation is needed."
	    << std::endl;
  
  error();
  return 0.;
}



template <>
Real FE<2,SZABAB>::shape(const Elem* elem,
			 const Order order,
			 const unsigned int i,
			 const Point& p)
{
  assert (elem != NULL);
  
  const ElemType type = elem->type();

  // Declare that we are using our own special power function
  // from the Utility namespace.  This saves typing later.
  using Utility::pow;
  
  switch (order)
    {      
      // 1st & 2nd-order Szabo-Babuska.
    case FIRST:
    case SECOND:
      {
	switch (type)
	  {

	    // Szabo-Babuska shape functions on the triangle.
	  case TRI6:
	    {
	      const Real l1 = 1-p(0)-p(1);
	      const Real l2 = p(0);
	      const Real l3 = p(1);	      
	      
	      assert (i<6);
	      
	      switch (i)
		{
		case 0: return l1;		
		case 1: return l2;		  
		case 2: return l3;
		  
		case 3: return l1*l2*(-4.*sqrt6);			  
		case 4: return l2*l3*(-4.*sqrt6);		  
		case 5: return l3*l1*(-4.*sqrt6);	
		  
		default:
		  error();
		}
	    }
    

	    // Szabo-Babuska shape functions on the quadrilateral.
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
	      
	      return (FE<1,SZABAB>::shape(EDGE3, order, i0[i], xi)*
		      FE<1,SZABAB>::shape(EDGE3, order, i1[i], eta));
	      
	    }

	    
	  default:
	    error();
	  }
      }
	   

      // 3rd-order Szabo-Babuska.
    case THIRD:
      {
	switch (type)
	  {
	 
	    // Szabo-Babuska shape functions on the triangle.
	  case TRI6:
	    {
	      Real l1 = 1-p(0)-p(1);
	      Real l2 = p(0);
	      Real l3 = p(1);
	      
	      Real f=1;
	      
	      assert (i<10);


	      if (i==4 && (elem->node(0) > elem->node(1)))f=-1;
	      if (i==6 && (elem->node(1) > elem->node(2)))f=-1;     
	      if (i==8 && (elem->node(2) > elem->node(0)))f=-1;


	      switch (i)
		{
		  //nodal modes
		case 0: return l1;		
		case 1: return l2;		  
		case 2: return l3;
		  
		  //side modes
		case 3: return   l1*l2*(-4.*sqrt6);			  
		case 4: return f*l1*l2*(-4.*sqrt10)*(l2-l1);

		case 5: return   l2*l3*(-4.*sqrt6);		  
		case 6: return f*l2*l3*(-4.*sqrt10)*(l3-l2);
	  
		case 7: return   l3*l1*(-4.*sqrt6);
		case 8: return f*l3*l1*(-4.*sqrt10)*(l1-l3);	
		
		  //internal modes
		case 9: return l1*l2*l3;
	
		default:
		  error();
		}
	    }
	    

	    // Szabo-Babuska shape functions on the quadrilateral.
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
	      	
	      Real f=1.;
	      	      
	      // take care of edge orientation, this is needed at
	      // edge shapes with (y=0)-asymmetric 1D shapes, these have
	      // one 1D shape index being 0 or 1, the other one being odd and >=3
	      	     
	      switch(i)
		{
		case  5: // edge 0 nodes    			
		  if (elem->node(0) > elem->node(1))f = -1.;
		  break;
	      	case  7: // edge 1 nodes
		  if (elem->node(1) > elem->node(2))f = -1.;
		  break;
	        case  9: // edge 2 nodes
		  if (elem->node(3) > elem->node(2))f = -1.;
		  break;
	        case 11: // edge 3 nodes
		  if (elem->node(0) > elem->node(3))f = -1.;
		  break;
		}	      
	      
	      return f*(FE<1,SZABAB>::shape(EDGE3, order, i0[i], xi)*
			FE<1,SZABAB>::shape(EDGE3, order, i1[i], eta));
	    }

	  default:
	    error();
	  }
      }
	   

      

      // 4th-order Szabo-Babuska.
    case FOURTH:
      {
	switch (type)
	  {
	    // Szabo-Babuska shape functions on the triangle.
	  case TRI6:
	    {
	      Real l1 = 1-p(0)-p(1);
	      Real l2 = p(0);
	      Real l3 = p(1);
	      
	      Real f=1;
	      
	      assert (i<15);
	      
	      
	      if (i== 4 && (elem->node(0) > elem->node(1)))f=-1;
	      if (i== 7 && (elem->node(1) > elem->node(2)))f=-1;     
	      if (i==10 && (elem->node(2) > elem->node(0)))f=-1;
	      

	      switch (i)
		{
		  //nodal modes
		case  0: return l1;		
		case  1: return l2;		  
		case  2: return l3;
		  
		  //side modes
		case  3: return   l1*l2*(-4.*sqrt6);			  
		case  4: return f*l1*l2*(-4.*sqrt10)*(l2-l1);
		case  5: return   l1*l2*(-sqrt14)*(5.*pow<2>(l2-l1)-1);		  
		  
		case  6: return   l2*l3*(-4.*sqrt6);	
		case  7: return f*l2*l3*(-4.*sqrt10)*(l3-l2);	  
		case  8: return   l2*l3*(-sqrt14)*(5.*pow<2>(l3-l2)-1);	  
		  
		case  9: return   l3*l1*(-4.*sqrt6);		  
		case 10: return f*l3*l1*(-4.*sqrt10)*(l1-l3);		
		case 11: return   l3*l1*(-sqrt14)*(5.*pow<2>(l1-l3)-1);
		  
		  //internal modes
		case 12: return l1*l2*l3;
		  
		case 13: return l1*l2*l3*(l2-l1);
		case 14: return l1*l2*l3*(2*l3-1);
	
		default:
		  error();
		}
	    }
	  

	    // Szabo-Babuska shape functions on the quadrilateral.
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
	      
	      Real f=1.;	      	      
	      
	      switch(i)
		{
		case  5: // edge 0 nodes    			
		  if (elem->node(0) > elem->node(1))f = -1.;
		  break;
	      	case  8: // edge 1 nodes
		  if (elem->node(1) > elem->node(2))f = -1.;
		  break;
	        case 11: // edge 2 nodes
		  if (elem->node(3) > elem->node(2))f = -1.;
		  break;
	        case 14: // edge 3 nodes
		  if (elem->node(0) > elem->node(3))f = -1.;
		  break;
		}	      
   	      
	      return f*(FE<1,SZABAB>::shape(EDGE3, order, i0[i], xi)*
			FE<1,SZABAB>::shape(EDGE3, order, i1[i], eta));
	    }
	    
	  default:
	    error();
	  }
      }
      

      

      // 5th-order Szabo-Babuska.
    case FIFTH:
      {
	switch (type)
	  {
	    // Szabo-Babuska shape functions on the triangle.
	  case TRI6:
	    {
	      Real l1 = 1-p(0)-p(1);
	      Real l2 = p(0);
	      Real l3 = p(1);

	      const Real x=l2-l1;
	      const Real y=2.*l3-1;	
	      
	      Real f=1;	      
	      
	      assert (i<21);      	      
	      

	      if ((i== 4||i== 6) && (elem->node(0) > elem->node(1)))f=-1;
	      if ((i== 8||i==10) && (elem->node(1) > elem->node(2)))f=-1;     
	      if ((i==12||i==14) && (elem->node(2) > elem->node(0)))f=-1;
	      

	      switch (i)
		{
		  //nodal modes
		case  0: return l1;		
		case  1: return l2;		  
		case  2: return l3;
		  
		  //side modes
		case  3: return   l1*l2*(-4.*sqrt6);			  
		case  4: return f*l1*l2*(-4.*sqrt10)*(l2-l1);	
		case  5: return   l1*l2*(-sqrt14)*(5.*pow<2>(l2-l1)-1.);		  
		case  6: return f*l1*l2*(-sqrt2)*(21.*pow<3>(l2-l1)-9.*(l2-l1));		  
	  
		case  7: return   l2*l3*(-4.*sqrt6);
		case  8: return f*l2*l3*(-4.*sqrt10)*(l3-l2);	  
		case  9: return   l2*l3*(-sqrt14)*(5.*pow<2>(l3-l2)-1.);	  
		case 10: return f*l2*l3*(-sqrt2)*(21.*pow<3>(l3-l2)-9.*(l3-l2));	  
		  
		case 11: return   l3*l1*(-4.*sqrt6);		  
		case 12: return f*l3*l1*(-4.*sqrt10)*(l1-l3);  
		case 13: return   l3*l1*(-sqrt14)*(5.*pow<2>(l1-l3)-1.);		  
		case 14: return f*l3*l1*(-sqrt2)*(21.*pow<3>(l1-l3)-9.*(l1-l3));
		  
		  //internal modes
		case 15: return l1*l2*l3;
		  
		case 16: return l1*l2*l3*x;	
		case 17: return l1*l2*l3*y;
		  
		case 18: return l1*l2*l3*(1.5*pow<2>(x)-0.5);
		case 19: return l1*l2*l3*x*y;	
		case 20: return l1*l2*l3*(1.5*pow<2>(y)-0.5);	
		  
		  
		default:
		  error();
		}
	    } // case TRI6
	    
	    // Szabo-Babuska shape functions on the quadrilateral.
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
	      
	      Real f=1.;
	      
	      switch(i)
		{
		case  5: // edge 0 nodes
		case  7:    			
		  if (elem->node(0) > elem->node(1))f = -1.;
		  break;
		case  9: // edge 1 nodes
	      	case 11:	      	
		  if (elem->node(1) > elem->node(2))f = -1.;
		  break;
	        case 13: // edge 2 nodes
	        case 15:
		  if (elem->node(3) > elem->node(2))f = -1.;
		  break;
	        case 14: // edge 3 nodes
	        case 19:
		  if (elem->node(0) > elem->node(3))f = -1.;
		  break;
		}	     
	      
	      return f*(FE<1,SZABAB>::shape(EDGE3, order, i0[i], xi)*
			FE<1,SZABAB>::shape(EDGE3, order, i1[i], eta));	      

	    } // case QUAD8/QUAD9

	  default:
	    error();

	  } // switch type

      } // case FIFTH
      
      // 6th-order Szabo-Babuska.
    case SIXTH:
      {
	switch (type)
	  {
	    // Szabo-Babuska shape functions on the triangle.
	  case TRI6:
	    {
	      Real l1 = 1-p(0)-p(1);
	      Real l2 = p(0);
	      Real l3 = p(1);

	      const Real x=l2-l1;
	      const Real y=2.*l3-1;	      
	      
	      Real f=1;

	      assert (i<28);
	      
	      
	      if ((i== 4||i== 6) && (elem->node(0) > elem->node(1)))f=-1;
	      if ((i== 9||i==11) && (elem->node(1) > elem->node(2)))f=-1;     
	      if ((i==14||i==16) && (elem->node(2) > elem->node(0)))f=-1;

	      
	      switch (i)
		{
		  //nodal modes
		case  0: return l1;		
		case  1: return l2;		  
		case  2: return l3;
		  
		  //side modes
		case  3: return   l1*l2*(-4.*sqrt6);
		case  4: return f*l1*l2*(-4.*sqrt10)*(l2-l1);		  
		case  5: return   l1*l2*(-sqrt14)*(5.*pow<2>(l2-l1)-1.);		  
		case  6: return f*l1*l2*(-sqrt2)*(21.*pow<3>(l2-l1)-9.*(l2-l1));		  
		case  7: return   l1*l2*(-sqrt22)*(10.5*pow<4>(l2-l1)-7.*pow<2>(l2-l1)+0.5);
			  
		case  8: return   l2*l3*(-4.*sqrt6);
		case  9: return f*l2*l3*(-4.*sqrt10)*(l3-l2);	  
		case 10: return   l2*l3*(-sqrt14)*(5.*pow<2>(l3-l2)-1.);	  
		case 11: return f*l2*l3*(-sqrt2)*(21.*pow<3>(l3-l2)-9.*(l3-l2));	  
		case 12: return   l2*l3*(-sqrt22)*(10.5*pow<4>(l3-l2)-7.*pow<2>(l3-l2)+0.5);
		  
		case 13: return   l3*l1*(-4.*sqrt6);
		case 14: return f*l3*l1*(-4.*sqrt10)*(l1-l3);
		case 15: return   l3*l1*(-sqrt14)*(5.*pow<2>(l1-l3)-1.);		  
		case 16: return f*l3*l1*(-sqrt2)*(21.*pow<3>(l1-l3)-9.*(l1-l3));		  
		case 17: return   l3*l1*(-sqrt22)*(10.5*pow<4>(l1-l3)-7.*pow<2>(l1-l3)+0.5);
		  
		  
		  //internal modes
		case 18: return l1*l2*l3;
		  
		case 19: return l1*l2*l3*x;	
		case 20: return l1*l2*l3*y;
		  
		case 21: return l1*l2*l3*0.5*(3.*pow<2>(x)-1);
		case 22: return l1*l2*l3*x*y;	
		case 23: return l1*l2*l3*0.5*(3.*pow<2>(y)-1);
		  
		case 24: return l1*l2*l3*0.5*(5.*pow<3>(x)-3.*x);
		case 25: return l1*l2*l3*0.5*(3.*pow<2>(x)-1.)*y;
		case 26: return l1*l2*l3*0.5*(3.*pow<2>(y)-1.)*x;
		case 27: return l1*l2*l3*0.5*(5.*pow<3>(y)-3.*y);
		  
		  
		default:
		  error();
		}
	    } // case TRI6

	    // Szabo-Babuska shape functions on the quadrilateral.
	  case QUAD8:
	  case QUAD9:
	    {
	      // Compute quad shape functions as a tensor-product
	      const Real xi  = p(0);
	      const Real eta = p(1);
	      
	      assert (i < 49);
	      
	      //                                0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48
	      static const unsigned int i0[] = {0, 1, 1, 0, 2, 3, 4, 5, 6, 1, 1, 1, 1, 1, 2, 3, 4, 5, 6, 0, 0, 0, 0, 0, 2, 3, 4, 5, 6, 2, 3, 4, 5, 6, 2, 3, 4, 5, 6, 2, 3, 4, 5, 6, 2, 3, 4, 5, 6};
	      static const unsigned int i1[] = {0, 0, 1, 1, 0, 0, 0, 0, 0, 2, 3, 4, 5, 6, 1, 1, 1, 1, 1, 2, 3, 4, 5, 6, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6};
	      
	      Real f=1.;
	      
	      switch(i)
		{
		case  5: // edge 0 nodes
		case  7:    			
		  if (elem->node(0) > elem->node(1))f = -1.;
		  break;
	      	case 10: // edge 1 nodes
	      	case 12:	      	
		  if (elem->node(1) > elem->node(2))f = -1.;
		  break;
	        case 15: // edge 2 nodes
	        case 17:
		  if (elem->node(3) > elem->node(2))f = -1.;
		  break;
	        case 20: // edge 3 nodes
	        case 22:
		  if (elem->node(0) > elem->node(3))f = -1.;
		  break;
		}	     
	      
	      return f*(FE<1,SZABAB>::shape(EDGE3, order, i0[i], xi)*
			FE<1,SZABAB>::shape(EDGE3, order, i1[i], eta));

	    } // case QUAD8/QUAD9

	  default:
	    error();

	  } // switch type

      } // case SIXTH


      // 7th-order Szabo-Babuska.
    case SEVENTH:
      {
	switch (type)
	  {
	    // Szabo-Babuska shape functions on the triangle.
	  case TRI6:
	    {
	      Real l1 = 1-p(0)-p(1);
	      Real l2 = p(0);
	      Real l3 = p(1);	      
	      
	      const Real x=l2-l1;
	      const Real y=2.*l3-1.;	      

	      Real f=1;

	      assert (i<36);
	      
	      
	      if ((i>= 4&&i<= 8) && (elem->node(0) > elem->node(1)))f=-1;
	      if ((i>=10&&i<=14) && (elem->node(1) > elem->node(2)))f=-1;     
	      if ((i>=16&&i<=20) && (elem->node(2) > elem->node(0)))f=-1;

	      
	      switch (i)
		{
		  //nodal modes
		case  0: return l1;		
		case  1: return l2;		  
		case  2: return l3;
		  
		  //side modes
		case  3: return   l1*l2*(-4.*sqrt6);
		case  4: return f*l1*l2*(-4.*sqrt10)*(l2-l1);		  
		case  5: return   l1*l2*(-sqrt14)*(5.*pow<2>(l2-l1)-1.);		  
		case  6: return f*l1*l2*(-sqrt2)*(21.*pow<3>(l2-l1)-9.*(l2-l1));		  
		case  7: return   l1*l2*(-sqrt22)*(10.5*pow<4>(l2-l1)-7.*pow<2>(l2-l1)+0.5);
		case  8: return f*l1*l2*(-sqrt26)*(16.5*pow<5>(l2-l1)-15.*pow<3>(l2-l1)+2.5*(l2-l1));
			  
		case  9: return   l2*l3*(-4.*sqrt6);
		case 10: return f*l2*l3*(-4.*sqrt10)*(l3-l2);	  
		case 11: return   l2*l3*(-sqrt14)*(5.*pow<2>(l3-l2)-1.);	  
		case 12: return f*l2*l3*(-sqrt2)*(21.*pow<3>(l3-l2)-9.*(l3-l2));	  
		case 13: return   l2*l3*(-sqrt22)*(10.5*pow<4>(l3-l2)-7.*pow<2>(l3-l2)+0.5);
		case 14: return f*l2*l3*(-sqrt26)*(16.5*pow<5>(l3-l2)-15.*pow<3>(l3-l2)+2.5*(l3-l2));
		  
		case 15: return   l3*l1*(-4.*sqrt6);
		case 16: return f*l3*l1*(-4.*sqrt10)*(l1-l3);
		case 17: return   l3*l1*(-sqrt14)*(5.*pow<2>(l1-l3)-1.);		  
		case 18: return f*l3*l1*(-sqrt2)*(21.*pow<3>(l1-l3)-9.*(l1-l3));		  
		case 19: return   l3*l1*(-sqrt22)*(10.5*pow<4>(l1-l3)-7.*pow<2>(l1-l3)+0.5);
		case 20: return f*l3*l1*(-sqrt26)*(16.5*pow<5>(l1-l3)-15.*pow<3>(l1-l3)+2.5*(l1-l3));
		  
		  //internal modes
		case 21: return l1*l2*l3;
		  
		case 22: return l1*l2*l3*x;	
		case 23: return l1*l2*l3*y;
		  
		case 24: return l1*l2*l3*0.5*(3.*pow<2>(x)-1.);
		case 25: return l1*l2*l3*x*y;	
		case 26: return l1*l2*l3*0.5*(3.*pow<2>(y)-1.);
		  
		case 27: return l1*l2*l3*0.5*(5.*pow<3>(x)-3.*x);	
		case 28: return l1*l2*l3*0.5*(3.*pow<2>(x)-1.)*y;
		case 29: return l1*l2*l3*0.5*(3.*pow<2>(y)-1.)*x;
		case 30: return l1*l2*l3*0.5*(5.*pow<3>(y)-3.*y);
		  
		case 31: return l1*l2*l3*0.125*(63.*pow<5>(x)-70.*pow<3>(x)+15.*x);
		case 32: return l1*l2*l3*0.5*(5.*pow<3>(x)-3*x)*y;
		case 33: return l1*l2*l3*0.25*(3.*pow<2>(x)-1)*(3*pow<2>(y)-1);
		case 34: return l1*l2*l3*0.5*(5.*pow<3>(y)-3*y)*x;
		case 35: return l1*l2*l3*0.125*(63.*pow<5>(y)-70.*pow<3>(y)+15.*y);
		  
		default:
		  error();
		}
	    } // case TRI6
	    
	    // Szabo-Babuska shape functions on the quadrilateral.
	  case QUAD8:
	  case QUAD9:
	    {
	      // Compute quad shape functions as a tensor-product
	      const Real xi  = p(0);
	      const Real eta = p(1);
	      
	      assert (i < 64);
	      
	      //                                0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63
	      static const unsigned int i0[] = {0, 1, 1, 0, 2, 3, 4, 5, 6, 7, 1, 1, 1, 1, 1, 1, 2, 3, 4, 5, 6, 7, 0, 0, 0, 0, 0, 0, 2, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7};
	      static const unsigned int i1[] = {0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 2, 3, 4, 5, 6, 7, 1, 1, 1, 1, 1, 1, 2, 3, 4, 5, 6, 7, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7};
	      
	      Real f=1.;
	      
	      switch(i)
		{
		case  5: // edge 0 nodes
		case  7:
		case  9:    			
		  if (elem->node(0) > elem->node(1))f = -1.;
		  break;
	      	case 11: // edge 1 nodes
	      	case 13:	      	
		case 15:
		  if (elem->node(1) > elem->node(2))f = -1.;
		  break;
	        case 17: // edge 2 nodes
	        case 19:
		case 21:
		  if (elem->node(3) > elem->node(2))f = -1.;
		  break;
	        case 23: // edge 3 nodes
	        case 25:
		case 27:
		  if (elem->node(0) > elem->node(3))f = -1.;
		  break;
		}	     
	      
	      return f*(FE<1,SZABAB>::shape(EDGE3, order, i0[i], xi)*
			FE<1,SZABAB>::shape(EDGE3, order, i1[i], eta));	      

	    } // case QUAD8/QUAD9
	    
	  default:
	    error();

	  } // switch type

      } // case SEVENTH


      // by default throw an error
    default:
      std::cerr << "ERROR: Unsupported polynomial order!" << std::endl;
      error();

    } // switch order
}

      



template <>
Real FE<2,SZABAB>::shape_deriv(const ElemType,
				   const Order,			    
				   const unsigned int,
				   const unsigned int,
				   const Point&)
{
  std::cerr << "Szabo-Babuska polynomials require the element type\n"
	    << "because edge orientation is needed."
	    << std::endl;

  error();
  return 0.;
}



template <>
Real FE<2,SZABAB>::shape_deriv(const Elem* elem,
			       const Order order,
			       const unsigned int i,
			       const unsigned int j,
			       const Point& p)
{
  assert (elem != NULL);

  const ElemType type = elem->type();
  
  switch (order)
    {

      // 1st & 2nd-order Szabo-Babuska.
    case FIRST:
    case SECOND:
      {
	switch (type)
	  {
	    
	    // Szabo-Babuska shape functions on the triangle.
	  case TRI6:
	    {
	      // Here we use finite differences to compute the derivatives!
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

		    return (FE<2,SZABAB>::shape(elem, order, i, pp) -
			    FE<2,SZABAB>::shape(elem, order, i, pm))/2./eps;
		  }

		  // d()/deta
		case 1:
		  {
		    const Point pp(p(0), p(1)+eps);
		    const Point pm(p(0), p(1)-eps);

		    return (FE<2,SZABAB>::shape(elem, order, i, pp) -
			    FE<2,SZABAB>::shape(elem, order, i, pm))/2./eps;
		  }
		  

		default:
		  error();
		}
	    }

	    

	    // Szabo-Babuska shape functions on the quadrilateral.
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
		  return (FE<1,SZABAB>::shape_deriv(EDGE3, order, i0[i], 0, xi)*
			  FE<1,SZABAB>::shape      (EDGE3, order, i1[i],    eta));

		  // d()/deta
		case 1:		      
		  return (FE<1,SZABAB>::shape      (EDGE3, order, i0[i],    xi)*
			  FE<1,SZABAB>::shape_deriv(EDGE3, order, i1[i], 0, eta));

		default:
		  error();
		}	      
	    }

	  default:
	    error();
	  }
      }

      

      // 3rd-order Szabo-Babuska.
    case THIRD:
      {
	switch (type)
	  {
	    // Szabo-Babuska shape functions on the triangle.
	  case TRI6:
	    {
	      // Here we use finite differences to compute the derivatives!
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

		    return (FE<2,SZABAB>::shape(elem, order, i, pp) -
			    FE<2,SZABAB>::shape(elem, order, i, pm))/2./eps;
		  }

		  // d()/deta
		case 1:
		  {
		    const Point pp(p(0), p(1)+eps);
		    const Point pm(p(0), p(1)-eps);

		    return (FE<2,SZABAB>::shape(elem, order, i, pp) -
			    FE<2,SZABAB>::shape(elem, order, i, pm))/2./eps;
		  }
		  

		default:
		  error();
		}
	    }
	  

	    // Szabo-Babuska shape functions on the quadrilateral.
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
	      
	      Real f=1.;

	      switch(i)
		{
	      	case  5: // edge 0 nodes    			
		  if (elem->node(0) > elem->node(1))f = -1.;
		  break;
	      	case  7: // edge 1 nodes
		  if (elem->node(1) > elem->node(2))f = -1.;
		  break;
	        case  9: // edge 2 nodes
		  if (elem->node(3) > elem->node(2))f = -1.;
		  break;
	        case 11: // edge 3 nodes
		  if (elem->node(0) > elem->node(3))f = -1.;
		  break;
		}	      

	      
	      switch (j)
		{
		  // d()/dxi
		case 0:		  		  
		  return f*(FE<1,SZABAB>::shape_deriv(EDGE3, order, i0[i], 0, xi)*
			    FE<1,SZABAB>::shape      (EDGE3, order, i1[i],    eta));
	      
		  // d()/deta
		case 1:		  		  
		  return f*(FE<1,SZABAB>::shape      (EDGE3, order, i0[i],    xi)*
			    FE<1,SZABAB>::shape_deriv(EDGE3, order, i1[i], 0, eta));

		default:
		  error();
		}
	    }

	  default:
	    error();
	  }
      }
	   

      

      // 4th-order Szabo-Babuska.
    case FOURTH:
      {
	switch (type)
	  {
	 
	    // Szabo-Babuska shape functions on the triangle.
	  case TRI6:
	    {
	      // Here we use finite differences to compute the derivatives!
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

		    return (FE<2,SZABAB>::shape(elem, order, i, pp) -
			    FE<2,SZABAB>::shape(elem, order, i, pm))/2./eps;
		  }

		  // d()/deta
		case 1:
		  {
		    const Point pp(p(0), p(1)+eps);
		    const Point pm(p(0), p(1)-eps);

		    return (FE<2,SZABAB>::shape(elem, order, i, pp) -
			    FE<2,SZABAB>::shape(elem, order, i, pm))/2./eps;
		  }
		  

		default:
		  error();
		}
	    }
	
	    

	    // Szabo-Babuska shape functions on the quadrilateral.
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
	      
	      Real f=1.;

	      switch(i)
		{
	      	case  5: // edge 0 nodes    			
		  if (elem->node(0) > elem->node(1))f = -1.;
		  break;
	      	case  8: // edge 1 nodes
		  if (elem->node(1) > elem->node(2))f = -1.;
		  break;
	        case 11: // edge 2 nodes
		  if (elem->node(3) > elem->node(2))f = -1.;
		  break;
	        case 14: // edge 3 nodes
		  if (elem->node(0) > elem->node(3))f = -1.;
		  break;
		}	      

	      
	      switch (j)
		{
		  // d()/dxi
		case 0:		  		  
		  return f*(FE<1,SZABAB>::shape_deriv(EDGE3, order, i0[i], 0, xi)*
			    FE<1,SZABAB>::shape      (EDGE3, order, i1[i],    eta));
	      
		  // d()/deta
		case 1:		  		  
		  return f*(FE<1,SZABAB>::shape      (EDGE3, order, i0[i],    xi)*
			    FE<1,SZABAB>::shape_deriv(EDGE3, order, i1[i], 0, eta));

		default:
		  error();
		}
	    }

	  default:
	    error();
	  }
      }
	   

      

      // 5th-order Szabo-Babuska.
    case FIFTH:
      {
	// Szabo-Babuska shape functions on the quadrilateral.
	switch (type)
	  {
	 
	    // Szabo-Babuska shape functions on the triangle.
	  case TRI6:
	    {
	      // Here we use finite differences to compute the derivatives!
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

		    return (FE<2,SZABAB>::shape(elem, order, i, pp) -
			    FE<2,SZABAB>::shape(elem, order, i, pm))/2./eps;
		  }

		  // d()/deta
		case 1:
		  {
		    const Point pp(p(0), p(1)+eps);
		    const Point pm(p(0), p(1)-eps);

		    return (FE<2,SZABAB>::shape(elem, order, i, pp) -
			    FE<2,SZABAB>::shape(elem, order, i, pm))/2./eps;
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
	      
	      Real f=1.;
	     
	      switch(i)
		{
		case  5: // edge 0 nodes
	      	case  7:    			
		  if (elem->node(0) > elem->node(1))f = -1.;
		  break;
	      	case  9: // edge 1 nodes
	      	case 11:	      	
		  if (elem->node(1) > elem->node(2))f = -1.;
		  break;
	        case 13: // edge 2 nodes
	        case 15:
		  if (elem->node(3) > elem->node(2))f = -1.;
		  break;
	        case 14: // edge 3 nodes
	        case 19:
		  if (elem->node(0) > elem->node(3))f = -1.;
		  break;
		}	     
	      
	      
	      switch (j)
		{
		  // d()/dxi
		case 0:		  		  
		  return f*(FE<1,SZABAB>::shape_deriv(EDGE3, order, i0[i], 0, xi)*
			    FE<1,SZABAB>::shape      (EDGE3, order, i1[i],    eta));
	      
		  // d()/deta
		case 1:		  		  
		  return f*(FE<1,SZABAB>::shape      (EDGE3, order, i0[i],    xi)*
			    FE<1,SZABAB>::shape_deriv(EDGE3, order, i1[i], 0, eta));

		default:
		  error();
		}
	    }

	  default:
	    error();
	  }
      }


    // 6th-order Szabo-Babuska.
    case SIXTH:
      {
	// Szabo-Babuska shape functions on the quadrilateral.
	switch (type)
	  {
	 
	    // Szabo-Babuska shape functions on the triangle.
	  case TRI6:
	    {
	      // Here we use finite differences to compute the derivatives!
	      const Real eps = 1.e-6;
	      
	      assert (i < 28);
	      assert (j < 2);
	      
	      switch (j)
		{
		  //  d()/dxi
		case 0:
		  {
		    const Point pp(p(0)+eps, p(1));
		    const Point pm(p(0)-eps, p(1));

		    return (FE<2,SZABAB>::shape(elem, order, i, pp) -
			    FE<2,SZABAB>::shape(elem, order, i, pm))/2./eps;
		  }

		  // d()/deta
		case 1:
		  {
		    const Point pp(p(0), p(1)+eps);
		    const Point pm(p(0), p(1)-eps);

		    return (FE<2,SZABAB>::shape(elem, order, i, pp) -
			    FE<2,SZABAB>::shape(elem, order, i, pm))/2./eps;
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
	      
	      assert (i < 49);

	      //                                0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48
	      static const unsigned int i0[] = {0, 1, 1, 0, 2, 3, 4, 5, 6, 1, 1, 1, 1, 1, 2, 3, 4, 5, 6, 0, 0, 0, 0, 0, 2, 3, 4, 5, 6, 2, 3, 4, 5, 6, 2, 3, 4, 5, 6, 2, 3, 4, 5, 6, 2, 3, 4, 5, 6};
	      static const unsigned int i1[] = {0, 0, 1, 1, 0, 0, 0, 0, 0, 2, 3, 4, 5, 6, 1, 1, 1, 1, 1, 2, 3, 4, 5, 6, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6};
	      
	      Real f=1.;
	     
	      switch(i)
		{
	      	case  5: // edge 0 nodes
	      	case  7:    			
		  if (elem->node(0) > elem->node(1))f = -1.;
		  break;
	      	case 10: // edge 1 nodes
	      	case 12:	      	
		  if (elem->node(1) > elem->node(2))f = -1.;
		  break;
	        case 15: // edge 2 nodes
	        case 17:
		  if (elem->node(3) > elem->node(2))f = -1.;
		  break;
	        case 20: // edge 3 nodes
	        case 22:
		  if (elem->node(0) > elem->node(3))f = -1.;
		  break;
		}	       

	      
	      switch (j)
		{
		  // d()/dxi
		case 0:		  		  
		  return f*(FE<1,SZABAB>::shape_deriv(EDGE3, order, i0[i], 0, xi)*
			    FE<1,SZABAB>::shape      (EDGE3, order, i1[i],    eta));
	      
		  // d()/deta
		case 1:		  		  
		  return f*(FE<1,SZABAB>::shape      (EDGE3, order, i0[i],    xi)*
			    FE<1,SZABAB>::shape_deriv(EDGE3, order, i1[i], 0, eta));

		default:
		  error();
		}
	    }

	  default:
	    error();
	  }
      }


    // 7th-order Szabo-Babuska.
    case SEVENTH:
      {
	// Szabo-Babuska shape functions on the quadrilateral.
	switch (type)
	  {
	 
	    // Szabo-Babuska shape functions on the triangle.
	  case TRI6:
	    {
	      // Here we use finite differences to compute the derivatives!
	      const Real eps = 1.e-6;
	      
	      assert (i < 36);
	      assert (j < 2);
	      
	      switch (j)
		{
		  //  d()/dxi
		case 0:
		  {
		    const Point pp(p(0)+eps, p(1));
		    const Point pm(p(0)-eps, p(1));

		    return (FE<2,SZABAB>::shape(elem, order, i, pp) -
			    FE<2,SZABAB>::shape(elem, order, i, pm))/2./eps;
		  }

		  // d()/deta
		case 1:
		  {
		    const Point pp(p(0), p(1)+eps);
		    const Point pm(p(0), p(1)-eps);

		    return (FE<2,SZABAB>::shape(elem, order, i, pp) -
			    FE<2,SZABAB>::shape(elem, order, i, pm))/2./eps;
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
	      
	      assert (i < 64);

	      //                                0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63
	      static const unsigned int i0[] = {0, 1, 1, 0, 2, 3, 4, 5, 6, 7, 1, 1, 1, 1, 1, 1, 2, 3, 4, 5, 6, 7, 0, 0, 0, 0, 0, 0, 2, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7};
	      static const unsigned int i1[] = {0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 2, 3, 4, 5, 6, 7, 1, 1, 1, 1, 1, 1, 2, 3, 4, 5, 6, 7, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7};
	      
	      Real f=1.;
	     
	      switch(i)
		{
	      	case  5: // edge 0 nodes
	      	case  7:
		case  9:    			
		  if (elem->node(0) > elem->node(1))f = -1.;
		  break;
	      	case 11: // edge 1 nodes
	      	case 13:	      	
		case 15:
		  if (elem->node(1) > elem->node(2))f = -1.;
		  break;
	        case 17: // edge 2 nodes
	        case 19:
		case 21:
		  if (elem->node(3) > elem->node(2))f = -1.;
		  break;
	        case 23: // edge 3 nodes
	        case 25:
		case 27:
		  if (elem->node(0) > elem->node(3))f = -1.;
		  break;
		}	     	       
	      
	      
	      switch (j)
		{
		  // d()/dxi
		case 0:		  		  
		  return f*(FE<1,SZABAB>::shape_deriv(EDGE3, order, i0[i], 0, xi)*
			    FE<1,SZABAB>::shape      (EDGE3, order, i1[i],    eta));
		  
		  // d()/deta
		case 1:		  		  
		  return f*(FE<1,SZABAB>::shape      (EDGE3, order, i0[i],    xi)*
			    FE<1,SZABAB>::shape_deriv(EDGE3, order, i1[i], 0, eta));
		  
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


#endif	// ENABLE_HIGHER_ORDER_SHAPES


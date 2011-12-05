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
#include "elem.h"
#include "number_lookups.h"
#include "utility.h"


namespace libMesh
{


template <>
Real FE<2,BERNSTEIN>::shape(const ElemType,
			    const Order,
			    const unsigned int,
			    const Point&)
{
  libMesh::err << "Bernstein polynomials require the element type\n"
	        << "because edge orientation is needed."
	        << std::endl;

  libmesh_error();
  return 0.;
}



template <>
Real FE<2,BERNSTEIN>::shape(const Elem* elem,
			    const Order order,
			    const unsigned int i,
			    const Point& p)
{
  libmesh_assert (elem != NULL);

  const ElemType type = elem->type();

  const Order totalorder = static_cast<Order>(order + elem->p_level());

  // Declare that we are using our own special power function
  // from the Utility namespace.  This saves typing later.
  using Utility::pow;

  switch (type)
    {
    // Hierarchic shape functions on the quadrilateral.
    case QUAD4:
    case QUAD9:
      {
        // Compute quad shape functions as a tensor-product
        const Real xi  = p(0);
        const Real eta = p(1);

        libmesh_assert (i < (totalorder+1u)*(totalorder+1u));

// Example i, i0, i1 values for totalorder = 5:
//                                    0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35
//  static const unsigned int i0[] = {0, 1, 1, 0, 2, 3, 4, 5, 1, 1, 1, 1, 2, 3, 4, 5, 0, 0, 0, 0, 2, 3, 3, 2, 4, 4, 4, 3, 2, 5, 5, 5, 5, 4, 3, 2};
//  static const unsigned int i1[] = {0, 0, 1, 1, 0, 0, 0, 0, 2, 3, 4, 5, 1, 1, 1, 1, 2, 3, 4, 5, 2, 2, 3, 3, 2, 3, 4, 4, 4, 2, 3, 4, 5, 5, 5, 5};

        unsigned int i0, i1;

        // Vertex DoFs
        if (i == 0)
          { i0 = 0; i1 = 0; }
        else if (i == 1)
          { i0 = 1; i1 = 0; }
        else if (i == 2)
          { i0 = 1; i1 = 1; }
        else if (i == 3)
          { i0 = 0; i1 = 1; }


        // Edge DoFs
        else if (i < totalorder + 3u)
          { i0 = i - 2; i1 = 0; }
        else if (i < 2u*totalorder + 2)
          { i0 = 1; i1 = i - totalorder - 1; }
        else if (i < 3u*totalorder + 1)
          { i0 = i - 2u*totalorder; i1 = 1; }
        else if (i < 4u*totalorder)
          { i0 = 0; i1 = i - 3u*totalorder + 1; }
        // Interior DoFs. Use Roy's number look up
        else
          {
	    unsigned int basisnum = i - 4*totalorder;
	    i0 = square_number_column[basisnum] + 2;
	    i1 = square_number_row[basisnum] + 2;
          }


	// Flip odd degree of freedom values if necessary
	// to keep continuity on sides.
	if     ((i>= 4                 && i<= 4+  totalorder-2u) && elem->point(0) > elem->point(1)) i0=totalorder+2-i0;	//
	else if((i>= 4+  totalorder-1u && i<= 4+2*totalorder-3u) && elem->point(1) > elem->point(2)) i1=totalorder+2-i1;
	else if((i>= 4+2*totalorder-2u && i<= 4+3*totalorder-4u) && elem->point(3) > elem->point(2)) i0=totalorder+2-i0;
	else if((i>= 4+3*totalorder-3u && i<= 4+4*totalorder-5u) && elem->point(0) > elem->point(3)) i1=totalorder+2-i1;


        return (FE<1,BERNSTEIN>::shape(EDGE3, totalorder, i0, xi)*
		FE<1,BERNSTEIN>::shape(EDGE3, totalorder, i1, eta));
      }
      // handle serendipity QUAD8 element separately
    case QUAD8:
      {
	libmesh_assert (totalorder < 3);

	const Real xi  = p(0);
	const Real eta = p(1);

	libmesh_assert (i < 8);

	//                                0  1  2  3  4  5  6  7  8
	static const unsigned int i0[] = {0, 1, 1, 0, 2, 1, 2, 0, 2};
	static const unsigned int i1[] = {0, 0, 1, 1, 0, 2, 1, 2, 2};
	static const Real scal[] = {-0.25, -0.25, -0.25, -0.25, 0.5, 0.5, 0.5, 0.5};

	//B_t,i0(i)|xi * B_s,i1(i)|eta + scal(i) * B_t,2|xi * B_t,2|eta
	return (FE<1,BERNSTEIN>::shape(EDGE3, totalorder, i0[i], xi)*
		FE<1,BERNSTEIN>::shape(EDGE3, totalorder, i1[i], eta)
		+scal[i]*
		FE<1,BERNSTEIN>::shape(EDGE3, totalorder, i0[8], xi)*
		FE<1,BERNSTEIN>::shape(EDGE3, totalorder, i1[8], eta));

      }

    case TRI3:
      libmesh_assert (totalorder<2);
    case TRI6:
      switch (totalorder)
	{
	case FIRST:
	  {
	    const Real x=p(0);
	    const Real y=p(1);
	    const Real r=1.-x-y;

	    libmesh_assert(i<3);

	    switch(i)
	      {
	      case 0: return r;  //f0,0,1
	      case 1: return x;  //f0,1,1
	      case 2: return y;  //f1,0,1

	      default: libmesh_error(); return 0;
	      }
	  }
	case SECOND:
	  {
	    const Real x=p(0);
	    const Real y=p(1);
	    const Real r=1.-x-y;

	    libmesh_assert(i<6);

	    switch(i)
	      {
	      case 0: return r*r;
	      case 1: return x*x;
	      case 2: return y*y;

	      case 3: return 2.*x*r;
	      case 4: return 2.*x*y;
	      case 5: return 2.*r*y;

	      default: libmesh_error(); return 0;
	      }
	  }
	case THIRD:
	  {
	    const Real x=p(0);
	    const Real y=p(1);
	    const Real r=1.-x-y;
	    libmesh_assert(i<10);

	    unsigned int shape=i;


	    if((i==3||i==4) && elem->point(0) > elem->point(1)) shape=7-i;
	    if((i==5||i==6) && elem->point(1) > elem->point(2)) shape=11-i;
	    if((i==7||i==8) && elem->point(0) > elem->point(2)) shape=15-i;

	    switch(shape)
	      {
	      case 0: return r*r*r;
	      case 1: return x*x*x;
	      case 2: return y*y*y;

	      case 3: return 3.*x*r*r;
	      case 4: return 3.*x*x*r;

	      case 5: return 3.*y*x*x;
	      case 6: return 3.*y*y*x;

	      case 7: return 3.*y*r*r;
	      case 8: return 3.*y*y*r;

	      case 9: return 6.*x*y*r;

	      default: libmesh_error(); return 0;
	      }
	  }
	case FOURTH:
	  {
	    const Real x=p(0);
	    const Real y=p(1);
	    const Real r=1-x-y;
	    unsigned int shape=i;

	    libmesh_assert(i<15);

	    if((i==3||i== 5) && elem->point(0) > elem->point(1))shape=8-i;
	    if((i==6||i== 8) && elem->point(1) > elem->point(2))shape=14-i;
	    if((i==9||i==11) && elem->point(0) > elem->point(2))shape=20-i;


	    switch(shape)
	      {
		// point functions
	      case  0: return r*r*r*r;
	      case  1: return x*x*x*x;
	      case  2: return y*y*y*y;

		// edge functions
	      case  3: return 4.*x*r*r*r;
	      case  4: return 6.*x*x*r*r;
	      case  5: return 4.*x*x*x*r;

	      case  6: return 4.*y*x*x*x;
	      case  7: return 6.*y*y*x*x;
	      case  8: return 4.*y*y*y*x;

	      case  9: return 4.*y*r*r*r;
	      case 10: return 6.*y*y*r*r;
	      case 11: return 4.*y*y*y*r;

		// inner functions
	      case 12: return 12.*x*y*r*r;
	      case 13: return 12.*x*x*y*r;
	      case 14: return 12.*x*y*y*r;

	      default: libmesh_error(); return 0;
	      }
	  }
	case FIFTH:
	  {
	    const Real x=p(0);
	    const Real y=p(1);
	    const Real r=1-x-y;
	    unsigned int shape=i;

	    libmesh_assert(i<21);

	    if((i>= 3&&i<= 6) && elem->point(0) > elem->point(1))shape=9-i;
	    if((i>= 7&&i<=10) && elem->point(1) > elem->point(2))shape=17-i;
	    if((i>=11&&i<=14) && elem->point(0) > elem->point(2))shape=25-i;

	    switch(shape)
	      {
		//point functions
	      case  0: return pow<5>(r);
	      case  1: return pow<5>(x);
	      case  2: return pow<5>(y);

		//edge functions
	      case  3: return  5.*x        *pow<4>(r);
	      case  4: return 10.*pow<2>(x)*pow<3>(r);
	      case  5: return 10.*pow<3>(x)*pow<2>(r);
	      case  6: return  5.*pow<4>(x)*r;

	      case  7: return  5.*y	   *pow<4>(x);
	      case  8: return 10.*pow<2>(y)*pow<3>(x);
	      case  9: return 10.*pow<3>(y)*pow<2>(x);
	      case 10: return  5.*pow<4>(y)*x;

	      case 11: return  5.*y	   *pow<4>(r);
	      case 12: return 10.*pow<2>(y)*pow<3>(r);
	      case 13: return 10.*pow<3>(y)*pow<2>(r);
	      case 14: return  5.*pow<4>(y)*r;

		//inner functions
	      case 15: return 20.*x*y*pow<3>(r);
	      case 16: return 30.*x*pow<2>(y)*pow<2>(r);
	      case 17: return 30.*pow<2>(x)*y*pow<2>(r);
	      case 18: return 20.*x*pow<3>(y)*r;
	      case 19: return 20.*pow<3>(x)*y*r;
	      case 20: return 30.*pow<2>(x)*pow<2>(y)*r;

	      default: libmesh_error(); return 0;
	      }
	  }
	case SIXTH:
	  {
	    const Real x=p(0);
	    const Real y=p(1);
	    const Real r=1-x-y;
	    unsigned int shape=i;

	    libmesh_assert(i<28);

	    if((i>= 3&&i<= 7) && elem->point(0) > elem->point(1))shape=10-i;
	    if((i>= 8&&i<=12) && elem->point(1) > elem->point(2))shape=20-i;
	    if((i>=13&&i<=17) && elem->point(0) > elem->point(2))shape=30-i;

	    switch(shape)
	      {
		//point functions
	      case  0: return pow<6>(r);
	      case  1: return pow<6>(x);
	      case  2: return pow<6>(y);

		//edge functions
	      case  3: return  6.*x        *pow<5>(r);
	      case  4: return 15.*pow<2>(x)*pow<4>(r);
	      case  5: return 20.*pow<3>(x)*pow<3>(r);
	      case  6: return 15.*pow<4>(x)*pow<2>(r);
	      case  7: return  6.*pow<5>(x)*r;

	      case  8: return  6.*y        *pow<5>(x);
	      case  9: return 15.*pow<2>(y)*pow<4>(x);
	      case 10: return 20.*pow<3>(y)*pow<3>(x);
	      case 11: return 15.*pow<4>(y)*pow<2>(x);
	      case 12: return  6.*pow<5>(y)*x;

	      case 13: return  6.*y        *pow<5>(r);
	      case 14: return 15.*pow<2>(y)*pow<4>(r);
	      case 15: return 20.*pow<3>(y)*pow<3>(r);
	      case 16: return 15.*pow<4>(y)*pow<2>(r);
	      case 17: return  6.*pow<5>(y)*r;

		//inner functions
	      case 18: return 30.*x*y*pow<4>(r);
	      case 19: return 60.*x*pow<2>(y)*pow<3>(r);
	      case 20: return 60.*  pow<2>(x)*y*pow<3>(r);
	      case 21: return 60.*x*pow<3>(y)*pow<2>(r);
	      case 22: return 60.*pow<3>(x)*y*pow<2>(r);
	      case 23: return 90.*pow<2>(x)*pow<2>(y)*pow<2>(r);
	      case 24: return 30.*x*pow<4>(y)*r;
	      case 25: return 60.*pow<2>(x)*pow<3>(y)*r;
	      case 26: return 60.*pow<3>(x)*pow<2>(y)*r;
	      case 27: return 30.*pow<4>(x)*y*r;

	      default: libmesh_error(); return 0;
	      } // switch shape
	  } // case TRI6
	default:
	  {
	    libMesh::err << "ERROR: element order!" << std::endl;
	    libmesh_error();
	  }


	} // switch order


    default:
      {
	libMesh::err << "ERROR: Unsupported element type!" << std::endl;
	libmesh_error();
      }

    } // switch type


  // old code
//   switch (totalorder)
//     {

//     case FIRST:

//       switch(type)
// 	{

// 	case TRI6:
// 	  {
// 	    const Real x=p(0);
// 	    const Real y=p(1);
// 	    const Real r=1.-x-y;

// 	    libmesh_assert(i<3);

// 	    switch(i)
// 	      {
// 	      case 0: return r;  //f0,0,1
// 	      case 1: return x;      //f0,1,1
// 	      case 2: return y;      //f1,0,1
// 	      }
// 	  }

// 	case QUAD8:
// 	case QUAD9:
// 	  {
// 	    // Compute quad shape functions as a tensor-product
// 	    const Real xi  = p(0);
// 	    const Real eta = p(1);

// 	    libmesh_assert (i < 4);

// 	    //                                0  1  2  3
// 	    static const unsigned int i0[] = {0, 1, 1, 0};
// 	    static const unsigned int i1[] = {0, 0, 1, 1};

// 	    return (FE<1,BERNSTEIN>::shape(EDGE3, totalorder, i0[i], xi)*
// 		    FE<1,BERNSTEIN>::shape(EDGE3, totalorder, i1[i], eta));

// 	  }
// 	default:
// 	  libmesh_error();
// 	}

//     case SECOND:
//       switch(type)
// 	{
// 	case TRI6:
// 	  {
// 	    const Real x=p(0);
// 	    const Real y=p(1);
// 	    const Real r=1.-x-y;

// 	    libmesh_assert(i<6);

// 	    switch(i)
// 	      {
// 	      case 0: return r*r;
// 	      case 1: return x*x;
// 	      case 2: return y*y;

// 	      case 3: return 2.*x*r;
// 	      case 4: return 2.*x*y;
// 	      case 5: return 2.*r*y;
// 	      }
// 	  }

// 	  // Bernstein shape functions on the 8-noded quadrilateral.
// 	case QUAD8:
// 	  {
// 	    const Real xi  = p(0);
// 	    const Real eta = p(1);

// 	    libmesh_assert (i < 8);

// 	    //                                0  1  2  3  4  5  6  7  8
// 	    static const unsigned int i0[] = {0, 1, 1, 0, 2, 1, 2, 0, 2};
// 	    static const unsigned int i1[] = {0, 0, 1, 1, 0, 2, 1, 2, 2};
// 	    static const Real scal[] = {-0.25, -0.25, -0.25, -0.25, 0.5, 0.5, 0.5, 0.5};

// 	    return (FE<1,BERNSTEIN>::shape(EDGE3, totalorder, i0[i], xi)*
// 		    FE<1,BERNSTEIN>::shape(EDGE3, totalorder, i1[i], eta)
// 		    +scal[i]*
// 		    FE<1,BERNSTEIN>::shape(EDGE3, totalorder, i0[8], xi)*
// 		    FE<1,BERNSTEIN>::shape(EDGE3, totalorder, i1[8], eta));

// 	    //B_t,i0(i)|xi * B_s,i1(i)|eta + scal(i) * B_t,2|xi * B_t,2|eta
// 	  }

// 	  // Bernstein shape functions on the 9-noded quadrilateral.
// 	case QUAD9:
// 	  {
// 	    // Compute quad shape functions as a tensor-product
// 	    const Real xi  = p(0);
// 	    const Real eta = p(1);

// 	    libmesh_assert (i < 9);

// 	    //                                0  1  2  3  4  5  6  7  8
// 	    static const unsigned int i0[] = {0, 1, 1, 0, 2, 1, 2, 0, 2};
// 	    static const unsigned int i1[] = {0, 0, 1, 1, 0, 2, 1, 2, 2};

// 	    return (FE<1,BERNSTEIN>::shape(EDGE3, totalorder, i0[i], xi)*
// 		    FE<1,BERNSTEIN>::shape(EDGE3, totalorder, i1[i], eta));

// 	  }
// 	default:
// 	  libmesh_error();
// 	}

//     case THIRD:
//       switch(type)
// 	{
// 	case TRI6:
// 	  {
// 	    const Real x=p(0);
// 	    const Real y=p(1);
// 	    const Real r=1.-x-y;
// 	    libmesh_assert(i<10);

// 	    unsigned int shape=i;


// 	    if((i==3||i==4) && elem->node(0) > elem->node(1))shape=7-i;
// 	    if((i==5||i==6) && elem->node(1) > elem->node(2))shape=11-i;
// 	    if((i==7||i==8) && elem->node(0) > elem->node(2))shape=15-i;


// 	    switch(shape)
// 	      {
// 	      case 0: return r*r*r;
// 	      case 1: return x*x*x;
// 	      case 2: return y*y*y;

// 	      case 3: return 3.*x*r*r;
// 	      case 4: return 3.*x*x*r;

// 	      case 5: return 3.*y*x*x;
// 	      case 6: return 3.*y*y*x;

// 	      case 7: return 3.*y*r*r;
// 	      case 8: return 3.*y*y*r;

// 	      case 9: return 6.*x*y*r;
// 	      }
// 	  }

// 	  // Bernstein shape functions on the quadrilateral.
// 	case QUAD8:
// 	case QUAD9:
// 	  {
// 	    // Compute quad shape functions as a tensor-product
// 	    Real xi  = p(0);
// 	    Real eta = p(1);

// 	    libmesh_assert (i < 16);

// 	    //                                    0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15
// 	    static const unsigned int i0_reg[] = {0,  1,  1,  0,  2,  3,  1,  1,  2,  3,  0,  0,  2,  3,  2,  3};
// 	    static const unsigned int i1_reg[] = {0,  0,  1,  1,  0,  0,  2,  3,  1,  1,  2,  3,  2,  2,  3,  3};

// 	    unsigned int i0=i0_reg[i];
// 	    unsigned int i1=i1_reg[i];

// 	    if((i== 4||i== 5) && elem->node(0) > elem->node(1)) i0=5-i0;  // 2->3   3->2
// 	    if((i== 6||i== 7) && elem->node(1) > elem->node(2)) i1=5-i1;
// 	    if((i== 8||i== 9) && elem->node(3) > elem->node(2)) i0=5-i0;
// 	    if((i==10||i==11) && elem->node(0) > elem->node(3)) i1=5-i1;

// 	    // element dof orientation is needed when used with ifems
// // 	    if(i > 11)
// // 	      {
// // 		const unsigned int min_node = std::min(elem->node(1),
// // 						       std::min(elem->node(2),
// // 								std::min(elem->node(0),
// // 									 elem->node(3))));
// // 		if (elem->node(0) == min_node)
// // 		  if (elem->node(1) == std::min(elem->node(1), elem->node(3)))
// // 		    {
// // 		      // Case 1
// // 		      xi  = xi;
// // 		      eta = eta;
// // 		    }
// // 		  else
// // 		    {
// // 		      // Case 2
// // 		      xi  = eta;
// // 		      eta = xi;
// // 		    }

// // 		else if (elem->node(3) == min_node)
// // 		  if (elem->node(0) == std::min(elem->node(0), elem->node(2)))
// // 		    {
// // 		      // Case 3
// // 		      xi  = -eta;
// // 		      eta = xi;
// // 		    }
// // 		  else
// // 		    {
// // 		      // Case 4
// // 		      xi  = xi;
// // 		      eta = -eta;
// // 		    }

// // 		else if (elem->node(2) == min_node)
// // 		  if (elem->node(3) == std::min(elem->node(3), elem->node(1)))
// // 		    {
// // 		      // Case 5
// // 		      xi  = -xi;
// // 		      eta = -eta;
// // 		    }
// // 		  else
// // 		    {
// // 		      // Case 6
// // 		      xi  = -eta;
// // 		      eta = -xi;
// // 		    }

// // 		else if (elem->node(1) == min_node)
// // 		  if (elem->node(2) == std::min(elem->node(2), elem->node(0)))
// // 		    {
// // 		      // Case 7
// // 		      xi  = eta;
// // 		      eta = -xi;
// // 		    }
// // 		  else
// // 		    {
// // 		      // Case 8
// // 		      xi  = -xi;
// // 		      eta = eta;
// // 		    }
// // 	      }


// 	    return (FE<1,BERNSTEIN>::shape(EDGE3, totalorder, i0, xi)*
// 		    FE<1,BERNSTEIN>::shape(EDGE3, totalorder, i1, eta));

// 	  }

// 	default:
// 	  libmesh_error();

// 	}

//     case FOURTH:
//       switch(type)
// 	{
// 	case TRI6:
// 	  {
// 	    const Real x=p(0);
// 	    const Real y=p(1);
// 	    const Real r=1-x-y;
// 	    unsigned int shape=i;

// 	    libmesh_assert(i<15);

// 	    if((i==3||i== 5) && elem->node(0) > elem->node(1))shape=8-i;
// 	    if((i==6||i== 8) && elem->node(1) > elem->node(2))shape=14-i;
// 	    if((i==9||i==11) && elem->node(0) > elem->node(2))shape=20-i;


// 	    switch(shape)
// 	      {
// 		// point functions
// 	      case  0: return r*r*r*r;
// 	      case  1: return x*x*x*x;
// 	      case  2: return y*y*y*y;

// 		// edge functions
// 	      case  3: return 4.*x*r*r*r;
// 	      case  4: return 6.*x*x*r*r;
// 	      case  5: return 4.*x*x*x*r;

// 	      case  6: return 4.*y*x*x*x;
// 	      case  7: return 6.*y*y*x*x;
// 	      case  8: return 4.*y*y*y*x;

// 	      case  9: return 4.*y*r*r*r;
// 	      case 10: return 6.*y*y*r*r;
// 	      case 11: return 4.*y*y*y*r;

// 		// inner functions
// 	      case 12: return 12.*x*y*r*r;
// 	      case 13: return 12.*x*x*y*r;
// 	      case 14: return 12.*x*y*y*r;

// 	      }
// 	  }


// 	  // Bernstein shape functions on the quadrilateral.
// 	case QUAD8:
// 	case QUAD9:
// 	  {
// 	    // Compute quad shape functions as a tensor-product
// 	    const Real xi  = p(0);
// 	    const Real eta = p(1);

// 	    libmesh_assert (i < 25);

// 	    //                                    0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24
// 	    static const unsigned int i0_reg[] = {0, 1, 1, 0, 2, 3, 4, 1, 1, 1, 2, 3, 4, 0, 0, 0, 2, 3, 4, 2, 3, 4, 2, 3, 4};
// 	    static const unsigned int i1_reg[] = {0, 0, 1, 1, 0, 0, 0, 2, 3, 4, 1, 1, 1, 2, 3, 4, 2, 2, 2, 3, 3, 3, 4, 4, 4};

// 	    unsigned int i0=i0_reg[i];
// 	    unsigned int i1=i1_reg[i];

// 	    if((i>= 4&&i<= 6) && elem->node(0) > elem->node(1)) i0=6-i0;	// 2->4,  4->2
// 	    if((i>= 7&&i<= 9) && elem->node(1) > elem->node(2)) i1=6-i1;
// 	    if((i>=10&&i<=12) && elem->node(3) > elem->node(2)) i0=6-i0;
// 	    if((i>=13&&i<=15) && elem->node(0) > elem->node(3)) i1=6-i1;

// 	    return (FE<1,BERNSTEIN>::shape(EDGE3, totalorder, i0, xi)*
// 		    FE<1,BERNSTEIN>::shape(EDGE3, totalorder, i1, eta));

// 	  }

// 	default:
// 	  libmesh_error();

// 	}

//     case FIFTH:
//       switch(type)
// 	{
// 	case TRI6:
// 	  {
// 	    const Real x=p(0);
// 	    const Real y=p(1);
// 	    const Real r=1-x-y;
// 	    unsigned int shape=i;

// 	    libmesh_assert(i<21);

// 	    if((i>= 3&&i<= 6) && elem->node(0) > elem->node(1))shape=9-i;
// 	    if((i>= 7&&i<=10) && elem->node(1) > elem->node(2))shape=17-i;
// 	    if((i>=11&&i<=14) && elem->node(0) > elem->node(2))shape=25-i;

// 	    switch(shape)
// 	      {
// 		//point functions
// 	      case  0: return pow<5>(r);
// 	      case  1: return pow<5>(x);
// 	      case  2: return pow<5>(y);

// 		//edge functions
// 	      case  3: return  5.*x        *pow<4>(r);
// 	      case  4: return 10.*pow<2>(x)*pow<3>(r);
// 	      case  5: return 10.*pow<3>(x)*pow<2>(r);
// 	      case  6: return  5.*pow<4>(x)*r;

// 	      case  7: return  5.*y	   *pow<4>(x);
// 	      case  8: return 10.*pow<2>(y)*pow<3>(x);
// 	      case  9: return 10.*pow<3>(y)*pow<2>(x);
// 	      case 10: return  5.*pow<4>(y)*x;

// 	      case 11: return  5.*y	   *pow<4>(r);
// 	      case 12: return 10.*pow<2>(y)*pow<3>(r);
// 	      case 13: return 10.*pow<3>(y)*pow<2>(r);
// 	      case 14: return  5.*pow<4>(y)*r;

// 		//inner functions
// 	      case 15: return 20.*x*y*pow<3>(r);
// 	      case 16: return 30.*x*pow<2>(y)*pow<2>(r);
// 	      case 17: return 30.*pow<2>(x)*y*pow<2>(r);
// 	      case 18: return 20.*x*pow<3>(y)*r;
// 	      case 19: return 20.*pow<3>(x)*y*r;
// 	      case 20: return 30.*pow<2>(x)*pow<2>(y)*r;

// 	      }
// 	  }


// 	  // Bernstein shape functions on the quadrilateral.
// 	case QUAD8:
// 	case QUAD9:
// 	  {
// 	    // Compute quad shape functions as a tensor-product
// 	    const Real xi  = p(0);
// 	    const Real eta = p(1);

// 	    libmesh_assert (i < 36);

// 	    //                                    0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35
// 	    static const unsigned int i0_reg[] = {0, 1, 1, 0, 2, 3, 4, 5, 1, 1, 1, 1, 2, 3, 4, 5, 0, 0, 0, 0, 2, 3, 4, 5, 2, 3, 4, 5, 2, 3, 4, 5, 2, 3, 4, 5};
// 	    static const unsigned int i1_reg[] = {0, 0, 1, 1, 0, 0, 0, 0, 2, 3, 4, 5, 1, 1, 1, 1, 2, 3, 4, 5, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5};

// 	    unsigned int i0=i0_reg[i];
// 	    unsigned int i1=i1_reg[i];

// 	    if((i>= 4&&i<= 7) && elem->node(0) > elem->node(1)) i0=7-i0;
// 	    if((i>= 8&&i<=11) && elem->node(1) > elem->node(2)) i1=7-i1;
// 	    if((i>=12&&i<=15) && elem->node(3) > elem->node(2)) i0=7-i0;
// 	    if((i>=16&&i<=19) && elem->node(0) > elem->node(3)) i1=7-i1;

// 	    return (FE<1,BERNSTEIN>::shape(EDGE3, totalorder, i0, xi)*
// 		    FE<1,BERNSTEIN>::shape(EDGE3, totalorder, i1, eta));

// 	  }
// 	default:
// 	  libmesh_error();
// 	}

//     case SIXTH:
//       switch(type)
// 	{
// 	case TRI6:
// 	  {
// 	    const Real x=p(0);
// 	    const Real y=p(1);
// 	    const Real r=1-x-y;
// 	    unsigned int shape=i;

// 	    libmesh_assert(i<28);

// 	    if((i>= 3&&i<= 7) && elem->node(0) > elem->node(1))shape=10-i;
// 	    if((i>= 8&&i<=12) && elem->node(1) > elem->node(2))shape=20-i;
// 	    if((i>=13&&i<=17) && elem->node(0) > elem->node(2))shape=30-i;

// 	    switch(shape)
// 	      {
// 		//point functions
// 	      case  0: return pow<6>(r);
// 	      case  1: return pow<6>(x);
// 	      case  2: return pow<6>(y);

// 		//edge functions
// 	      case  3: return  6.*x        *pow<5>(r);
// 	      case  4: return 15.*pow<2>(x)*pow<4>(r);
// 	      case  5: return 20.*pow<3>(x)*pow<3>(r);
// 	      case  6: return 15.*pow<4>(x)*pow<2>(r);
// 	      case  7: return  6.*pow<5>(x)*r;

// 	      case  8: return  6.*y        *pow<5>(x);
// 	      case  9: return 15.*pow<2>(y)*pow<4>(x);
// 	      case 10: return 20.*pow<3>(y)*pow<3>(x);
// 	      case 11: return 15.*pow<4>(y)*pow<2>(x);
// 	      case 12: return  6.*pow<5>(y)*x;

// 	      case 13: return  6.*y        *pow<5>(r);
// 	      case 14: return 15.*pow<2>(y)*pow<4>(r);
// 	      case 15: return 20.*pow<3>(y)*pow<3>(r);
// 	      case 16: return 15.*pow<4>(y)*pow<2>(r);
// 	      case 17: return  6.*pow<5>(y)*r;

// 		//inner functions
// 	      case 18: return 30.*x*y*pow<4>(r);
// 	      case 19: return 60.*x*pow<2>(y)*pow<3>(r);
// 	      case 20: return 60.*  pow<2>(x)*y*pow<3>(r);
// 	      case 21: return 60.*x*pow<3>(y)*pow<2>(r);
// 	      case 22: return 60.*pow<3>(x)*y*pow<2>(r);
// 	      case 23: return 90.*pow<2>(x)*pow<2>(y)*pow<2>(r);
// 	      case 24: return 30.*x*pow<4>(y)*r;
// 	      case 25: return 60.*pow<2>(x)*pow<3>(y)*r;
// 	      case 26: return 60.*pow<3>(x)*pow<2>(y)*r;
// 	      case 27: return 30.*pow<4>(x)*y*r;

// 	      } // switch shape
// 	  } // case TRI6



// 	  // Bernstein shape functions on the quadrilateral.
// 	case QUAD8:
// 	case QUAD9:
// 	  {
// 	    // Compute quad shape functions as a tensor-product
// 	    const Real xi  = p(0);
// 	    const Real eta = p(1);

// 	    libmesh_assert (i < 49);

// 	    //                                    0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48
// 	    static const unsigned int i0_reg[] = {0, 1, 1, 0, 2, 3, 4, 5, 6, 1, 1, 1, 1, 1, 2, 3, 4, 5, 6, 0, 0, 0, 0, 0, 2, 3, 4, 5, 6, 2, 3, 4, 5, 6, 2, 3, 4, 5, 6, 2, 3, 4, 5, 6, 2, 3, 4, 5, 6};
// 	    static const unsigned int i1_reg[] = {0, 0, 1, 1, 0, 0, 0, 0, 0, 2, 3, 4, 5, 6, 1, 1, 1, 1, 1, 2, 3, 4, 5, 6, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6};

// 	    unsigned int i0=i0_reg[i];
// 	    unsigned int i1=i1_reg[i];

// 	    if((i>= 4&&i<= 8) && elem->node(0) > elem->node(1)) i0=8-i0;
// 	    if((i>= 9&&i<=13) && elem->node(1) > elem->node(2)) i1=8-i1;
// 	    if((i>=14&&i<=18) && elem->node(3) > elem->node(2)) i0=8-i0;
// 	    if((i>=19&&i<=23) && elem->node(0) > elem->node(3)) i1=8-i1;

// 	    return (FE<1,BERNSTEIN>::shape(EDGE3, totalorder, i0, xi)*
// 		    FE<1,BERNSTEIN>::shape(EDGE3, totalorder, i1, eta));

// 	  } // case QUAD8/9

// 	default:
// 	  libmesh_error();

// 	}
//       // 7th-order Bernstein.
//     case SEVENTH:
//       {
// 	switch (type)
// 	  {

// 	    // Szabo-Babuska shape functions on the quadrilateral.
// 	  case QUAD9:
// 	    {
// 	      // Compute quad shape functions as a tensor-product
// 	      const Real xi  = p(0);
// 	      const Real eta = p(1);

// 	      libmesh_assert (i < 64);

// 	      //                                0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63
// 	      static const unsigned int i0[] = {0, 1, 1, 0, 2, 3, 4, 5, 6, 7, 1, 1, 1, 1, 1, 1, 2, 3, 4, 5, 6, 7, 0, 0, 0, 0, 0, 0, 2, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7};
// 	      static const unsigned int i1[] = {0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 2, 3, 4, 5, 6, 7, 1, 1, 1, 1, 1, 1, 2, 3, 4, 5, 6, 7, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7};

// 	      Real f=1.;

// 	      switch(i)
// 		{
// 		case  5: // edge 0 nodes
// 		case  7:
// 		case  9:
// 		  if (elem->node(0) > elem->node(1))f = -1.;
// 		  break;
// 	      	case 11: // edge 1 nodes
// 	      	case 13:
// 		case 15:
// 		  if (elem->node(1) > elem->node(2))f = -1.;
// 		  break;
// 	        case 17: // edge 2 nodes
// 	        case 19:
// 		case 21:
// 		  if (elem->node(3) > elem->node(2))f = -1.;
// 		  break;
// 	        case 23: // edge 3 nodes
// 	        case 25:
// 		case 27:
// 		  if (elem->node(0) > elem->node(3))f = -1.;
// 		  break;
// 		}

// 	      return f*(FE<1,BERNSTEIN>::shape(EDGE3, totalorder, i0[i], xi)*
// 			FE<1,BERNSTEIN>::shape(EDGE3, totalorder, i1[i], eta));

// 	    } // case QUAD8/QUAD9

// 	  default:
// 	    libmesh_error();

// 	  } // switch type

//       } // case SEVENTH


//       // by default throw an error
//     default:
//       libMesh::err << "ERROR: Unsupported polynomial order!" << std::endl;
//       libmesh_error();
//     }

  libmesh_error();
  return 0.;
}



template <>
Real FE<2,BERNSTEIN>::shape_deriv(const ElemType,
				  const Order,
				  const unsigned int,
				  const unsigned int,
				  const Point&)
{
  libMesh::err << "Bernstein polynomials require the element type\n"
	        << "because edge orientation is needed."
	        << std::endl;

  libmesh_error();
  return 0.;
}



template <>
Real FE<2,BERNSTEIN>::shape_deriv(const Elem* elem,
				  const Order order,
				  const unsigned int i,
				  const unsigned int j,
				  const Point& p)
{
  libmesh_assert (elem != NULL);

  const ElemType type = elem->type();

  const Order totalorder = static_cast<Order>(order + elem->p_level());

  switch (type)
    {
      // Hierarchic shape functions on the quadrilateral.
    case QUAD4:
    case QUAD9:
      {
        // Compute quad shape functions as a tensor-product
        const Real xi  = p(0);
        const Real eta = p(1);

        libmesh_assert (i < (totalorder+1u)*(totalorder+1u));

        unsigned int i0, i1;

        // Vertex DoFs
        if (i == 0)
          { i0 = 0; i1 = 0; }
        else if (i == 1)
          { i0 = 1; i1 = 0; }
        else if (i == 2)
          { i0 = 1; i1 = 1; }
        else if (i == 3)
          { i0 = 0; i1 = 1; }


        // Edge DoFs
        else if (i < totalorder + 3u)
          { i0 = i - 2; i1 = 0; }
        else if (i < 2u*totalorder + 2)
          { i0 = 1; i1 = i - totalorder - 1; }
        else if (i < 3u*totalorder + 1)
          { i0 = i - 2u*totalorder; i1 = 1; }
        else if (i < 4u*totalorder)
          { i0 = 0; i1 = i - 3u*totalorder + 1; }
        // Interior DoFs
        else
          {
	    unsigned int basisnum = i - 4*totalorder;
	    i0 = square_number_column[basisnum] + 2;
	    i1 = square_number_row[basisnum] + 2;
          }


	// Flip odd degree of freedom values if necessary
	// to keep continuity on sides
	if     ((i>= 4                 && i<= 4+  totalorder-2u) && elem->point(0) > elem->point(1)) i0=totalorder+2-i0;	//
	else if((i>= 4+  totalorder-1u && i<= 4+2*totalorder-3u) && elem->point(1) > elem->point(2)) i1=totalorder+2-i1;
	else if((i>= 4+2*totalorder-2u && i<= 4+3*totalorder-4u) && elem->point(3) > elem->point(2)) i0=totalorder+2-i0;
	else if((i>= 4+3*totalorder-3u && i<= 4+4*totalorder-5u) && elem->point(0) > elem->point(3)) i1=totalorder+2-i1;

	switch (j)
	  {
	    // d()/dxi
	  case 0:
	    return (FE<1,BERNSTEIN>::shape_deriv(EDGE3, totalorder, i0, 0, xi)*
		    FE<1,BERNSTEIN>::shape      (EDGE3, totalorder, i1,    eta));

	    // d()/deta
	  case 1:
	    return (FE<1,BERNSTEIN>::shape      (EDGE3, totalorder, i0,    xi)*
		    FE<1,BERNSTEIN>::shape_deriv(EDGE3, totalorder, i1, 0, eta));

	  }
      }

      // Bernstein shape functions on the 8-noded quadrilateral
      // is handled separately.
    case QUAD8:
      {
	libmesh_assert (totalorder < 3);

	const Real xi  = p(0);
	const Real eta = p(1);

	libmesh_assert (i < 8);

	//                                0  1  2  3  4  5  6  7  8
	static const unsigned int i0[] = {0, 1, 1, 0, 2, 1, 2, 0, 2};
	static const unsigned int i1[] = {0, 0, 1, 1, 0, 2, 1, 2, 2};
	static const Real scal[] = {-0.25, -0.25, -0.25, -0.25, 0.5, 0.5, 0.5, 0.5};
	switch (j)
	  {
	    // d()/dxi
	  case 0:
	    return (FE<1,BERNSTEIN>::shape_deriv(EDGE3, totalorder, i0[i], 0, xi)*
		    FE<1,BERNSTEIN>::shape      (EDGE3, totalorder, i1[i],    eta)
		    +scal[i]*
		    FE<1,BERNSTEIN>::shape_deriv(EDGE3, totalorder, i0[8], 0, xi)*
		    FE<1,BERNSTEIN>::shape      (EDGE3, totalorder, i1[8],    eta));

	    // d()/deta
	  case 1:
	    return (FE<1,BERNSTEIN>::shape      (EDGE3, totalorder, i0[i],    xi)*
		    FE<1,BERNSTEIN>::shape_deriv(EDGE3, totalorder, i1[i], 0, eta)
		    +scal[i]*
		    FE<1,BERNSTEIN>::shape      (EDGE3, totalorder, i0[8],    xi)*
		    FE<1,BERNSTEIN>::shape_deriv(EDGE3, totalorder, i1[8], 0, eta));

	  default:
	    libmesh_error();
	  }
      }

    case TRI3:
      libmesh_assert (totalorder<2);
    case TRI6:
      {
	// I have been lazy here and am using finite differences
	// to compute the derivatives!
	const Real eps = 1.e-6;

	switch (j)
	  {
	    //  d()/dxi
	  case 0:
	    {
	      const Point pp(p(0)+eps, p(1));
	      const Point pm(p(0)-eps, p(1));

	      return (FE<2,BERNSTEIN>::shape(elem, totalorder, i, pp) -
		      FE<2,BERNSTEIN>::shape(elem, totalorder, i, pm))/2./eps;
	    }

	    // d()/deta
	  case 1:
	    {
	      const Point pp(p(0), p(1)+eps);
	      const Point pm(p(0), p(1)-eps);

	      return (FE<2,BERNSTEIN>::shape(elem, totalorder, i, pp) -
		      FE<2,BERNSTEIN>::shape(elem, totalorder, i, pm))/2./eps;
	    }


	  default:
	    libmesh_error();
	  }
      }

    default:
      {
	libMesh::err << "ERROR: Unsupported element type!" << std::endl;
	libmesh_error();
      }

    }

  // old code
//   switch (totalorder)
//     {

//       // 1st order Bernsteins are aquivalent to Lagrange elements.
//     case FIRST:
//       {
// 	switch (type)
// 	  {
// 	    // Bernstein shape functions on the triangle.
// 	  case TRI6:
// 	    {
// 	      // I have been lazy here and am using finite differences
// 	      // to compute the derivatives!
// 	      const Real eps = 1.e-6;

// 	      libmesh_assert (i < 3);
// 	      libmesh_assert (j < 2);

// 	      switch (j)
// 		{
// 		  //  d()/dxi
// 		case 0:
// 		  {
// 		    const Point pp(p(0)+eps, p(1));
// 		    const Point pm(p(0)-eps, p(1));

// 		    return (FE<2,BERNSTEIN>::shape(elem, order, i, pp) -
// 			    FE<2,BERNSTEIN>::shape(elem, order, i, pm))/2./eps;
// 		  }

// 		  // d()/deta
// 		case 1:
// 		  {
// 		    const Point pp(p(0), p(1)+eps);
// 		    const Point pm(p(0), p(1)-eps);

// 		    return (FE<2,BERNSTEIN>::shape(elem, order, i, pp) -
// 			    FE<2,BERNSTEIN>::shape(elem, order, i, pm))/2./eps;
// 		  }


// 		default:
// 		  libmesh_error();
// 		}
// 	    }

// 	    // Bernstein shape functions on the quadrilateral.
// 	  case QUAD8:
// 	  case QUAD9:
// 	    {
// 	      // Compute quad shape functions as a tensor-product
// 	      const Real xi  = p(0);
// 	      const Real eta = p(1);

// 	      libmesh_assert (i < 4);

// 	      //                                0  1  2  3
// 	      static const unsigned int i0[] = {0, 1, 1, 0};
// 	      static const unsigned int i1[] = {0, 0, 1, 1};

// 	      switch (j)
// 		{
// 		  // d()/dxi
// 		case 0:
// 		  return (FE<1,BERNSTEIN>::shape_deriv(EDGE3, totalorder, i0[i], 0, xi)*
// 			  FE<1,BERNSTEIN>::shape      (EDGE3, totalorder, i1[i],    eta));

// 		  // d()/deta
// 		case 1:
// 		  return (FE<1,BERNSTEIN>::shape      (EDGE3, totalorder, i0[i],    xi)*
// 			  FE<1,BERNSTEIN>::shape_deriv(EDGE3, totalorder, i1[i], 0, eta));

// 		default:
// 		  libmesh_error();
// 		}
// 	    }

// 	  default:
// 	    libmesh_error();
// 	  }
//       }

//       // second order Bernsteins.
//     case SECOND:
//       {
// 	switch (type)
// 	  {
// 	    // Bernstein shape functions on the triangle.
// 	  case TRI6:
// 	    {
// 	      // I have been lazy here and am using finite differences
// 	      // to compute the derivatives!
// 	      const Real eps = 1.e-6;

// 	      libmesh_assert (i < 6);
// 	      libmesh_assert (j < 2);

// 	      switch (j)
// 		{
// 		  //  d()/dxi
// 		case 0:
// 		  {
// 		    const Point pp(p(0)+eps, p(1));
// 		    const Point pm(p(0)-eps, p(1));

// 		    return (FE<2,BERNSTEIN>::shape(elem, order, i, pp) -
// 			    FE<2,BERNSTEIN>::shape(elem, order, i, pm))/2./eps;
// 		  }

// 		  // d()/deta
// 		case 1:
// 		  {
// 		    const Point pp(p(0), p(1)+eps);
// 		    const Point pm(p(0), p(1)-eps);

// 		    return (FE<2,BERNSTEIN>::shape(elem, order, i, pp) -
// 			    FE<2,BERNSTEIN>::shape(elem, order, i, pm))/2./eps;
// 		  }


// 		default:
// 		  libmesh_error();
// 		}
// 	    }

// 	    // Bernstein shape functions on the 8-noded quadrilateral.
// 	  case QUAD8:
// 	    {
// 	      const Real xi  = p(0);
// 	      const Real eta = p(1);

// 	      libmesh_assert (i < 8);

// 	      //                                0  1  2  3  4  5  6  7  8
// 	      static const unsigned int i0[] = {0, 1, 1, 0, 2, 1, 2, 0, 2};
// 	      static const unsigned int i1[] = {0, 0, 1, 1, 0, 2, 1, 2, 2};
// 	      static const Real scal[] = {-0.25, -0.25, -0.25, -0.25, 0.5, 0.5, 0.5, 0.5};
// 	      switch (j)
// 		{
// 		  // d()/dxi
// 		case 0:
// 		  return (FE<1,BERNSTEIN>::shape_deriv(EDGE3, totalorder, i0[i], 0, xi)*
// 			  FE<1,BERNSTEIN>::shape      (EDGE3, totalorder, i1[i],    eta)
// 			  +scal[i]*
// 			  FE<1,BERNSTEIN>::shape_deriv(EDGE3, totalorder, i0[8], 0, xi)*
// 			  FE<1,BERNSTEIN>::shape      (EDGE3, totalorder, i1[8],    eta));

// 		  // d()/deta
// 		case 1:
// 		  return (FE<1,BERNSTEIN>::shape      (EDGE3, totalorder, i0[i],    xi)*
// 			  FE<1,BERNSTEIN>::shape_deriv(EDGE3, totalorder, i1[i], 0, eta)
// 			  +scal[i]*
// 			  FE<1,BERNSTEIN>::shape      (EDGE3, totalorder, i0[8],    xi)*
// 			  FE<1,BERNSTEIN>::shape_deriv(EDGE3, totalorder, i1[8], 0, eta));

// 		default:
// 		  libmesh_error();
// 		}
// 	    }


// 	    // Bernstein shape functions on the 9-noded quadrilateral.
// 	  case QUAD9:
// 	    {
// 	      // Compute quad shape functions as a tensor-product
// 	      const Real xi  = p(0);
// 	      const Real eta = p(1);

// 	      libmesh_assert (i < 9);

// 	      //                                0  1  2  3  4  5  6  7  8
// 	      static const unsigned int i0[] = {0, 1, 1, 0, 2, 1, 2, 0, 2};
// 	      static const unsigned int i1[] = {0, 0, 1, 1, 0, 2, 1, 2, 2};

// 	      switch (j)
// 		{
// 		  // d()/dxi
// 		case 0:
// 		  return (FE<1,BERNSTEIN>::shape_deriv(EDGE3, totalorder, i0[i], 0, xi)*
// 			  FE<1,BERNSTEIN>::shape      (EDGE3, totalorder, i1[i],    eta));

// 		  // d()/deta
// 		case 1:
// 		  return (FE<1,BERNSTEIN>::shape      (EDGE3, totalorder, i0[i],    xi)*
// 			  FE<1,BERNSTEIN>::shape_deriv(EDGE3, totalorder, i1[i], 0, eta));

// 		default:
// 		  libmesh_error();
// 		}
// 	    }

// 	  default:
// 	    libmesh_error();
// 	  }
//       }


//       // 3rd-order Bernsteins.
//     case THIRD:
//       {
// 	switch (type)
// 	  {
// 	    // Bernstein shape functions on the triangle.
// 	  case TRI6:
// 	    {
// 	      // I have been lazy here and am using finite differences
// 	      // to compute the derivatives!
// 	      const Real eps = 1.e-6;

// 	      libmesh_assert (i < 10);
// 	      libmesh_assert (j < 2);


// 	      unsigned int shape=i;

// 	      /**
// 	       * Note, since the derivatives are computed using FE::shape
// 	       * the continuity along neighboring elements (shape_flip)
// 	       * is already considered.
// 	       */


// 	      switch (j)
// 		{
// 		  //  d()/dxi
// 		case 0:
// 		  {
// 		    const Point pp(p(0)+eps, p(1));
// 		    const Point pm(p(0)-eps, p(1));

// 		    return (FE<2,BERNSTEIN>::shape(elem, totalorder, shape, pp) -
// 			    FE<2,BERNSTEIN>::shape(elem, totalorder, shape, pm))/2./eps;
// 		  }

// 		  // d()/deta
// 		case 1:
// 		  {
// 		    const Point pp(p(0), p(1)+eps);
// 		    const Point pm(p(0), p(1)-eps);

// 		    return (FE<2,BERNSTEIN>::shape(elem, totalorder, shape, pp) -
// 			    FE<2,BERNSTEIN>::shape(elem, totalorder, shape, pm))/2./eps;
// 		  }


// 		default:
// 		  libmesh_error();
// 		}
// 	    }



// 	    // Bernstein shape functions on the quadrilateral.
// 	  case QUAD8:
// 	  case QUAD9:
// 	    {
// 	      // Compute quad shape functions as a tensor-product
// 	      const Real xi  = p(0);
// 	      const Real eta = p(1);

// 	      libmesh_assert (i < 16);

// 		  //                                    0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15
// 	      static const unsigned int i0_reg[] = {0,  1,  1,  0,  2,  3,  1,  1,  2,  3,  0,  0,  2,  3,  2,  3};
// 	      static const unsigned int i1_reg[] = {0,  0,  1,  1,  0,  0,  2,  3,  1,  1,  2,  3,  2,  2,  3,  3};

// 	      unsigned int i0=i0_reg[i];
// 	      unsigned int i1=i1_reg[i];

// 	      if((i== 4||i== 5) && elem->node(0) > elem->node(1)) i0=5-i0;
// 	      if((i== 6||i== 7) && elem->node(1) > elem->node(2)) i1=5-i1;
// 	      if((i== 8||i== 9) && elem->node(3) > elem->node(2)) i0=5-i0;
// 	      if((i==10||i==11) && elem->node(0) > elem->node(3)) i1=5-i1;

// 	      switch (j)
// 		{
// 			  // d()/dxi
// 		case 0:
// 		  return (FE<1,BERNSTEIN>::shape_deriv(EDGE3, totalorder, i0, 0, xi)*
// 			  FE<1,BERNSTEIN>::shape      (EDGE3, totalorder, i1,    eta));

// 		  // d()/deta
// 		case 1:
// 		  return (FE<1,BERNSTEIN>::shape      (EDGE3, totalorder, i0,    xi)*
// 			  FE<1,BERNSTEIN>::shape_deriv(EDGE3, totalorder, i1, 0, eta));

// 		default:
// 		  libmesh_error();
// 		}
// 	    }

// 	  default:
// 	    libmesh_error();
// 	  }
//       }




//       // 4th-order Bernsteins.
//     case FOURTH:
//       {
// 	switch (type)
// 	  {
// 	    // Bernstein shape functions on the triangle.
// 	  case TRI6:
// 	    {
// 	      // I have been lazy here and am using finite differences
// 	      // to compute the derivatives!
// 	      const Real eps = 1.e-6;

// 	      libmesh_assert (i < 15);
// 	      libmesh_assert (j < 2);

// 	      unsigned int shape=i;

// 	      switch (j)
// 		{
// 		  //  d()/dxi
// 		case 0:
// 		  {
// 		    const Point pp(p(0)+eps, p(1));
// 		    const Point pm(p(0)-eps, p(1));

// 		    return (FE<2,BERNSTEIN>::shape(elem, order, shape, pp) -
// 			    FE<2,BERNSTEIN>::shape(elem, order, shape, pm))/2./eps;
// 		  }

// 		  // d()/deta
// 		case 1:
// 		  {
// 		    const Point pp(p(0), p(1)+eps);
// 		    const Point pm(p(0), p(1)-eps);

// 		    return (FE<2,BERNSTEIN>::shape(elem, order, shape, pp) -
// 			    FE<2,BERNSTEIN>::shape(elem, order, shape, pm))/2./eps;
// 		  }


// 		default:
// 		  libmesh_error();
// 		}
// 	    }



// 	    // Bernstein shape functions on the quadrilateral.
// 	  case QUAD8:
// 	  case QUAD9:
// 	    {
// 	      // Compute quad shape functions as a tensor-product
// 	      const Real xi  = p(0);
// 	      const Real eta = p(1);

// 	      libmesh_assert (i < 25);

// 	      //                                      0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24
// 	      	static const unsigned int i0_reg[] = {0, 1, 1, 0, 2, 3, 4, 1, 1, 1, 2, 3, 4, 0, 0, 0, 2, 3, 4, 2, 3, 4, 2, 3, 4};
// 	      	static const unsigned int i1_reg[] = {0, 0, 1, 1, 0, 0, 0, 2, 3, 4, 1, 1, 1, 2, 3, 4, 2, 2, 2, 3, 3, 3, 4, 4, 4};

// 	      	unsigned int i0=i0_reg[i];
// 	      	unsigned int i1=i1_reg[i];

// 		if((i== 4||i== 6) && elem->node(0) > elem->node(1)) i0=6-i0;
// 		if((i== 7||i== 9) && elem->node(1) > elem->node(2)) i1=6-i1;
// 		if((i==10||i==12) && elem->node(3) > elem->node(2)) i0=6-i0;
// 	        if((i==13||i==15) && elem->node(0) > elem->node(3)) i1=6-i1;


// 		switch (j)
// 		  {
// 		    // d()/dxi
// 		  case 0:
// 		    return (FE<1,BERNSTEIN>::shape_deriv(EDGE3, totalorder, i0, 0, xi)*
// 			    FE<1,BERNSTEIN>::shape      (EDGE3, totalorder, i1,    eta));

// 		  // d()/deta
// 		  case 1:
// 		    return (FE<1,BERNSTEIN>::shape      (EDGE3, totalorder, i0,    xi)*
// 			    FE<1,BERNSTEIN>::shape_deriv(EDGE3, totalorder, i1, 0, eta));

// 		default:
// 		  libmesh_error();
// 		}
// 	    }

// 	  default:
// 	    libmesh_error();
// 	  }
//       }




//       // 5th-order Bernsteins.
//     case FIFTH:
//       {
// 	// Bernstein shape functions on the quadrilateral.
// 	switch (type)
// 	  {
// 	    // Bernstein shape functions on the triangle.
// 	  case TRI6:
// 	    {
// 	      // I have been lazy here and am using finite differences
// 	      // to compute the derivatives!
// 	      const Real eps = 1.e-6;

// 	      libmesh_assert (i < 21);
// 	      libmesh_assert (j < 2);

// 	      unsigned int shape=i;

// 	      switch (j)
// 		{
// 		  //  d()/dxi
// 		case 0:
// 		  {
// 		    const Point pp(p(0)+eps, p(1));
// 		    const Point pm(p(0)-eps, p(1));

// 		    return (FE<2,BERNSTEIN>::shape(elem, order, shape, pp) -
// 			    FE<2,BERNSTEIN>::shape(elem, order, shape, pm))/2./eps;
// 		  }

// 		  // d()/deta
// 		case 1:
// 		  {
// 		    const Point pp(p(0), p(1)+eps);
// 		    const Point pm(p(0), p(1)-eps);

// 		    return (FE<2,BERNSTEIN>::shape(elem, order, shape, pp) -
// 			    FE<2,BERNSTEIN>::shape(elem, order, shape, pm))/2./eps;
// 		  }


// 		default:
// 		  libmesh_error();
// 		}
// 	    } // case TRI6



// 	  case QUAD8:
// 	  case QUAD9:
// 	    {
// 	      // Compute quad shape functions as a tensor-product
// 	      const Real xi  = p(0);
// 	      const Real eta = p(1);

// 	      libmesh_assert (i < 36);

// 	      //                                          0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35
// 		    static const unsigned int i0_reg[] = {0, 1, 1, 0, 2, 3, 4, 5, 1, 1, 1, 1, 2, 3, 4, 5, 0, 0, 0, 0, 2, 3, 4, 5, 2, 3, 4, 5, 2, 3, 4, 5, 2, 3, 4, 5};
// 		    static const unsigned int i1_reg[] = {0, 0, 1, 1, 0, 0, 0, 0, 2, 3, 4, 5, 1, 1, 1, 1, 2, 3, 4, 5, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5};

// 		    unsigned int i0=i0_reg[i];
// 		    unsigned int i1=i1_reg[i];

// 		    if((i>= 4&&i<= 7) && elem->node(0) > elem->node(1)) i0=7-i0;
// 		    if((i>= 8&&i<=11) && elem->node(1) > elem->node(2)) i1=7-i1;
// 		    if((i>=12&&i<=15) && elem->node(3) > elem->node(2)) i0=7-i0;
// 		    if((i>=16&&i<=19) && elem->node(0) > elem->node(3)) i1=7-i1;


// 		    switch (j)
// 		      {
// 			// d()/dxi
// 		      case 0:
// 			return (FE<1,BERNSTEIN>::shape_deriv(EDGE3, totalorder, i0, 0, xi)*
// 				FE<1,BERNSTEIN>::shape      (EDGE3, totalorder, i1,    eta));

// 			// d()/deta
// 		      case 1:
// 			return (FE<1,BERNSTEIN>::shape      (EDGE3, totalorder, i0,    xi)*
// 				FE<1,BERNSTEIN>::shape_deriv(EDGE3, totalorder, i1, 0, eta));

// 		      default:
// 			libmesh_error();
// 		      }
// 	    } // case QUAD8/9

// 	  default:
// 	    libmesh_error();
// 	  }

//       }


//       // 6th-order Bernsteins.
//     case SIXTH:
//       {
// 	// Bernstein shape functions on the quadrilateral.
// 	switch (type)
// 	  {
// 	    // Bernstein shape functions on the triangle.
// 	  case TRI6:
// 	    {
// 	      // I have been lazy here and am using finite differences
// 	      // to compute the derivatives!
// 	      const Real eps = 1.e-6;

// 	      libmesh_assert (i < 28);
// 	      libmesh_assert (j < 2);

// 	      unsigned int shape=i;

// 	      switch (j)
// 		{
// 		  //  d()/dxi
// 		case 0:
// 		  {
// 		    const Point pp(p(0)+eps, p(1));
// 		    const Point pm(p(0)-eps, p(1));

// 		    return (FE<2,BERNSTEIN>::shape(elem, order, shape, pp) -
// 			    FE<2,BERNSTEIN>::shape(elem, order, shape, pm))/2./eps;
// 		  }

// 		  // d()/deta
// 		case 1:
// 		  {
// 		    const Point pp(p(0), p(1)+eps);
// 		    const Point pm(p(0), p(1)-eps);

// 		    return (FE<2,BERNSTEIN>::shape(elem, order, shape, pp) -
// 			    FE<2,BERNSTEIN>::shape(elem, order, shape, pm))/2./eps;
// 		  }


// 		default:
// 		  libmesh_error();
// 		}
// 	    } // case TRI6



// 	  case QUAD8:
// 	  case QUAD9:
// 	    {
// 	      // Compute quad shape functions as a tensor-product
// 	      const Real xi  = p(0);
// 	      const Real eta = p(1);

// 	      libmesh_assert (i < 49);

// 	      //                                    0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48
// 	      static const unsigned int i0_reg[] = {0, 1, 1, 0, 2, 3, 4, 5, 6, 1, 1, 1, 1, 1, 2, 3, 4, 5, 6, 0, 0, 0, 0, 0, 2, 3, 4, 5, 6, 2, 3, 4, 5, 6, 2, 3, 4, 5, 6, 2, 3, 4, 5, 6, 2, 3, 4, 5, 6};
// 	      static const unsigned int i1_reg[] = {0, 0, 1, 1, 0, 0, 0, 0, 0, 2, 3, 4, 5, 6, 1, 1, 1, 1, 1, 2, 3, 4, 5, 6, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6};

// 	      unsigned int i0=i0_reg[i];
// 	      unsigned int i1=i1_reg[i];

// 	      if((i>= 4&&i<= 8) && elem->node(0) > elem->node(1)) i0=8-i0;
// 	      if((i>= 9&&i<=13) && elem->node(1) > elem->node(2)) i1=8-i1;
// 	      if((i>=14&&i<=18) && elem->node(3) > elem->node(2)) i0=8-i0;
// 	      if((i>=19&&i<=23) && elem->node(0) > elem->node(3)) i1=8-i1;


// 	      switch (j)
// 		{
// 		  // d()/dxi
// 		case 0:
// 		  return (FE<1,BERNSTEIN>::shape_deriv(EDGE3, totalorder, i0, 0, xi)*
// 			  FE<1,BERNSTEIN>::shape      (EDGE3, totalorder, i1,    eta));

// 		  // d()/deta
// 		case 1:
// 		  return (FE<1,BERNSTEIN>::shape      (EDGE3, totalorder, i0,    xi)*
// 			  FE<1,BERNSTEIN>::shape_deriv(EDGE3, totalorder, i1, 0, eta));

// 		default:
// 		  libmesh_error();
// 		}
// 	    } // case QUAD8/9

// 	  default:
// 	    libmesh_error();
// 	  }

// 	// 7th-order Bernstein.
//       case SEVENTH:
// 	{
// 	  // Szabo-Babuska shape functions on the quadrilateral.
// 	  switch (type)
// 	    {


// 	    case QUAD9:
// 	      {
// 		// Compute quad shape functions as a tensor-product
// 		const Real xi  = p(0);
// 		const Real eta = p(1);

// 		libmesh_assert (i < 64);

// 		//                                0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63
// 		static const unsigned int i0[] = {0, 1, 1, 0, 2, 3, 4, 5, 6, 7, 1, 1, 1, 1, 1, 1, 2, 3, 4, 5, 6, 7, 0, 0, 0, 0, 0, 0, 2, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7};
// 		static const unsigned int i1[] = {0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 2, 3, 4, 5, 6, 7, 1, 1, 1, 1, 1, 1, 2, 3, 4, 5, 6, 7, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7};

// 		Real f=1.;

// 		switch(i)
// 		  {
// 		  case  5: // edge 0 nodes
// 		  case  7:
// 		  case  9:
// 		    if (elem->node(0) > elem->node(1))f = -1.;
// 		    break;
// 		  case 11: // edge 1 nodes
// 		  case 13:
// 		  case 15:
// 		    if (elem->node(1) > elem->node(2))f = -1.;
// 		    break;
// 		  case 17: // edge 2 nodes
// 		  case 19:
// 		  case 21:
// 		    if (elem->node(3) > elem->node(2))f = -1.;
// 		    break;
// 		  case 23: // edge 3 nodes
// 		  case 25:
// 		  case 27:
// 		    if (elem->node(0) > elem->node(3))f = -1.;
// 		  break;
// 		  }


// 		switch (j)
// 		  {
// 		    // d()/dxi
// 		  case 0:
// 		    return f*(FE<1,BERNSTEIN>::shape_deriv(EDGE3, totalorder, i0[i], 0, xi)*
// 			      FE<1,BERNSTEIN>::shape      (EDGE3, totalorder, i1[i],    eta));

// 		    // d()/deta
// 		  case 1:
// 		    return f*(FE<1,BERNSTEIN>::shape      (EDGE3, totalorder, i0[i],    xi)*
// 			      FE<1,BERNSTEIN>::shape_deriv(EDGE3, totalorder, i1[i], 0, eta));

// 		  default:
// 		    libmesh_error();
// 		  }
// 	      }

// 	    default:
// 	      libmesh_error();
// 	    }
// 	}

//       }
//       // by default throw an error;call the orientation-independent shape functions
//     default:
//       libMesh::err << "ERROR: Unsupported polynomial order!" << std::endl;
//       libmesh_error();
//     }


  libmesh_error();
  return 0.;
}



template <>
Real FE<2,BERNSTEIN>::shape_second_deriv(const ElemType,
					 const Order,
					 const unsigned int,
					 const unsigned int,
					 const Point&)
{
  static bool warning_given = false;

  if (!warning_given)
  libMesh::err << "Second derivatives for Bernstein elements "
                << "are not yet implemented!"
                << std::endl;

  warning_given = true;
  return 0.;
}



template <>
Real FE<2,BERNSTEIN>::shape_second_deriv(const Elem*,
					 const Order,
					 const unsigned int,
					 const unsigned int,
					 const Point&)
{
  static bool warning_given = false;

  if (!warning_given)
  libMesh::err << "Second derivatives for Bernstein elements "
                << "are not yet implemented!"
                << std::endl;

  warning_given = true;
  return 0.;
}

} // namespace libMesh


#endif	// LIBMESH_ENABLE_HIGHER_ORDER_SHAPES

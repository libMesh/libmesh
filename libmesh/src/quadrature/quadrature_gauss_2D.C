// $Id: quadrature_gauss_2D.C,v 1.4 2003-01-24 17:24:45 jwpeterson Exp $

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



// Local includes
#include "quadrature_gauss.h"



void QGauss::init_2D(const ElemType _type)
{
#if DIM > 1
  
  //-----------------------------------------------------------------------
  // 2D quadrature rules
  switch (_type)
    {


      //---------------------------------------------
      // Quadrilateral quadrature rules
    case QUAD4:
    case QUAD8:
    case QUAD9:
      {
	// We compute the 2D quadrature rule as a tensor
	// product of the 1D quadrature rule.
	QGauss q1D(1,_order);
	q1D.init(EDGE2);
	tensor_product_quad( &q1D );
	return;
      }

	    
      //---------------------------------------------
      // Triangle quadrature rules
    case TRI3:
    case TRI6:
      {
	switch(_order)
	  {
	  case CONST:
	  case FIRST:
	    {
	      // Exact for linears
	      _points.resize(1);
	      _weights.resize(1);
		  
	      _points[0](0) = .33333333333333333333333333333333;
	      _points[0](1) = .33333333333333333333333333333333;

	      _weights[0] = .5;

	      return;
	    }
	  case SECOND:
	    {
	      // Exact for quadratics
	      _points.resize(3);
	      _weights.resize(3);

	      _points[0](0) = .5;
	      _points[0](1) = .5;

	      _points[1](0) = 0.;
	      _points[1](1) = .5;

	      _points[2](0) = .5;
	      _points[2](1) = .0;


	      _weights[0] = 1./6.;
	      _weights[1] = 1./6.;
	      _weights[2] = 1./6.;

	      return;
	    }
	  case THIRD:
	    {
	      // Exact for cubics
	      _points.resize(4);
	      _weights.resize(4);
		  
	      _points[0](0) = .33333333333333333333333333333333;
	      _points[0](1) = .33333333333333333333333333333333;

	      _points[1](0) = .2;
	      _points[1](1) = .6;

	      _points[2](0) = .2;
	      _points[2](1) = .2;

	      _points[3](0) = .6;
	      _points[3](1) = .2;


	      _weights[0] = -27./96.;
	      _weights[1] =  25./96.;
	      _weights[2] =  25./96.;
	      _weights[3] =  25./96.;

	      return;
	    }
	  case FOURTH:
	  case FIFTH:
	    {
	      // Exact for quintics
	      // Taken from pg. 222 of "The finite element method," vol. 1
	      // ed. 5 by Zienkiewicz & Taylor
	      _points.resize(7);
	      _weights.resize(7);
		  
	      const real a1 = .0597158717;
	      const real b1 = .4701420641;
	      const real a2 = .7974269853;
	      const real b2 = .1012865073;
		  
	      _points[0](0) = .33333333333333333333333333333333;
	      _points[0](1) = .33333333333333333333333333333333;

	      _points[1](0) = a1;
	      _points[1](1) = b1;

	      _points[2](0) = b1;
	      _points[2](1) = a1;

	      _points[3](0) = b1;
	      _points[3](1) = b1;

	      _points[4](0) = a2;
	      _points[4](1) = b2;

	      _points[5](0) = b2;
	      _points[5](1) = a2;

	      _points[6](0) = b2;
	      _points[6](1) = b2;


	      _weights[0] = .1125;
	      _weights[1] = .06619707635;
	      _weights[2] = _weights[1];
	      _weights[3] = _weights[1];
	      _weights[4] = .06296959025;
	      _weights[5] = _weights[4];
	      _weights[6] = _weights[4];

	      return;
	    }

	  case SIXTH:
	  case SEVENTH:
	    {
	      // This quadrature rule combines two
	      // four-point rules on the unit interval
	      // to create a 16-point, 7th-order accurate
	      // quadrature rule.  For information, see
	      // Approximate Calculation of Multiple Integrals,
	      // Stroud, A. H. p 314
	      std::vector<real> r(4);
	      std::vector<real> s(4);
	      std::vector<real> A(4);
	      std::vector<real> B(4);

	      // Interval Quadrature points
	      r[0] = 0.0694318422; s[0] = 0.0571041961;
	      r[1] = 0.3300094782; s[1] = 0.2768430136;
	      r[2] = 0.6699905218; s[2] = 0.5835904324;
	      r[3] = 0.9305681558; s[3] = 0.8602401357;

	      // Interval Quadrature Weights
	      A[0] = 0.1739274226; B[0] = 0.1355069134;
	      A[1] = 0.3260725774; B[1] = 0.2034645680;
	      A[2] = 0.3260725774; B[2] = 0.1298475476;
	      A[3] = 0.1739274226; B[3] = 0.0311809709;

	      // Compute the conical products
	      _points.resize(16);
	      _weights.resize(16);
	      unsigned int gp = 0;
	      for (unsigned int i=0; i<4; i++)
		for (unsigned int j=0; j<4; j++)
		  {
		    _points[gp](0) = s[j];
		    _points[gp](1) = r[i]*(1-s[j]);
		    _weights[gp]   = A[i]*B[j];
		    gp++;
		  }
	      return;
	    }
	    
	  default:
	    {
	      std::cout << "Quadrature rule not supported!" << std::endl;

	      error();
	    }
	  };
      }

	    
      //---------------------------------------------
      // Unsupported type
    default:
      {
	std::cerr << "Element type not supported!:" << _type << std::endl;
	error();
      }
    };

  error();

  return;

#endif
};



void QGauss::init_2D(const ElemType _type,
		     const unsigned int side)
{
#if DIM > 1
  
  //-----------------------------------------------------------------------
  // 2D quadrature rules
  switch (_type)
    {            
      //---------------------------------------------
      // Quadrilateral quadrature rules
    case QUAD4:
    case QUAD8:
    case QUAD9:
      {
	QGauss q1D(1,_order);
	q1D.init(EDGE2);
	side_rule_quad(&q1D, side);
	return;
      }

	    
      //---------------------------------------------
      // Triangle quadrature rules
    case TRI3:
    case TRI6:
      {
	QGauss q1D(1,_order);
	q1D.init(EDGE2);
	side_rule_tri(&q1D, side);
	return;
      }
	    
    default:
      {
	std::cerr << "Element type not supported!:" << _type << std::endl;
	error();
      }
    };
  
  error();
  
  return;

#endif
};




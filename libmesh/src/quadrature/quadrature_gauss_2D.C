// $Id: quadrature_gauss_2D.C,v 1.11 2003-05-22 17:06:23 jwpeterson Exp $

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
#include "quadrature_jacobi.h"


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
	tensor_product_quad( q1D );
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
		  
	      const Real a1 = .0597158717;
	      const Real b1 = .4701420641;
	      const Real a2 = .7974269853;
	      const Real b2 = .1012865073;
		  
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
	  case EIGHTH:
	  case NINTH:     
	  case TENTH:        
	  case ELEVENTH:     
	  case TWELFTH:      
	  case THIRTEENTH:   
	  case FOURTEENTH:   
	  case FIFTEENTH:    
	  case SIXTEENTH:    
	  case SEVENTEENTH:  
	  case EIGHTTEENTH:  
	  case NINTEENTH:    
	  case TWENTIETH:    
	  case TWENTYFIRST:  
	  case TWENTYSECOND: 
	  case TWENTYTHIRD:  
	    {
	      // The following quadrature rules are
	      // generated as conical products.  These
	      // tend to be non-optimal (use too many
	      // points, cluster points in certain
	      // regions of the domain) but they are
	      // quite easy to automatically generate
	      // using a 1D Gauss rule on [0,1] and a
	      // 1D Jacobi-Gauss rule on [0,1].

	      // Define the quadrature rules...
	      QGauss  gauss1D(1,_order);
	      QJacobi jac1D(1,_order,1,0);
	      
	      // The Gauss rule needs to be scaled to [0,1]
	      std::pair<Real, Real> old_range(-1,1);
	      std::pair<Real, Real> new_range(0,1);
	      gauss1D.scale(old_range,
			    new_range);

	      // Compute the tensor product
	      tensor_product_tri(gauss1D, jac1D);
	      return;
	    }
	    
	  default:
	    {
	      std::cout << "Quadrature rule not supported!" << std::endl;

	      error();
	    }
	  }
      }

	    
      //---------------------------------------------
      // Unsupported type
    default:
      {
	std::cerr << "Element type not supported!:" << _type << std::endl;
	error();
      }
    }

  error();

  return;

#endif
}

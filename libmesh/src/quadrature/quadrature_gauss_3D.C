// $Id: quadrature_gauss_3D.C,v 1.11 2003-05-22 17:06:23 jwpeterson Exp $

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


void QGauss::init_3D(const ElemType _type)
{
#if DIM == 3
  
  //-----------------------------------------------------------------------
  // 3D quadrature rules
  switch (_type)
    {
      //---------------------------------------------
      // Hex quadrature rules
    case HEX8:
    case HEX20:
    case HEX27:
      {
	// We compute the 3D quadrature rule as a tensor
	// product of the 1D quadrature rule.
	QGauss q1D(1,_order);
	q1D.init(EDGE2);

	tensor_product_hex( q1D );
	
	return;
      }


      
      //---------------------------------------------
      // Tetrahedral quadrature rules
    case TET4:
    case TET10:
      {
	// Taken from pg. 222 of "The finite element method," vol. 1
	// ed. 5 by Zienkiewicz & Taylor
	      
	switch(_order)
	  {
	  case CONST:
	  case FIRST:
	    {
	      // Exact for linears
	      _points.resize(1);
	      _weights.resize(1);
		    
		    
	      _points[0](0) = .25;
	      _points[0](1) = .25;
	      _points[0](2) = .25;
		    
	      _weights[0] = .1666666666666666666666666666666666666666666667;
		    
	      return;
	    }
	  case SECOND:
	    {
	      // Exact for quadratics
	      _points.resize(4);
	      _weights.resize(4);
		    
		    
	      const Real a = .585410196624969;
	      const Real b = .138196601125011;
		    
	      _points[0](0) = a;
	      _points[0](1) = b;
	      _points[0](2) = b;
		    
	      _points[1](0) = b;
	      _points[1](1) = a;
	      _points[1](2) = b;
		    
	      _points[2](0) = b;
	      _points[2](1) = b;
	      _points[2](2) = a;
		    
	      _points[3](0) = b;
	      _points[3](1) = b;
	      _points[3](2) = b;
		    
		    
		    
	      _weights[0] = .0416666666666666666666666666666666666666666667;
	      _weights[1] = _weights[0];
	      _weights[2] = _weights[0];
	      _weights[3] = _weights[0];
		    
	      return;
	    }
	  case THIRD:
	    {
	      // Exact for cubics
	      _points.resize(5);
	      _weights.resize(5);
		    
		    
	      _points[0](0) = .25;
	      _points[0](1) = .25;
	      _points[0](2) = .25;
		    
	      _points[1](0) = .5;
	      _points[1](1) = .16666666666666666666666666666666666666666667;
	      _points[1](2) = .16666666666666666666666666666666666666666667;
		    
	      _points[2](0) = .16666666666666666666666666666666666666666667;
	      _points[2](1) = .5;
	      _points[2](2) = .16666666666666666666666666666666666666666667;
		    
	      _points[3](0) = .16666666666666666666666666666666666666666667;
	      _points[3](1) = .16666666666666666666666666666666666666666667;
	      _points[3](2) = .5;
		    
	      _points[4](0) = .16666666666666666666666666666666666666666667;
	      _points[4](1) = .16666666666666666666666666666666666666666667;
	      _points[4](2) = .16666666666666666666666666666666666666666667;
		    
		    
	      _weights[0] = -.133333333333333333333333333333333333333333333;
	      _weights[1] = .075;
	      _weights[2] = _weights[1];
	      _weights[3] = _weights[1];
	      _weights[4] = _weights[1];
		    
	      return;
	    }	    
	  case FOURTH:
	    {
	      _points.resize(11);
	      _weights.resize(11);
		    
	      _points[0](0) = 0.25;
	      _points[0](1) = 0.25;
	      _points[0](2) = 0.25;

	      {
		const Real a = 0.785714285714286;
		const Real b = 0.071428571428571;
		
		_points[1](0) = a;
		_points[1](1) = b;
		_points[1](2) = b;
		
		_points[2](0) = b;
		_points[2](1) = a;
		_points[2](2) = b;
		
		_points[3](0) = b;
		_points[3](1) = b;
		_points[3](2) = a;
		
		_points[4](0) = b;
		_points[4](1) = b;
		_points[4](2) = b;
	      }
	      {
		const Real a = 0.399403576166799;
		const Real b = 0.100596423833201;
		
		_points[5](0) = a;
		_points[5](1) = a;
		_points[5](2) = a;
		
		_points[6](0) = a;
		_points[6](1) = a;
		_points[6](2) = b;
		
		_points[7](0) = a;
		_points[7](1) = b;
		_points[7](2) = b;
		
		_points[8](0) = b;
		_points[8](1) = a;
		_points[8](2) = b;
		
		_points[9](0) = b;
		_points[9](1) = b;
		_points[9](2) = a;
		
		_points[10](0) = b;
		_points[10](1) = b;
		_points[10](2) = b;
	      }
	      
	      _weights[0]  = -0.013155555555555555555555555555555555555555555555555555555555555555555556;
	      _weights[1]  =  0.007622222222222222222222222222222222222222222222222222222222222222222222;
	      _weights[2]  = _weights[1];
	      _weights[3]  = _weights[1];
	      _weights[4]  = _weights[1];
	      _weights[5]  =  0.024888888888888888888888888888888888888888888888888888888888888888888889;
	      _weights[6]  = _weights[5];
	      _weights[7]  = _weights[5];
	      _weights[8]  = _weights[5];
	      _weights[9]  = _weights[5];
	      _weights[10] = _weights[5];
		    
	      return;
	    }
	  case FIFTH:
	    {
	      _points.resize(17);
	      _weights.resize(17);
		    
	      _points[0](0) = 0.25;
	      _points[0](1) = 0.25;
	      _points[0](2) = 0.25;

	      {
		const Real a = 0.;
		const Real b = 0.333333333333333333333333333333333333333;
		
		_points[1](0) = a;
		_points[1](1) = b;
		_points[1](2) = b;
		
		_points[2](0) = b;
		_points[2](1) = a;
		_points[2](2) = b;
		
		_points[3](0) = b;
		_points[3](1) = b;
		_points[3](2) = a;
		
		_points[4](0) = b;
		_points[4](1) = b;
		_points[4](2) = b;
	      }
	      {
		const Real a = 0.7272727272727272727272727272727272727272727272727272727;
		const Real b = 0.0909090909090909090909090909090909090909090909090909091;
		
		_points[5](0) = a;
		_points[5](1) = b;
		_points[5](2) = b;
		
		_points[6](0) = b;
		_points[6](1) = a;
		_points[6](2) = b;
		
		_points[7](0) = b;
		_points[7](1) = b;
		_points[7](2) = a;
		
		_points[8](0) = b;
		_points[8](1) = b;
		_points[8](2) = b;
	      }
	      {
		const Real a = 0.066550153573664;
		const Real b = 0.433449846426336;
		
		_points[9](0) = a;
		_points[9](1) = a;
		_points[9](2) = a;
		
		_points[10](0) = a;
		_points[10](1) = a;
		_points[10](2) = b;
		
		_points[11](0) = a;
		_points[11](1) = b;
		_points[11](2) = b;
		
		_points[12](0) = b;
		_points[12](1) = a;
		_points[12](2) = b;
		
		_points[13](0) = b;
		_points[13](1) = b;
		_points[13](2) = a;
		
		_points[14](0) = b;
		_points[14](1) = b;
		_points[14](2) = b;		
	      }
	      
	      _weights[0]  = 0.030283678097089;
	      _weights[1]  = 0.006026785714286;
	      _weights[2]  = _weights[1];
	      _weights[3]  = _weights[1];
	      _weights[4]  = _weights[1];
	      _weights[5]  = 0.011645249086029;
	      _weights[6]  = _weights[5];
	      _weights[7]  = _weights[5];
	      _weights[8]  = _weights[5];
	      _weights[9]  = 0.010949141561386;
	      _weights[10] = _weights[9];
	      _weights[11] = _weights[9];
	      _weights[12] = _weights[9];
	      _weights[13] = _weights[9];
	      _weights[14] = _weights[9];
		    
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
	      // using a 1D Gauss rule on [0,1] and two 
	      // 1D Jacobi-Gauss rules on [0,1].

	      // Points (1D)            // Weights
	      // std::vector<Real> r(4);   std::vector<Real> A(4);
// 	      std::vector<Real> s(4);   std::vector<Real> B(4);
// 	      std::vector<Real> t(4);   std::vector<Real> C(4);

// 	      // Fill in the interval points
// 	      r[0] = 0.0694318422; s[0] = 0.0571041961; t[0] = 0.0485005494;
// 	      r[1] = 0.3300094782; s[1] = 0.2768430136; t[1] = 0.2386007376;
// 	      r[2] = 0.6699905218; s[2] = 0.5835904324; t[2] = 0.5170472951;
// 	      r[3] = 0.9305681558; s[3] = 0.8602401357; t[3] = 0.7958514179;

// 	      // Fill in the weights
// 	      A[0] = 0.1739274226; B[0] = 0.1355069134; C[0] = 0.1108884156;
// 	      A[1] = 0.3260725774; B[1] = 0.2034645680; C[1] = 0.1434587898;
// 	      A[2] = 0.3260725774; B[2] = 0.1298475476; C[2] = 0.0686338872;
// 	      A[3] = 0.1739274226; B[3] = 0.0311809709; C[3] = 0.0103522407;

// 	      // Make room for the points and weights
// 	      _points.resize(64);
// 	      _weights.resize(64);

// 	      // Compute the conical products
// 	      unsigned int gp = 0;
// 	      for (unsigned int i=0; i<4; i++)
// 		for (unsigned int j=0; j<4; j++)
// 		  for (unsigned int k=0; k<4; k++)
// 		    {
// 		      _points[gp](0) = t[k];
// 		      _points[gp](1) = s[j]*(1.-t[k]);
// 		      _points[gp](2) = r[i]*(1.-s[j])*(1.-t[k]);
// 		      _weights[gp]   = A[i]*B[j]*C[k];
// 		      gp++;
// 		    }
	      
	      // Define the quadrature rules...
	      QGauss  gauss1D(1,_order);
	      QJacobi jacA1D(1,_order,1,0);
	      QJacobi jacB1D(1,_order,2,0);

	      // The Gauss rule needs to be scaled to [0,1]
	      std::pair<Real, Real> old_range(-1,1);
	      std::pair<Real, Real> new_range(0,1);
	      gauss1D.scale(old_range,
			    new_range);

	      // Compute the tensor product
	      tensor_product_tet(gauss1D, jacA1D, jacB1D);
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
      // Prism quadrature rules
    case PRISM6:
    case PRISM15:
    case PRISM18:
      {
	// We compute the 3D quadrature rule as a tensor
	// product of the 1D quadrature rule and a 2D
	// triangle quadrature rule
	    
	QGauss q1D(1,_order);
	QGauss q2D(2,_order);

	// Initialize 
	q1D.init(EDGE2);
	q2D.init(TRI3);

	tensor_product_prism(q1D, q2D);
	
	return;
      }
      

      
      //---------------------------------------------
      // Pyramid
    case PYRAMID5:
      {
	// We compute the Pyramid rule as a conical
	// product of the interval [0,1] and the
	// reference square [-1,1] x [-1,1] as per
	// Stroud, A.H. "Approximate Calculation of
	// Multiple Integrals."
	// This should be exact for quadratics, (Stroud, 32)

	// Get a rule for the reference quad
	QGauss q2D(2,_order); 
	q2D.init(QUAD4);
	
	// Get a rule for the interval [-1,1] 
	QGauss q1D(1,_order);
	q1D.init(EDGE2);

	// Storage for our temporary rule
	std::vector<Real> pts(q1D.n_points());
	std::vector<Real> wts(q1D.n_points());

	// Modify the 1D rule to fit in [0,1] instead
	for (unsigned int i=0; i<q1D.n_points(); i++)
	  {
	    pts[i] = 0.5*q1D.qp(i)(0) + 0.5;
	    wts[i] = 0.5*q1D.w(i); 
	  }

	// Allocate space for the new rule
	_points.resize(q1D.n_points() * q2D.n_points());
	_weights.resize(q1D.n_points() * q2D.n_points());
	      
	// Compute the conical product of the 1D and 2D rules
	unsigned int qp = 0;
	for (unsigned int i=0; i<q1D.n_points(); i++)
	  for (unsigned int j=0; j<q2D.n_points(); j++)
	    {
	      _points[qp](0) = q2D.qp(j)(0)*(1.0-pts[i]);
	      _points[qp](1) = q2D.qp(j)(1)*(1.0-pts[i]);
	      _points[qp](2) = pts[i];

	      _weights[qp] = 0.33333333333333 * // Scale factor!
		wts[i] * q2D.w(j);
		    
	      qp++;
	    }

	return;

      }


      
      //---------------------------------------------
      // Unsupported type
    default:
      {
	std::cerr << "ERROR: Unsupported type: " << _type << std::endl;
	error();
      }
    }

  error();

  return;
  
#endif
}

// $Id$

// The libMesh Finite Element Library.
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
#include "quadrature_monomial.h"
#include "quadrature_gauss.h"


void QMonomial::init_3D(const ElemType _type,
			unsigned int p)
{

  switch (_type)
    {
      //---------------------------------------------
      // Hex quadrature rules
    case HEX8:
    case HEX20:
    case HEX27:
      {
	switch(_order + 2*p)
	  {

	  case SEVENTH:
	    {
	      // A degree 7, 38-point, "rotationally-symmetric" rule by
	      // Kim and Song, Comm. Korean Math. Soc vol. 13, no. 4, 1998, pp. 913-931.
	      //
	      // A SEVENTH-order Gauss product rule (which integrates tri-7th order polynomials)
	      // would have 4^3=64 points.
	      const Real data[3][4] =
		{
		  {9.01687807821291289082811566285950e-01L, 0.00000000000000000000000000000000e+00L, 0.00000000000000000000000000000000e+00L, 2.95189738262622903181631100062774e-01L}, 
		  {4.08372221499474674069588900002128e-01L, 4.08372221499474674069588900002128e-01L, 4.08372221499474674069588900002128e-01L, 4.04055417266200582425904380777126e-01L}, 
		  {8.59523090201054193116477875786220e-01L, 8.59523090201054193116477875786220e-01L, 4.14735913727987720499709244748633e-01L, 1.24850759678944080062624098058597e-01L}  
		};

	      const unsigned int rule_id[3] = {
		1, // (x,0,0) -> 6 permutations
		4, // (x,x,x) -> 8 permutations
		5  // (x,x,z) -> 24 permutations
	      };

	      _points.resize(38);
	      _weights.resize(38);

	      kim_rule(data, rule_id, 3);
	      return;
	    } // end case SEVENTH

	  case EIGHTH:
	    {
	      // A degree 8, 47-point, "rotationally-symmetric" rule by
	      // Kim and Song, Comm. Korean Math. Soc vol. 13, no. 4, 1998, pp. 913-931.
	      //
	      // A EIGHTH-order Gauss product rule (which integrates tri-8th order polynomials)
	      // would have 5^3=125 points.
	      const Real data[5][4] =
		{
		  {0.00000000000000000000000000000000e+00L, 0.00000000000000000000000000000000e+00L, 0.00000000000000000000000000000000e+00L, 4.51903714875199690490763818699555e-01L}, 
		  {7.82460796435951590652813975429717e-01L, 0.00000000000000000000000000000000e+00L, 0.00000000000000000000000000000000e+00L, 2.99379177352338919703385618576171e-01L}, 
		  {4.88094669706366480526729301468686e-01L, 4.88094669706366480526729301468686e-01L, 4.88094669706366480526729301468686e-01L, 3.00876159371240019939698689791164e-01L}, 
		  {8.62218927661481188856422891110042e-01L, 8.62218927661481188856422891110042e-01L, 8.62218927661481188856422891110042e-01L, 4.94843255877038125738173175714853e-02L},  
		  {2.81113909408341856058098281846420e-01L, 9.44196578292008195318687494773744e-01L, 6.97574833707236996779391729948984e-01L, 1.22872389222467338799199767122592e-01L}  
		};

	      const unsigned int rule_id[5] = {
		0, // (0,0,0) -> 1 permutation
		1, // (x,0,0) -> 6 permutations
		4, // (x,x,x) -> 8 permutations
		4, // (x,x,x) -> 8 permutations
		6  // (x,y,z) -> 24 permutations
	      };

	      _points.resize(47);
	      _weights.resize(47);

	      kim_rule(data, rule_id, 5);
	      return;
	    } // end case EIGHTH

	    
	    // By default: construct and use a Gauss quadrature rule
	  default:
	    {
	      // Break out and fall down into the default: case for the
	      // outer switch statement.
	      break;
	    }
	    
	  } // end switch(_order + 2*p)
      } // end case HEX8/20/27

      
      // By default: construct and use a Gauss quadrature rule
    default:
      {
	QGauss gauss_rule(3, _order);
	gauss_rule.init(_type, p);

	// Swap points and weights with the about-to-be destroyed rule.
	_points.swap (gauss_rule.get_points() );
	_weights.swap(gauss_rule.get_weights());

	return;
      }
    } // end switch (_type)
}

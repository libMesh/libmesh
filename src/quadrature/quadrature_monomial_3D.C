// The libMesh Finite Element Library.
// Copyright (C) 2002-2012 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/quadrature_monomial.h"
#include "libmesh/quadrature_gauss.h"

namespace libMesh
{


void QMonomial::init_3D(const ElemType type_in,
			unsigned int p)
{

  switch (type_in)
    {
      //---------------------------------------------
      // Hex quadrature rules
    case HEX8:
    case HEX20:
    case HEX27:
      {
	switch(_order + 2*p)
	  {

	    // The CONSTANT/FIRST rule is the 1-point Gauss "product" rule...we fall
	    // through to the default case for this rule.

	  case SECOND:
	  case THIRD:
	    {
	      // A degree 3, 6-point, "rotationally-symmetric" rule by
	      // Kim and Song, Comm. Korean Math. Soc vol. 13, no. 4, 1998, pp. 913-931.
	      //
	      // Warning: this rule contains points on the boundary of the reference
	      // element, and therefore may be unsuitable for some problems.  The alternative
	      // would be a 2x2x2 Gauss product rule.
	      const Real data[1][4] =
		{
		  {1.0L, 0.0L, 0.0L, static_cast<Real>(4.0L/3.0L)}
		};

	      const unsigned int rule_id[1] = {
		1 // (x,0,0) -> 6 permutations
	      };

	      _points.resize(6);
	      _weights.resize(6);

	      kim_rule(data, rule_id, 1);
	      return;
	    } // end case SECOND,THIRD

	  case FOURTH:
	  case FIFTH:
	    {
	      // A degree 5, 13-point rule by Stroud,
	      // AH Stroud, "Some Fifth Degree Integration Formulas for Symmetric Regions II.",
	      // Numerische Mathematik 9, pp. 460-468 (1967).
	      //
	      // This rule is provably minimal in the number of points.  The equations given for
	      // the n-cube on pg. 466 of the paper for mu/gamma and gamma are wrong, at least for
	      // the n=3 case.  The analytical values given here were computed by me [JWP] in Maple.

	      // Convenient intermediate values.
	      const Real sqrt19 = std::sqrt(19.L);
	      const Real tp     = std::sqrt(71440.L + 6802.L*sqrt19);

	      // Point data for permutations.
	      const Real eta    =  0.00000000000000000000000000000000e+00L;

	      const Real lambda =  std::sqrt(1919.L/3285.L - 148.L*sqrt19/3285.L + 4.L*tp/3285.L);
	      // 8.8030440669930978047737818209860e-01L;

	      const Real xi     = -std::sqrt(1121.L/3285.L +  74.L*sqrt19/3285.L - 2.L*tp/3285.L);
	      // -4.9584817142571115281421242364290e-01L;

	      const Real mu     =  std::sqrt(1121.L/3285.L +  74.L*sqrt19/3285.L + 2.L*tp/3285.L);
	      // 7.9562142216409541542982482567580e-01L;

	      const Real gamma  =  std::sqrt(1919.L/3285.L - 148.L*sqrt19/3285.L - 4.L*tp/3285.L);
	      // 2.5293711744842581347389255929324e-02L;

	      // Weights: the centroid weight is given analytically.  Weight B (resp C) goes
	      // with the {lambda,xi} (resp {gamma,mu}) permutation.  The single-precision
	      // results reported by Stroud are given for reference.

	      const Real A      = 32.0L / 19.0L;
	      // Stroud: 0.21052632  * 8.0 = 1.684210560;

	      const Real B      = 1.L / ( 260072.L/133225.L  - 1520*sqrt19/133225.L + (133.L - 37.L*sqrt19)*tp/133225.L );
	      // 5.4498735127757671684690782180890e-01L; // Stroud: 0.068123420 * 8.0 = 0.544987360;

	      const Real C      = 1.L / ( 260072.L/133225.L  - 1520*sqrt19/133225.L - (133.L - 37.L*sqrt19)*tp/133225.L );
	      // 5.0764422766979170420572375713840e-01L; // Stroud: 0.063455527 * 8.0 = 0.507644216;

 	      _points.resize(13);
 	      _weights.resize(13);

	      unsigned int c=0;

	      // Point with weight A (origin)
	      _points[c] = Point(eta, eta, eta);
	      _weights[c++] = A;

	      // Points with weight B
	      _points[c] = Point(lambda, xi, xi);
	      _weights[c++] = B;
	      _points[c] = -_points[c-1];
	      _weights[c++] = B;

	      _points[c] = Point(xi, lambda, xi);
	      _weights[c++] = B;
	      _points[c] = -_points[c-1];
	      _weights[c++] = B;

	      _points[c] = Point(xi, xi, lambda);
	      _weights[c++] = B;
	      _points[c] = -_points[c-1];
	      _weights[c++] = B;

	      // Points with weight C
	      _points[c] = Point(mu, mu, gamma);
	      _weights[c++] = C;
	      _points[c] = -_points[c-1];
	      _weights[c++] = C;

	      _points[c] = Point(mu, gamma, mu);
	      _weights[c++] = C;
	      _points[c] = -_points[c-1];
	      _weights[c++] = C;

	      _points[c] = Point(gamma, mu, mu);
	      _weights[c++] = C;
	      _points[c] = -_points[c-1];
	      _weights[c++] = C;

	      return;


// 	      // A degree 5, 14-point, "rotationally-symmetric" rule by
// 	      // Kim and Song, Comm. Korean Math. Soc vol. 13, no. 4, 1998, pp. 913-931.
// 	      // Was also reported in Stroud's 1971 book.
// 	      const Real data[2][4] =
// 		{
// 		  {7.95822425754221463264548820476135e-01L, 0.00000000000000000000000000000000e+00L, 0.00000000000000000000000000000000e+00L, 8.86426592797783933518005540166204e-01L},
// 		  {7.58786910639328146269034278112267e-01L, 7.58786910639328146269034278112267e-01L, 7.58786910639328146269034278112267e-01L, 3.35180055401662049861495844875346e-01L}
// 		};

// 	      const unsigned int rule_id[2] = {
// 		1, // (x,0,0) -> 6 permutations
// 		4  // (x,x,x) -> 8 permutations
// 	      };

// 	      _points.resize(14);
// 	      _weights.resize(14);

// 	      kim_rule(data, rule_id, 2);
// 	      return;
	    } // end case FOURTH,FIFTH

	  case SIXTH:
	  case SEVENTH:
	    {
	      if (allow_rules_with_negative_weights)
		{
		  // A degree 7, 31-point, "rotationally-symmetric" rule by
		  // Kim and Song, Comm. Korean Math. Soc vol. 13, no. 4, 1998, pp. 913-931.
		  // This rule contains a negative weight, so only use it if such type of
		  // rules are allowed.
		  const Real data[3][4] =
		    {
		      {0.00000000000000000000000000000000e+00L, 0.00000000000000000000000000000000e+00L, 0.00000000000000000000000000000000e+00L, -1.27536231884057971014492753623188e+00L},
		      {5.85540043769119907612630781744060e-01L, 0.00000000000000000000000000000000e+00L, 0.00000000000000000000000000000000e+00L,  8.71111111111111111111111111111111e-01L},
		      {6.94470135991704766602025803883310e-01L, 9.37161638568208038511047377665396e-01L, 4.15659267604065126239606672567031e-01L,  1.68695652173913043478260869565217e-01L}
		    };

		  const unsigned int rule_id[3] = {
		    0, // (0,0,0) -> 1 permutation
		    1, // (x,0,0) -> 6 permutations
		    6  // (x,y,z) -> 24 permutations
		  };

		  _points.resize(31);
		  _weights.resize(31);

		  kim_rule(data, rule_id, 3);
		  return;
		} // end if (allow_rules_with_negative_weights)


	      // A degree 7, 34-point, "fully-symmetric" rule, first published in
	      // P.C. Hammer and A.H. Stroud, "Numerical Evaluation of Multiple Integrals II",
	      // Mathmatical Tables and Other Aids to Computation, vol 12., no 64, 1958, pp. 272-280
	      //
	      // This rule happens to fall under the same general
	      // construction as the Kim rules, so we've re-used
	      // that code here.  Stroud gives 16 digits for his rule,
	      // and this is the most accurate version I've found.
	      //
	      // For comparison, a SEVENTH-order Gauss product rule
	      // (which integrates tri-7th order polynomials) would
	      // have 4^3=64 points.
	      const Real
		r  = std::sqrt(6.L/7.L),
		s  = std::sqrt((960.L - 3.L*std::sqrt(28798.L)) / 2726.L),
		t  = std::sqrt((960.L + 3.L*std::sqrt(28798.L)) / 2726.L),
		B1 = 8624.L / 29160.L,
		B2 = 2744.L / 29160.L,
		B3 = 8.L*(774.L*t*t - 230.L)/(9720.L*(t*t-s*s)),
		B4 = 8.L*(230.L - 774.L*s*s)/(9720.L*(t*t-s*s));

	      const Real data[4][4] =
		{
		  {r, 0.L, 0.L, B1},
		  {r,   r, 0.L, B2},
		  {s,   s,   s, B3},
		  {t,   t,   t, B4}
		};

	      const unsigned int rule_id[4] = {
		1, // (x,0,0) -> 6 permutations
		2, // (x,x,0) -> 12 permutations
		4, // (x,x,x) -> 8 permutations
		4  // (x,x,x) -> 8 permutations
		  };

	      _points.resize(34);
	      _weights.resize(34);

	      kim_rule(data, rule_id, 4);
	      return;


//	      // A degree 7, 38-point, "rotationally-symmetric" rule by
//	      // Kim and Song, Comm. Korean Math. Soc vol. 13, no. 4, 1998, pp. 913-931.
//	      //
//	      // This rule is obviously inferior to the 34-point rule above...
//	      const Real data[3][4] =
//		{
//		  {9.01687807821291289082811566285950e-01L, 0.00000000000000000000000000000000e+00L, 0.00000000000000000000000000000000e+00L, 2.95189738262622903181631100062774e-01L},
//		  {4.08372221499474674069588900002128e-01L, 4.08372221499474674069588900002128e-01L, 4.08372221499474674069588900002128e-01L, 4.04055417266200582425904380777126e-01L},
//		  {8.59523090201054193116477875786220e-01L, 8.59523090201054193116477875786220e-01L, 4.14735913727987720499709244748633e-01L, 1.24850759678944080062624098058597e-01L}
//		};
//
//	      const unsigned int rule_id[3] = {
//		1, // (x,0,0) -> 6 permutations
//		4, // (x,x,x) -> 8 permutations
//		5  // (x,x,z) -> 24 permutations
//	      };
//
//	      _points.resize(38);
//	      _weights.resize(38);
//
//	      kim_rule(data, rule_id, 3);
//	      return;
	    } // end case SIXTH,SEVENTH

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
	gauss_rule.init(type_in, p);

	// Swap points and weights with the about-to-be destroyed rule.
	_points.swap (gauss_rule.get_points() );
	_weights.swap(gauss_rule.get_weights());

	return;
      }
    } // end switch (type_in)
}

} // namespace libMesh

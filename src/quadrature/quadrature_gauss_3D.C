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
#include "libmesh/quadrature_gauss.h"
#include "libmesh/quadrature_conical.h"
#include "libmesh/quadrature_gm.h"

namespace libMesh
{

void QGauss::init_3D(const ElemType type_in,
                     unsigned int p)
{
#if LIBMESH_DIM == 3

  //-----------------------------------------------------------------------
  // 3D quadrature rules
  switch (type_in)
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
	q1D.init(EDGE2,p);

	tensor_product_hex( q1D );

	return;
      }



      //---------------------------------------------
      // Tetrahedral quadrature rules
    case TET4:
    case TET10:
      {
	switch(_order + 2*p)
	  {
	    // Taken from pg. 222 of "The finite element method," vol. 1
	    // ed. 5 by Zienkiewicz & Taylor
	  case CONSTANT:
	  case FIRST:
	    {
	      // Exact for linears
	      _points.resize(1);
	      _weights.resize(1);


	      _points[0](0) = .25;
	      _points[0](1) = .25;
	      _points[0](2) = .25;

	      _weights[0] = .1666666666666666666666666666666666666666666667L;

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



	    // Can be found in the class notes
	    // http://www.cs.rpi.edu/~flaherje/FEM/fem6.ps
	    // by Flaherty.
	    //
	    // Caution: this rule has a negative weight and may be
	    // unsuitable for some problems.
	    // Exact for cubics.
	    //
	    // Note: Keast (see ref. elsewhere in this file) also gives
	    // a third-order rule with positive weights, but it contains points
	    // on the ref. elt. boundary, making it less suitable for FEM calculations.
	  case THIRD:
	    {
	      if (allow_rules_with_negative_weights)
		{
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
		} // end if (allow_rules_with_negative_weights)
	      else
		{
		  // If a rule with postive weights is required, the 2x2x2 conical
		  // product rule is third-order accurate and has less points than
		  // the next-available positive-weight rule at FIFTH order.
		  QConical conical_rule(3, _order);
		  conical_rule.init(type_in, p);

		  // Swap points and weights with the about-to-be destroyed rule.
		  _points.swap (conical_rule.get_points() );
		  _weights.swap(conical_rule.get_weights());

		  return;
		}
	      // Note: if !allow_rules_with_negative_weights, fall through to next case.
	    }



	    // Originally a Keast rule,
	    //    Patrick Keast,
	    //    Moderate Degree Tetrahedral Quadrature Formulas,
	    //    Computer Methods in Applied Mechanics and Engineering,
	    //    Volume 55, Number 3, May 1986, pages 339-348.
	    //
	    // Can also be found the class notes
	    // http://www.cs.rpi.edu/~flaherje/FEM/fem6.ps
	    // by Flaherty.
	    //
	    // Caution: this rule has a negative weight and may be
	    // unsuitable for some problems.
	  case FOURTH:
	    {
	      if (allow_rules_with_negative_weights)
		{
		  _points.resize(11);
		  _weights.resize(11);

		  // The raw data for the quadrature rule.
		  const Real rule_data[3][4] = {
		    {0.250000000000000000e+00,                         0.,                            0.,  -0.131555555555555556e-01},  // 1
		    {0.785714285714285714e+00,   0.714285714285714285e-01,                            0.,   0.762222222222222222e-02},  // 4
		    {0.399403576166799219e+00,                         0.,      0.100596423833200785e+00,   0.248888888888888889e-01}   // 6
		  };


		  // Now call the keast routine to generate _points and _weights
		  keast_rule(rule_data, 3);

		  return;
		} // end if (allow_rules_with_negative_weights)
	      // Note: if !allow_rules_with_negative_weights, fall through to next case.
	    }




	    // Walkington's fifth-order 14-point rule from
	    // "Quadrature on Simplices of Arbitrary Dimension"
	    //
	    // We originally had a Keast rule here, but this rule had
	    // more points than an equivalent rule by Walkington and
	    // also contained points on the boundary of the ref. elt,
	    // making it less suitable for FEM calculations.
	  case FIFTH:
	    {
	      _points.resize(14);
	      _weights.resize(14);

	      // permutations of these points and suitably-modified versions of
	      // these points are the quadrature point locations
	      const Real a[3] = {0.31088591926330060980,    // a1 from the paper
				 0.092735250310891226402,   // a2 from the paper
				 0.045503704125649649492};  // a3 from the paper

	      // weights.  a[] and wt[] are the only floating-point inputs required
	      // for this rule.
	      const Real wt[3] = {0.018781320953002641800,    // w1 from the paper
                                  0.012248840519393658257,    // w2 from the paper
                                  0.0070910034628469110730};  // w3 from the paper

	      // The first two sets of 4 points are formed in a similar manner
	      for (unsigned int i=0; i<2; ++i)
		{
		  // Where we will insert values into _points and _weights
		  const unsigned int offset=4*i;

		  // Stuff points and weights values into their arrays
		  const Real b = 1. - 3.*a[i];

		  // Here are the permutations.  Order of these is not important,
		  // all have the same weight
		  _points[offset + 0] = Point(a[i], a[i], a[i]);
		  _points[offset + 1] = Point(a[i],    b, a[i]);
		  _points[offset + 2] = Point(   b, a[i], a[i]);
		  _points[offset + 3] = Point(a[i], a[i],    b);

		  // These 4 points all have the same weights
		  for (unsigned int j=0; j<4; ++j)
		    _weights[offset + j] = wt[i];
		} // end for


	      {
		// The third set contains 6 points and is formed a little differently
		const unsigned int offset = 8;
		const Real b = 0.5*(1. - 2.*a[2]);

		// Here are the permutations.  Order of these is not important,
		// all have the same weight
		_points[offset + 0] = Point(b   ,    b, a[2]);
		_points[offset + 1] = Point(b   , a[2], a[2]);
		_points[offset + 2] = Point(a[2], a[2],    b);
		_points[offset + 3] = Point(a[2],    b, a[2]);
		_points[offset + 4] = Point(   b, a[2],    b);
		_points[offset + 5] = Point(a[2],    b,    b);

		// These 6 points all have the same weights
		for (unsigned int j=0; j<6; ++j)
		  _weights[offset + j] = wt[2];
	      }


	      // Original rule by Keast, unsuitable because it has points on the
	      // reference element boundary!
	      // 	      _points.resize(15);
	      // 	      _weights.resize(15);

	      // 	      _points[0](0) = 0.25;
	      // 	      _points[0](1) = 0.25;
	      // 	      _points[0](2) = 0.25;

	      // 	      {
	      // 		const Real a = 0.;
	      // 		const Real b = 0.333333333333333333333333333333333333333;

	      // 		_points[1](0) = a;
	      // 		_points[1](1) = b;
	      // 		_points[1](2) = b;

	      // 		_points[2](0) = b;
	      // 		_points[2](1) = a;
	      // 		_points[2](2) = b;

	      // 		_points[3](0) = b;
	      // 		_points[3](1) = b;
	      // 		_points[3](2) = a;

	      // 		_points[4](0) = b;
	      // 		_points[4](1) = b;
	      // 		_points[4](2) = b;
	      // 	      }
	      // 	      {
	      // 		const Real a = 0.7272727272727272727272727272727272727272727272727272727;
	      // 		const Real b = 0.0909090909090909090909090909090909090909090909090909091;

	      // 		_points[5](0) = a;
	      // 		_points[5](1) = b;
	      // 		_points[5](2) = b;

	      // 		_points[6](0) = b;
	      // 		_points[6](1) = a;
	      // 		_points[6](2) = b;

	      // 		_points[7](0) = b;
	      // 		_points[7](1) = b;
	      // 		_points[7](2) = a;

	      // 		_points[8](0) = b;
	      // 		_points[8](1) = b;
	      // 		_points[8](2) = b;
	      // 	      }
	      // 	      {
	      // 		const Real a = 0.066550153573664;
	      // 		const Real b = 0.433449846426336;

	      // 		_points[9](0) = b;
	      // 		_points[9](1) = a;
	      // 		_points[9](2) = a;

	      // 		_points[10](0) = a;
	      // 		_points[10](1) = a;
	      // 		_points[10](2) = b;

	      // 		_points[11](0) = a;
	      // 		_points[11](1) = b;
	      // 		_points[11](2) = b;

	      // 		_points[12](0) = b;
	      // 		_points[12](1) = a;
	      // 		_points[12](2) = b;

	      // 		_points[13](0) = b;
	      // 		_points[13](1) = b;
	      // 		_points[13](2) = a;

	      // 		_points[14](0) = a;
	      // 		_points[14](1) = b;
	      // 		_points[14](2) = a;
	      // 	      }

	      // 	      _weights[0]  = 0.030283678097089;
	      // 	      _weights[1]  = 0.006026785714286;
	      // 	      _weights[2]  = _weights[1];
	      // 	      _weights[3]  = _weights[1];
	      // 	      _weights[4]  = _weights[1];
	      // 	      _weights[5]  = 0.011645249086029;
	      // 	      _weights[6]  = _weights[5];
	      // 	      _weights[7]  = _weights[5];
	      // 	      _weights[8]  = _weights[5];
	      // 	      _weights[9]  = 0.010949141561386;
	      // 	      _weights[10] = _weights[9];
	      // 	      _weights[11] = _weights[9];
	      // 	      _weights[12] = _weights[9];
	      // 	      _weights[13] = _weights[9];
	      // 	      _weights[14] = _weights[9];

	      return;
	    }




	    // This rule is originally from Keast:
	    //    Patrick Keast,
	    //    Moderate Degree Tetrahedral Quadrature Formulas,
	    //    Computer Methods in Applied Mechanics and Engineering,
	    //    Volume 55, Number 3, May 1986, pages 339-348.
	    //
	    // It is accurate on 6th-degree polynomials and has 24 points
	    // vs. 64 for the comparable conical product rule.
	    //
	    // Values copied 24th June 2008 from:
	    // http://people.scs.fsu.edu/~burkardt/f_src/keast/keast.f90
	  case SIXTH:
	    {
	      _points.resize (24);
	      _weights.resize(24);

	      // The raw data for the quadrature rule.
	      const Real rule_data[4][4] = {
		{0.356191386222544953e+00 , 0.214602871259151684e+00 ,                       0., 0.00665379170969464506e+00}, // 4
		{0.877978124396165982e+00 , 0.0406739585346113397e+00,                       0., 0.00167953517588677620e+00}, // 4
		{0.0329863295731730594e+00, 0.322337890142275646e+00 ,                       0., 0.00922619692394239843e+00}, // 4
		{0.0636610018750175299e+00, 0.269672331458315867e+00 , 0.603005664791649076e+00, 0.00803571428571428248e+00}  // 12
	      };


	      // Now call the keast routine to generate _points and _weights
	      keast_rule(rule_data, 4);

	      return;
	    }


	    // Keast's 31 point, 7th-order rule contains points on the reference
	    // element boundary, so we've decided not to include it here.
	    //
	    // Keast's 8th-order rule has 45 points.  and a negative
	    // weight, so if you've explicitly disallowed such rules
	    // you will fall through to the conical product rules
	    // below.
	  case SEVENTH:
	  case EIGHTH:
	    {
	      if (allow_rules_with_negative_weights)
		{
		  _points.resize (45);
		  _weights.resize(45);

		  // The raw data for the quadrature rule.
		  const Real rule_data[7][4] = {
		    {0.250000000000000000e+00,                        0.,                        0.,   -0.393270066412926145e-01},  // 1
		    {0.617587190300082967e+00,  0.127470936566639015e+00,                        0.,    0.408131605934270525e-02},  // 4
		    {0.903763508822103123e+00,  0.320788303926322960e-01,                        0.,    0.658086773304341943e-03},  // 4
		    {0.497770956432810185e-01,                        0.,  0.450222904356718978e+00,    0.438425882512284693e-02},  // 6
		    {0.183730447398549945e+00,                        0.,  0.316269552601450060e+00,    0.138300638425098166e-01},  // 6
		    {0.231901089397150906e+00,  0.229177878448171174e-01,  0.513280033360881072e+00,    0.424043742468372453e-02},  // 12
		    {0.379700484718286102e-01,  0.730313427807538396e+00,  0.193746475248804382e+00,    0.223873973961420164e-02}   // 12
		  };


		  // Now call the keast routine to generate _points and _weights
		  keast_rule(rule_data, 7);

		  return;
		} // end if (allow_rules_with_negative_weights)
	      // Note: if !allow_rules_with_negative_weights, fall through to next case.
	    }

	    // Fall back on Grundmann-Moller or Conical Product rules at high orders.
	  default:
	    {
	      if ((allow_rules_with_negative_weights) && (_order + 2*p < 34))
		{
		  // The Grundmann-Moller rules are defined to arbitrary order and
		  // can have significantly fewer evaluation points than concial product
		  // rules.  If you allow rules with negative weights, the GM rules
		  // will be more efficient for degree up to 33 (for degree 34 and
		  // higher, CP is more efficient!) but may be more susceptible
		  // to round-off error.  Safest is to disallow rules with negative
		  // weights, but this decision should be made on a case-by-case basis.
		  QGrundmann_Moller gm_rule(3, _order);
		  gm_rule.init(type_in, p);

		  // Swap points and weights with the about-to-be destroyed rule.
		  _points.swap (gm_rule.get_points() );
		  _weights.swap(gm_rule.get_weights());

		  return;
		}

	      else
		{
		  // The following quadrature rules are generated as
		  // conical products.  These tend to be non-optimal
		  // (use too many points, cluster points in certain
		  // regions of the domain) but they are quite easy to
		  // automatically generate using a 1D Gauss rule on
		  // [0,1] and two 1D Jacobi-Gauss rules on [0,1].
		  QConical conical_rule(3, _order);
		  conical_rule.init(type_in, p);

		  // Swap points and weights with the about-to-be destroyed rule.
		  _points.swap (conical_rule.get_points() );
		  _weights.swap(conical_rule.get_weights());

		  return;
		}
	    }
	  }
      } // end case TET4,TET10



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
	q1D.init(EDGE2,p);
	q2D.init(TRI3,p);

	tensor_product_prism(q1D, q2D);

	return;
      }



      //---------------------------------------------
      // Pyramid
    case PYRAMID5:
    case PYRAMID14:
      {
	// We compute the Pyramid rule as a conical product of a
	// Jacobi rule with alpha==2 on the interval [0,1] two 1D
	// Gauss rules with suitably modified points.  The idea comes
	// from: Stroud, A.H. "Approximate Calculation of Multiple
	// Integrals."
	QConical conical_rule(3, _order);
	conical_rule.init(type_in, p);

	// Swap points and weights with the about-to-be destroyed rule.
	_points.swap (conical_rule.get_points() );
	_weights.swap(conical_rule.get_weights());

	return;

      }



      //---------------------------------------------
      // Unsupported type
    default:
      {
	libMesh::err << "ERROR: Unsupported type: " << type_in << std::endl;
	libmesh_error();
      }
    }

  libmesh_error();

  return;

#endif
}

} // namespace libMesh

// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2007  Benjamin S. Kirk, John W. Peterson
  
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


void QGauss::init_2D(const ElemType _type,
                     unsigned int p)
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
	//
	// For QUADs, a quadrature rule of order 'p' must be able to integrate
	// bilinear (p=1), biquadratic (p=2), bicubic (p=3), etc. polynomials of the form
	//
	// (x^p + x^{p-1} + ... + 1) * (y^p + y^{p-1} + ... + 1)
	//
	// These polynomials have terms *up to* degree 2p but they are *not* complete
	// polynomials of degree 2p. For example, when p=2 we have
	//        1
	//     x      y
	// x^2    xy     y^2
	//    yx^2   xy^2
	//       x^2y^2
	QGauss q1D(1,_order);
	q1D.init(EDGE2,p);
	tensor_product_quad( q1D );
	return;
      }

	    
      //---------------------------------------------
      // Triangle quadrature rules
    case TRI3:
    case TRI6:
      {
	switch(_order + 2*p)
	  {
	  case CONSTANT:
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

	      // Alternate rule with points on ref. elt. boundaries.
	      // Not ideal for problems with material coefficient discontinuities
	      // aligned along element boundaries.
	      // _points[0](0) = .5;
	      // _points[0](1) = .5;
	      // _points[1](0) = 0.;
	      // _points[1](1) = .5;
	      // _points[2](0) = .5;
	      // _points[2](1) = .0;
	      
	      _points[0](0) = 2./3.;
	      _points[0](1) = 1./6.;

	      _points[1](0) = 1./6.;
	      _points[1](1) = 2./3.;

	      _points[2](0) = 1./6.;
	      _points[2](1) = 1./6.;


	      _weights[0] = 1./6.;
	      _weights[1] = 1./6.;
	      _weights[2] = 1./6.;

	      return;
	    }
	  case THIRD:
	    {
	      if (allow_rules_with_negative_weights)
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
		} // end if (allow_rules_with_negative_weights)
	      // Note: if !allow_rules_with_negative_weights, fall through to next case.
	    }


	    
	  case FOURTH:
	    {
	      // A degree 4 rule with six points.  This rule can be found in many places
	      // including:
	      //
	      // J.N. Lyness and D. Jespersen, Moderate degree symmetric
	      // quadrature rules for the triangle, J. Inst. Math. Appl.  15 (1975),
	      // 19--32.

	      _points.resize(6);
	      _weights.resize(6);
	      
	      // The points are arranged symmetrically in two sets ('a' and 'b') of three.
	      const Real a = 9.1576213509770743e-02;  const Real wa = 5.4975871827660933e-02;
	      const Real b = 4.4594849091596488e-01;  const Real wb = 1.1169079483900573e-01;

	      _points[0] = Point(a      ,       a); _weights[0] = wa;
	      _points[1] = Point(a      , 1.-2.*a); _weights[1] = wa;
	      _points[2] = Point(1.-2.*a,       a); _weights[2] = wa;

	      _points[3] = Point(b      ,      b); _weights[3] = wb;
	      _points[4] = Point(b      ,1.-2.*b); _weights[4] = wb;
	      _points[5] = Point(1.-2.*b,      b); _weights[5] = wb;
	      
	      return;
	    }


	    
	  case FIFTH:
	    {
	      // Exact for quintics
	      // Taken from "Quadrature on Simplices of Arbitrary
	      // Dimension" by Walkington
	      _points.resize(7);
	      _weights.resize(7);
		  
	      const Real b1 = 2./7. + std::sqrt(15.)/21.;
	      const Real a1 = 1. - 2.*b1;
	      const Real b2 = 2./7. - std::sqrt(15.)/21.;
	      const Real a2 = 1. - 2.*b2;
		  
	      _points[0](0) = 1./3.;
	      _points[0](1) = 1./3.;

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


	      _weights[0] = 9./80.;
	      _weights[1] = 31./480. + std::sqrt(15.)/2400.;
	      _weights[2] = _weights[1];
	      _weights[3] = _weights[1];
	      _weights[4] = 31./480. - std::sqrt(15.)/2400.;
	      _weights[5] = _weights[4];
	      _weights[6] = _weights[4];

	      return;
	    }




	    
	    // A degree 7 rule with 12 points.  This rule can be found in:
	    //
	    // K. Gatermann, The construction of symmetric cubature
	    // formulas for the square and the triangle, Computing 40
	    // (1988), 229--240.
	    //
	    // This rule, which is provably minimal in the number of
	    // integration points, is said to be 'Ro3 invariant' which
	    // means that a given set of barycentric coordinates
	    // (z1,z2,z3) implies the quadrature points (z1,z2),
	    // (z3,z1), (z2,z3) which are formed by taking the first
	    // two entries in cyclic permutations of the barycentric
	    // point.  Barycentric coordinates are related in the
	    // sense that: z3 = 1 - z1 - z2.
	    //
	    // The 12-point sixth-order rule for triangles given in
	    // Flaherty's (http://www.cs.rpi.edu/~flaherje/FEM/fem6.ps)
	    // lecture notes has been removed in favor of this rule
	    // which is higher-order (for the same number of
	    // quadrature points points) and has a few more digits of
	    // precision in the points and weights.  Some 10-point
	    // degree 6 rules exist for the triangle but they have
	    // quadrature points outside the region of integration.
	  case SIXTH:
	  case SEVENTH:
	    {
	      _points.resize (12);
	      _weights.resize(12);
	      
	      const unsigned int nrows=4;
	      
	      // In each of the rows below, the first two entries are (z1, z2) which imply
	      // z3.  The third entry is the weight for each of the points in the cyclic permutation.
	      const Real p[nrows][3] = {
		{6.2382265094402118e-02, 6.7517867073916085e-02, 2.6517028157436251e-02}, // group A
		{5.5225456656926611e-02, 3.2150249385198182e-01, 4.3881408714446055e-02}, // group B
		{3.4324302945097146e-02, 6.6094919618673565e-01, 2.8775042784981585e-02}, // group C
		{5.1584233435359177e-01, 2.7771616697639178e-01, 6.7493187009802774e-02}  // group D
	      };

	      for (unsigned int i=0, offset=0; i<nrows; ++i)
		{
		  _points[offset + 0] = Point(p[i][0],            p[i][1]); // (z1,z2)
		  _points[offset + 1] = Point(1.-p[i][0]-p[i][1], p[i][0]); // (z3,z1)
		  _points[offset + 2] = Point(p[i][1], 1.-p[i][0]-p[i][1]); // (z2,z3)

		  // All these points get the same weight
		  _weights[offset + 0] = p[i][2];
		  _weights[offset + 1] = p[i][2];
		  _weights[offset + 2] = p[i][2];

		  // Increment offset
		  offset += 3;
		}

	      return;
	      
	      // // This rule, which is originally due to:
	      // // Dunavant, "High degree efficient symmetrical Gaussian quadrature rules for
	      // // the triangle", IJNME 21 p. 1129--1148, 1985.
	      // //
	      // // It was copied 23rd June 2008 from:
	      // // http://people.scs.fsu.edu/~burkardt/f_src/dunavant/dunavant.f90
	      // //
	      // // This rule contains a negative weight and will fall through to the next case if you have
	      // // chosen to disallow such rules.
	      // if (allow_rules_with_negative_weights)
	      // 	{
	      // 	  _points.resize(13);
	      // 	  _weights.resize(13);
	      // 
	      // 	  // The raw data for the quadrature rule.
	      // 	  const Real p[4][4] = {
	      // 	    {                1./3.,                    0.,                    0., -0.149570044467682e+00 / 2.0}, // 1-perm
	      // 	    {0.479308067841920e+00, 0.260345966079040e+00,                    0., 0.175615257433208e+00  / 2.0}, // 3-perm
	      // 	    {0.869739794195568e+00, 0.065130102902216e+00,                    0., 0.053347235608838e+00  / 2.0}, // 3-perm
	      // 	    {0.048690315425316e+00, 0.312865496004874e+00, 0.638444188569810e+00, 0.077113760890257e+00  / 2.0}  // 6-perm
	      // 	  };
	      // 
	      // 
	      // 	  // Now call the dunavant routine to generate _points and _weights
	      // 	  dunavant_rule(p, 4);
	      // 	  
 	      // 	  return;
	      // 	} // end if (allow_rules_with_negative_weights)
	      // // Note: if !allow_rules_with_negative_weights, fall through to next case.
	    }



	    
	    // Another Dunavant rule.  This one has all positive weights.  This rule has
	    // 16 points while a comparable conical product rule would have 5*5=25.
	    //
	    // It was copied 23rd June 2008 from:
	    // http://people.scs.fsu.edu/~burkardt/f_src/dunavant/dunavant.f90
	  case EIGHTH:
	    {
	      _points.resize(16);
	      _weights.resize(16);

	      // The raw data for the quadrature rule.
	      const Real p[5][4] = {
		{                1./3.,                    0.,                    0., 0.144315607677787e+00 / 2.0}, // 1-perm
		{0.081414823414554e+00, 0.459292588292723e+00,                    0., 0.095091634267285e+00 / 2.0}, // 3-perm
		{0.658861384496480e+00, 0.170569307751760e+00,                    0., 0.103217370534718e+00 / 2.0}, // 3-perm
		{0.898905543365938e+00, 0.050547228317031e+00,                    0., 0.032458497623198e+00 / 2.0}, // 3-perm
		{0.008394777409958e+00, 0.263112829634638e+00, 0.728492392955404e+00, 0.027230314174435e+00 / 2.0}  // 6-perm
	      };

	      
	      // Now call the dunavant routine to generate _points and _weights
	      dunavant_rule(p, 5);

	      return;
	    }


	    
	    // Another Dunavant rule.  This one has all positive weights.  This rule has 19
	    // points. The comparable conical product rule would have 25.
	    // It was copied 23rd June 2008 from:
	    // http://people.scs.fsu.edu/~burkardt/f_src/dunavant/dunavant.f90
	  case NINTH:
	    {
	      _points.resize(19);
	      _weights.resize(19);


	      // The raw data for the quadrature rule.
	      const Real p[6][4] = {
		{                1./3.,                    0.,                    0., 0.097135796282799e+00 / 2.0}, // 1-perm
		{0.020634961602525e+00, 0.489682519198738e+00,                    0., 0.031334700227139e+00 / 2.0}, // 3-perm
		{0.125820817014127e+00, 0.437089591492937e+00,                    0., 0.077827541004774e+00 / 2.0}, // 3-perm
		{0.623592928761935e+00, 0.188203535619033e+00,                    0., 0.079647738927210e+00 / 2.0}, // 3-perm
		{0.910540973211095e+00, 0.044729513394453e+00,                    0., 0.025577675658698e+00 / 2.0}, // 3-perm 
		{0.036838412054736e+00, 0.221962989160766e+00, 0.741198598784498e+00, 0.043283539377289e+00 / 2.0}  // 6-perm
	      };

	      
	      // Now call the dunavant routine to generate _points and _weights
	      dunavant_rule(p, 6);

	      return;
	    }

	    
	    // Another Dunavant rule with all positive weights.  This rule has 25
	    // points. The comparable conical product rule would have 36.
	    // It was copied 23rd June 2008 from:
	    // http://people.scs.fsu.edu/~burkardt/f_src/dunavant/dunavant.f90
	  case TENTH:
	    {
	      _points.resize (25);
	      _weights.resize(25);

	      // The raw data for the quadrature rule.
	      const Real p[6][4] = {
		{                1./3.,                    0.,                    0., 0.090817990382754e+00 / 2.0}, // 1-perm
		{0.028844733232685e+00, 0.485577633383657e+00,                    0., 0.036725957756467e+00 / 2.0}, // 3-perm
		{0.781036849029926e+00, 0.109481575485037e+00,                    0., 0.045321059435528e+00 / 2.0}, // 3-perm
		{0.141707219414880e+00, 0.307939838764121e+00, 0.550352941820999e+00, 0.072757916845420e+00 / 2.0}, // 6-perm
		{0.025003534762686e+00, 0.246672560639903e+00, 0.728323904597411e+00, 0.028327242531057e+00 / 2.0}, // 6-perm 
		{0.009540815400299e+00, 0.066803251012200e+00, 0.923655933587500e+00, 0.009421666963733e+00 / 2.0}  // 6-perm
	      };

	      
	      // Now call the dunavant routine to generate _points and _weights
	      dunavant_rule(p, 6);

	      return;
	    }


	    
	    // Another Dunavant rule with all positive weights.  This rule has 33
	    // points. The comparable conical product rule would have 36 (ELEVENTH) or 49 (TWELFTH).
	    //
	    // Dunavant's 11th-order rule contains points outside the region of
	    // integration, and is thus unacceptable for our FEM calculations.
	    // 
	    // It was copied 23rd June 2008 from:
	    // http://people.scs.fsu.edu/~burkardt/f_src/dunavant/dunavant.f90
	  case ELEVENTH:     
	  case TWELFTH:
	    {
	      _points.resize (33);
	      _weights.resize(33);

	      // The raw data for the quadrature rule.
	      const Real p[8][4] = {
		{0.023565220452390e+00, 0.488217389773805e+00,                    0., 0.025731066440455e+00 / 2.0 }, // 3-perm
		{0.120551215411079e+00, 0.439724392294460e+00,                    0., 0.043692544538038e+00 / 2.0 }, // 3-perm
		{0.457579229975768e+00, 0.271210385012116e+00,                    0., 0.062858224217885e+00 / 2.0 }, // 3-perm
		{0.744847708916828e+00, 0.127576145541586e+00,                    0., 0.034796112930709e+00 / 2.0 }, // 3-perm
		{0.957365299093579e+00, 0.021317350453210e+00,                    0., 0.006166261051559e+00 / 2.0 }, // 3-perm
		{0.115343494534698e+00, 0.275713269685514e+00, 0.608943235779788e+00, 0.040371557766381e+00 / 2.0 }, // 6-perm
		{0.022838332222257e+00, 0.281325580989940e+00, 0.695836086787803e+00, 0.022356773202303e+00 / 2.0 }, // 6-perm 
		{0.025734050548330e+00, 0.116251915907597e+00, 0.858014033544073e+00, 0.017316231108659e+00 / 2.0 }  // 6-perm
	      };

	      
	      // Now call the dunavant routine to generate _points and _weights
	      dunavant_rule(p, 8);

	      return;
	    }

	    
	    // Another Dunavant rule with all positive weights.  This rule has 37
	    // points. The comparable conical product rule would have 49 points.
	    //
	    // It was copied 23rd June 2008 from:
	    // http://people.scs.fsu.edu/~burkardt/f_src/dunavant/dunavant.f90
	  case THIRTEENTH:
	    {
	      _points.resize (37);
	      _weights.resize(37);

	      // The raw data for the quadrature rule.
	      const Real p[10][4] = {
		{                1./3.,                    0.,                    0., 0.052520923400802e+00 / 2.0}, // 1-perm
		{0.009903630120591e+00, 0.495048184939705e+00,                    0., 0.011280145209330e+00 / 2.0}, // 3-perm
		{0.062566729780852e+00, 0.468716635109574e+00,                    0., 0.031423518362454e+00 / 2.0}, // 3-perm
		{0.170957326397447e+00, 0.414521336801277e+00,                    0., 0.047072502504194e+00 / 2.0}, // 3-perm
		{0.541200855914337e+00, 0.229399572042831e+00,                    0., 0.047363586536355e+00 / 2.0}, // 3-perm
		{0.771151009607340e+00, 0.114424495196330e+00,                    0., 0.031167529045794e+00 / 2.0}, // 3-perm
		{0.950377217273082e+00, 0.024811391363459e+00,                    0., 0.007975771465074e+00 / 2.0}, // 3-perm
		{0.094853828379579e+00, 0.268794997058761e+00, 0.636351174561660e+00, 0.036848402728732e+00 / 2.0}, // 6-perm
		{0.018100773278807e+00, 0.291730066734288e+00, 0.690169159986905e+00, 0.017401463303822e+00 / 2.0}, // 6-perm 
		{0.022233076674090e+00, 0.126357385491669e+00, 0.851409537834241e+00, 0.015521786839045e+00 / 2.0}  // 6-perm
	      };

	      
	      // Now call the dunavant routine to generate _points and _weights
	      dunavant_rule(p, 10);

	      return;
	    }

	    
	    // Another Dunavant rule.  This rule has 42 points, while
	    // a comparable conical product rule would have 64.
	    //
	    // It was copied 23rd June 2008 from:
	    // http://people.scs.fsu.edu/~burkardt/f_src/dunavant/dunavant.f90
	  case FOURTEENTH:
	    {
	      _points.resize (42);
	      _weights.resize(42);

	      // The raw data for the quadrature rule.
	      const Real p[10][4] = {
		{0.022072179275643e+00, 0.488963910362179e+00,                    0., 0.021883581369429e+00 / 2.0}, // 3-perm
		{0.164710561319092e+00, 0.417644719340454e+00,                    0., 0.032788353544125e+00 / 2.0}, // 3-perm
		{0.453044943382323e+00, 0.273477528308839e+00,                    0., 0.051774104507292e+00 / 2.0}, // 3-perm
		{0.645588935174913e+00, 0.177205532412543e+00,                    0., 0.042162588736993e+00 / 2.0}, // 3-perm
		{0.876400233818255e+00, 0.061799883090873e+00,                    0., 0.014433699669777e+00 / 2.0}, // 3-perm
		{0.961218077502598e+00, 0.019390961248701e+00,                    0., 0.004923403602400e+00 / 2.0}, // 3-perm
		{0.057124757403648e+00, 0.172266687821356e+00, 0.770608554774996e+00, 0.024665753212564e+00 / 2.0}, // 6-perm
		{0.092916249356972e+00, 0.336861459796345e+00, 0.570222290846683e+00, 0.038571510787061e+00 / 2.0}, // 6-perm
		{0.014646950055654e+00, 0.298372882136258e+00, 0.686980167808088e+00, 0.014436308113534e+00 / 2.0}, // 6-perm 
		{0.001268330932872e+00, 0.118974497696957e+00, 0.879757171370171e+00, 0.005010228838501e+00 / 2.0}  // 6-perm
	      };

	      
	      // Now call the dunavant routine to generate _points and _weights
	      dunavant_rule(p, 10);

	      return;
	    }

	    
	    // 15th-order rule by Wandzura.  
	    //
	    // Stephen Wandzura, Hong Xiao,
	    // Symmetric Quadrature Rules on a Triangle,
	    // Computers and Mathematics with Applications,
	    // Volume 45, Number 12, June 2003, pages 1829-1840.
	    //
	    // Wandzura's work extends the work of Dunavant by providing degree
	    // 5,10,15,20,25, and 30 rules with positive weights for the triangle.
	    //
	    // Copied on 3rd July 2008 from:
	    // http://people.scs.fsu.edu/~burkardt/f_src/wandzura/wandzura.f90
	  case FIFTEENTH:
	    {
	      _points.resize (54);
	      _weights.resize(54);

	      // The raw data for the quadrature rule.
	      const Real p[12][4] = {
		{0.08343840726175e+00, 0.45828079636912e+00,                  0.0, 0.3266181884880529e-01 / 2.0}, // 3-perm 
		{0.19277907084174e+00, 0.40361046457913e+00,                  0.0, 0.2741281803136436e-01 / 2.0}, // 3-perm 
		{0.41360566417395e+00, 0.29319716791303e+00,                  0.0, 0.2651003659870330e-01 / 2.0}, // 3-perm 
		{0.70706442611445e+00, 0.14646778694277e+00,                  0.0, 0.2921596213648611e-01 / 2.0}, // 3-perm 
		{0.88727426466879e+00, 0.05636286766560e+00,                  0.0, 0.1058460806624399e-01 / 2.0}, // 3-perm 
		{0.96684974628326e+00, 0.01657512685837e+00,                  0.0, 0.3614643064092035e-02 / 2.0}, // 3-perm 
		{0.00991220330923e+00, 0.23953455415479e+00, 0.75055324253598e+00, 0.8527748101709436e-02 / 2.0}, // 6-perm 
		{0.01580377063023e+00, 0.40487880731834e+00, 0.57931742205143e+00, 0.1391617651669193e-01 / 2.0}, // 6-perm 
		{0.00514360881697e+00, 0.09500211311304e+00, 0.89985427806998e+00, 0.4291932940734835e-02 / 2.0}, // 6-perm 
		{0.04892232575299e+00, 0.14975310732227e+00, 0.80132456692474e+00, 0.1623532928177489e-01 / 2.0}, // 6-perm 
		{0.06876874863252e+00, 0.28691961244133e+00, 0.64431163892615e+00, 0.2560734092126239e-01 / 2.0}, // 6-perm 
		{0.16840441812470e+00, 0.28183566809908e+00, 0.54975991377622e+00, 0.3308819553164567e-01 / 2.0}  // 6-perm 
	      };

	      
	      // Now call the dunavant routine to generate _points and _weights
	      dunavant_rule(p, 12);

	      return;
	    }
	    
	    // Dunavant's 15th and 16th-order rules contain points
	    // outside the domain of integration and thus are not
	    // suitable for FEM calculations.  Dunavant's 17th-order
	    // rule has 61 points, while a comparable conical product
	    // rule would have 64 (15th-order) or 81 (16th and 17th
	    // orders).
	    //
	    // It was copied 23rd June 2008 from:
	    // http://people.scs.fsu.edu/~burkardt/f_src/dunavant/dunavant.f90
	  case SIXTEENTH:    
	  case SEVENTEENTH:
	    {
	      _points.resize (61);
	      _weights.resize(61);

	      // The raw data for the quadrature rule.
	      const Real p[15][4] = {
		{                1./3.,                    0.,                    0., 0.033437199290803e+00 / 2.0}, // 1-perm
		{0.005658918886452e+00, 0.497170540556774e+00,                    0., 0.005093415440507e+00 / 2.0}, // 3-perm
		{0.035647354750751e+00, 0.482176322624625e+00,                    0., 0.014670864527638e+00 / 2.0}, // 3-perm
		{0.099520061958437e+00, 0.450239969020782e+00,                    0., 0.024350878353672e+00 / 2.0}, // 3-perm
		{0.199467521245206e+00, 0.400266239377397e+00,                    0., 0.031107550868969e+00 / 2.0}, // 3-perm
		{0.495717464058095e+00, 0.252141267970953e+00,                    0., 0.031257111218620e+00 / 2.0}, // 3-perm
		{0.675905990683077e+00, 0.162047004658461e+00,                    0., 0.024815654339665e+00 / 2.0}, // 3-perm
		{0.848248235478508e+00, 0.075875882260746e+00,                    0., 0.014056073070557e+00 / 2.0}, // 3-perm
		{0.968690546064356e+00, 0.015654726967822e+00,                    0., 0.003194676173779e+00 / 2.0}, // 3-perm
		{0.010186928826919e+00, 0.334319867363658e+00, 0.655493203809423e+00, 0.008119655318993e+00 / 2.0}, // 6-perm
		{0.135440871671036e+00, 0.292221537796944e+00, 0.572337590532020e+00, 0.026805742283163e+00 / 2.0}, // 6-perm
		{0.054423924290583e+00, 0.319574885423190e+00, 0.626001190286228e+00, 0.018459993210822e+00 / 2.0}, // 6-perm 
		{0.012868560833637e+00, 0.190704224192292e+00, 0.796427214974071e+00, 0.008476868534328e+00 / 2.0}, // 6-perm
		{0.067165782413524e+00, 0.180483211648746e+00, 0.752351005937729e+00, 0.018292796770025e+00 / 2.0}, // 6-perm 
		{0.014663182224828e+00, 0.080711313679564e+00, 0.904625504095608e+00, 0.006665632004165e+00 / 2.0}  // 6-perm
	      };

	      
	      // Now call the dunavant routine to generate _points and _weights
	      dunavant_rule(p, 15);

	      return;
	    }


	    
	    // Dunavant's 18th-order rule contains points outside the region and is therefore unsuitable
	    // for our FEM calculations.  His 19th-order rule has 73 points, compared with 100 points for
	    // a comparable-order conical product rule.
	    //
	    // It was copied 23rd June 2008 from:
	    // http://people.scs.fsu.edu/~burkardt/f_src/dunavant/dunavant.f90
	  case EIGHTTEENTH:  
	  case NINTEENTH:
	    {
	      _points.resize (73);
	      _weights.resize(73);

	      // The raw data for the quadrature rule.
	      const Real p[17][4] = {
		{                1./3.,                    0.,                    0., 0.032906331388919e+00 / 2.0}, // 1-perm
		{0.020780025853987e+00, 0.489609987073006e+00,                    0., 0.010330731891272e+00 / 2.0}, // 3-perm
		{0.090926214604215e+00, 0.454536892697893e+00,                    0., 0.022387247263016e+00 / 2.0}, // 3-perm
		{0.197166638701138e+00, 0.401416680649431e+00,                    0., 0.030266125869468e+00 / 2.0}, // 3-perm
		{0.488896691193805e+00, 0.255551654403098e+00,                    0., 0.030490967802198e+00 / 2.0}, // 3-perm
		{0.645844115695741e+00, 0.177077942152130e+00,                    0., 0.024159212741641e+00 / 2.0}, // 3-perm
		{0.779877893544096e+00, 0.110061053227952e+00,                    0., 0.016050803586801e+00 / 2.0}, // 3-perm
		{0.888942751496321e+00, 0.055528624251840e+00,                    0., 0.008084580261784e+00 / 2.0}, // 3-perm
		{0.974756272445543e+00, 0.012621863777229e+00,                    0., 0.002079362027485e+00 / 2.0}, // 3-perm
		{0.003611417848412e+00, 0.395754787356943e+00, 0.600633794794645e+00, 0.003884876904981e+00 / 2.0}, // 6-perm
		{0.134466754530780e+00, 0.307929983880436e+00, 0.557603261588784e+00, 0.025574160612022e+00 / 2.0}, // 6-perm
		{0.014446025776115e+00, 0.264566948406520e+00, 0.720987025817365e+00, 0.008880903573338e+00 / 2.0}, // 6-perm 
		{0.046933578838178e+00, 0.358539352205951e+00, 0.594527068955871e+00, 0.016124546761731e+00 / 2.0}, // 6-perm
		{0.002861120350567e+00, 0.157807405968595e+00, 0.839331473680839e+00, 0.002491941817491e+00 / 2.0}, // 6-perm 
		{0.223861424097916e+00, 0.075050596975911e+00, 0.701087978926173e+00, 0.018242840118951e+00 / 2.0}, // 6-perm
		{0.034647074816760e+00, 0.142421601113383e+00, 0.822931324069857e+00, 0.010258563736199e+00 / 2.0}, // 6-perm 
		{0.010161119296278e+00, 0.065494628082938e+00, 0.924344252620784e+00, 0.003799928855302e+00 / 2.0}  // 6-perm
	      };

	      
	      // Now call the dunavant routine to generate _points and _weights
	      dunavant_rule(p, 17);

	      return;
	    }

	    
	    // 20th-order rule by Wandzura.  
	    //
	    // Stephen Wandzura, Hong Xiao,
	    // Symmetric Quadrature Rules on a Triangle,
	    // Computers and Mathematics with Applications,
	    // Volume 45, Number 12, June 2003, pages 1829-1840.
	    //
	    // Wandzura's work extends the work of Dunavant by providing degree
	    // 5,10,15,20,25, and 30 rules with positive weights for the triangle.
	    //
	    // Copied on 3rd July 2008 from:
	    // http://people.scs.fsu.edu/~burkardt/f_src/wandzura/wandzura.f90
	  case TWENTIETH:
	    {
	      // The equivalent concial product rule would have 121 points
	      _points.resize (85);
	      _weights.resize(85);

	      // The raw data for the quadrature rule.
	      const Real p[19][4] = {
		{0.33333333333333e+00,                  0.0,                  0.0, 0.2761042699769952e-01 / 2.0}, // 1-perm
		{0.00150064932443e+00, 0.49924967533779e+00,                  0.0, 0.1779029547326740e-02 / 2.0}, // 3-perm 
		{0.09413975193895e+00, 0.45293012403052e+00,                  0.0, 0.2011239811396117e-01 / 2.0}, // 3-perm
		{0.20447212408953e+00, 0.39776393795524e+00,                  0.0, 0.2681784725933157e-01 / 2.0}, // 3-perm
		{0.47099959493443e+00, 0.26450020253279e+00,                  0.0, 0.2452313380150201e-01 / 2.0}, // 3-perm
		{0.57796207181585e+00, 0.21101896409208e+00,                  0.0, 0.1639457841069539e-01 / 2.0}, // 3-perm
		{0.78452878565746e+00, 0.10773560717127e+00,                  0.0, 0.1479590739864960e-01 / 2.0}, // 3-perm
		{0.92186182432439e+00, 0.03906908783780e+00,                  0.0, 0.4579282277704251e-02 / 2.0}, // 3-perm
		{0.97765124054134e+00, 0.01117437972933e+00,                  0.0, 0.1651826515576217e-02 / 2.0}, // 3-perm
		{0.00534961818734e+00, 0.06354966590835e+00, 0.93110071590431e+00, 0.2349170908575584e-02 / 2.0}, // 6-perm
		{0.00795481706620e+00, 0.15710691894071e+00, 0.83493826399309e+00, 0.4465925754181793e-02 / 2.0}, // 6-perm
		{0.01042239828126e+00, 0.39564211436437e+00, 0.59393548735436e+00, 0.6099566807907972e-02 / 2.0}, // 6-perm
		{0.01096441479612e+00, 0.27316757071291e+00, 0.71586801449097e+00, 0.6891081327188203e-02 / 2.0}, // 6-perm
		{0.03856671208546e+00, 0.10178538248502e+00, 0.85964790542952e+00, 0.7997475072478163e-02 / 2.0}, // 6-perm
		{0.03558050781722e+00, 0.44665854917641e+00, 0.51776094300637e+00, 0.7386134285336024e-02 / 2.0}, // 6-perm
		{0.04967081636276e+00, 0.19901079414950e+00, 0.75131838948773e+00, 0.1279933187864826e-01 / 2.0}, // 6-perm
		{0.05851972508433e+00, 0.32426118369228e+00, 0.61721909122339e+00, 0.1725807117569655e-01 / 2.0}, // 6-perm
		{0.12149778700439e+00, 0.20853136321013e+00, 0.66997084978547e+00, 0.1867294590293547e-01 / 2.0}, // 6-perm
		{0.14071084494394e+00, 0.32317056653626e+00, 0.53611858851980e+00, 0.2281822405839526e-01 / 2.0}  // 6-perm
	      };

	      
	      // Now call the dunavant routine to generate _points and _weights
	      dunavant_rule(p, 19);

	      return;
	    }



	    // 25th-order rule by Wandzura.  
	    //
	    // Stephen Wandzura, Hong Xiao,
	    // Symmetric Quadrature Rules on a Triangle,
	    // Computers and Mathematics with Applications,
	    // Volume 45, Number 12, June 2003, pages 1829-1840.
	    //
	    // Wandzura's work extends the work of Dunavant by providing degree
	    // 5,10,15,20,25, and 30 rules with positive weights for the triangle.
	    //
	    // Copied on 3rd July 2008 from:
	    // http://people.scs.fsu.edu/~burkardt/f_src/wandzura/wandzura.f90
	    // case TWENTYFIRST: // fall through to 121 point conical product rule below
	  case TWENTYSECOND:  
	  case TWENTYTHIRD:   
	  case TWENTYFOURTH:  	    
	  case TWENTYFIFTH:
	    {
	      // The equivalent concial product rule would have 169 points
	      _points.resize (126);
	      _weights.resize(126);

	      // The raw data for the quadrature rule.
	      const Real p[26][4] = {
		{0.02794648307317e+00, 0.48602675846341e+00,                  0.0, 0.8005581880020417e-02 / 2.0},  // 3-perm
		{0.13117860132765e+00, 0.43441069933617e+00,                  0.0, 0.1594707683239050e-01 / 2.0},  // 3-perm
		{0.22022172951207e+00, 0.38988913524396e+00,                  0.0, 0.1310914123079553e-01 / 2.0},  // 3-perm
		{0.40311353196039e+00, 0.29844323401980e+00,                  0.0, 0.1958300096563562e-01 / 2.0},  // 3-perm
		{0.53191165532526e+00, 0.23404417233737e+00,                  0.0, 0.1647088544153727e-01 / 2.0},  // 3-perm
		{0.69706333078196e+00, 0.15146833460902e+00,                  0.0, 0.8547279074092100e-02 / 2.0},  // 3-perm
		{0.77453221290801e+00, 0.11273389354599e+00,                  0.0, 0.8161885857226492e-02 / 2.0},  // 3-perm
		{0.84456861581695e+00, 0.07771569209153e+00,                  0.0, 0.6121146539983779e-02 / 2.0},  // 3-perm
		{0.93021381277141e+00, 0.03489309361430e+00,                  0.0, 0.2908498264936665e-02 / 2.0},  // 3-perm
		{0.98548363075813e+00, 0.00725818462093e+00,                  0.0, 0.6922752456619963e-03 / 2.0},  // 3-perm
		{0.00129235270444e+00, 0.22721445215336e+00, 0.77149319514219e+00, 0.1248289199277397e-02 / 2.0},  // 6-perm
		{0.00539970127212e+00, 0.43501055485357e+00, 0.55958974387431e+00, 0.3404752908803022e-02 / 2.0},  // 6-perm
		{0.00638400303398e+00, 0.32030959927220e+00, 0.67330639769382e+00, 0.3359654326064051e-02 / 2.0},  // 6-perm
		{0.00502821150199e+00, 0.09175032228001e+00, 0.90322146621800e+00, 0.1716156539496754e-02 / 2.0},  // 6-perm
		{0.00682675862178e+00, 0.03801083585872e+00, 0.95516240551949e+00, 0.1480856316715606e-02 / 2.0},  // 6-perm
		{0.01001619963993e+00, 0.15742521848531e+00, 0.83255858187476e+00, 0.3511312610728685e-02 / 2.0},  // 6-perm
		{0.02575781317339e+00, 0.23988965977853e+00, 0.73435252704808e+00, 0.7393550149706484e-02 / 2.0},  // 6-perm
		{0.03022789811992e+00, 0.36194311812606e+00, 0.60782898375402e+00, 0.7983087477376558e-02 / 2.0},  // 6-perm
		{0.03050499010716e+00, 0.08355196095483e+00, 0.88594304893801e+00, 0.4355962613158041e-02 / 2.0},  // 6-perm
		{0.04595654736257e+00, 0.14844322073242e+00, 0.80560023190501e+00, 0.7365056701417832e-02 / 2.0},  // 6-perm
		{0.06744280054028e+00, 0.28373970872753e+00, 0.64881749073219e+00, 0.1096357284641955e-01 / 2.0},  // 6-perm
		{0.07004509141591e+00, 0.40689937511879e+00, 0.52305553346530e+00, 0.1174996174354112e-01 / 2.0},  // 6-perm
		{0.08391152464012e+00, 0.19411398702489e+00, 0.72197448833499e+00, 0.1001560071379857e-01 / 2.0},  // 6-perm
		{0.12037553567715e+00, 0.32413434700070e+00, 0.55549011732214e+00, 0.1330964078762868e-01 / 2.0},  // 6-perm
		{0.14806689915737e+00, 0.22927748355598e+00, 0.62265561728665e+00, 0.1415444650522614e-01 / 2.0},  // 6-perm
		{0.19177186586733e+00, 0.32561812259598e+00, 0.48261001153669e+00, 0.1488137956116801e-01 / 2.0}   // 6-perm
	      };

	      
	      // Now call the dunavant routine to generate _points and _weights
	      dunavant_rule(p, 26);

	      return;
	    }



	    // 30th-order rule by Wandzura.  
	    //
	    // Stephen Wandzura, Hong Xiao,
	    // Symmetric Quadrature Rules on a Triangle,
	    // Computers and Mathematics with Applications,
	    // Volume 45, Number 12, June 2003, pages 1829-1840.
	    //
	    // Wandzura's work extends the work of Dunavant by providing degree
	    // 5,10,15,20,25, and 30 rules with positive weights for the triangle.
	    //
	    // Copied on 3rd July 2008 from:
	    // http://people.scs.fsu.edu/~burkardt/f_src/wandzura/wandzura.f90
	  case TWENTYSIXTH:   
	  case TWENTYSEVENTH: 
	  case TWENTYEIGHTH:  
	  case TWENTYNINTH:   	    
	  case THIRTIETH:
	    {
	      // The equivalent concial product rule would have 256 points
	      _points.resize (175);
	      _weights.resize(175);

	      // The raw data for the quadrature rule.
	      const Real p[36][4] = {
		{0.33333333333333e+00,                  0.0,                  0.0, 0.1557996020289920e-01 / 2.0}, // 1-perm 
		{0.00733011643277e+00, 0.49633494178362e+00,                  0.0, 0.3177233700534134e-02 / 2.0}, // 3-perm  
		{0.08299567580296e+00, 0.45850216209852e+00,                  0.0, 0.1048342663573077e-01 / 2.0}, // 3-perm  
		{0.15098095612541e+00, 0.42450952193729e+00,                  0.0, 0.1320945957774363e-01 / 2.0}, // 3-perm  
		{0.23590585989217e+00, 0.38204707005392e+00,                  0.0, 0.1497500696627150e-01 / 2.0}, // 3-perm  
		{0.43802430840785e+00, 0.28098784579608e+00,                  0.0, 0.1498790444338419e-01 / 2.0}, // 3-perm  
		{0.54530204829193e+00, 0.22734897585403e+00,                  0.0, 0.1333886474102166e-01 / 2.0}, // 3-perm  
		{0.65088177698254e+00, 0.17455911150873e+00,                  0.0, 0.1088917111390201e-01 / 2.0}, // 3-perm  
		{0.75348314559713e+00, 0.12325842720144e+00,                  0.0, 0.8189440660893461e-02 / 2.0}, // 3-perm  
		{0.83983154221561e+00, 0.08008422889220e+00,                  0.0, 0.5575387588607785e-02 / 2.0}, // 3-perm  
		{0.90445106518420e+00, 0.04777446740790e+00,                  0.0, 0.3191216473411976e-02 / 2.0}, // 3-perm  
		{0.95655897063972e+00, 0.02172051468014e+00,                  0.0, 0.1296715144327045e-02 / 2.0}, // 3-perm  
		{0.99047064476913e+00, 0.00476467761544e+00,                  0.0, 0.2982628261349172e-03 / 2.0}, // 3-perm  
		{0.00092537119335e+00, 0.41529527091331e+00, 0.58377935789334e+00, 0.9989056850788964e-03 / 2.0}, // 6-perm  
		{0.00138592585556e+00, 0.06118990978535e+00, 0.93742416435909e+00, 0.4628508491732533e-03 / 2.0}, // 6-perm  
		{0.00368241545591e+00, 0.16490869013691e+00, 0.83140889440718e+00, 0.1234451336382413e-02 / 2.0}, // 6-perm  
		{0.00390322342416e+00, 0.02503506223200e+00, 0.97106171434384e+00, 0.5707198522432062e-03 / 2.0}, // 6-perm  
		{0.00323324815501e+00, 0.30606446515110e+00, 0.69070228669389e+00, 0.1126946125877624e-02 / 2.0}, // 6-perm  
		{0.00646743211224e+00, 0.10707328373022e+00, 0.88645928415754e+00, 0.1747866949407337e-02 / 2.0}, // 6-perm  
		{0.00324747549133e+00, 0.22995754934558e+00, 0.76679497516308e+00, 0.1182818815031657e-02 / 2.0}, // 6-perm  
		{0.00867509080675e+00, 0.33703663330578e+00, 0.65428827588746e+00, 0.1990839294675034e-02 / 2.0}, // 6-perm  
		{0.01559702646731e+00, 0.05625657618206e+00, 0.92814639735063e+00, 0.1900412795035980e-02 / 2.0}, // 6-perm  
		{0.01797672125369e+00, 0.40245137521240e+00, 0.57957190353391e+00, 0.4498365808817451e-02 / 2.0}, // 6-perm  
		{0.01712424535389e+00, 0.24365470201083e+00, 0.73922105263528e+00, 0.3478719460274719e-02 / 2.0}, // 6-perm  
		{0.02288340534658e+00, 0.16538958561453e+00, 0.81172700903888e+00, 0.4102399036723953e-02 / 2.0}, // 6-perm  
		{0.03273759728777e+00, 0.09930187449585e+00, 0.86796052821639e+00, 0.4021761549744162e-02 / 2.0}, // 6-perm  
		{0.03382101234234e+00, 0.30847833306905e+00, 0.65770065458860e+00, 0.6033164660795066e-02 / 2.0}, // 6-perm  
		{0.03554761446002e+00, 0.46066831859211e+00, 0.50378406694787e+00, 0.3946290302129598e-02 / 2.0}, // 6-perm  
		{0.05053979030687e+00, 0.21881529945393e+00, 0.73064491023920e+00, 0.6644044537680268e-02 / 2.0}, // 6-perm  
		{0.05701471491573e+00, 0.37920955156027e+00, 0.56377573352399e+00, 0.8254305856078458e-02 / 2.0}, // 6-perm  
		{0.06415280642120e+00, 0.14296081941819e+00, 0.79288637416061e+00, 0.6496056633406411e-02 / 2.0}, // 6-perm  
		{0.08050114828763e+00, 0.28373128210592e+00, 0.63576756960645e+00, 0.9252778144146602e-02 / 2.0}, // 6-perm  
		{0.10436706813453e+00, 0.19673744100444e+00, 0.69889549086103e+00, 0.9164920726294280e-02 / 2.0}, // 6-perm  
		{0.11384489442875e+00, 0.35588914121166e+00, 0.53026596435959e+00, 0.1156952462809767e-01 / 2.0}, // 6-perm  
		{0.14536348771552e+00, 0.25981868535191e+00, 0.59481782693256e+00, 0.1176111646760917e-01 / 2.0}, // 6-perm  
		{0.18994565282198e+00, 0.32192318123130e+00, 0.48813116594672e+00, 0.1382470218216540e-01 / 2.0}  // 6-perm   
	      };

	      
	      // Now call the dunavant routine to generate _points and _weights
	      dunavant_rule(p, 36);

	      return;
	    }
	    
	    
	    // By default, we fall back on the conical product rules.  If the user
	    // requests an order higher than what is currently available in the 1D
	    // rules, an error will be thrown from the respective 1D code.
	  default:
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
              // FIXME - why can't we init() these explicitly? [RHS]
	      QGauss  gauss1D(1,static_cast<Order>(_order+2*p));
	      QJacobi jac1D(1,static_cast<Order>(_order+2*p),1,0);
	      
	      // The Gauss rule needs to be scaled to [0,1]
	      std::pair<Real, Real> old_range(-1,1);
	      std::pair<Real, Real> new_range(0,1);
	      gauss1D.scale(old_range,
			    new_range);

	      // Compute the tensor product
	      tensor_product_tri(gauss1D, jac1D);
	      return;
	    }
	    
	  }
      }

	    
      //---------------------------------------------
      // Unsupported type
    default:
      {
	std::cerr << "Element type not supported!:" << _type << std::endl;
	libmesh_error();
      }
    }

  libmesh_error();

  return;

#endif
}

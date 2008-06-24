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

	  case SIXTH:
	    {
	      // Exact for sixth degree polynomials
	      // Taken from http://www.cs.rpi.edu/~flaherje/FEM/fem6.ps
	      // by Flaherty
	      _points.resize(12);
	      _weights.resize(12);
	      const Real w1 = 0.050844906370207 / 2.0;
	      const Real w2 = 0.116786275726379 / 2.0;
	      const Real w3 = 0.082851075618374 / 2.0;
	      const Real a1 = 0.873821971016996;
	      const Real a2 = 0.063089014491502;
	      const Real b1 = 0.501426509658179;
	      const Real b2 = 0.249286745170910;
	      const Real c1 = 0.636502499121399;
	      const Real c2 = 0.310352451033785;
	      const Real c3 = 0.053145049844816;

	      _points[0](0) = a1;
	      _points[0](1) = a2;

	      _points[1](0) = a2;
	      _points[1](1) = a1;

	      _points[2](0) = a2;
	      _points[2](1) = a2;

	      _points[3](0) = b1;
	      _points[3](1) = b2;

	      _points[4](0) = b2;
	      _points[4](1) = b1;

	      _points[5](0) = b2;
	      _points[5](1) = b2;

	      _points[6](0) = c1;
	      _points[6](1) = c2;

	      _points[7](0) = c1;
	      _points[7](1) = c3;

	      _points[8](0) = c2;
	      _points[8](1) = c1;

	      _points[9](0) = c2;
	      _points[9](1) = c3;

	      _points[10](0) = c3;
	      _points[10](1) = c1;

	      _points[11](0) = c3;
	      _points[11](1) = c2;

	      _weights[0] = w1;
	      _weights[1] = w1;
	      _weights[2] = w1;
	      _weights[3] = w2;
	      _weights[4] = w2;
	      _weights[5] = w2;
	      _weights[6] = w3;
	      _weights[7] = w3;
	      _weights[8] = w3;
	      _weights[9] = w3;
	      _weights[10] = w3;
	      _weights[11] = w3;

	      return;
	    }

	    
	    // This rule, which is originally due to:
	    // Dunavant, "High degree efficient symmetrical Gaussian quadrature rules for
	    // the triangle", IJNME 21 p. 1129--1148, 1985.
	    //
	    // It was copied 23rd June 2008 from:
	    // http://people.scs.fsu.edu/~burkardt/f_src/dunavant/dunavant.f90
	    //
	    // This rule contains a negative weight and will fall through to EIGHTH if you have
	    // chosen to disallow such rules.
	  case SEVENTH:
	    {
	      if (allow_rules_with_negative_weights)
		{
		  _points.resize(13);
		  _weights.resize(13);

		  // The raw data for the quadrature rule.
		  const Real p[4][4] = {
		    {                1./3.,                    0.,                    0., -0.149570044467682e+00 / 2.0}, // 1-perm
		    {0.479308067841920e+00, 0.260345966079040e+00,                    0., 0.175615257433208e+00  / 2.0}, // 3-perm
		    {0.869739794195568e+00, 0.065130102902216e+00,                    0., 0.053347235608838e+00  / 2.0}, // 3-perm
		    {0.048690315425316e+00, 0.312865496004874e+00, 0.638444188569810e+00, 0.077113760890257e+00  / 2.0}  // 6-perm
		  };

	      
		  // Now call the dunavant routine to generate _points and _weights
		  dunavant_rule(p, 4);
		  
 		  return;
		} // end if (allow_rules_with_negative_weights)
	      // Note: if !allow_rules_with_negative_weights, fall through to next case.
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

	    
	    // Dunavant's 15th and 16th-order rules contain points
	    // outside the domain of integration and thus are not
	    // suitable for FEM calculations.  Dunavant's 17th-order
	    // rule has 61 points, while a comparable conical product
	    // rule would have 64 (15th-order) or 81 (16th and 17th
	    // orders).
	    //
	    // It was copied 23rd June 2008 from:
	    // http://people.scs.fsu.edu/~burkardt/f_src/dunavant/dunavant.f90
	  case FIFTEENTH:    
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
	    
	  default:
	    {
	      std::cout << "Quadrature rule not supported!" << std::endl;

	      libmesh_error();
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

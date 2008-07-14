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
#include "quadrature_monomial.h"
#include "quadrature_gauss.h"


void QMonomial::init_2D(const ElemType _type,
			unsigned int p)
{

  switch (_type)
    {
      //---------------------------------------------
      // Quadrilateral quadrature rules
    case QUAD4:
    case QUAD8:
    case QUAD9:
      {
	switch(_order + 2*p)
	  {
	  case SECOND:
	    {
	      // A degree=2 rule for the QUAD with 3 points.
	      // A tensor product degree-2 Gauss would have 4 points.
	      // This rule (or a variation on it) is probably available in
	      //
	      // A.H. Stroud, Approximate calculation of multiple integrals,
	      // Prentice-Hall, Englewood Cliffs, N.J., 1971.
	      //
	      // though I have never actually seen a reference for it.
	      // Luckily it's fairly easy to derive, which is what I've done
	      // here [JWP].
	      const Real
		s=std::sqrt(1./3.),
		t=std::sqrt(2./3.);
	      
	      const Real data[2][3] =
		{
		  {0.0,  s,  2.0},
		  {  t, -s,  1.0}
		};

	      _points.resize(3);
	      _weights.resize(3);

	      wissmann_rule(data, 2);
	      
	      return;
	    } // end case SECOND


	    
	  case FOURTH:
	    {
	      // A pair of degree=4 rules for the QUAD "C2" due to
	      // Wissmann and Becker. These rules both have six points.
	      // A tensor product degree-4 Gauss would have 9 points.
	      //
	      // J. W. Wissmann and T. Becker, Partially symmetric cubature 
	      // formulas for even degrees of exactness, SIAM J. Numer. Anal.  23
	      // (1986), 676--685.                                           
	      const Real data[4][3] =
		{
		  // First of 2 degree-4 rules given by Wissmann
		  {0.0000000000000000e+00,  0.0000000000000000e+00,  1.1428571428571428e+00},
		  {0.0000000000000000e+00,  9.6609178307929590e-01,  4.3956043956043956e-01},
		  {8.5191465330460049e-01,  4.5560372783619284e-01,  5.6607220700753210e-01},
		  {6.3091278897675402e-01, -7.3162995157313452e-01,  6.4271900178367668e-01}
		  //
		  // Second of 2 degree-4 rules given by Wissmann.  These both
		  // yield 4th-order accurate rules, I just chose the one that
		  // happened to contain the origin.
		  // {0.000000000000000, -0.356822089773090,  1.286412084888852}, 
		  // {0.000000000000000,  0.934172358962716,  0.491365692888926},
		  // {0.774596669241483,  0.390885162530071,  0.761883709085613},
		  // {0.774596669241483, -0.852765377881771,  0.349227402025498} 
		};

	      _points.resize(6);
	      _weights.resize(6);

	      wissmann_rule(data, 4);
	      
	      return;
	    } // end case FOURTH



	    
	  case SIXTH:
	    {
	      // A pair of degree=6 rules for the QUAD "C2" due to
	      // Wissmann and Becker. These rules both have 10 points.
	      // A tensor product degree-6 Gauss would have 16 points.
	      //
	      // J. W. Wissmann and T. Becker, Partially symmetric cubature 
	      // formulas for even degrees of exactness, SIAM J. Numer. Anal.  23
	      // (1986), 676--685.                                           
	      const Real data[6][3] =
		{
		  // First of 2 degree-6, 10 point rules given by Wissmann
		  // {0.000000000000000,  0.836405633697626,  0.455343245714174},
		  // {0.000000000000000, -0.357460165391307,  0.827395973202966},
		  // {0.888764014654765,  0.872101531193131,  0.144000884599645},
		  // {0.604857639464685,  0.305985162155427,  0.668259104262665},
		  // {0.955447506641064, -0.410270899466658,  0.225474004890679},
		  // {0.565459993438754, -0.872869311156879,  0.320896396788441} 
		  //
		  // Second of 2 degree-6, 10 point rules given by Wissmann.
		  // Either of these will work, I just chose the one with points
		  // slightly further into the element interior.
		  {0.0000000000000000e+00,  8.6983337525005900e-01,  3.9275059096434794e-01},
		  {0.0000000000000000e+00, -4.7940635161211124e-01,  7.5476288124261053e-01},
		  {8.6374282634615388e-01,  8.0283751620765670e-01,  2.0616605058827902e-01},
		  {5.1869052139258234e-01,  2.6214366550805818e-01,  6.8999213848986375e-01},
		  {9.3397254497284950e-01, -3.6309658314806653e-01,  2.6051748873231697e-01},
		  {6.0897753601635630e-01, -8.9660863276245265e-01,  2.6956758608606100e-01}
		};

	      _points.resize(10);
	      _weights.resize(10);

	      wissmann_rule(data, 6);

	      return; 
	    } // end case SIXTH

	  case EIGHTH:
	    {
	      // A pair of degree=8 rules for the QUAD "C2" due to
	      // Wissmann and Becker. These rules both have 16 points.
	      // A tensor product degree-6 Gauss would have 25 points.
	      //
	      // J. W. Wissmann and T. Becker, Partially symmetric cubature 
	      // formulas for even degrees of exactness, SIAM J. Numer. Anal.  23
	      // (1986), 676--685.                                           
	      const Real data[10][3] =
		{
		  // First of 2 degree-8, 16 point rules given by Wissmann
		  // {0.000000000000000,  0.000000000000000,  0.055364705621440},
		  // {0.000000000000000,  0.757629177660505,  0.404389368726076},
		  // {0.000000000000000, -0.236871842255702,  0.533546604952635},
		  // {0.000000000000000, -0.989717929044527,  0.117054188786739},
		  // {0.639091304900370,  0.950520955645667,  0.125614417613747},
		  // {0.937069076924990,  0.663882736885633,  0.136544584733588},
		  // {0.537083530541494,  0.304210681724104,  0.483408479211257},
		  // {0.887188506449625, -0.236496718536120,  0.252528506429544},
		  // {0.494698820670197, -0.698953476086564,  0.361262323882172},
		  // {0.897495818279768, -0.900390774211580,  0.085464254086247}
		  //
		  // Second of 2 degree-8, 16 point rules given by Wissmann.
		  // Either of these will work, I just chose the one with points
		  // further into the element interior.
		  {0.0000000000000000e+00,  6.5956013196034176e-01,  4.5027677630559029e-01},
		  {0.0000000000000000e+00, -9.4914292304312538e-01,  1.6657042677781274e-01},
		  {9.5250946607156228e-01,  7.6505181955768362e-01,  9.8869459933431422e-02},
		  {5.3232745407420624e-01,  9.3697598108841598e-01,  1.5369674714081197e-01},
		  {6.8473629795173504e-01,  3.3365671773574759e-01,  3.9668697607290278e-01},
		  {2.3314324080140552e-01, -7.9583272377396852e-02,  3.5201436794569501e-01},
		  {9.2768331930611748e-01, -2.7224008061253425e-01,  1.8958905457779799e-01},
		  {4.5312068740374942e-01, -6.1373535339802760e-01,  3.7510100114758727e-01},
		  {8.3750364042281223e-01, -8.8847765053597136e-01,  1.2561879164007201e-01}
		};

	      _points.resize(16);
	      _weights.resize(16);

	      wissmann_rule(data, /*10*/ 9);

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
      } // end case QUAD4/8/9

      
      // By default: construct and use a Gauss quadrature rule
    default:
      {
	QGauss gauss_rule(2, _order);
	gauss_rule.init(_type, p);

	// Swap points and weights with the about-to-be destroyed rule.
	_points.swap (gauss_rule.get_points() );
	_weights.swap(gauss_rule.get_weights());

	return;
      }
    } // end switch (_type)
}

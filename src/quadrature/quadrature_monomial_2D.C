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


void QMonomial::init_2D(const ElemType type_in,
			unsigned int p)
{

  switch (type_in)
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



	  // For third-order, fall through to default case, use 2x2 Gauss product rule.
	  // case THIRD:
	  //   {
	  //   }  // end case THIRD

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




	  case FIFTH:
	    {
	      // A degree 5, 7-point rule due to Stroud.
	      //
	      // A.H. Stroud, Approximate calculation of multiple integrals,
	      // Prentice-Hall, Englewood Cliffs, N.J., 1971.
	      //
	      // This rule is provably minimal in the number of points.
	      // A tensor-product rule accurate for "bi-quintic" polynomials would have 9 points.
	      const Real data[3][3] =
		{
		  {                                  0.L,                                     0.L, static_cast<Real>(8.L  /  7.L)}, // 1
		  {                                  0.L, static_cast<Real>(std::sqrt(14.L/15.L)), static_cast<Real>(20.L / 63.L)}, // 2
		  {static_cast<Real>(std::sqrt(3.L/5.L)),   static_cast<Real>(std::sqrt(1.L/3.L)), static_cast<Real>(20.L / 36.L)}  // 4
		};

	      const unsigned int symmetry[3] = {
		0, // Origin
		7, // Central Symmetry
		6  // Rectangular
	      };

	      _points.resize (7);
	      _weights.resize(7);

	      stroud_rule(data, symmetry, 3);

	      return;
	    } // end case FIFTH




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




	  case SEVENTH:
	    {
	      // A degree 7, 12-point rule due to Tyler, can be found in Stroud's book
	      //
	      // A.H. Stroud, Approximate calculation of multiple integrals,
	      // Prentice-Hall, Englewood Cliffs, N.J., 1971.
	      //
	      // This rule is fully-symmetric and provably minimal in the number of points.
	      // A tensor-product rule accurate for "bi-septic" polynomials would have 16 points.
	      const Real
		r  = std::sqrt(6.L/7.L),
		s  = std::sqrt( (114.L - 3.L*std::sqrt(583.L)) / 287.L ),
		t  = std::sqrt( (114.L + 3.L*std::sqrt(583.L)) / 287.L ),
		B1 = 196.L / 810.L,
		B2 = 4.L * (178981.L + 2769.L*std::sqrt(583.L)) / 1888920.L,
		B3 = 4.L * (178981.L - 2769.L*std::sqrt(583.L)) / 1888920.L;

	      const Real data[3][3] =
		{
		  {r, 0.0, B1}, // 4
		  {s, 0.0, B2}, // 4
		  {t, 0.0, B3}  // 4
		};

	      const unsigned int symmetry[3] = {
		3, // Full Symmetry, (x,0)
		2, // Full Symmetry, (x,x)
		2  // Full Symmetry, (x,x)
	      };

	      _points.resize (12);
	      _weights.resize(12);

	      stroud_rule(data, symmetry, 3);

	      return;
	    } // end case SEVENTH




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




	  case NINTH:
	    {
	      // A degree 9, 17-point rule due to Moller.
	      //
	      // H.M. Moller,  Kubaturformeln mit minimaler Knotenzahl,
	      // Numer. Math.  25 (1976), 185--200.
	      //
	      // This rule is provably minimal in the number of points.
	      // A tensor-product rule accurate for "bi-ninth" degree polynomials would have 25 points.
	      const Real data[5][3] =
		{
		  {0.0000000000000000e+00, 0.0000000000000000e+00, 5.2674897119341563e-01}, // 1
		  {6.3068011973166885e-01, 9.6884996636197772e-01, 8.8879378170198706e-02}, // 4
		  {9.2796164595956966e-01, 7.5027709997890053e-01, 1.1209960212959648e-01}, // 4
		  {4.5333982113564719e-01, 5.2373582021442933e-01, 3.9828243926207009e-01}, // 4
		  {8.5261572933366230e-01, 7.6208328192617173e-02, 2.6905133763978080e-01}  // 4
		};

	      const unsigned int symmetry[5] = {
		0, // Single point
		4, // Rotational Invariant
		4, // Rotational Invariant
		4, // Rotational Invariant
		4  // Rotational Invariant
	      };

	      _points.resize (17);
	      _weights.resize(17);

	      stroud_rule(data, symmetry, 5);

	      return;
	    } // end case NINTH




	  case TENTH:
	  case ELEVENTH:
	    {
	      // A degree 11, 24-point rule due to Cools and Haegemans.
	      //
	      // R. Cools and A. Haegemans, Another step forward in searching for
	      // cubature formulae with a minimal number of knots for the square,
	      // Computing 40 (1988), 139--146.
	      //
	      // P. Verlinden and R. Cools, The algebraic construction of a minimal
	      // cubature formula of degree 11 for the square, Cubature Formulas
	      // and their Applications (Russian) (Krasnoyarsk) (M.V. Noskov, ed.),
	      // 1994, pp. 13--23.
	      //
	      // This rule is provably minimal in the number of points.
	      // A tensor-product rule accurate for "bi-tenth" or "bi-eleventh" degree polynomials would have 36 points.
	      const Real data[6][3] =
		{
		  {6.9807610454956756e-01, 9.8263922354085547e-01, 4.8020763350723814e-02}, // 4
		  {9.3948638281673690e-01, 8.2577583590296393e-01, 6.6071329164550595e-02}, // 4
		  {9.5353952820153201e-01, 1.8858613871864195e-01, 9.7386777358668164e-02}, // 4
		  {3.1562343291525419e-01, 8.1252054830481310e-01, 2.1173634999894860e-01}, // 4
		  {7.1200191307533630e-01, 5.2532025036454776e-01, 2.2562606172886338e-01}, // 4
		  {4.2484724884866925e-01, 4.1658071912022368e-02, 3.5115871839824543e-01}  // 4
		};

	      const unsigned int symmetry[6] = {
		4, // Rotational Invariant
		4, // Rotational Invariant
		4, // Rotational Invariant
		4, // Rotational Invariant
		4, // Rotational Invariant
		4  // Rotational Invariant
	      };

	      _points.resize (24);
	      _weights.resize(24);

	      stroud_rule(data, symmetry, 6);

	      return;
	    } // end case TENTH,ELEVENTH




	  case TWELFTH:
	  case THIRTEENTH:
	    {
	      // A degree 13, 33-point rule due to Cools and Haegemans.
	      //
	      // R. Cools and A. Haegemans, Another step forward in searching for
	      // cubature formulae with a minimal number of knots for the square,
	      // Computing 40 (1988), 139--146.
	      //
	      // A tensor-product rule accurate for "bi-12" or "bi-13" degree polynomials would have 49 points.
	      const Real data[9][3] =
		{
		  {0.0000000000000000e+00, 0.0000000000000000e+00, 3.0038211543122536e-01}, // 1
		  {9.8348668243987226e-01, 7.7880971155441942e-01, 2.9991838864499131e-02}, // 4
		  {8.5955600564163892e-01, 9.5729769978630736e-01, 3.8174421317083669e-02}, // 4
		  {9.5892517028753485e-01, 1.3818345986246535e-01, 6.0424923817749980e-02}, // 4
		  {3.9073621612946100e-01, 9.4132722587292523e-01, 7.7492738533105339e-02}, // 4
		  {8.5007667369974857e-01, 4.7580862521827590e-01, 1.1884466730059560e-01}, // 4
		  {6.4782163718701073e-01, 7.5580535657208143e-01, 1.2976355037000271e-01}, // 4
		  {7.0741508996444936e-02, 6.9625007849174941e-01, 2.1334158145718938e-01}, // 4
		  {4.0930456169403884e-01, 3.4271655604040678e-01, 2.5687074948196783e-01}  // 4
		};

	      const unsigned int symmetry[9] = {
		0, // Single point
		4, // Rotational Invariant
		4, // Rotational Invariant
		4, // Rotational Invariant
		4, // Rotational Invariant
		4, // Rotational Invariant
		4, // Rotational Invariant
		4, // Rotational Invariant
		4  // Rotational Invariant
	      };

	      _points.resize (33);
	      _weights.resize(33);

	      stroud_rule(data, symmetry, 9);

	      return;
	    } // end case TWELFTH,THIRTEENTH




	  case FOURTEENTH:
	  case FIFTEENTH:
	    {
	      // A degree-15, 48 point rule originally due to Rabinowitz and Richter,
	      // can be found in Cools' 1971 book.
	      //
	      // A.H. Stroud, Approximate calculation of multiple integrals,
	      // Prentice-Hall, Englewood Cliffs, N.J., 1971.
	      //
	      // The product Gauss rule for this order has 8^2=64 points.
	      const Real data[9][3] =
		{
		  {9.915377816777667e-01L, 0.0000000000000000e+00,  3.01245207981210e-02L}, // 4
		  {8.020163879230440e-01L, 0.0000000000000000e+00,  8.71146840209092e-02L}, // 4
		  {5.648674875232742e-01L, 0.0000000000000000e+00, 1.250080294351494e-01L}, // 4
		  {9.354392392539896e-01L, 0.0000000000000000e+00,  2.67651407861666e-02L}, // 4
		  {7.624563338825799e-01L, 0.0000000000000000e+00,  9.59651863624437e-02L}, // 4
		  {2.156164241427213e-01L, 0.0000000000000000e+00, 1.750832998343375e-01L}, // 4
		  {9.769662659711761e-01L, 6.684480048977932e-01L,  2.83136372033274e-02L}, // 4
		  {8.937128379503403e-01L, 3.735205277617582e-01L,  8.66414716025093e-02L}, // 4
		  {6.122485619312083e-01L, 4.078983303613935e-01L, 1.150144605755996e-01L}  // 4
		};

	      const unsigned int symmetry[9] = {
		3, // Full Symmetry, (x,0)
		3, // Full Symmetry, (x,0)
		3, // Full Symmetry, (x,0)
		2, // Full Symmetry, (x,x)
		2, // Full Symmetry, (x,x)
		2, // Full Symmetry, (x,x)
		1, // Full Symmetry, (x,y)
		1, // Full Symmetry, (x,y)
		1, // Full Symmetry, (x,y)
	      };

	      _points.resize (48);
	      _weights.resize(48);

	      stroud_rule(data, symmetry, 9);

	      return;
	    } // 	  case FOURTEENTH, FIFTEENTH:




	  case SIXTEENTH:
	  case SEVENTEENTH:
	    {
	      // A degree 17, 60-point rule due to Cools and Haegemans.
	      //
	      // R. Cools and A. Haegemans, Another step forward in searching for
	      // cubature formulae with a minimal number of knots for the square,
	      // Computing 40 (1988), 139--146.
	      //
	      // A tensor-product rule accurate for "bi-14" or "bi-15" degree polynomials would have 64 points.
	      // A tensor-product rule accurate for "bi-16" or "bi-17" degree polynomials would have 81 points.
	      const Real data[10][3] =
		{
		  {9.8935307451260049e-01, 0.0000000000000000e+00, 2.0614915919990959e-02}, // 4
		  {3.7628520715797329e-01, 0.0000000000000000e+00, 1.2802571617990983e-01}, // 4
		  {9.7884827926223311e-01, 0.0000000000000000e+00, 5.5117395340318905e-03}, // 4
		  {8.8579472916411612e-01, 0.0000000000000000e+00, 3.9207712457141880e-02}, // 4
		  {1.7175612383834817e-01, 0.0000000000000000e+00, 7.6396945079863302e-02}, // 4
		  {5.9049927380600241e-01, 3.1950503663457394e-01, 1.4151372994997245e-01}, // 8
		  {7.9907913191686325e-01, 5.9797245192945738e-01, 8.3903279363797602e-02}, // 8
		  {8.0374396295874471e-01, 5.8344481776550529e-02, 6.0394163649684546e-02}, // 8
		  {9.3650627612749478e-01, 3.4738631616620267e-01, 5.7387752969212695e-02}, // 8
		  {9.8132117980545229e-01, 7.0600028779864611e-01, 2.1922559481863763e-02}, // 8
		};

	      const unsigned int symmetry[10] = {
		3, // Fully symmetric (x,0)
		3, // Fully symmetric (x,0)
		2, // Fully symmetric (x,x)
		2, // Fully symmetric (x,x)
		2, // Fully symmetric (x,x)
		1, // Fully symmetric (x,y)
		1, // Fully symmetric (x,y)
		1, // Fully symmetric (x,y)
		1, // Fully symmetric (x,y)
		1  // Fully symmetric (x,y)
	      };

	      _points.resize (60);
	      _weights.resize(60);

	      stroud_rule(data, symmetry, 10);

	      return;
	    } // end case FOURTEENTH through SEVENTEENTH



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
	gauss_rule.init(type_in, p);

	// Swap points and weights with the about-to-be destroyed rule.
	_points.swap (gauss_rule.get_points() );
	_weights.swap(gauss_rule.get_weights());

	return;
      }
    } // end switch (type_in)
}

} // namespace libMesh

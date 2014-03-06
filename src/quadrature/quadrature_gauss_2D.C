// The libMesh Finite Element Library.
// Copyright (C) 2002-2014 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

namespace libMesh
{


void QGauss::init_2D(const ElemType type_in,
                     unsigned int p)
{
#if LIBMESH_DIM > 1

  //-----------------------------------------------------------------------
  // 2D quadrature rules
  switch (type_in)
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
    case TRI3SD:
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

              _points[0](0) = 1.0L/3.0L;
              _points[0](1) = 1.0L/3.0L;

              _weights[0] = 0.5;

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

              _points[0](0) = 2.0L/3.0L;
              _points[0](1) = 1.0L/6.0L;

              _points[1](0) = 1.0L/6.0L;
              _points[1](1) = 2.0L/3.0L;

              _points[2](0) = 1.0L/6.0L;
              _points[2](1) = 1.0L/6.0L;


              _weights[0] = 1.0L/6.0L;
              _weights[1] = 1.0L/6.0L;
              _weights[2] = 1.0L/6.0L;

              return;
            }
          case THIRD:
            {
              // Exact for cubics
              _points.resize(4);
              _weights.resize(4);

              // This rule is formed from a tensor product of
              // appropriately-scaled Gauss and Jacobi rules.  (See
              // also the QConical quadrature class, this is a
              // hard-coded version of one of those rules.)  For high
              // orders these rules generally have too many points,
              // but at extremely low order they are competitive and
              // have the additional benefit of having all positive
              // weights.
              _points[0](0) = 1.5505102572168219018027159252941e-01L;
              _points[0](1) = 1.7855872826361642311703513337422e-01L;
              _points[1](0) = 6.4494897427831780981972840747059e-01L;
              _points[1](1) = 7.5031110222608118177475598324603e-02L;
              _points[2](0) = 1.5505102572168219018027159252941e-01L;
              _points[2](1) = 6.6639024601470138670269327409637e-01L;
              _points[3](0) = 6.4494897427831780981972840747059e-01L;
              _points[3](1) = 2.8001991549907407200279599420481e-01L;

              _weights[0] = 1.5902069087198858469718450103758e-01L;
              _weights[1] = 9.0979309128011415302815498962418e-02L;
              _weights[2] = 1.5902069087198858469718450103758e-01L;
              _weights[3] = 9.0979309128011415302815498962418e-02L;

              return;


              // The following third-order rule is quite commonly cited
              // in the literature and most likely works fine.  However,
              // we generally prefer a rule with all positive weights
              // and an equal number of points, when available.
              //
              //  (allow_rules_with_negative_weights)
              // {
              //   // Exact for cubics
              //   _points.resize(4);
              //   _weights.resize(4);
              //
              //   _points[0](0) = .33333333333333333333333333333333;
              //   _points[0](1) = .33333333333333333333333333333333;
              //
              //   _points[1](0) = .2;
              //   _points[1](1) = .6;
              //
              //   _points[2](0) = .2;
              //   _points[2](1) = .2;
              //
              //   _points[3](0) = .6;
              //   _points[3](1) = .2;
              //
              //
              //   _weights[0] = -27./96.;
              //   _weights[1] =  25./96.;
              //   _weights[2] =  25./96.;
              //   _weights[3] =  25./96.;
              //
              //   return;
              // } // end if (allow_rules_with_negative_weights)
              // Note: if !allow_rules_with_negative_weights, fall through to next case.
            }



            // A degree 4 rule with six points.  This rule can be found in many places
            // including:
            //
            // J.N. Lyness and D. Jespersen, Moderate degree symmetric
            // quadrature rules for the triangle, J. Inst. Math. Appl.  15 (1975),
            // 19--32.
            //
            // We used the code in:
            // L. Zhang, T. Cui, and H. Liu. "A set of symmetric quadrature rules
            // on triangles and tetrahedra"  Journal of Computational Mathematics,
            // v. 27, no. 1, 2009, pp. 89-96.
            // to generate additional precision.
          case FOURTH:
            {
              const unsigned int n_wts = 2;
              const Real wts[n_wts] =
                {
                  1.1169079483900573284750350421656140e-01L,
                  5.4975871827660933819163162450105264e-02L
                };

              const Real a[n_wts] =
                {
                  4.4594849091596488631832925388305199e-01L,
                  9.1576213509770743459571463402201508e-02L
                };

              const Real b[n_wts] = {0., 0.}; // not used
              const unsigned int permutation_ids[n_wts] = {3, 3};

              dunavant_rule2(wts, a, b, permutation_ids, n_wts); // 6 total points

              return;
            }



            // Exact for quintics
            // Can be found in "Quadrature on Simplices of Arbitrary
            // Dimension" by Walkington.
          case FIFTH:
            {
              const unsigned int n_wts = 3;
              const Real wts[n_wts] =
                {
                  static_cast<Real>(9.0L/80.0L),
                  static_cast<Real>(31.0L/480.0L + std::sqrt(15.0L)/2400.0L),
                  static_cast<Real>(31.0L/480.0L - std::sqrt(15.0L)/2400.0L)
                };

              const Real a[n_wts] =
                {
                  0., // 'a' parameter not used for origin permutation
                  static_cast<Real>(2.0L/7.0L + std::sqrt(15.0L)/21.0L),
                  static_cast<Real>(2.0L/7.0L - std::sqrt(15.0L)/21.0L)
                };

              const Real b[n_wts] = {0., 0., 0.}; // not used
              const unsigned int permutation_ids[n_wts] = {1, 3, 3};

              dunavant_rule2(wts, a, b, permutation_ids, n_wts); // 7 total points

              return;
            }



            // A degree 6 rule with 12 points.  This rule can be found in many places
            // including:
            //
            // J.N. Lyness and D. Jespersen, Moderate degree symmetric
            // quadrature rules for the triangle, J. Inst. Math. Appl.  15 (1975),
            // 19--32.
            //
            // We used the code in:
            // L. Zhang, T. Cui, and H. Liu. "A set of symmetric quadrature rules
            // on triangles and tetrahedra"  Journal of Computational Mathematics,
            // v. 27, no. 1, 2009, pp. 89-96.
            // to generate additional precision.
            //
            // Note that the following 7th-order Ro3-invariant rule also has only 12 points,
            // which technically makes it the superior rule.  This one is here for completeness.
          case SIXTH:
            {
              const unsigned int n_wts = 3;
              const Real wts[n_wts] =
                {
                  5.8393137863189683012644805692789721e-02L,
                  2.5422453185103408460468404553434492e-02L,
                  4.1425537809186787596776728210221227e-02L
                };

              const Real a[n_wts] =
                {
                  2.4928674517091042129163855310701908e-01L,
                  6.3089014491502228340331602870819157e-02L,
                  3.1035245103378440541660773395655215e-01L
                };

              const Real b[n_wts] =
                {
                  0.,
                  0.,
                  6.3650249912139864723014259441204970e-01L
                };

              const unsigned int permutation_ids[n_wts] = {3, 3, 6}; // 12 total points

              dunavant_rule2(wts, a, b, permutation_ids, n_wts);

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
            // quadrature points) and has a few more digits of
            // precision in the points and weights.  Some 10-point
            // degree 6 rules exist for the triangle but they have
            // quadrature points outside the region of integration.
          case SEVENTH:
            {
              _points.resize (12);
              _weights.resize(12);

              const unsigned int nrows=4;

              // In each of the rows below, the first two entries are (z1, z2) which imply
              // z3.  The third entry is the weight for each of the points in the cyclic permutation.
              const Real rule_data[nrows][3] = {
                {6.2382265094402118e-02, 6.7517867073916085e-02, 2.6517028157436251e-02}, // group A
                {5.5225456656926611e-02, 3.2150249385198182e-01, 4.3881408714446055e-02}, // group B
                {3.4324302945097146e-02, 6.6094919618673565e-01, 2.8775042784981585e-02}, // group C
                {5.1584233435359177e-01, 2.7771616697639178e-01, 6.7493187009802774e-02}  // group D
              };

              for (unsigned int i=0, offset=0; i<nrows; ++i)
                {
                  _points[offset + 0] = Point(rule_data[i][0],                    rule_data[i][1]); // (z1,z2)
                  _points[offset + 1] = Point(1.-rule_data[i][0]-rule_data[i][1], rule_data[i][0]); // (z3,z1)
                  _points[offset + 2] = Point(rule_data[i][1], 1.-rule_data[i][0]-rule_data[i][1]); // (z2,z3)

                  // All these points get the same weight
                  _weights[offset + 0] = rule_data[i][2];
                  _weights[offset + 1] = rule_data[i][2];
                  _weights[offset + 2] = rule_data[i][2];

                  // Increment offset
                  offset += 3;
                }

              return;


              //       // The following is an inferior 7th-order Lyness-style rule with 15 points.
              //       // It's here only for completeness and the Ro3-invariant rule above should
              //       // be used instead!
              //       const unsigned int n_wts = 3;
              //       const Real wts[n_wts] =
              // {
              //   2.6538900895116205835977487499847719e-02L,
              //   3.5426541846066783659206291623201826e-02L,
              //   3.4637341039708446756138297960207647e-02L
              // };
              //
              //       const Real a[n_wts] =
              // {
              //   6.4930513159164863078379776030396538e-02L,
              //   2.8457558424917033519741605734978046e-01L,
              //   3.1355918438493150795585190219862865e-01L
              // };
              //
              //       const Real b[n_wts] =
              // {
              //   0.,
              //   1.9838447668150671917987659863332941e-01L,
              //   4.3863471792372471511798695971295936e-02L
              // };
              //
              //       const unsigned int permutation_ids[n_wts] = {3, 6, 6}; // 15 total points
              //
              //       dunavant_rule2(wts, a, b, permutation_ids, n_wts);
              //
              //       return;
            }




            // Another Dunavant rule.  This one has all positive weights.  This rule has
            // 16 points while a comparable conical product rule would have 5*5=25.
            //
            // It was copied 23rd June 2008 from:
            // http://people.scs.fsu.edu/~burkardt/f_src/dunavant/dunavant.f90
            //
            // Additional precision obtained from the code in:
            // L. Zhang, T. Cui, and H. Liu. "A set of symmetric quadrature rules
            // on triangles and tetrahedra"  Journal of Computational Mathematics,
            // v. 27, no. 1, 2009, pp. 89-96.
          case EIGHTH:
            {
              const unsigned int n_wts = 5;
              const Real wts[n_wts] =
                {
                  7.2157803838893584125545555244532310e-02L,
                  4.7545817133642312396948052194292159e-02L,
                  5.1608685267359125140895775146064515e-02L,
                  1.6229248811599040155462964170890299e-02L,
                  1.3615157087217497132422345036954462e-02L
                };

              const Real a[n_wts] =
                {
                  0.0, // 'a' parameter not used for origin permutation
                  4.5929258829272315602881551449416932e-01L,
                  1.7056930775176020662229350149146450e-01L,
                  5.0547228317030975458423550596598947e-02L,
                  2.6311282963463811342178578628464359e-01L,
                };

              const Real b[n_wts] =
                {
                  0.,
                  0.,
                  0.,
                  0.,
                  7.2849239295540428124100037917606196e-01L
                };

              const unsigned int permutation_ids[n_wts] = {1, 3, 3, 3, 6}; // 16 total points

              dunavant_rule2(wts, a, b, permutation_ids, n_wts);

              return;
            }



            // Another Dunavant rule.  This one has all positive weights.  This rule has 19
            // points. The comparable conical product rule would have 25.
            // It was copied 23rd June 2008 from:
            // http://people.scs.fsu.edu/~burkardt/f_src/dunavant/dunavant.f90
            //
            // Additional precision obtained from the code in:
            // L. Zhang, T. Cui, and H. Liu. "A set of symmetric quadrature rules
            // on triangles and tetrahedra"  Journal of Computational Mathematics,
            // v. 27, no. 1, 2009, pp. 89-96.
          case NINTH:
            {
              const unsigned int n_wts = 6;
              const Real wts[n_wts] =
                {
                  4.8567898141399416909620991253644315e-02L,
                  1.5667350113569535268427415643604658e-02L,
                  1.2788837829349015630839399279499912e-02L,
                  3.8913770502387139658369678149701978e-02L,
                  3.9823869463605126516445887132022637e-02L,
                  2.1641769688644688644688644688644689e-02L
                };

              const Real a[n_wts] =
                {
                  0.0, // 'a' parameter not used for origin permutation
                  4.8968251919873762778370692483619280e-01L,
                  4.4729513394452709865106589966276365e-02L,
                  4.3708959149293663726993036443535497e-01L,
                  1.8820353561903273024096128046733557e-01L,
                  2.2196298916076569567510252769319107e-01L
                };

              const Real b[n_wts] =
                {
                  0.,
                  0.,
                  0.,
                  0.,
                  0.,
                  7.4119859878449802069007987352342383e-01L
                };

              const unsigned int permutation_ids[n_wts] = {1, 3, 3, 3, 3, 6}; // 19 total points

              dunavant_rule2(wts, a, b, permutation_ids, n_wts);

              return;
            }


            // Another Dunavant rule with all positive weights.  This rule has 25
            // points. The comparable conical product rule would have 36.
            // It was copied 23rd June 2008 from:
            // http://people.scs.fsu.edu/~burkardt/f_src/dunavant/dunavant.f90
            //
            // Additional precision obtained from the code in:
            // L. Zhang, T. Cui, and H. Liu. "A set of symmetric quadrature rules
            // on triangles and tetrahedra"  Journal of Computational Mathematics,
            // v. 27, no. 1, 2009, pp. 89-96.
          case TENTH:
            {
              const unsigned int n_wts = 6;
              const Real wts[n_wts] =
                {
                  4.5408995191376790047643297550014267e-02L,
                  1.8362978878233352358503035945683300e-02L,
                  2.2660529717763967391302822369298659e-02L,
                  3.6378958422710054302157588309680344e-02L,
                  1.4163621265528742418368530791049552e-02L,
                  4.7108334818664117299637354834434138e-03L
                };

              const Real a[n_wts] =
                {
                  0.0, // 'a' parameter not used for origin permutation
                  4.8557763338365737736750753220812615e-01L,
                  1.0948157548503705479545863134052284e-01L,
                  3.0793983876412095016515502293063162e-01L,
                  2.4667256063990269391727646541117681e-01L,
                  6.6803251012200265773540212762024737e-02L
                };

              const Real b[n_wts] =
                {
                  0.,
                  0.,
                  0.,
                  5.5035294182099909507816172659300821e-01L,
                  7.2832390459741092000873505358107866e-01L,
                  9.2365593358750027664630697761508843e-01L
                };

              const unsigned int permutation_ids[n_wts] = {1, 3, 3, 6, 6, 6}; // 25 total points

              dunavant_rule2(wts, a, b, permutation_ids, n_wts);

              return;
            }


            // Dunavant's 11th-order rule contains points outside the region of
            // integration, and is thus unacceptable for our FEM calculations.
            //
            // This 30-point, 11th-order rule was obtained by me [JWP] using the code in
            //
            // Additional precision obtained from the code in:
            // L. Zhang, T. Cui, and H. Liu. "A set of symmetric quadrature rules
            // on triangles and tetrahedra"  Journal of Computational Mathematics,
            // v. 27, no. 1, 2009, pp. 89-96.
            //
            // Note: the 28-point 11th-order rule obtained by Zhang in the paper above
            // does not appear to be unique.  It is a solution in the sense that it
            // minimizes the error in the least-squares minimization problem, but
            // it involves too many unknowns and the Jacobian is therefore singular
            // when attempting to improve the solution via Newton's method.
          case ELEVENTH:
            {
              const unsigned int n_wts = 6;
              const Real wts[n_wts] =
                {
                  3.6089021198604635216985338480426484e-02L,
                  2.1607717807680420303346736867931050e-02L,
                  3.1144524293927978774861144478241807e-03L,
                  2.9086855161081509446654185084988077e-02L,
                  8.4879241614917017182977532679947624e-03L,
                  1.3795732078224796530729242858347546e-02L
                };

              const Real a[n_wts] =
                {
                  3.9355079629947969884346551941969960e-01L,
                  4.7979065808897448654107733982929214e-01L,
                  5.1003445645828061436081405648347852e-03L,
                  2.6597620190330158952732822450744488e-01L,
                  2.8536418538696461608233522814483715e-01L,
                  1.3723536747817085036455583801851025e-01L
                };

              const Real b[n_wts] =
                {
                  0.,
                  0.,
                  5.6817155788572446538150614865768991e-02L,
                  1.2539956353662088473247489775203396e-01L,
                  1.2409970153698532116262152247041742e-02L,
                  5.2792057988217708934207928630851643e-02L
                };

              const unsigned int permutation_ids[n_wts] = {3, 3, 6, 6, 6, 6}; // 30 total points

              dunavant_rule2(wts, a, b, permutation_ids, n_wts);

              return;
            }




            // Another Dunavant rule with all positive weights.  This rule has 33
            // points. The comparable conical product rule would have 36 (ELEVENTH) or 49 (TWELFTH).
            //
            // It was copied 23rd June 2008 from:
            // http://people.scs.fsu.edu/~burkardt/f_src/dunavant/dunavant.f90
            //
            // Additional precision obtained from the code in:
            // L. Zhang, T. Cui, and H. Liu. "A set of symmetric quadrature rules
            // on triangles and tetrahedra"  Journal of Computational Mathematics,
            // v. 27, no. 1, 2009, pp. 89-96.
          case TWELFTH:
            {
              const unsigned int n_wts = 8;
              const Real wts[n_wts] =
                {
                  3.0831305257795086169332418926151771e-03L,
                  3.1429112108942550177135256546441273e-02L,
                  1.7398056465354471494664198647499687e-02L,
                  2.1846272269019201067728631278737487e-02L,
                  1.2865533220227667708895461535782215e-02L,
                  1.1178386601151722855919538351159995e-02L,
                  8.6581155543294461858210504055170332e-03L,
                  2.0185778883190464758914349626118386e-02L
                };

              const Real a[n_wts] =
                {
                  2.1317350453210370246856975515728246e-02L,
                  2.7121038501211592234595134039689474e-01L,
                  1.2757614554158592467389632515428357e-01L,
                  4.3972439229446027297973662348436108e-01L,
                  4.8821738977380488256466206525881104e-01L,
                  2.8132558098993954824813069297455275e-01L,
                  1.1625191590759714124135414784260182e-01L,
                  2.7571326968551419397479634607976398e-01L
                };

              const Real b[n_wts] =
                {
                  0.,
                  0.,
                  0.,
                  0.,
                  0.,
                  6.9583608678780342214163552323607254e-01L,
                  8.5801403354407263059053661662617818e-01L,
                  6.0894323577978780685619243776371007e-01L
                };

              const unsigned int permutation_ids[n_wts] = {3, 3, 3, 3, 3, 6, 6, 6}; // 33 total points

              dunavant_rule2(wts, a, b, permutation_ids, n_wts);

              return;
            }


            // Another Dunavant rule with all positive weights.  This rule has 37
            // points. The comparable conical product rule would have 49 points.
            //
            // It was copied 23rd June 2008 from:
            // http://people.scs.fsu.edu/~burkardt/f_src/dunavant/dunavant.f90
            //
            // A second rule with additional precision obtained from the code in:
            // L. Zhang, T. Cui, and H. Liu. "A set of symmetric quadrature rules
            // on triangles and tetrahedra"  Journal of Computational Mathematics,
            // v. 27, no. 1, 2009, pp. 89-96.
          case THIRTEENTH:
            {
              const unsigned int n_wts = 9;
              const Real wts[n_wts] =
                {
                  3.3980018293415822140887212340442440e-02L,
                  2.7800983765226664353628733005230734e-02L,
                  2.9139242559599990702383541756669905e-02L,
                  3.0261685517695859208964000161454122e-03L,
                  1.1997200964447365386855399725479827e-02L,
                  1.7320638070424185232993414255459110e-02L,
                  7.4827005525828336316229285664517190e-03L,
                  1.2089519905796909568722872786530380e-02L,
                  4.7953405017716313612975450830554457e-03L
                };

              const Real a[n_wts] =
                {
                  0., // 'a' parameter not used for origin permutation
                  4.2694141425980040602081253503137421e-01L,
                  2.2137228629183290065481255470507908e-01L,
                  2.1509681108843183869291313534052083e-02L,
                  4.8907694645253934990068971909020439e-01L,
                  3.0844176089211777465847185254124531e-01L,
                  1.1092204280346339541286954522167452e-01L,
                  1.6359740106785048023388790171095725e-01L,
                  2.7251581777342966618005046435408685e-01L
                };

              const Real b[n_wts] =
                {
                  0.,
                  0.,
                  0.,
                  0.,
                  0.,
                  6.2354599555367557081585435318623659e-01L,
                  8.6470777029544277530254595089569318e-01L,
                  7.4850711589995219517301859578870965e-01L,
                  7.2235779312418796526062013230478405e-01L
                };

              const unsigned int permutation_ids[n_wts] = {1, 3, 3, 3, 3, 6, 6, 6, 6}; // 37 total points

              dunavant_rule2(wts, a, b, permutation_ids, n_wts);

              return;
            }


            // Another Dunavant rule.  This rule has 42 points, while
            // a comparable conical product rule would have 64.
            //
            // It was copied 23rd June 2008 from:
            // http://people.scs.fsu.edu/~burkardt/f_src/dunavant/dunavant.f90
            //
            // Additional precision obtained from the code in:
            // L. Zhang, T. Cui, and H. Liu. "A set of symmetric quadrature rules
            // on triangles and tetrahedra"  Journal of Computational Mathematics,
            // v. 27, no. 1, 2009, pp. 89-96.
          case FOURTEENTH:
            {
              const unsigned int n_wts = 10;
              const Real wts[n_wts] =
                {
                  1.0941790684714445320422472981662986e-02L,
                  1.6394176772062675320655489369312672e-02L,
                  2.5887052253645793157392455083198201e-02L,
                  2.1081294368496508769115218662093065e-02L,
                  7.2168498348883338008549607403266583e-03L,
                  2.4617018012000408409130117545210774e-03L,
                  1.2332876606281836981437622591818114e-02L,
                  1.9285755393530341614244513905205430e-02L,
                  7.2181540567669202480443459995079017e-03L,
                  2.5051144192503358849300465412445582e-03L
                };

              const Real a[n_wts] =
                {
                  4.8896391036217863867737602045239024e-01L,
                  4.1764471934045392250944082218564344e-01L,
                  2.7347752830883865975494428326269856e-01L,
                  1.7720553241254343695661069046505908e-01L,
                  6.1799883090872601267478828436935788e-02L,
                  1.9390961248701048178250095054529511e-02L,
                  1.7226668782135557837528960161365733e-01L,
                  3.3686145979634500174405519708892539e-01L,
                  2.9837288213625775297083151805961273e-01L,
                  1.1897449769695684539818196192990548e-01L
                };

              const Real b[n_wts] =
                {
                  0.,
                  0.,
                  0.,
                  0.,
                  0.,
                  0.,
                  7.7060855477499648258903327416742796e-01L,
                  5.7022229084668317349769621336235426e-01L,
                  6.8698016780808783735862715402031306e-01L,
                  8.7975717137017112951457163697460183e-01L
                };

              const unsigned int permutation_ids[n_wts]
                = {3, 3, 3, 3, 3, 3, 6, 6, 6, 6}; // 42 total points

              dunavant_rule2(wts, a, b, permutation_ids, n_wts);

              return;
            }


            // This 49-point rule was found by me [JWP] using the code in:
            //
            // L. Zhang, T. Cui, and H. Liu. "A set of symmetric quadrature rules
            // on triangles and tetrahedra"  Journal of Computational Mathematics,
            // v. 27, no. 1, 2009, pp. 89-96.
            //
            // A 54-point, 15th-order rule is reported by
            //
            // Stephen Wandzura, Hong Xiao,
            // Symmetric Quadrature Rules on a Triangle,
            // Computers and Mathematics with Applications,
            // Volume 45, Number 12, June 2003, pages 1829-1840.
            //
            // can be found here:
            // http://people.scs.fsu.edu/~burkardt/f_src/wandzura/wandzura.f90
            //
            // but this 49-point rule is superior.
          case FIFTEENTH:
            {
              const unsigned int n_wts = 11;
              const Real wts[n_wts] =
                {
                  2.4777380743035579804788826970198951e-02L,
                  9.2433943023307730591540642828347660e-03L,
                  2.2485768962175402793245929133296627e-03L,
                  6.7052581900064143760518398833360903e-03L,
                  1.9011381726930579256700190357527956e-02L,
                  1.4605445387471889398286155981802858e-02L,
                  1.5087322572773133722829435011138258e-02L,
                  1.5630213780078803020711746273129099e-02L,
                  6.1808086085778203192616856133701233e-03L,
                  3.2209366452594664857296985751120513e-03L,
                  5.8747373242569702667677969985668817e-03L
                };

              const Real a[n_wts] =
                {
                  0.0, // 'a' parameter not used for origin
                  7.9031013655541635005816956762252155e-02L,
                  1.8789501810770077611247984432284226e-02L,
                  4.9250168823249670532514526605352905e-01L,
                  4.0886316907744105975059040108092775e-01L,
                  5.3877851064220142445952549348423733e-01L,
                  2.0250549804829997692885033941362673e-01L,
                  5.5349674918711643207148086558288110e-01L,
                  7.8345022567320812359258882143250181e-01L,
                  8.9514624528794883409864566727625002e-01L,
                  3.2515745241110782862789881780746490e-01L
                };

              const Real b[n_wts] =
                {
                  0.,
                  0.,
                  0.,
                  0.,
                  0.,
                  1.9412620368774630292701241080996842e-01L,
                  9.8765911355712115933807754318089099e-02L,
                  7.7663767064308164090246588765178087e-02L,
                  2.1594628433980258573654682690950798e-02L,
                  1.2563596287784997705599005477153617e-02L,
                  1.5082654870922784345283124845552190e-02L
                };

              const unsigned int permutation_ids[n_wts]
                = {1, 3, 3, 3, 3, 6, 6, 6, 6, 6, 6}; // 49 total points

              dunavant_rule2(wts, a, b, permutation_ids, n_wts);

              return;
            }




            // Dunavant's 16th-order rule contains points outside the region of
            // integration, and is thus unacceptable for our FEM calculations.
            //
            // This 55-point, 16th-order rule was obtained by me [JWP] using the code in
            //
            // Additional precision obtained from the code in:
            // L. Zhang, T. Cui, and H. Liu. "A set of symmetric quadrature rules
            // on triangles and tetrahedra"  Journal of Computational Mathematics,
            // v. 27, no. 1, 2009, pp. 89-96.
            //
            // Note: the 55-point 16th-order rule obtained by Zhang in the paper above
            // does not appear to be unique.  It is a solution in the sense that it
            // minimizes the error in the least-squares minimization problem, but
            // it involves too many unknowns and the Jacobian is therefore singular
            // when attempting to improve the solution via Newton's method.
          case SIXTEENTH:
            {
              const unsigned int n_wts = 12;
              const Real wts[n_wts] =
                {
                  2.2668082505910087151996321171534230e-02L,
                  8.4043060714818596159798961899306135e-03L,
                  1.0850949634049747713966288634484161e-03L,
                  7.2252773375423638869298219383808751e-03L,
                  1.2997715227338366024036316182572871e-02L,
                  2.0054466616677715883228810959112227e-02L,
                  9.7299841600417010281624372720122710e-03L,
                  1.1651974438298104227427176444311766e-02L,
                  9.1291185550484450744725847363097389e-03L,
                  3.5568614040947150231712567900113671e-03L,
                  5.8355861686234326181790822005304303e-03L,
                  4.7411314396804228041879331486234396e-03L
                };

              const Real a[n_wts] =
                {
                  0.0, // 'a' parameter not used for centroid weight
                  8.5402539407933203673769900926355911e-02L,
                  1.2425572001444092841183633409631260e-02L,
                  4.9174838341891594024701017768490960e-01L,
                  4.5669426695387464162068900231444462e-01L,
                  4.8506759880447437974189793537259677e-01L,
                  2.0622099278664205707909858461264083e-01L,
                  3.2374950270039093446805340265853956e-01L,
                  7.3834330556606586255186213302750029e-01L,
                  9.1210673061680792565673823935174611e-01L,
                  6.6129919222598721544966837350891531e-01L,
                  1.7807138906021476039088828811346122e-01L
                };

              const Real b[n_wts] =
                {
                  0.0,
                  0.0,
                  0.0,
                  0.0,
                  0.0,
                  3.2315912848634384647700266402091638e-01L,
                  1.5341553679414688425981898952416987e-01L,
                  7.4295478991330687632977899141707872e-02L,
                  7.1278762832147862035977841733532020e-02L,
                  1.6623223223705792825395256602140459e-02L,
                  1.4160772533794791868984026749196156e-02L,
                  1.4539694958941854654807449467759690e-02L
                };

              const unsigned int permutation_ids[n_wts]
                = {1, 3, 3, 3, 3, 6, 6, 6, 6, 6, 6, 6}; // 55 total points

              dunavant_rule2(wts, a, b, permutation_ids, n_wts);

              return;
            }


            // Dunavant's 17th-order rule has 61 points, while a
            // comparable conical product rule would have 81 (16th and 17th orders).
            //
            // It can be found here:
            // http://people.scs.fsu.edu/~burkardt/f_src/dunavant/dunavant.f90
            //
            // Zhang reports an identical rule in:
            // L. Zhang, T. Cui, and H. Liu. "A set of symmetric quadrature rules
            // on triangles and tetrahedra"  Journal of Computational Mathematics,
            // v. 27, no. 1, 2009, pp. 89-96.
            //
            // Note: the 61-point 17th-order rule obtained by Dunavant and Zhang
            // does not appear to be unique.  It is a solution in the sense that it
            // minimizes the error in the least-squares minimization problem, but
            // it involves too many unknowns and the Jacobian is therefore singular
            // when attempting to improve the solution via Newton's method.
            //
            // Therefore, we prefer the following 63-point rule which
            // I [JWP] found.  It appears to be more accurate than the
            // rule reported by Dunavant and Zhang, even though it has
            // a few more points.
          case SEVENTEENTH:
            {
              const unsigned int n_wts = 12;
              const Real wts[n_wts] =
                {
                  1.7464603792572004485690588092246146e-02L,
                  5.9429003555801725246549713984660076e-03L,
                  1.2490753345169579649319736639588729e-02L,
                  1.5386987188875607593083456905596468e-02L,
                  1.1185807311917706362674684312990270e-02L,
                  1.0301845740670206831327304917180007e-02L,
                  1.1767783072977049696840016810370464e-02L,
                  3.8045312849431209558329128678945240e-03L,
                  4.5139302178876351271037137230354382e-03L,
                  2.2178812517580586419412547665472893e-03L,
                  5.2216271537483672304731416553063103e-03L,
                  9.8381136389470256422419930926212114e-04L
                };

              const Real a[n_wts] =
                {
                  2.8796825754667362165337965123570514e-01L,
                  4.9216175986208465345536805750663939e-01L,
                  4.6252866763171173685916780827044612e-01L,
                  1.6730292951631792248498303276090273e-01L,
                  1.5816335500814652972296428532213019e-01L,
                  1.6352252138387564873002458959679529e-01L,
                  6.2447680488959768233910286168417367e-01L,
                  8.7317249935244454285263604347964179e-01L,
                  3.4428164322282694677972239461699271e-01L,
                  9.1584484467813674010523309855340209e-02L,
                  2.0172088013378989086826623852040632e-01L,
                  9.6538762758254643474731509845084691e-01L
                };

              const Real b[n_wts] =
                {
                  0.0,
                  0.0,
                  0.0,
                  3.4429160695501713926320695771253348e-01L,
                  2.2541623431550639817203145525444726e-01L,
                  8.0670083153531811694942222940484991e-02L,
                  6.5967451375050925655738829747288190e-02L,
                  4.5677879890996762665044366994439565e-02L,
                  1.1528411723154215812386518751976084e-02L,
                  9.3057714323900610398389176844165892e-03L,
                  1.5916814107619812717966560404970160e-02L,
                  1.0734733163764032541125434215228937e-02L
                };

              const unsigned int permutation_ids[n_wts]
                = {3, 3, 3, 6, 6, 6, 6, 6, 6, 6, 6, 6}; // 63 total points

              dunavant_rule2(wts, a, b, permutation_ids, n_wts);

              return;

              //       _points.resize (61);
              //       _weights.resize(61);

              //       // The raw data for the quadrature rule.
              //       const Real p[15][4] = {
              // {                1./3.,                    0.,                    0., 0.033437199290803e+00 / 2.0}, // 1-perm
              // {0.005658918886452e+00, 0.497170540556774e+00,                    0., 0.005093415440507e+00 / 2.0}, // 3-perm
              // {0.035647354750751e+00, 0.482176322624625e+00,                    0., 0.014670864527638e+00 / 2.0}, // 3-perm
              // {0.099520061958437e+00, 0.450239969020782e+00,                    0., 0.024350878353672e+00 / 2.0}, // 3-perm
              // {0.199467521245206e+00, 0.400266239377397e+00,                    0., 0.031107550868969e+00 / 2.0}, // 3-perm
              // {0.495717464058095e+00, 0.252141267970953e+00,                    0., 0.031257111218620e+00 / 2.0}, // 3-perm
              // {0.675905990683077e+00, 0.162047004658461e+00,                    0., 0.024815654339665e+00 / 2.0}, // 3-perm
              // {0.848248235478508e+00, 0.075875882260746e+00,                    0., 0.014056073070557e+00 / 2.0}, // 3-perm
              // {0.968690546064356e+00, 0.015654726967822e+00,                    0., 0.003194676173779e+00 / 2.0}, // 3-perm
              // {0.010186928826919e+00, 0.334319867363658e+00, 0.655493203809423e+00, 0.008119655318993e+00 / 2.0}, // 6-perm
              // {0.135440871671036e+00, 0.292221537796944e+00, 0.572337590532020e+00, 0.026805742283163e+00 / 2.0}, // 6-perm
              // {0.054423924290583e+00, 0.319574885423190e+00, 0.626001190286228e+00, 0.018459993210822e+00 / 2.0}, // 6-perm
              // {0.012868560833637e+00, 0.190704224192292e+00, 0.796427214974071e+00, 0.008476868534328e+00 / 2.0}, // 6-perm
              // {0.067165782413524e+00, 0.180483211648746e+00, 0.752351005937729e+00, 0.018292796770025e+00 / 2.0}, // 6-perm
              // {0.014663182224828e+00, 0.080711313679564e+00, 0.904625504095608e+00, 0.006665632004165e+00 / 2.0}  // 6-perm
              //       };


              //       // Now call the dunavant routine to generate _points and _weights
              //       dunavant_rule(p, 15);

              //       return;
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
              const Real rule_data[17][4] = {
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
              dunavant_rule(rule_data, 17);

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
              const Real rule_data[19][4] = {
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
              dunavant_rule(rule_data, 19);

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
              const Real rule_data[26][4] = {
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
              dunavant_rule(rule_data, 26);

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
              const Real rule_data[36][4] = {
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
              dunavant_rule(rule_data, 36);

              return;
            }


            // By default, we fall back on the conical product rules.  If the user
            // requests an order higher than what is currently available in the 1D
            // rules, an error will be thrown from the respective 1D code.
          default:
            {
              // The following quadrature rules are generated as
              // conical products.  These tend to be non-optimal
              // (use too many points, cluster points in certain
              // regions of the domain) but they are quite easy to
              // automatically generate using a 1D Gauss rule on
              // [0,1] and two 1D Jacobi-Gauss rules on [0,1].
              QConical conical_rule(2, _order);
              conical_rule.init(type_in, p);

              // Swap points and weights with the about-to-be destroyed rule.
              _points.swap (conical_rule.get_points() );
              _weights.swap(conical_rule.get_weights());

              return;
            }
          }
      }


      //---------------------------------------------
      // Unsupported type
    default:
      {
        libMesh::err << "Element type not supported!:" << type_in << std::endl;
        libmesh_error();
      }
    }

  libmesh_error();

  return;

#endif
}

} // namespace libMesh

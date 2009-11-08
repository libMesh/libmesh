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



#ifndef __quadrature_gm_h__
#define __quadrature_gm_h__

// Local includes
#include "quadrature.h"




/**
 * This class implements the Grundmann-Moller quadrature rules for
 * tetrahedra.  The GM rules are well-defined for simplices of
 * arbitrary dimension and to any order, but the rules by Dunavant for
 * two-dimensional simplices are in general superior.  This is primarily
 * due to the fact that the GM rules contain a significant proportion 
 * of negative weights, making them susceptible to round-off error
 * at high-order.
 *
 * The GM rules are interesting in 3D because they overlap with the
 * conical product rules at higher order while having significantly
 * fewer evaluation points, making them potentially much more
 * efficient.  The table below gives a comparison between the number
 * of points in a conical product (CP) rule and the GM rule of
 * equivalent order.  The GM rules are defined to be exact for
 * polynomials of degree d=2*s+1, s=0,1,2,3,... The table also gives
 * the percentage of each GM rule's weights which are negative.
 * Although the percentage of negative weights does not grow
 * particularly quickly, the amplification factor (a measure of the
 * effect of round-off) defined as
 *
 * amp. factor = \f$ \frac{1}{V} \sum \|w_i\|, \f$
 *
 * where V is the volume of the reference element, does grow quickly.
 * (A rule with all positive has has an amplification factor of 1.0 by
 * definition.)
\verbatim
  s  | d     | N. CP        | N. GM   | % neg wts | amp. factor
-----------------------------------------------------------------
  0  | 1     |              | 1       |           | 
  1  | 2-3   |              | 5       |           |
  2  | 4-5   |              | 15      |           |
  3  | 6-7   |              | 35      | 31.43     |   11.94
  4  | 8-9   |  5^3=125     | 70      | 34.29     |   25.35 
  5  | 10-11 |  6^3=216     | 126     | 36.51     |   54.14
  6  | 12-13 |  7^3=343     | 210     | 38.10     |  116.30
  7  | 14-15 |  8^3=512     | 330     | 39.39     |  251.10
  8  | 16-17 |  9^3=729     | 495     | 40.40     |  544.68
  9  | 18-19 | 10^3=1,000   | 715     | 41.26     | 1186.16
 10  | 20-21 | 11^3=1,331   | 1,001   | 41.96     | 2591.97  
 11  | 22-23 | 12^3=1,728   | 1,365   | 42.56     | 5680.75
 ...
 16  | 32-33 | 17^3=4,913   | 4,845   |
 17  | 34-35 | 18^3=5,832   | 5,985   | <= Cross-over point, CP has fewer points for d >= 34
 18  | 36-37 | 19^3=6,859   | 7,315   | 
 ...
 21  | 42-43 | 22^3=10,648  | 12,650  | 
\endverbatim
 *
 * 
 * Reference:
 *    Axel Grundmann and Michael M\"{o}ller,
 *    "Invariant Integration Formulas for the N-Simplex 
 *    by Combinatorial Methods,"
 *    SIAM Journal on Numerical Analysis,
 *    Volume 15, Number 2, April 1978, pages 282-290.
 *
 * Reference LGPL Fortran90 code by John Burkardt can be found here: 
 * http://people.scs.fsu.edu/~burkardt/f_src/gm_rules/gm_rules.html
 *
 * @author John W. Peterson, 2008
 */

// ------------------------------------------------------------
// QGrundmann_Moller class definition

class QGrundmann_Moller : public QBase
{
 public:

  /**
   * Constructor.  Declares the order of the quadrature rule.
   */
  QGrundmann_Moller (const unsigned int _dim,
		     const Order _order=INVALID_ORDER);

  /**
   * Destructor.
   */
  ~QGrundmann_Moller();

  /**
   * @returns \p QGRUNDMANN_MOLLER
   */
  QuadratureType type() const { return QGRUNDMANN_MOLLER; }

 
 private:

  void init_1D (const ElemType,
		unsigned int =0)
  {
    // See about making this non-pure virtual in the base class
    libmesh_error();
  }

  /**
   * The GM rules are only defined for 3D since better 2D rules
   * for simplexes are available.
   */
  void init_3D (const ElemType _type=INVALID_ELEM,
		unsigned int p_level=0);

  /**
   * This routine is called from the different cases of init_3D().
   * It actually fills the _points and _weights vectors for a given rule index, s.
   */
  void gm_rule(unsigned int s);
     
  
  /**
   * Routine which generates p-compositions of a given order, s,
   * as well as permutations thereof.  This routine is called internally by
   * the gm_rule() routine, you should not call this yourself!
   */
  void compose_all(unsigned int s, // number to be compositioned
		   unsigned int p, // # of partitions
		   std::vector<std::vector<unsigned int> >& result);

};







#endif

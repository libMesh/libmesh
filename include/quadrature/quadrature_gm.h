// The libMesh Finite Element Library.
// Copyright (C) 2002-2016 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_QUADRATURE_GM_H
#define LIBMESH_QUADRATURE_GM_H

// Local includes
#include "libmesh/quadrature.h"

// C++ includes

namespace libMesh
{

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
 * The percentage of negative weights appears to approach 50, and the
 * amplification factor (a measure of the effect of round-off) defined
 * as
 *
 * amp. factor = \f$ \frac{1}{V} \sum \|w_i\|, \f$
 *
 * where V is the volume of the reference element, grows like exp(C*s).
 * (A rule with all positive weights has an amplification factor of
 * 1.0 by definition.)
 *
 * \verbatim
 * s   degree  n_pts(conical)  n_pts(GM)   % neg wts  amp. factor
 * ------------------------------------------------------------------------
 * 0   1       1               1            0.00      1.00e+00
 * 1   3       8               5           20.00      2.60e+00
 * 2   5       27              15          26.67      5.63e+00
 * 3   7       64              35          31.43      1.19e+01
 * 4   9       125             70          34.29      2.54e+01
 * 5   11      216             126         36.51      5.41e+01
 * 6   13      343             210         38.10      1.16e+02
 * 7   15      512             330         39.39      2.51e+02
 * 8   17      729             495         40.40      5.45e+02
 * 9   19      1000            715         41.26      1.19e+03
 * 10  21      1331            1001        41.96      2.59e+03
 * 11  23      1728            1365        42.56      5.68e+03
 * 12  25      2197            1820        43.08      1.25e+04
 * 13  27      2744            2380        43.53      2.75e+04
 * 14  29      3375            3060        43.92      6.07e+04
 * 15  31      4096            3876        44.27      1.34e+05
 * 16  33      4913            4845        44.58      2.97e+05
 * 17  35      5832            5985        44.86      6.59e+05 <= Conical rule has fewer points for degree >= 34
 * 18  37      6859            7315        45.11      1.46e+06
 * 19  39      8000            8855        45.34      3.25e+06
 * 20  41      9261            10626       45.55      7.23e+06
 * 21  43      10648           12650       45.74      1.61e+07
 * \endverbatim
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
 * \author John W. Peterson
 * \date 2008
 */
class QGrundmann_Moller libmesh_final : public QBase
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
  virtual QuadratureType type() const libmesh_override { return QGRUNDMANN_MOLLER; }


private:

  virtual void init_1D (const ElemType,
                        unsigned int =0) libmesh_override
  {
    // See about making this non-pure virtual in the base class
    libmesh_not_implemented();
  }

  /**
   * Initialize a 3D GM rule.  Only makes sense for Tets.
   */
  virtual void init_3D (const ElemType _type=INVALID_ELEM,
                        unsigned int p_level=0) libmesh_override;

  /**
   * Initialize a 2D GM rule.  Only makes sense for Tris.
   */
  virtual void init_2D (const ElemType _type=INVALID_ELEM,
                        unsigned int p_level=0) libmesh_override;

  /**
   * This routine is called from init_2D() and init_3D().  It actually
   * fills the _points and _weights vectors for a given rule index, s
   * and dimension, dim.
   */
  void gm_rule(unsigned int s, unsigned int dim);

  /**
   * Routine which generates p-compositions of a given order, s,
   * as well as permutations thereof.  This routine is called internally by
   * the gm_rule() routine, you should not call this yourself!
   */
  void compose_all(unsigned int s, // number to be compositioned
                   unsigned int p, // # of partitions
                   std::vector<std::vector<unsigned int> > & result);
};



} // namespace libMesh





#endif // LIBMESH_QUADRATURE_GM_H

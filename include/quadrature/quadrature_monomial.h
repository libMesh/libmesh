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



#ifndef LIBMESH_QUADRATURE_MONOMIAL_H
#define LIBMESH_QUADRATURE_MONOMIAL_H

// Local includes
#include "libmesh/quadrature.h"

// C++ includes

namespace libMesh
{

/**
 * This class defines alternate quadrature rules on
 * "tensor-product" elements (QUADs and HEXes) which can be
 * useful when integrating monomial finite element bases.
 *
 * While tensor product rules are ideal for integrating
 * bi/tri-linear, bi/tri-quadratic, etc. (i.e. tensor product)
 * bases (which consist of incomplete polynomials up to degree=
 * dim*p) they are not optimal for the MONOMIAL or FEXYZ bases,
 * which consist of complete polynomials of degree=p.
 *
 * This class is implemented to provide quadrature rules which are
 * more efficient than tensor product rules when they are available,
 * and fall back on Gaussian quadrature rules when necessary.
 *
 * A number of these rules have been helpfully collected in electronic form by:
 *
 * Prof. Ronald Cools
 * Katholieke Universiteit Leuven,
 * Dept. Computerwetenschappen
 * http://www.cs.kuleuven.ac.be/~nines/research/ecf/ecf.html
 *
 * (A username and password to access the tables is available by request.)
 *
 * We also provide the original reference for each rule, as available,
 * in the source code file.
 *
 * \author John W. Peterson
 * \date 2008
 */
class QMonomial libmesh_final : public QBase
{
public:

  /**
   * Constructor.  Declares the order of the quadrature rule.
   */
  QMonomial (const unsigned int _dim,
             const Order _order=INVALID_ORDER);

  /**
   * Destructor.
   */
  ~QMonomial();

  /**
   * @returns \p QMONOMIAL
   */
  virtual QuadratureType type() const libmesh_override { return QMONOMIAL; }


private:

  /**
   * Just uses a Gauss rule in 1D.
   */
  virtual void init_1D (const ElemType,
                        unsigned int =0) libmesh_override;

  /**
   * More efficient rules for QUADs
   */
  virtual void init_2D (const ElemType _type=INVALID_ELEM,
                        unsigned int p_level=0) libmesh_override;

  /**
   * More efficient rules for HEXes
   */
  virtual void init_3D (const ElemType _type=INVALID_ELEM,
                        unsigned int p_level=0) libmesh_override;



  /**
   * Wissmann published three interesting "partially symmetric" rules
   * for integrating degree 4, 6, and 8 polynomials exactly on QUADs.
   * These rules have all positive weights, all points inside the
   * reference element, and have fewer points than tensor-product
   * rules of equivalent order, making them superior to those rules
   * for monomial bases.
   *
   * J. W. Wissman and T. Becker, Partially symmetric cubature
   * formulas for even degrees of exactness, SIAM J. Numer. Anal.  23
   * (1986), 676--685.
   */
  void wissmann_rule(const Real rule_data[][3],
                     const unsigned int n_pts);

  /**
   * Stroud's rules for QUADs and HEXes can have one of several
   * different types of symmetry.  The rule_symmetry array describes
   * how the different lines of the rule_data array are to be
   * applied.  The different rule_symmetry possibilities are:
   * 0)  Origin or single-point: (x,y)
   * Fully-symmetric, 3 cases:
   *   1) (x,y) -> (x,y), (-x,y), (x,-y), (-x,-y)
   *               (y,x), (-y,x), (y,-x), (-y,-x)
   *   2) (x,x) -> (x,x), (-x,x), (x,-x), (-x,-x)
   *   3) (x,0) -> (x,0), (-x,0), (0, x), ( 0,-x)
   * 4) Rotational Invariant, (x,y) -> (x,y), (-x,-y), (-y, x), (y,-x)
   * 5) Partial Symmetry,     (x,y) -> (x,y), (-x, y) [x!=0]
   * 6) Rectangular Symmetry, (x,y) -> (x,y), (-x, y), (-x,-y), (x,-y)
   * 7) Central Symmetry,     (0,y) -> (0,y), ( 0,-y)
   *
   * Not all rules with these symmetries are due to Stroud, however,
   * his book is probably the most frequently-cited compendium of
   * quadrature rules and later authors certainly built upon his work.
   */
  void stroud_rule(const Real rule_data[][3],
                   const unsigned int * rule_symmetry,
                   const unsigned int n_pts);

  /**
   * Rules from Kim and Song, Comm. Korean Math. Soc vol. 13, no. 4, 1998, pp. 913-931.
   * The rules are obtained by considering the group G^{rot} of rotations of the reference
   * hex, and the invariant polynomials of this group.
   *
   * In Kim and Song's rules, quadrauture points are described by the following points
   * and their unique permutations under the G^{rot} group:
   *
   * 0.) (0,0,0) ( 1 perm ) -> [0, 0, 0]
   * 1.) (x,0,0) ( 6 perms) -> [x, 0, 0], [0, -x, 0], [-x, 0, 0], [0, x, 0], [0, 0, -x], [0, 0, x]
   * 2.) (x,x,0) (12 perms) -> [x, x, 0], [x, -x, 0], [-x, -x, 0], [-x, x, 0], [x, 0, -x], [x, 0, x],
   *                           [0, x, -x], [0, x, x], [0, -x, -x], [-x, 0, -x], [0, -x, x], [-x, 0, x]
   * 3.) (x,y,0) (24 perms) -> [x, y, 0], [y, -x, 0], [-x, -y, 0], [-y, x, 0], [x, 0, -y], [x, -y, 0],
   *                           [x, 0, y], [0, y, -x], [-x, y, 0],
   *                           [0, y, x], [y, 0, -x], [0, -y, -x], [-y, 0, -x], [y, x, 0], [-y, -x, 0],
   *                           [y, 0, x], [0, -y, x], [-y, 0, x],
   *                           [-x, 0, y], [0, -x, -y], [0, -x, y], [-x, 0, -y], [0, x, y], [0, x, -y]
   * 4.) (x,x,x) ( 8 perms) -> [x, x, x], [x, -x, x], [-x, -x, x], [-x, x, x], [x, x, -x], [x, -x, -x],
   *                           [-x, x, -x], [-x, -x, -x]
   * 5.) (x,x,z) (24 perms) -> [x, x, z], [x, -x, z], [-x, -x, z], [-x, x, z], [x, z, -x], [x, -x, -z],
   *                           [x, -z, x], [z, x, -x], [-x, x, -z],
   *                           [-z, x, x], [x, -z, -x], [-z, -x, -x], [-x, z, -x], [x, x, -z],
   *                           [-x, -x, -z], [x, z, x], [z, -x, x],
   *                           [-x, -z, x], [-x, z, x], [z, -x, -x], [-z, -x, x], [-x, -z, -x],
   *                           [z, x, x], [-z, x, -x]
   * 6.) (x,y,z) (24 perms) -> [x, y, z], [y, -x, z], [-x, -y, z], [-y, x, z], [x, z, -y], [x, -y, -z],
   *                           [x, -z, y], [z, y, -x], [-x, y, -z],
   *                           [-z, y, x], [y, -z, -x], [-z, -y, -x], [-y, z, -x], [y, x, -z], [-y, -x, -z],
   *                           [y, z, x], [z, -y, x], [-y, -z, x], [-x, z, y], [z, -x, -y],
   *                           [-z, -x, y], [-x, -z, -y], [z, x, y], [-z, x, -y]
   *
   * Only two of Kim and Song's rules are particularly useful for FEM calculations: the degree 7,
   * 38-point rule and their degree 8, 47-point rule.  The others either contain negative weights or points
   * outside the reference interval.  The points and weights, to 32 digits, were obtained from:
   * Ronald Cools' website (http://www.cs.kuleuven.ac.be/~nines/research/ecf/ecf.html) and
   * the unique permutations of G^{rot} were computed by me [JWP] using Maple.
   */
  void kim_rule(const Real rule_data[][4],
                const unsigned int * rule_id,
                const unsigned int n_pts);
};

} // namespace libMesh

#endif // LIBMESH_QUADRATURE_MONOMIAL_H

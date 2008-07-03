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



#ifndef __quadrature_gauss_h__
#define __quadrature_gauss_h__

// Local includes
#include "quadrature.h"




/**
 * This class implemenets specific orders of Gauss quadrature.
 * Gauss quadrature rules of Order \p p have the property of
 * integrating polynomials of degree \p p exactly.
 */

// ------------------------------------------------------------
// QGauss class definition

class QGauss : public QBase
{
 public:

  /**
   * Constructor.  Declares the order of the quadrature rule.
   */
  QGauss (const unsigned int _dim,
	  const Order _order=INVALID_ORDER);

  /**
   * Destructor.
   */
  ~QGauss();

  /**
   * @returns \p QGAUSS
   */
  QuadratureType type() const { return QGAUSS; }

 
 private:

  void init_1D (const ElemType _type=INVALID_ELEM,
		unsigned int p_level=0);
  void init_2D (const ElemType _type=INVALID_ELEM,
		unsigned int p_level=0);
  void init_3D (const ElemType _type=INVALID_ELEM,
		unsigned int p_level=0);


  /**
   * The Dunavant rule is for triangles.  It takes permutation points and
   * weights in a specific format as input and fills the pre-sized _points and _weights
   * vectors.  This function is only used internally by the TRI geometric
   * elements.
   */
  void dunavant_rule(const Real rule_data[][4],
		     const unsigned int n_pts);

  /**
   * The Keast rule is for tets.  It takes permutation points and weights
   * in a specific format as input and fills the pre-sized _points and _weights
   * vectors.  This function is only used internally by the TET geometric
   * elements.
   */
  void keast_rule(const Real rule_data[][4],
		  const unsigned int n_pts);

  /**
   * Wissmann published three interesting "partially symmetric" rules
   * for integrating degree 4, 6, and 8 polynomials exactly on QUADs.
   * These rules have all positive weights, all points inside the
   * reference element, and have fewer points than tensor-product
   * rules of the same order, making them superior to those rules.
   * 
   * J. W. Wissman and T. Becker, Partially symmetric cubature 
   * formulas for even degrees of exactness, SIAM J. Numer. Anal.  23
   * (1986), 676--685.
  */
  void wissmann_rule(const Real rule_data[][3],
		     const unsigned int n_pts);
};



// ------------------------------------------------------------
// QGauss class members
inline
QGauss::QGauss(const unsigned int d,
	       const Order o) : QBase(d,o)
{
  // explicitly call the init function in 1D since the
  // other tensor-product rules require this one.
  // note that EDGE will not be used internally, however
  // if we called the function with INVALID_ELEM it would try to
  // be smart and return, thinking it had already done the work.
  if (_dim == 1)
    init(EDGE2);
}




inline
QGauss::~QGauss()
{
}




#endif

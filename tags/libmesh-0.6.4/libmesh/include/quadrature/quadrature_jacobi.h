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



#ifndef __quadrature_jacobi_h__
#define __quadrature_jacobi_h__

// C++ includes

// Local includes
#include "quadrature.h"


// ------------------------------------------------------------
// QJacobi class definition


/**
 * This class implements two (for now) Jacobi-Gauss quadrature
 * rules.  These rules have the same order of accuracy as the
 * normal Gauss quadrature rules, but instead of integrating
 * a weight function w(x)=1, they integrate w(x)=(1-x)^alpha * (1+x)^beta.
 * The reason that they are useful is that they allow you to
 * implement conical product rules for triangles and tetrahedra.
 * Although these product rules are non-optimal (use more points
 * than necessary) they are automatically constructable for high
 * orders of accuracy where other formulae may not exist. 
 *
 * There is not much sense in using this class directly,
 * since it only provides 1D rules, weighted, as described before.
 * Still, this class is particularly helpful: check \p QGauss
 * for triangles and tetrahedra, with orders beyond \p FIFTH.
 */
class QJacobi : public QBase
{
 public:

  /**
   * Constructor.  Currently, only one-dimensional rules provided.
   * Check \p QGauss for versions of Jacobi quadrature rule for
   * higher dimensions.
   */
  QJacobi (const unsigned int _dim,
	   const Order _order=INVALID_ORDER,
	   const unsigned int _alpha=1,
	   const unsigned int _beta=0);

  /**
   * Destructor. Empty.
   */
  ~QJacobi() {}

  /**
   * @returns the \p QuadratureType, either 
   * \p QJACOBI_1_0 or \p QJACOBI_2_0.
   */
  QuadratureType type() const;


 private:
  const unsigned int _alpha;
  const unsigned int _beta;
  
  void init_1D (const ElemType _type=INVALID_ELEM,
		unsigned int p_level=0);

};




// ------------------------------------------------------------
// QJacobi class members
inline
QJacobi::QJacobi(const unsigned int d,
		 const Order o,
		 const unsigned int a,
		 const unsigned int b) : QBase(d,o), _alpha(a), _beta(b)
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
QuadratureType QJacobi::type() const
{
  if ((_alpha == 1) && (_beta == 0))
    return QJACOBI_1_0;

  else if ((_alpha == 2) && (_beta == 0))
    return QJACOBI_2_0;

  else
    { 
      libmesh_error();
      return INVALID_Q_RULE;
    }
}



#endif

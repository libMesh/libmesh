// $Id: quadrature_gauss.h,v 1.3 2003-01-21 19:24:35 benkirk Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2002  Benjamin S. Kirk, John W. Peterson
  
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
 * Gauss quadrature rules of order \p p have the property of
 * integrating polynomials of degree \p 2p-1 exactly.
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

  
 private:

  void init_1D (const ElemType _type=INVALID_ELEM);
  void init_2D (const ElemType _type=INVALID_ELEM);
  void init_3D (const ElemType _type=INVALID_ELEM);
  
  void init_2D (const ElemType _type,
		const unsigned int side);
  void init_3D (const ElemType _type,
		const unsigned int side);

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
};




inline
QGauss::~QGauss()
{
};






#endif

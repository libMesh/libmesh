// $Id: quadrature_simpson.h,v 1.2 2004-01-03 15:37:42 benkirk Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2004  Benjamin S. Kirk, John W. Peterson
  
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



#ifndef __quadrature_simpson_h__
#define __quadrature_simpson_h__

// Local includes
#include "quadrature.h"




/**
 * This class implemenets Simpson quadrature.
 * This is the same thing as Newton-Cotes quadrature with three points.
 * Simpson's rule can integrate polynomials of degree three exactly.
 *
 * @author John W. Peterson, 2003
 */

// ------------------------------------------------------------
// QSimpson class definition

class QSimpson : public QBase
{
 public:

  /**
   * Constructor.  Declares the order of the quadrature rule.
   */
  QSimpson (const unsigned int _dim);

  /**
   * Destructor. Empty.
   */
  ~QSimpson() {}

  /**
   * @returns \p QSIMPSON
   */
  QuadratureType type() const { return QSIMPSON; }


 private:

  void init_1D (const ElemType _type=INVALID_ELEM);
  void init_2D (const ElemType _type=INVALID_ELEM);
  void init_3D (const ElemType _type=INVALID_ELEM);
  
};



// ------------------------------------------------------------
// QSimpson class members
inline
QSimpson::QSimpson(const unsigned int d) : QBase(d,THIRD)
{
  // explicitly call the init function in 1D since the
  // other tensor-product rules require this one.
  // note that EDGE will not be used internally, however
  // if we called the function with INVALID_ELEM it would try to
  // be smart and return, thinking it had already done the work.
  if (_dim == 1)
    init(EDGE2);
}




#endif

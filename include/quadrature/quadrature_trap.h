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



#ifndef LIBMESH_QUADRATURE_TRAP_H
#define LIBMESH_QUADRATURE_TRAP_H

// Local includes
#include "libmesh/quadrature.h"

// C++ includes

namespace libMesh
{


// ------------------------------------------------------------
// QTrap class definition


/**
 * This class implemenets trapezoidal quadratue.  Sometimes
 * also known as Newton-Cotes quadrature with two points.  These rules
 * sample at the corners and will integrate linears exactly.
 */
class QTrap : public QBase
{
public:

  /**
   * Constructor.  Declares the dimension of the quadrature rule.
   */
  explicit
  QTrap (const unsigned int _dim,
         const Order o=FIRST);

  /**
   * Destructor. Empty.
   */
  ~QTrap() {}

  /**
   * @returns \p QTRAP
   */
  QuadratureType type() const { return QTRAP; }

private:

  void init_1D (const ElemType _type=INVALID_ELEM,
                unsigned int p_level=0);
  void init_2D (const ElemType _type=INVALID_ELEM,
                unsigned int p_level=0);
  void init_3D (const ElemType _type=INVALID_ELEM,
                unsigned int p_level=0);

};

// ------------------------------------------------------------
// QTrap class members
inline
QTrap::QTrap(const unsigned int d,
             const Order) : QBase(d,FIRST)
{
  // explicitly call the init function in 1D since the
  // other tensor-product rules require this one.
  // note that EDGE will not be used internally, however
  // if we called the function with INVALID_ELEM it would try to
  // be smart and return, thinking it had already done the work.
  if (_dim == 1)
    init(EDGE2);
}


} // namespace libMesh



#endif // LIBMESH_QUADRATURE_TRAP_H

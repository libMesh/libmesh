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



#ifndef LIBMESH_QUADRATURE_SIMPSON_H
#define LIBMESH_QUADRATURE_SIMPSON_H

// Local includes
#include "libmesh/quadrature.h"

// C++ includes

namespace libMesh
{

/**
 * This class implemenets Simpson quadrature.
 * This is the same thing as Newton-Cotes quadrature with three points.
 * Simpson's rule can integrate polynomials of degree three exactly.
 *
 * \author John W. Peterson
 * \date 2003
 */
class QSimpson : public QBase
{
public:

  /**
   * Constructor.  Declares the order of the quadrature rule.
   */
  explicit
  QSimpson (const unsigned int _dim,
            const Order o=THIRD) :
    QBase(_dim, o)
  {
    // explicitly call the init function in 1D since the
    // other tensor-product rules require this one.
    // note that EDGE will not be used internally, however
    // if we called the function with INVALID_ELEM it would try to
    // be smart and return, thinking it had already done the work.
    if (_dim == 1)
      init(EDGE2);
  }

  /**
   * Destructor. Empty.
   */
  ~QSimpson() {}

  /**
   * @returns \p QSIMPSON
   */
  virtual QuadratureType type() const libmesh_override { return QSIMPSON; }


private:

  virtual void init_1D (const ElemType _type=INVALID_ELEM,
                        unsigned int p_level=0) libmesh_override;
  virtual void init_2D (const ElemType _type=INVALID_ELEM,
                        unsigned int p_level=0) libmesh_override;
  virtual void init_3D (const ElemType _type=INVALID_ELEM,
                        unsigned int p_level=0) libmesh_override;
};


} // namespace libMesh



#endif // LIBMESH_QUADRATURE_SIMPSON_H

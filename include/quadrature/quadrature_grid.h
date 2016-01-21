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



#ifndef LIBMESH_QUADRATURE_GRID_H
#define LIBMESH_QUADRATURE_GRID_H

// Local includes
#include "libmesh/quadrature.h"

// C++ includes

namespace libMesh
{

/**
 * This class creates quadrature points on a uniform grid, with
 * order+1 points on an edge.
 *
 * Unlike most libMesh quadrature rules, QGrid does *not* reduce
 * error exponentially on smooth functions as you increase the
 * quadrature order.  Instead, it reduces the error quadratically.
 * However, this error reduction is more reliable on non-smooth
 * functions.
 *
 * This quadrature type may be useful iff you are integrating
 * functions which have discontinuities or discontinuous derivatives
 * on scales smaller than your element size.
 *
 * \author Roy Stogner
 * \date 2005
 */
class QGrid libmesh_final : public QBase
{
public:

  /**
   * Constructor.  Declares the order of the quadrature rule.
   */
  QGrid (const unsigned int _dim,
         const Order _order=INVALID_ORDER) :
    QBase(_dim, _order)
  {}

  /**
   * Destructor.
   */
  ~QGrid() {}

  /**
   * @returns \p QGRID
   */
  virtual QuadratureType type() const libmesh_override { return QGRID; }


private:

  virtual void init_1D (const ElemType _type=INVALID_ELEM,
                        unsigned int p_level=0) libmesh_override;
  virtual void init_2D (const ElemType _type=INVALID_ELEM,
                        unsigned int p_level=0) libmesh_override;
  virtual void init_3D (const ElemType _type=INVALID_ELEM,
                        unsigned int p_level=0) libmesh_override;
};

} // namespace libMesh



#endif // LIBMESH_QUADRATURE_GRID_H

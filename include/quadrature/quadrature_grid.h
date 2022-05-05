// The libMesh Finite Element Library.
// Copyright (C) 2002-2022 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

namespace libMesh
{

/**
 * This class creates quadrature points on a uniform grid, with
 * order+1 points on an edge.
 *
 * Unlike most libMesh quadrature rules, QGrid does *not* reduce the
 * integration error exponentially on smooth functions as you increase
 * the quadrature order.  Instead, it reduces the error quadratically.
 * However, this error reduction is more reliable on non-smooth
 * functions.
 *
 * This quadrature type may be useful if you are integrating functions
 * which have discontinuities or discontinuous derivatives on scales
 * smaller than your element size.
 *
 * \author Roy Stogner
 * \date 2005
 * \brief Implements grid-based quadrature rules suitable for non-smooth functions.
 */
class QGrid final : public QBase
{
public:

  /**
   * Constructor.  Declares the order of the quadrature rule.
   */
  QGrid (unsigned int dim,
         Order order=INVALID_ORDER) :
    QBase(dim, order)
  {}

  /**
   * Copy/move ctor, copy/move assignment operator, and destructor are
   * all explicitly defaulted for this simple class.
   */
  QGrid (const QGrid &) = default;
  QGrid (QGrid &&) = default;
  QGrid & operator= (const QGrid &) = default;
  QGrid & operator= (QGrid &&) = default;
  virtual ~QGrid() = default;

  /**
   * \returns \p QGRID.
   */
  virtual QuadratureType type() const override;


private:

  virtual void init_1D (const ElemType, unsigned int) override;
  virtual void init_2D (const ElemType, unsigned int) override;
  virtual void init_3D (const ElemType, unsigned int) override;
};

} // namespace libMesh

#endif // LIBMESH_QUADRATURE_GRID_H

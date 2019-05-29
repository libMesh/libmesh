// The libMesh Finite Element Library.
// Copyright (C) 2002-2019 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

namespace libMesh
{

/**
 * This class implements trapezoidal quadrature.  Sometimes also known
 * as Newton-Cotes quadrature with two points.  These rules sample at
 * the corners and will integrate linears exactly.
 *
 * \author Benjamin Kirk
 * \date 2003
 * \brief Implements trapezoidal rule, i.e. nodal quadrature for linear elements.
 */
class QTrap : public QBase
{
public:

  /**
   * Constructor.  Declares the order of the quadrature rule.  We
   * explicitly call the \p init function in 1D since the other
   * tensor-product rules require this one.
   *
   * \note The element type, EDGE2, will not be used internally,
   * however if we called the function with INVALID_ELEM it would try
   * to be smart and return, thinking it had already done the work.
   */
  explicit
  QTrap (unsigned int dim,
         Order order=FIRST) :
    QBase(dim,order)
  {
    if (_dim == 1)
      init(EDGE2);
  }

  /**
   * Copy/move ctor, copy/move assignment operator, and destructor are
   * all explicitly defaulted for this simple class.
   */
  QTrap (const QTrap &) = default;
  QTrap (QTrap &&) = default;
  QTrap & operator= (const QTrap &) = default;
  QTrap & operator= (QTrap &&) = default;
  virtual ~QTrap() = default;

  /**
   * \returns \p QTRAP.
   */
  virtual QuadratureType type() const override;

private:

  virtual void init_1D (const ElemType, unsigned int) override;
  virtual void init_2D (const ElemType, unsigned int) override;
  virtual void init_3D (const ElemType, unsigned int) override;
};

} // namespace libMesh

#endif // LIBMESH_QUADRATURE_TRAP_H

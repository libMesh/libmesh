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



#ifndef LIBMESH_QUADRATURE_CLOUGH_H
#define LIBMESH_QUADRATURE_CLOUGH_H

// Local includes
#include "libmesh/quadrature.h"

namespace libMesh
{

/**
 * This class creates a Gaussian quadrature rule duplicated for each
 * subelement of a Clough-Tocher divided macroelement.
 *
 * \author Roy Stogner
 * \date 2005
 * \brief Implements quadrature rules for Clough-Tocher macroelements.
 */
class QClough final : public QBase
{
public:

  /**
   * Constructor.  Declares the order of the quadrature rule.
   */
  QClough (unsigned int dim,
           Order order=INVALID_ORDER) :
    QBase(dim, order)
  {}

  /**
   * Copy/move ctor, copy/move assignment operator, and destructor are
   * all explicitly defaulted for this simple class.
   */
  QClough (const QClough &) = default;
  QClough (QClough &&) = default;
  QClough & operator= (const QClough &) = default;
  QClough & operator= (QClough &&) = default;
  virtual ~QClough() = default;

  /**
   * \returns \p QCLOUGH.
   */
  virtual QuadratureType type() const override;


private:

  virtual void init_1D (const ElemType, unsigned int) override;
  virtual void init_2D (const ElemType, unsigned int) override;
  virtual void init_3D (const ElemType, unsigned int) override;
};

} // namespace libMesh

#endif // LIBMESH_QUADRATURE_CLOUGH_H

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



#ifndef LIBMESH_QUADRATURE_JACOBI_H
#define LIBMESH_QUADRATURE_JACOBI_H

// Local includes
#include "libmesh/quadrature.h"

// C++ includes

namespace libMesh
{

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
 * for triangles and tetrahedra, with orders beyond \p THIRTIETH.
 *
 * \author John W. Peterson
 * \date 2003
 */
class QJacobi libmesh_final : public QBase
{
public:

  /**
   * Constructor.  Currently, only one-dimensional rules provided.
   */
  QJacobi (const unsigned int _dim,
           const Order _order=INVALID_ORDER,
           const unsigned int a=1,
           const unsigned int b=0) :
    QBase(_dim, _order),
    _alpha(a),
    _beta(b)
  {
    if (_dim == 1)
      init(EDGE2);
  }

  /**
   * Destructor. Empty.
   */
  ~QJacobi() {}

  /**
   * @returns the \p QuadratureType, either
   * \p QJACOBI_1_0 or \p QJACOBI_2_0.
   */
  virtual QuadratureType type() const libmesh_override;

private:
  const unsigned int _alpha;
  const unsigned int _beta;

  virtual void init_1D (const ElemType _type=INVALID_ELEM,
                        unsigned int p_level=0) libmesh_override;
};

} // namespace libMesh

#endif // LIBMESH_QUADRATURE_JACOBI_H

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



#ifndef LIBMESH_QUADRATURE_GAUSS_H
#define LIBMESH_QUADRATURE_GAUSS_H

// Local includes
#include "libmesh/quadrature.h"

// C++ includes

namespace libMesh
{

/**
 * This class implemenets specific orders of Gauss quadrature.
 * Gauss quadrature rules of Order \p p have the property of
 * integrating polynomials of degree \p p exactly.
 *
 * \author Benjamin Kirk
 * \author John W. Peterson
 * \date 2003
 */
class QGauss : public QBase
{
public:

  /**
   * Constructor.  Declares the order of the quadrature rule.
   */
  QGauss (const unsigned int _dim,
          const Order _order=INVALID_ORDER) :
    QBase(_dim, _order)
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
   * Destructor.
   */
  ~QGauss() {}

  /**
   * @returns \p QGAUSS
   */
  virtual QuadratureType type() const libmesh_override { return QGAUSS; }


private:

  virtual void init_1D (const ElemType _type=INVALID_ELEM,
                        unsigned int p_level=0) libmesh_override;
  virtual void init_2D (const ElemType _type=INVALID_ELEM,
                        unsigned int p_level=0) libmesh_override;
  virtual void init_3D (const ElemType _type=INVALID_ELEM,
                        unsigned int p_level=0) libmesh_override;

  /**
   * The Dunavant rule is for triangles.  It takes permutation points and
   * weights in a specific format as input and fills the pre-sized _points and _weights
   * vectors.  This function is only used internally by the TRI geometric
   * elements.
   */
  void dunavant_rule(const Real rule_data[][4],
                     const unsigned int n_pts);

  void dunavant_rule2(const Real * wts,
                      const Real * a,
                      const Real * b,
                      const unsigned int * permutation_ids,
                      const unsigned int n_wts);

  /**
   * The Keast rule is for tets.  It takes permutation points and weights
   * in a specific format as input and fills the pre-sized _points and _weights
   * vectors.  This function is only used internally by the TET geometric
   * elements.
   */
  void keast_rule(const Real rule_data[][4],
                  const unsigned int n_pts);
};

} // namespace libMesh

#endif // LIBMESH_QUADRATURE_GAUSS_H

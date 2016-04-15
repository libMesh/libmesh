// The libMesh Finite Element Library.
// Copyright (C) 2002-2012 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_QUADRATURE_COMPOSITE_H
#define LIBMESH_QUADRATURE_COMPOSITE_H

#include "libmesh/libmesh_config.h"

#if defined(LIBMESH_HAVE_TRIANGLE) && defined(LIBMESH_HAVE_TETGEN)

// Local includes
#include "libmesh/quadrature.h"
#include "libmesh/elem_cutter.h"
#include "libmesh/fe_base.h"
#include "libmesh/auto_ptr.h"

// C++ includes

namespace libMesh
{
/**
 * This class implements generic composite quadrature rules.
 * Composite quadrature rules are constructed from any of the
 * supported rules by breaking an element into subelements and
 * applying the base rule on each subelement.  This class uses the
 * ElemCutter, which is only available if libmesh is configured with
 * --disable-strict-lgpl.
 *
 * \author Benjamin Kirk
 * \date 2013
 */
template <class QSubCell>
class QComposite libmesh_final : public QSubCell
{
public:

  /**
   * Import base class data into namespace scope.
   */
  using QSubCell::_dim;
  using QSubCell::_points;
  using QSubCell::_weights;
  using QSubCell::init;

  /**
   * Constructor.  Declares the order of the quadrature rule.
   */
  QComposite (const unsigned int _dim,
              const Order _order=INVALID_ORDER);

  /**
   * Destructor.
   */
  ~QComposite();

  /**
   * @returns \p QCOMPOSITE
   */
  virtual QuadratureType type() const libmesh_override { return QCOMPOSITE; }

  /**
   * Initializes the data structures for a specific, potentially cut
   * element.  The array \p vertex_distance_func contains vertex
   * values of a signed distance function that cuts the element.  This
   * interface is indended to be extended by derived classes that can
   * cut the element into subelements, for example, and constuct a
   * composite quadrature rule for the cut element.
   */
  virtual void init (const Elem & elem,
                     const std::vector<Real> & vertex_distance_func,
                     unsigned int p_level=0) libmesh_override;

private:

  /**
   *
   */
  void add_subelem_values (const std::vector<Elem const *> & subelem);

  /**
   * Subcell quadrature object.
   */
  QSubCell _q_subcell;

  /**
   * ElemCutter object.
   */
  ElemCutter _elem_cutter;

  /**
   * Lagrange FE to use for subcell mapping.
   */
  UniquePtr<FEBase> _lagrange_fe;
};


} // namespace libMesh



#endif // LIBMESH_HAVE_TRIANGLE && LIBMESH_HAVE_TETGEN
#endif // LIBMESH_QUADRATURE_COMPOSITE_H

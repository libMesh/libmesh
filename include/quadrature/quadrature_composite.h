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



#ifndef LIBMESH_QUADRATURE_COMPOSITE_H
#define LIBMESH_QUADRATURE_COMPOSITE_H

#include "libmesh/libmesh_config.h"

#if defined(LIBMESH_HAVE_TRIANGLE) && defined(LIBMESH_HAVE_TETGEN)

// Local includes
#include "libmesh/quadrature.h"
#include "libmesh/elem_cutter.h"
#include "libmesh/fe_base.h"

// C++ includes
#include <memory>

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
 * \brief A quadrature rule for subdivided elements.
 */
template <class QSubCell>
class QComposite final : public QSubCell
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
  QComposite (unsigned int dim,
              Order order=INVALID_ORDER);

  /**
   * This class contains a unique_ptr member, so it can't be default
   * copy constructed or assigned.
   */
  QComposite (const QComposite &) = delete;
  QComposite & operator= (const QComposite &) = delete;

  /**
   * Copy/move ctor, copy/move assignment operator, and destructor are
   * all explicitly defaulted for this simple class.
   */
  QComposite (QComposite &&) = default;
  QComposite & operator= (QComposite &&) = default;
  virtual ~QComposite() = default;

  /**
   * \returns \p QCOMPOSITE.
   */
  virtual QuadratureType type() const override;

  /**
   * Overrides the base class init() function, and uses the ElemCutter to
   * subdivide the element into "inside" and "outside" subelements.
   */
  virtual void init (const Elem & elem,
                     const std::vector<Real> & vertex_distance_func,
                     unsigned int p_level=0) override;

private:

  /**
   * Helper function called from init() to collect all the points and
   * weights of the subelement quadrature rules.
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
  std::unique_ptr<FEBase> _lagrange_fe;
};

} // namespace libMesh

#endif // LIBMESH_HAVE_TRIANGLE && LIBMESH_HAVE_TETGEN
#endif // LIBMESH_QUADRATURE_COMPOSITE_H

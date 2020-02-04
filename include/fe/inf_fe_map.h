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



#ifndef LIBMESH_INF_FE_MAP_H
#define LIBMESH_INF_FE_MAP_H

#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

// Local includes
#include "libmesh/fe_type.h"
#include "libmesh/point.h"

// C++ includes
#include <vector>

namespace libMesh
{

// forward declarations
template <typename> class ElemTempl;
typedef ElemTempl<Real> Elem;
class FEComputeData;



/**
 * Class that encapsulates mapping (i.e. from physical space to
 * reference space and vice-versa) operations on infinite elements.
 */
class InfFEMap
{
public:
  /**
   * \returns The location (in physical space) of the point
   * \p p located on the reference element.
   */
  static Point map (const unsigned int dim,
                    const Elem * inf_elem,
                    const Point & reference_point);

  /**
   * \returns The location (on the reference element) of the
   * point \p p located in physical space.  First, the location
   * in the base face is computed. This requires inverting the
   * (possibly nonlinear) transformation map in the base, so it is
   * not trivial. The optional parameter \p tolerance defines
   * how close is "good enough."  The map inversion iteration
   * computes the sequence \f$ \{ p_n \} \f$, and the iteration is
   * terminated when \f$ \|p - p_n\| < \mbox{\texttt{tolerance}} \f$.
   * Once the base face point is determined, the radial local
   * coordinate is directly evaluated.
   * If \p interpolated is true, the interpolated distance from the
   * base element to the infinite element origin is used for the map
   * in radial direction.
   */
  static Point inverse_map (const unsigned int dim,
                            const Elem * elem,
                            const Point & p,
                            const Real tolerance = TOLERANCE,
                            const bool secure = true);


  /**
   * Takes a number points in physical space (in the \p physical_points
   * vector) and finds their location on the reference element for the
   * input element \p elem.  The values on the reference element are
   * returned in the vector \p reference_points
   */
  static void inverse_map (const unsigned int dim,
                           const Elem * elem,
                           const std::vector<Point> & physical_points,
                           std::vector<Point> &       reference_points,
                           const Real tolerance = TOLERANCE,
                           const bool secure = true);

  /**
   * \returns The value of the \f$ i^{th} \f$ polynomial evaluated
   * at \p v.  This method provides the approximation
   * in radial direction for the overall shape functions,
   * which is defined in \p InfFE::shape().
   * This method is allowed to be static, since it is independent
   * of dimension and base_family.  It is templated, though,
   * w.r.t. to radial \p FEFamily.
   *
   * \returns The value of the \f$ i^{th} \f$ mapping shape function
   * in radial direction evaluated at \p v when T_radial ==
   * INFINITE_MAP.  Currently, only one specific mapping shape is
   * used.  Namely the one by Marques JMMC, Owen DRJ: Infinite
   * elements in quasi-static materially nonlinear problems, Computers
   * and Structures, 1984.
   */
  static Real eval(Real v,
                   Order o,
                   unsigned int i);

  /**
   * \returns The value of the first derivative of the
   * \f$ i^{th} \f$ polynomial at coordinate \p v.
   * See \p eval for details.
   */
  static Real eval_deriv(Real v,
                         Order o,
                         unsigned int i);
};


} // namespace libMesh


#endif // ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS


#endif // LIBMESH_INF_FE_MAP_H

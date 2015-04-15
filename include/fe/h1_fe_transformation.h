// The libMesh Finite Element Library.
// Copyright (C) 2002-2014 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#ifndef LIBMESH_H1_FE_TRANSFORMATION_H
#define LIBMESH_H1_FE_TRANSFORMATION_H

#include "libmesh/fe_transformation_base.h"
#include "libmesh/compare_types.h"
#include "libmesh/elem.h"

namespace libMesh
{
/**
 * This class handles the computation of the shape functions in the physical domain
 * for H1 conforming elements. This class assumes the \p FEGenericBase object has been
 * initialized in the reference domain (i.e. \p init_shape_functions has been called).
 *
 * @author Paul T. Bauman, 2012
 */
template< typename OutputShape >
class H1FETransformation : public FETransformationBase<OutputShape>
{
public:

  H1FETransformation()
    : FETransformationBase<OutputShape>(){}

  virtual ~H1FETransformation(){}

  /**
   * Evaluates shape functions in physical coordinates for H1 conforming elements.
   * In this case \f$ \phi(x) = \phi(\xi) \f$
   */
  virtual void map_phi( const unsigned int,
                        const Elem* const,
                        const std::vector<Point>&,
                        const FEGenericBase<OutputShape>&,
                        std::vector<std::vector<OutputShape> >& ) const;

  /**
   * Evaluates shape function gradients in physical coordinates for H1 conforming
   * elements. dphi/dx = dphi/dxi * dxi/dx, etc.
   */
  virtual void map_dphi( const unsigned int dim,
                         const Elem* const elem,
                         const std::vector<Point>& qp,
                         const FEGenericBase<OutputShape>& fe,
                         std::vector<std::vector<typename FEGenericBase<OutputShape>::OutputGradient> >& dphi,
                         std::vector<std::vector<OutputShape> >& dphidx,
                         std::vector<std::vector<OutputShape> >& dphidy,
                         std::vector<std::vector<OutputShape> >& dphidz) const;

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
  /**
   * Evaluates shape function Hessians in physical coordinates based on H1 conforming
   * finite element transformation.
   * FIXME: These are currently not calculated correctly for non-affine elements.
   *        The second derivative terms of the FEMap are not implemented.
   */
  virtual void map_d2phi( const unsigned int dim,
                          const std::vector<Point>& qp,
                          const FEGenericBase<OutputShape>& fe,
                          std::vector<std::vector<typename FEGenericBase<OutputShape>::OutputTensor> >& d2phi,
                          std::vector<std::vector<OutputShape> >& d2phidx2,
                          std::vector<std::vector<OutputShape> >& d2phidxdy,
                          std::vector<std::vector<OutputShape> >& d2phidxdz,
                          std::vector<std::vector<OutputShape> >& d2phidy2,
                          std::vector<std::vector<OutputShape> >& d2phidydz,
                          std::vector<std::vector<OutputShape> >& d2phidz2  ) const;
#endif //LIBMESH_ENABLE_SECOND_DERIVATIVES

  /**
   * Evaluates the shape function curl in physical coordinates based on H1 conforming
   * finite element transformation.
   */
  virtual void map_curl( const unsigned int dim,
                         const Elem* const elem,
                         const std::vector<Point>& qp,
                         const FEGenericBase<OutputShape>& fe,
                         std::vector<std::vector<OutputShape> >& curl_phi ) const;

  /**
   * Evaluates the shape function divergence in physical coordinates based on H1 conforming
   * finite element transformation.
   */
  virtual void map_div( const unsigned int dim,
                        const Elem* const elem,
                        const std::vector<Point>& qp,
                        const FEGenericBase<OutputShape>& fe,
                        std::vector<std::vector<typename FEGenericBase<OutputShape>::OutputDivergence> >& div_phi ) const;

}; // class H1FETransformation

}

#endif // LIBMESH_H1_FE_TRANSFORMATION_H

// The libMesh Finite Element Library.
// Copyright (C) 2002-2025 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#ifndef LIBMESH_HDIV_FE_TRANSFORMATION_H
#define LIBMESH_HDIV_FE_TRANSFORMATION_H

#include "libmesh/fe_transformation_base.h"

namespace libMesh
{

/**
 * This class handles the computation of the shape functions in the
 * physical domain for HDiv conforming elements. This class assumes
 * the \p FEGenericBase object has been initialized in the reference
 * domain (i.e. \p init_shape_functions has been called).
 *
 * \author Nuno Nobre
 * \date 2023
 */
template<typename OutputShape>
class HDivFETransformation : public FETransformationBase<OutputShape>
{
public:

  HDivFETransformation()
    : FETransformationBase<OutputShape>() {}

  virtual ~HDivFETransformation() = default;

  /**
   * Pre-requests any necessary data from FEMap
   */
  virtual void init_map_phi(const FEGenericBase<OutputShape> & fe) const override;

  /**
   * Pre-requests any necessary data from FEMap
   */
  virtual void init_map_dphi(const FEGenericBase<OutputShape> & fe) const override;

  /**
   * Pre-requests any necessary data from FEMap
   */
  virtual void init_map_d2phi(const FEGenericBase<OutputShape> & fe) const override;

  /**
   * Evaluates shape functions in physical coordinates for \f$ H(div)
   * \f$ conforming elements.  In this case \f$ \phi = J^{-1} (dx/d\xi)
   * \hat{\phi} \f$, where \f$ (dx/d\xi) \f$ is the Jacobian matrix of the
   * element map and J = \det( dx/d\xi ).
   *
   * \note Here \f$ x, \xi \f$ are vectors.
   */
  virtual void map_phi(const unsigned int dim,
                       const Elem * const elem,
                       const std::vector<Point> & qp,
                       const FEGenericBase<OutputShape> & fe,
                       std::vector<std::vector<OutputShape>> & phi,
                       bool add_p_level = true) const override;

  /**
   * Evaluates shape function gradients in physical coordinates for
   * \f$ H(div) \f$ conforming elements.
   */
  virtual void map_dphi(const unsigned int /*dim*/,
                        const Elem * const /*elem*/,
                        const std::vector<Point> & /*qp*/,
                        const FEGenericBase<OutputShape> & /*fe*/,
                        std::vector<std::vector<typename FEGenericBase<OutputShape>::OutputGradient>> & /*dphi*/,
                        std::vector<std::vector<OutputShape>> & /*dphidx*/,
                        std::vector<std::vector<OutputShape>> & /*dphidy*/,
                        std::vector<std::vector<OutputShape>> & /*dphidz*/) const override
  {
    libmesh_warning("WARNING: Shape function gradients for HDiv elements are not currently being computed!");
  }

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
  /**
   * Evaluates shape function Hessians in physical coordinates based
   * on \f$ H(div) \f$ conforming finite element transformation.
   */
  virtual void map_d2phi(const unsigned int /*dim*/,
                         const std::vector<Point> & /*qp*/,
                         const FEGenericBase<OutputShape> & /*fe*/,
                         std::vector<std::vector<typename FEGenericBase<OutputShape>::OutputTensor>> & /*d2phi*/,
                         std::vector<std::vector<OutputShape>> & /*d2phidx2*/,
                         std::vector<std::vector<OutputShape>> & /*d2phidxdy*/,
                         std::vector<std::vector<OutputShape>> & /*d2phidxdz*/,
                         std::vector<std::vector<OutputShape>> & /*d2phidy2*/,
                         std::vector<std::vector<OutputShape>> & /*d2phidydz*/,
                         std::vector<std::vector<OutputShape>> & /*d2phidz2*/) const override
  {
    libmesh_warning("WARNING: Shape function Hessians for HDiv elements are not currently being computed!");
  }
#endif //LIBMESH_ENABLE_SECOND_DERIVATIVES

  /**
   * Evaluates the shape function curl in physical coordinates
   * based on \f$ H(div) \f$ conforming finite element transformation.
   */
  virtual void map_curl(const unsigned int /*dim*/,
                        const Elem * const /*elem*/,
                        const std::vector<Point> & /*qp*/,
                        const FEGenericBase<OutputShape> & /*fe*/,
                        std::vector<std::vector<OutputShape>> & /*curl_phi*/) const override
  {
    libmesh_warning("WARNING: Shape function curls for HDiv elements are not currently being computed!");
  }

  /**
   * Evaluates the shape function divergence in physical coordinates based on \f$ H(div) \f$ conforming
   * finite element transformation.
   * The transformation is \f$ \nabla \cdot \phi = J^{-1} \nabla \cdot \hat{\phi} \f$ where
   * \f$ J = \det( dx/d\xi ) \f$
   */
  virtual void map_div(const unsigned int dim,
                       const Elem * const elem,
                       const std::vector<Point> & qp,
                       const FEGenericBase<OutputShape> & fe,
                       std::vector<std::vector<typename FEGenericBase<OutputShape>::OutputDivergence>> & div_phi) const override;

}; // class HDivFETransformation

}

#endif // LIBMESH_HDIV_FE_TRANSFORMATION_H

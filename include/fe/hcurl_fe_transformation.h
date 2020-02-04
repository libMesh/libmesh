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

#ifndef LIBMESH_HCURL_FE_TRANSFORMATION_H
#define LIBMESH_HCURL_FE_TRANSFORMATION_H

#include "libmesh/fe_transformation_base.h"

namespace libMesh
{

/**
 * This class handles the computation of the shape functions in the
 * physical domain for HCurl conforming elements. This class assumes
 * the \p FEGenericBase object has been initialized in the reference
 * domain (i.e. \p init_shape_functions has been called).
 *
 * \author Paul T. Bauman
 * \date 2012
 */
template<typename MathType, typename RealType = Real>
class HCurlFETransformation : public FETransformationBase<MathType,RealType>
{
public:
  typedef ElemTempl<RealType> Elem;
  typedef PointTempl<RealType> Point;
  typedef typename FEOutputTwoType<MathType,RealType>::OutputShape OutputShape;

  HCurlFETransformation()
      : FETransformationBase<MathType,RealType>(){}

  virtual ~HCurlFETransformation(){}

  /**
   * Pre-requests any necessary data from FEMap
   */
  virtual void init_map_phi(const FEGenericBase<MathType,RealType> & fe) const override;

  /**
   * Pre-requests any necessary data from FEMap
   */
  virtual void init_map_dphi(const FEGenericBase<MathType,RealType> & fe) const override;

  /**
   * Pre-requests any necessary data from FEMap
   */
  virtual void init_map_d2phi(const FEGenericBase<MathType,RealType> & fe) const override;

  /**
   * Evaluates shape functions in physical coordinates for \f$ H(curl)
   * \f$ conforming elements.  In this case \f$ \phi = (dx/d\xi)^{-T}
   * \hat{\phi} \f$, where \f$ (dx/d\xi)^{-T} \f$ is the
   * inverse-transpose of the Jacobian matrix of the element map.
   *
   * \note Here \f$ x, \xi \f$ are vectors.
   */
  virtual void map_phi(const unsigned int dim,
                       const Elem * const elem,
                       const std::vector<Point> & qp,
                       const FEGenericBase<MathType,RealType> & fe,
                       std::vector<std::vector<OutputShape>> & phi) const override;

  /**
   * Evaluates shape function gradients in physical coordinates for
   * \f$ H(curl) \f$ conforming elements.
   */
  virtual void map_dphi(const unsigned int /*dim*/,
                        const Elem * const /*elem*/,
                        const std::vector<Point> & /*qp*/,
                        const FEGenericBase<MathType,RealType> & /*fe*/,
                        std::vector<std::vector<typename FEOutputTwoType<MathType,RealType>::OutputGradient>> & /*dphi*/,
                        std::vector<std::vector<OutputShape>> & /*dphidx*/,
                        std::vector<std::vector<OutputShape>> & /*dphidy*/,
                        std::vector<std::vector<OutputShape>> & /*dphidz*/) const override
  {
    libmesh_warning("WARNING: Shape function gradients for HCurl elements are not currently being computed!");
  }

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
  /**
   * Evaluates shape function Hessians in physical coordinates based
   * on \f$ H(curl) \f$ conforming finite element transformation.
   */
  virtual void map_d2phi(const unsigned int /*dim*/,
                         const std::vector<Point> & /*qp*/,
                         const FEGenericBase<MathType,RealType> & /*fe*/,
                         std::vector<std::vector<typename FEOutputTwoType<MathType,RealType>::OutputTensor>> & /*d2phi*/,
                         std::vector<std::vector<OutputShape>> & /*d2phidx2*/,
                         std::vector<std::vector<OutputShape>> & /*d2phidxdy*/,
                         std::vector<std::vector<OutputShape>> & /*d2phidxdz*/,
                         std::vector<std::vector<OutputShape>> & /*d2phidy2*/,
                         std::vector<std::vector<OutputShape>> & /*d2phidydz*/,
                         std::vector<std::vector<OutputShape>> & /*d2phidz2*/) const override
  {
    libmesh_warning("WARNING: Shape function Hessians for HCurl elements are not currently being computed!");
  }
#endif //LIBMESH_ENABLE_SECOND_DERIVATIVES

  /**
   * Evaluates the shape function curl in physical coordinates based on \f$ H(curl) \f$ conforming
   * finite element transformation.
   * In 2-D, the transformation is \f$ \nabla \times \phi = J^{-1} * \nabla \times \hat{\phi} \f$ where
   * \f$ J = \det( dx/d\xi ) \f$
   * In 3-D, the transformation is \f$ \nabla \times \phi = J^{-1} dx/d\xi \nabla \times \hat{\phi} \f$
   */
  virtual void map_curl(const unsigned int dim,
                        const Elem * const elem,
                        const std::vector<Point> & qp,
                        const FEGenericBase<MathType,RealType> & fe,
                        std::vector<std::vector<OutputShape>> & curl_phi) const override;

  /**
   * Evaluates the shape function divergence in physical coordinates
   * based on \f$ H(curl) \f$ conforming finite element transformation.
   */
  virtual void map_div(const unsigned int /*dim*/,
                       const Elem * const /*elem*/,
                       const std::vector<Point> & /*qp*/,
                       const FEGenericBase<MathType,RealType> & /*fe*/,
                       std::vector<std::vector<typename FEOutputTwoType<MathType,RealType>::OutputDivergence>> & /*div_phi*/) const override
  {
    libmesh_warning("WARNING: Shape function divergences for HCurl elements are not currently being computed!");
  }

}; // class HCurlFETransformation

}

#endif // LIBMESH_HCURL_FE_TRANSFORMATION_H

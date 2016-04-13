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

#ifndef LIBMESH_FE_TRANSFORMATION_BASE_H
#define LIBMESH_FE_TRANSFORMATION_BASE_H

#include "libmesh/fe_base.h"

namespace libMesh
{

//Forward Declarations
template< typename T > class FEGenericBase;
template< typename T > class H1FETransformation;
template< typename T > class HCurlFETransformation;
class FEType;

/**
 * This class handles the computation of the shape functions in the physical domain.
 * Derived classes implement the particular mapping: H1, HCurl, HDiv, or L2. This
 * class assumes the \p FEGenericBase object has been initialized in the reference
 * domain (i.e. init_shape_functions has been called).
 *
 * \author Paul T. Bauman
 * \date 2012
 */
template< typename OutputShape >
class FETransformationBase
{
public:

  FETransformationBase(){}
  virtual ~FETransformationBase(){}

  /**
   * Builds an FETransformation object based on the finite element type
   */
  static UniquePtr<FETransformationBase<OutputShape> > build(const FEType & type);

  /**
   * Pre-requests any necessary data from FEMap
   */
  virtual void init_map_phi(const FEGenericBase<OutputShape> & fe) const = 0;

  /**
   * Pre-requests any necessary data from FEMap
   */
  virtual void init_map_dphi(const FEGenericBase<OutputShape> & fe) const = 0;

  /**
   * Pre-requests any necessary data from FEMap
   */
  virtual void init_map_d2phi(const FEGenericBase<OutputShape> & fe) const = 0;

  /**
   * Evaluates shape functions in physical coordinates based on proper
   * finite element transformation.
   */
  virtual void map_phi(const unsigned int dim,
                       const Elem * const elem,
                       const std::vector<Point> & qp,
                       const FEGenericBase<OutputShape> & fe,
                       std::vector<std::vector<OutputShape> > & phi) const = 0;

  /**
   * Evaluates shape function gradients in physical coordinates based on proper
   * finite element transformation.
   */
  virtual void map_dphi(const unsigned int dim,
                        const Elem * const elem,
                        const std::vector<Point> & qp,
                        const FEGenericBase<OutputShape> & fe,
                        std::vector<std::vector<typename FEGenericBase<OutputShape>::OutputGradient> > & dphi,
                        std::vector<std::vector<OutputShape> > & dphidx,
                        std::vector<std::vector<OutputShape> > & dphidy,
                        std::vector<std::vector<OutputShape> > & dphidz) const = 0;

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
  /**
   * Evaluates shape function Hessians in physical coordinates based on proper
   * finite element transformation.
   */
  virtual void map_d2phi(const unsigned int dim,
                         const std::vector<Point> & qp,
                         const FEGenericBase<OutputShape> & fe,
                         std::vector<std::vector<typename FEGenericBase<OutputShape>::OutputTensor> > & d2phi,
                         std::vector<std::vector<OutputShape> > & d2phidx2,
                         std::vector<std::vector<OutputShape> > & d2phidxdy,
                         std::vector<std::vector<OutputShape> > & d2phidxdz,
                         std::vector<std::vector<OutputShape> > & d2phidy2,
                         std::vector<std::vector<OutputShape> > & d2phidydz,
                         std::vector<std::vector<OutputShape> > & d2phidz2) const = 0;
#endif //LIBMESH_ENABLE_SECOND_DERIVATIVES


  /**
   * Evaluates the shape function curl in physical coordinates based on proper
   * finite element transformation.
   */
  virtual void map_curl(const unsigned int dim,
                        const Elem * const elem,
                        const std::vector<Point> & qp,
                        const FEGenericBase<OutputShape> & fe,
                        std::vector<std::vector<OutputShape> > & curl_phi) const = 0;

  /**
   * Evaluates the shape function divergence in physical coordinates based on proper
   * finite element transformation.
   */
  virtual void map_div(const unsigned int dim,
                       const Elem * const elem,
                       const std::vector<Point> & qp,
                       const FEGenericBase<OutputShape> & fe,
                       std::vector<std::vector<typename FEGenericBase<OutputShape>::OutputDivergence> > & div_phi) const = 0;

}; // class FETransformationBase

}

#endif // LIBMESH_FE_TRANSFORMATION_BASE_H
